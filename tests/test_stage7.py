"""Stage 7 tests: one fidelity landscape, Loschmidt (projector) vs SWAP-test finite-shot gradients."""

import numpy as np
import pytest
from scipy.stats import binom

from qlo.stage5.theory import a_k, f_plus_minus, log_a_k, s_k
from qlo.stage5.theory import shot_variance as stage5_shot_variance
from qlo.stage5.theory import snr_squared as stage5_snr_squared
from qlo.stage6.directional import direction_probabilities
from qlo.stage6.exact_distribution import difference_probabilities, difference_probabilities_bruteforce, p_equal_zero_signal
from qlo.stage7 import analysis as AN
from qlo.stage7 import estimators as E
from qlo.stage7 import required_shots as R


def _F(theta):
    return float(np.prod(np.cos(np.asarray(theta) / 2.0) ** 2))


@pytest.fixture(scope="module")
def cases():
    rng = np.random.default_rng(77)
    out = []
    for _ in range(30):
        n = int(rng.integers(2, 9))
        out.append((n, rng.uniform(-np.pi, np.pi, n), int(rng.integers(0, n))))
    return out


def _moments(scheme, A, s, M):
    pp, pm = E.outcome_probabilities(scheme, A, s)
    r = np.arange(M + 1)
    W = binom.pmf(r, M, pp)[:, None] * binom.pmf(r, M, pm)[None, :]
    G = E.gradient_from_counts(scheme, r[:, None], r[None, :], M)
    mean = float(np.sum(W * G))
    return mean, float(np.sum(W * G * G)) - mean**2


# 1 / 13: the common landscape and the exact gradient are identical for both schemes
def test_common_landscape_identical(cases):
    for n, th, k in cases:
        A, s = a_k(th, k), s_k(th, k)
        g = A * s / 2
        for sc in E.SCHEMES:
            p_pos, p_neg, c = E.pair(sc, A, s)
            assert abs(c * (p_pos - p_neg) - g) < 1e-15
            assert abs(c * E.pair_delta(sc, A, s) - g) < 1e-17
    chk = AN.common_landscape_check(10, 5000, 0)
    assert max(v for k, v in chk.items() if k.endswith("minus_g")) < 1e-15


# 2: F_plus / F_minus against the circuit fidelity at the shifted points
def test_shifted_fidelities(cases):
    for n, th, k in cases:
        e = np.zeros(n); e[k] = np.pi / 2
        fp, fm = f_plus_minus(a_k(th, k), s_k(th, k))
        assert abs(_F(th + e) - fp) < 1e-15 and abs(_F(th - e) - fm) < 1e-15
        h = np.zeros(n); h[k] = 1e-6
        assert abs(-(_F(th + h) - _F(th - h)) / 2e-6 - a_k(th, k) * s_k(th, k) / 2) < 1e-8


# 3 / 7: Loschmidt unbiased with the Stage 5 variance (exact enumeration)
def test_loschmidt_unbiased_and_variance():
    for A, s, M in ((0.8, 0.3, 7), (0.05, -0.7, 12), (1e-3, 0.9, 25)):
        mean, var = _moments("loschmidt", A, s, M)
        assert abs(mean - A * s / 2) < 1e-15
        assert abs(var - float(E.gradient_variance("loschmidt", A, s, M))) < 1e-15
        assert abs(float(E.gradient_variance("loschmidt", A, s, M)) - float(stage5_shot_variance(A, s, M))) < 1e-17


# 4: SWAP estimator algebra and SIGN, built from the C_hat definitions
def test_swap_algebra_and_sign():
    rng = np.random.default_rng(4)
    for _ in range(50):
        M = int(rng.integers(1, 500)); kp, km = rng.integers(0, M + 1, 2)
        c_plus, c_minus = 1 - (2 * kp / M - 1), 1 - (2 * km / M - 1)
        assert abs(float(E.gradient_from_counts("swap", kp, km, M)) - 0.5 * (c_plus - c_minus)) < 1e-15
        assert abs(float(E.gradient_from_counts("swap", kp, km, M)) - (km - kp) / M) < 1e-15
    # the sign: a larger minus-shift count means a positive gradient, as dC/dtheta = (C_+ - C_-)/2 requires
    th = np.array([0.4, -1.1, 2.0]); A, s = a_k(th, 0), s_k(th, 0)
    q_plus, q_minus = E.outcome_probabilities("swap", A, s)
    assert np.sign(q_minus - q_plus) == np.sign(A * s / 2)


# 5 / 6: SWAP unbiased, analytic variance, Var(F_hat) = (1 - F^2)/M
def test_swap_unbiased_and_variance():
    for A, s, M in ((0.8, 0.3, 7), (0.05, -0.7, 12), (0.5, 0.0, 9), (1e-3, 0.9, 25)):
        mean, var = _moments("swap", A, s, M)
        assert abs(mean - A * s / 2) < 1e-15
        closed = (2 - A**2 * (1 + s**2) / 2) / (4 * M)
        assert abs(var - closed) < 1e-15 and abs(float(E.gradient_variance("swap", A, s, M)) - closed) < 1e-17
    for F in (0.0, 0.3, 0.99):
        q = (1 + F) / 2; M = 11
        r = np.arange(M + 1); w = binom.pmf(r, M, q); fh = 2 * r / M - 1
        assert abs(np.sum(w * fh) - F) < 1e-14
        assert abs(np.sum(w * fh**2) - F**2 - (1 - F**2) / M) < 1e-14


# 8 / 9: exact SNR for both schemes (and Stage 5 agreement for Loschmidt)
def test_snr_formulas():
    for A, s in ((0.7, 0.4), (1e-3, -0.9), (1e-8, 0.5)):
        for M in (1.0, 64.0, 1e6):
            for sc in E.SCHEMES:
                snr2 = float(E.snr_squared(sc, A, s, M))
                assert abs(snr2 - (A * s / 2) ** 2 / float(E.gradient_variance(sc, A, s, M))) <= 1e-12 * snr2
            assert abs(float(E.snr_squared("loschmidt", A, s, M)) - float(stage5_snr_squared(A, s, M))) <= 1e-12 * float(stage5_snr_squared(A, s, M))
    A, s, M = 1e-8, 0.6, 1e3
    assert abs(float(E.snr_squared("swap", A, s, M)) / (M * A**2 * s**2 / 2) - 1) < 1e-12
    assert abs(float(E.snr_squared("loschmidt", A, s, M)) / (M * A * s**2) - 1) < 1e-7


# 10: SWAP P_zero vs brute force
def test_swap_p_zero_bruteforce():
    rng = np.random.default_rng(10)
    for _ in range(40):
        A, s, M = rng.uniform(0, 1), rng.uniform(-1, 1), int(rng.integers(1, 60))
        pp, pn, _ = E.pair("swap", A, s)
        got = difference_probabilities(pp, pn, M); ref = difference_probabilities_bruteforce(float(pp), float(pn), M)
        assert max(abs(float(a) - b) for a, b in zip(got, ref)) < 1e-12


# 11: deep limit -> central binomial for SWAP, -> 1 for Loschmidt
def test_deep_limit_zero_probability():
    from math import comb

    for M in (1, 7, 40, 200):
        pp, pn, _ = E.pair("swap", 1e-12, 0.5)
        assert abs(float(difference_probabilities(pp, pn, M)[1]) - comb(2 * M, M) / 4**M) < 1e-10
    pp, pn, _ = E.pair("loschmidt", 1e-9, 0.5)
    assert float(difference_probabilities(pp, pn, 1000)[1]) > 0.999998
    assert abs(float(p_equal_zero_signal(1e4)) - 1 / np.sqrt(np.pi * 1e4)) < 1e-6


# 12: probabilities sum to one, both schemes
def test_direction_probabilities_sum():
    th = AN.sample_theta(6, 300, 0, 3)
    A, s = np.exp(log_a_k(th, 0)), s_k(th, 0)
    for sc in E.SCHEMES:
        pp, pn, _ = E.pair(sc, A, s)
        for M in (5.0, 300.0):
            d = direction_probabilities(pp, pn, M)
            assert np.max(np.abs(d["p_correct"] + d["p_wrong"] + d["p_zero"] - 1)) < 1e-10


# 14: paired theta handling - both schemes see the identical theta, A, s and exact gradient
def test_paired_theta():
    a = AN.required_shots_for_n(4, 300, 0, keep=300)
    assert {r["scheme"] for r in a["rows"]} == set(E.SCHEMES)
    r0, r1 = a["rows"]
    assert r0["empirical_mean_logA"] == r1["empirical_mean_logA"] and r0["log10_abs_g_median"] == r1["log10_abs_g_median"]
    th1, th2 = AN.sample_theta(5, 10, 0, 20), AN.sample_theta(5, 10, 0, 20)
    assert np.array_equal(th1, th2)
    runs, _ = AN.optimization_diagnostic((4,), [("exact", None, 0), ("LE", "loschmidt", 16), ("SWAP", "swap", 16)], 3, 2, 0.3, 0)
    f0 = runs.pivot_table(index="seed", columns="variant", values="F0")
    assert np.allclose(f0["exact"], f0["LE"]) and np.allclose(f0["exact"], f0["SWAP"])


# 15: required-shot inversion reaches the target (and M - 1 does not)
def test_required_shot_inversion():
    th = AN.sample_theta(4, 25, 0, 99)
    logA, s = log_a_k(th, 0), s_k(th, 0)
    A = np.exp(logA)
    for sc in E.SCHEMES:
        l, m = R.direction_shots(sc, logA, s, 0.75)
        pp, pn, _ = E.pair(sc, A, s)
        ex = m == R.METHOD["exact"]
        M = np.round(10.0 ** l[ex])
        assert np.all(direction_probabilities(pp[ex], pn[ex], M)["p_correct"] >= 0.75 - 1e-12)
        below = M > 1
        assert np.all(direction_probabilities(pp[ex][below], pn[ex][below], M[below] - 1)["p_correct"] < 0.75)
        l, _ = R.nonzero_shots(sc, logA, s, 0.9)
        M = np.round(10.0 ** l)
        assert np.all(1 - difference_probabilities(pp, pn, M)[1] >= 0.9 - 1e-12)
        for rho in (1.0, 2.0):
            M = 10.0 ** R.snr_shots(sc, logA, s, rho)
            assert np.all(E.snr_squared(sc, A, s, np.round(M)) >= rho**2 * (1 - 1e-9))


# 16 / 17: PennyLane SWAP test: mean = fidelity; finite-shot gradient consistent with the binomial model
def test_pennylane_swap():
    from qlo.stage7.pennylane_check import make_swap, swap_expval

    th = np.array([0.5, -0.8, 1.7])
    assert abs(swap_expval(3, th) - _F(th)) < 1e-12
    M, reps = 64, 150
    count = make_swap(3, M, np.random.default_rng(16))
    e = np.array([np.pi / 2, 0, 0])
    kp = np.array([count(th + e) for _ in range(reps)]); km = np.array([count(th - e) for _ in range(reps)])
    g = E.gradient_from_counts("swap", kp, km, M)
    A, s = a_k(th, 0), s_k(th, 0)
    var = float(E.gradient_variance("swap", A, s, M))
    assert abs(g.mean() - A * s / 2) < 5 * np.sqrt(var / reps)
    assert 0.7 < g.var(ddof=1) / var < 1.35
    assert np.allclose(g * M, np.round(g * M))  # lattice j/M


# 18: resource accounting separates sample count from state-copy cost
def test_resource_accounting():
    le, sw = E.resource_accounting("loschmidt", 8, 1024), E.resource_accounting("swap", 8, 1024)
    assert le["measurement_shots"] == sw["measurement_shots"] == 2048
    assert sw["total_state_copies"] == 2 * le["total_state_copies"]
    assert sw["qubits_per_shot"] == 17 and le["qubits_per_shot"] == 8 and sw["controlled_swaps_per_shot"] == 8


# information distances: cancellation-free for SWAP-like probabilities
def test_bernoulli_distances_stable():
    d = 1e-13
    out = AN.bernoulli_distances(0.5 + d / 2, 0.5 - d / 2, d)
    assert abs(out["hellinger2"] / (d**2 / 2) - 1) < 1e-6   # H^2 ~ delta^2 / (8 p(1-p)) at p = 1/2
    assert abs(out["kl"] / (2 * d**2) - 1) < 1e-6
    assert out["tv"] == d
