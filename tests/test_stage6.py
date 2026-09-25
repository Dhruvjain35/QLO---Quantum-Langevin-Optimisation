"""Stage 6 tests: global parity benchmark vs the frozen projector benchmark at finite shots."""

import numpy as np
import pennylane as qml
import pytest
from scipy.integrate import quad

from qlo.stage5.theory import a_k, f_plus_minus
from qlo.stage5.theory import s_k as proj_s_k
from qlo.stage6.analysis import (
    information_distances,
    matched_signal_table,
    parity_pairs,
    projector_pairs,
    wilson_ci,
)
from qlo.stage6.directional import (
    direction_probabilities,
    normal_approximation,
    shots_for_direction_exact,
    shots_for_direction_normal,
)
from qlo.stage6.exact_distribution import (
    binom_cdf,
    difference_probabilities,
    difference_probabilities_bruteforce,
    p_equal_zero_signal,
    p_equal_zero_signal_asymptotic,
)
from qlo.stage6.finite_shot import (
    estimator_from_counts,
    sample_gradient,
    shot_variance,
    shot_variance_from_probabilities,
    shots_for_snr,
    snr_squared,
)
from qlo.stage6.parity_benchmark import (
    E_LOG_ABS_COS,
    E_LOG_ABS_SIN,
    VAR_LOG_ABS_COS,
    b_k,
    cost,
    exact_gradient,
    exact_gradient_k,
    mu_z,
    s_k,
    shifted_mu,
    shifted_probabilities,
    theoretical_gradient_variance,
)


@pytest.fixture(scope="module")
def thetas():
    rng = np.random.default_rng(66)
    out = []
    for n in rng.integers(2, 9, 30):
        n = int(n)
        out.append((n, rng.uniform(-np.pi, np.pi, n), int(rng.integers(0, n))))
    return out


# ---- 1 parity cost closed form -------------------------------------------------------------------
def test_parity_cost_matches_statevector():
    rng = np.random.default_rng(1)
    for n in (1, 2, 3, 5):
        th = rng.uniform(-np.pi, np.pi, n)
        psi = np.array([1.0 + 0j])
        for t in th:
            psi = np.kron(psi, np.array([np.cos(t / 2), -1j * np.sin(t / 2)]))
        z = np.array([(-1) ** bin(i).count("1") for i in range(2**n)], dtype=float)
        mu = float(np.sum(np.abs(psi) ** 2 * z))
        assert abs(mu_z(th) - mu) < 1e-14
        assert abs(cost(th) - (1 - mu) / 2) < 1e-14
        assert 0.0 <= cost(th) <= 1.0


# ---- 2 exact gradient ----------------------------------------------------------------------------
def test_parity_exact_gradient(thetas):
    h = 1e-6
    for n, th, k in thetas:
        e = np.zeros(n); e[k] = 1.0
        fd = (cost(th + h * e) - cost(th - h * e)) / (2 * h)
        g = exact_gradient(th)
        assert abs(g[k] - fd) < 1e-8
        assert abs(g[k] - exact_gradient_k(s_k(th, k), b_k(th, k))) < 1e-15
    # prefix/suffix product handles an exact cos = 0 without division
    th = np.array([np.pi / 2, 0.3, -1.1])
    ref = np.array([0.5 * np.sin(th[k]) * np.prod(np.cos(np.delete(th, k))) for k in range(3)])
    assert np.allclose(exact_gradient(th), ref, atol=1e-16)


# ---- 3 theoretical gradient variance -------------------------------------------------------------
def test_parity_gradient_variance_theory():
    e_sin2 = quad(lambda t: np.sin(t) ** 2, -np.pi, np.pi)[0] / (2 * np.pi)
    e_cos2 = quad(lambda t: np.cos(t) ** 2, -np.pi, np.pi)[0] / (2 * np.pi)
    for n in (1, 2, 4, 7):
        assert abs(0.25 * e_sin2 * e_cos2 ** (n - 1) - theoretical_gradient_variance(n)) < 1e-15
        assert theoretical_gradient_variance(n) == 2.0 ** (-(n + 2))
    rng = np.random.default_rng(np.random.SeedSequence([6, 1, 4]))
    TH = rng.uniform(-np.pi, np.pi, (200_000, 4))
    g = 0.5 * s_k(TH, 0) * b_k(TH, 0)
    v = g.var(ddof=1)
    # E[g^4] = (1/16) E[sin^4] E[cos^4]^3 = (1/16)(3/8)^4 -> rel SE of the sample variance
    kurt = (1 / 16) * (3 / 8) ** 4 / theoretical_gradient_variance(4) ** 2
    se = theoretical_gradient_variance(4) * np.sqrt((kurt - 1) / TH.shape[0])
    assert abs(v - theoretical_gradient_variance(4)) < 5 * se
    assert abs(g.mean()) < 5 * np.sqrt(v / g.size)


# ---- 4 no structural-zero family -----------------------------------------------------------------
def test_no_structural_zero_family():
    for n in (2, 3, 4, 6, 8):
        for k in range(n):
            th = np.zeros(n); th[k] = np.pi / 2
            assert exact_gradient(th)[k] == 0.5
        TH = np.random.default_rng(n).uniform(-np.pi, np.pi, (500, n))
        G = np.stack([0.5 * s_k(TH, k) * b_k(TH, k) for k in range(n)], axis=1)
        assert np.all(np.max(np.abs(G), axis=0) > 1e-3)
        assert not np.any(G == 0.0)


# ---- 5 shifted probabilities ---------------------------------------------------------------------
def test_shifted_mu_and_probabilities(thetas):
    for n, th, k in thetas:
        s, B = s_k(th, k), b_k(th, k)
        e = np.zeros(n); e[k] = np.pi / 2
        mp, mm = shifted_mu(s, B)
        assert abs(mu_z(th + e) - mp) < 1e-14 and abs(mu_z(th - e) - mm) < 1e-14
        pp, pm = shifted_probabilities(s, B)
        assert abs(cost(th + e) - pp) < 1e-14 and abs(cost(th - e) - pm) < 1e-14
        assert abs(pp + pm - 1.0) < 1e-15
        assert abs((pp - pm) / 2 - exact_gradient(th)[k]) < 1e-15   # parameter-shift identity


# ---- 6 unbiased finite-shot estimator ------------------------------------------------------------
def test_parity_estimator_unbiased_exactly():
    from scipy.stats import binom

    for sB, M in ((0.3, 7), (-0.8, 12), (0.0, 5), (1e-3, 20)):
        pp, pm = shifted_probabilities(sB, 1.0)
        r = np.arange(M + 1)
        fp, fm = binom.pmf(r, M, pp), binom.pmf(r, M, pm)
        ghat = estimator_from_counts(r[:, None], r[None, :], M)
        mean = float(np.sum(fp[:, None] * fm[None, :] * ghat))
        var = float(np.sum(fp[:, None] * fm[None, :] * ghat**2)) - mean**2
        assert abs(mean - sB / 2) < 1e-14
        assert abs(var - shot_variance(sB, M)) < 1e-14   # 7 (exact, by enumeration)
    rng = np.random.default_rng(3)
    x = np.array([sample_gradient(0.4, 0.5, 64, rng) for _ in range(20_000)])
    assert abs(x.mean() - 0.1) < 5 * np.sqrt(shot_variance(0.2, 64) / x.size)


# ---- 7 variance formula --------------------------------------------------------------------------
def test_parity_variance_formula():
    for sB in (-1.0, -0.4, 0.0, 1e-8, 0.9, 1.0):
        pp, pm = shifted_probabilities(sB, 1.0)
        for M in (1, 16, 1000):
            assert abs(shot_variance(sB, M) - shot_variance_from_probabilities(pp, pm, M)) < 1e-16
            g = sB / 2
            assert abs(shot_variance(sB, M) - (1 - 4 * g**2) / (8 * M)) < 1e-16


# ---- 8 SNR formula and required shots ------------------------------------------------------------
def test_parity_snr_and_required_shots():
    for sB in (0.5, -0.1, 1e-4):
        for M in (1, 64, 1e6):
            var = shot_variance(sB, M)
            assert abs(snr_squared(sB, M) - (sB / 2) ** 2 / var) / snr_squared(sB, M) < 1e-12
        for rho in (1.0, 2.0):
            M = shots_for_snr(sB, rho)
            assert abs(snr_squared(sB, M) - rho**2) < 1e-9 * rho**2
    assert snr_squared(0.0, 100) == 0.0 and shots_for_snr(0.0, 1.0) == np.inf
    assert snr_squared(1.0, 100) == np.inf and shots_for_snr(1.0, 1.0) == 0.0


# ---- 9 exact P_zero vs brute force at small M ----------------------------------------------------
def test_difference_probabilities_vs_bruteforce():
    rng = np.random.default_rng(9)
    for _ in range(60):
        pa, pb = rng.uniform(0, 1, 2) * rng.choice([1.0, 1e-2])
        M = int(rng.integers(1, 80))
        got = difference_probabilities(pa, pb, M)
        ref = difference_probabilities_bruteforce(pa, pb, M)
        assert max(abs(float(g) - r) for g, r in zip(got, ref)) < 1e-12
    # the CDF helper
    from scipy.stats import binom

    for M, p, r in ((10, 0.3, 4), (50, 0.9, 49), (7, 0.5, -1), (7, 0.5, 7)):
        assert abs(float(binom_cdf(M, p, r)) - binom.cdf(r, M, p)) < 1e-13


# ---- 10 central-binomial identity at zero signal -------------------------------------------------
def test_zero_signal_central_binomial():
    from math import comb

    for M in (1, 2, 5, 17, 60):
        _, eq, _ = difference_probabilities(0.5, 0.5, M)
        assert abs(float(eq) - comb(2 * M, M) / 4**M) < 1e-13
        assert abs(float(p_equal_zero_signal(M)) - comb(2 * M, M) / 4**M) < 1e-12


# ---- 11 large-M 1/sqrt(pi M) ---------------------------------------------------------------------
def test_zero_signal_asymptotic():
    for M in (1e3, 1e4, 1e6):
        _, eq, _ = difference_probabilities(0.5, 0.5, M)
        rel = abs(float(eq) / float(p_equal_zero_signal_asymptotic(M)) - 1)
        assert rel < 1 / (4 * M)   # C(2M,M)/4^M = (1 - 1/(8M) + ...)/sqrt(pi M)
        assert rel > 0


# ---- 12 probabilities sum to one -----------------------------------------------------------------
def test_direction_probabilities_sum_to_one():
    rng = np.random.default_rng(12)
    pa, pb = rng.uniform(0, 1, 200), rng.uniform(0, 1, 200)
    M = rng.integers(1, 3000, 200).astype(float)
    d = direction_probabilities(pa, pb, M)
    s = d["p_correct"] + d["p_wrong"] + d["p_zero"]
    assert np.max(np.abs(s - 1)) < 1e-10
    d0 = direction_probabilities(0.5, 0.5, 100.0)
    assert bool(d0["zero_signal"]) and abs(float(d0["p_correct"] - d0["p_wrong"])) < 1e-15


# ---- 13 sign swap symmetry -----------------------------------------------------------------------
def test_sign_swap_symmetry():
    for sB in (0.3, 0.01, 1e-4):
        for M in (16, 1000):
            a = direction_probabilities(*shifted_probabilities(sB, 1.0), M)
            b = direction_probabilities(*shifted_probabilities(-sB, 1.0), M)
            for key in ("p_correct", "p_wrong", "p_zero"):
                assert abs(float(a[key]) - float(b[key])) < 1e-10


# ---- 14 normal approximation improves with M ----------------------------------------------------
def test_normal_approximation_improves():
    # fixed SNR ~ 1: the discreteness (tie) error of Phi(SNR) shrinks like 1/sqrt(M)
    errs = []
    for M in (16, 256, 4096, 65536):
        sB = float(np.sqrt(1.0 / (2 * M)))  # SNR^2 = 2 M sB^2/(1-sB^2) ~ 1
        pp, pm = shifted_probabilities(sB, 1.0)
        exact = float(direction_probabilities(pp, pm, M)["p_correct"])
        approx = float(normal_approximation(pp, pm, M, continuity=True)["p_correct_normal"])
        errs.append(abs(exact - approx))
    assert all(errs[i + 1] < errs[i] for i in range(len(errs) - 1))
    assert errs[-1] < 1e-3
    # exact inversion vs normal inversion agree to within a few percent at large M
    sB = 1e-3
    pp, pm = shifted_probabilities(sB, 1.0)
    m_ex, ok = shots_for_direction_exact(pp, pm, 0.75)
    assert bool(ok)
    assert abs(np.log10(float(m_ex)) - np.log10(float(shots_for_direction_normal(sB, 0.75)))) < 0.05


# ---- 15/16 log moments ---------------------------------------------------------------------------
def test_log_moment_constants():
    e_cos = quad(lambda t: np.log(abs(np.cos(t))), -np.pi, np.pi, limit=500, points=[-np.pi / 2, np.pi / 2])[0] / (2 * np.pi)
    e_sin = quad(lambda t: np.log(abs(np.sin(t))), -np.pi, np.pi, limit=500, points=[0.0])[0] / (2 * np.pi)
    v_cos = quad(lambda t: (np.log(abs(np.cos(t))) + np.log(2)) ** 2, -np.pi, np.pi, limit=500,
                 points=[-np.pi / 2, np.pi / 2])[0] / (2 * np.pi)
    assert abs(e_cos - E_LOG_ABS_COS) < 1e-9 and abs(E_LOG_ABS_COS + np.log(2)) < 1e-15
    assert abs(e_sin - E_LOG_ABS_SIN) < 1e-9
    assert abs(v_cos - VAR_LOG_ABS_COS) < 1e-7


# ---- 17 matched-signal bins are never mixed -----------------------------------------------------
def test_matched_signal_bins_consistent():
    ms = matched_signal_table((4, 6), 3000, (64,), seed=1, min_count=50)
    assert len(ms) > 0
    for (lo, hi), grp in ms.groupby(["bin_low", "bin_high"]):
        assert set(grp.benchmark) == {"parity", "projector"}    # a bin is reported for both or neither
        assert np.all((np.log10(grp.mean_abs_g) >= lo - 0.5) & (np.log10(grp.mean_abs_g) <= hi + 0.5))
    assert np.all(ms.n_in_bin >= 50)
    # the pair helpers carry the right sign conventions
    th = np.random.default_rng(0).uniform(-np.pi, np.pi, (100, 5))
    for fn in (parity_pairs, projector_pairs):
        pp, pn, g = fn(th, 0)
        assert np.allclose((pp - pn) / 2, g, atol=1e-16)
    A, s = a_k(th, 0), proj_s_k(th, 0)
    fp, fm = f_plus_minus(A, s)
    assert np.allclose(projector_pairs(th, 0)[0], fm)


# ---- 18 PennyLane parity statistics vs binomial model -------------------------------------------
def test_pennylane_parity_matches_binomial():
    n, M, reps = 3, 32, 300
    th = np.array([0.7, -0.4, 1.9])
    dev = qml.device("default.qubit", wires=n, seed=np.random.default_rng(18))
    obs = qml.prod(*[qml.PauliZ(j) for j in range(n)])

    @qml.set_shots(shots=M)
    @qml.qnode(dev, diff_method=None)
    def ev(t):
        for j in range(n):
            qml.RX(t[j], wires=j)
        return qml.expval(obs)

    e = np.zeros(n); e[0] = np.pi / 2
    samples = np.array([0.5 * ((1 - float(ev(th + e))) / 2 - (1 - float(ev(th - e))) / 2) for _ in range(reps)])
    s, B = s_k(th, 0), b_k(th, 0)
    assert np.allclose(samples * 2 * M, np.round(samples * 2 * M))    # lattice j/(2M)
    assert abs(samples.mean() - s * B / 2) < 5 * np.sqrt(shot_variance(s * B, M) / reps)
    ratio = samples.var(ddof=1) / shot_variance(s * B, M)
    assert 0.7 < ratio < 1.3
    d = direction_probabilities(*shifted_probabilities(s, B), M)
    for key, emp in (("p_zero", np.mean(samples == 0)), ("p_correct", np.mean(np.sign(samples) == np.sign(s * B)))):
        lo, hi = wilson_ci(int(round(emp * reps)), reps)
        assert lo - 0.02 <= float(d[key]) <= hi + 0.02


# ---- information distances ----------------------------------------------------------------------
def test_information_distances():
    d = information_distances(np.array([0.5, 0.2]), np.array([0.5, 0.1]))
    assert d["tv"][0] == 0 and d["hellinger2"][0] == 0 and d["kl"][0] == 0
    kl = 0.2 * np.log(2) + 0.8 * np.log(0.8 / 0.9)
    assert abs(d["kl"][1] - kl) < 1e-15
    assert 0 < d["hellinger2"][1] <= d["tv"][1]


# ---- regression: binomial CDF must not form 1 - p (p ~ 1e-18 rounds 1 - p to 1) ----------------
def test_binom_cdf_tiny_p_huge_M():
    from scipy.stats import poisson

    for M, p, r in ((1e20, 5e-18, 300.0), (1e24, 3e-22, 250.0), (1e18, 1e-16, 90.0)):
        ref = poisson.cdf(r, M * p)   # Poisson limit is exact to O(p) here
        assert abs(float(binom_cdf(M, p, r)) / ref - 1) < 1e-9
    gt, eq, lt = difference_probabilities(5e-17 * 1.065 / 2, 5e-17 * 0.935 / 2, 2e19)
    assert 0.97 < float(gt) < 0.99 and abs(float(gt + eq + lt) - 1) < 1e-10


# ---- the exact parity direction table is an upper bound that tracks per-sample exact bisection ----
def test_parity_direction_table_upper_bound():
    from qlo.stage6.scaling import parity_direction_table

    g, gl = parity_direction_table(0.75, -1.5, n_grid=300)
    assert np.all(np.diff(gl) <= 1e-12)                        # non-increasing in |sB|
    l = np.random.default_rng(4).uniform(-1.49, -0.01, 60)
    j = np.searchsorted(g, l, side="right") - 1
    pp, pm = shifted_probabilities(10.0**l, 1.0)
    m_ex, ok = shots_for_direction_exact(pp, pm, 0.75, lo_log10=0.0, hi_log10=5.0, n_iter=30)
    assert np.all(ok)
    d = gl[j] - np.log10(m_ex)
    assert np.all(d >= -1e-12) and np.all(d < 0.2)
