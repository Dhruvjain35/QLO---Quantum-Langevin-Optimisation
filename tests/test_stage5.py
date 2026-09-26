"""Stage 5 tests: exact finite-shot estimator distribution on the projector benchmark."""

import numpy as np
import pennylane as qml
import pytest
from scipy.integrate import quad

from qlo.stage4.landscape import exact_gradient, shifted_fidelities
from qlo.stage4.seeds import initial_theta
from qlo.stage5.analysis import wilson_ci
from qlo.stage5.theory import (
    E_LOG_COS2_HALF,
    a_k,
    arithmetic_mean_a,
    e_log_a,
    estimator_from_counts,
    exact_gradient_k,
    f_plus_minus,
    geometric_typical_a,
    local_estimator_from_counts,
    local_gradient_k,
    local_p_plus_minus,
    local_shots_for_snr,
    local_snr_squared,
    local_variance,
    p_both_zero,
    p_zero_bruteforce,
    p_zero_exact,
    p_zero_poisson_limit,
    s_k,
    shot_variance,
    shot_variance_from_F,
    shots_for_nonzero_exact,
    shots_for_snr,
    snr_squared,
)


@pytest.fixture(scope="module")
def thetas():
    rng = np.random.default_rng(55)
    return [(int(n), rng.uniform(-np.pi, np.pi, int(n)), int(k)) for n, k in zip(rng.integers(2, 9, 30), rng.integers(0, 2, 30))]


# ---- 1-4 algebraic identities against the Stage 4 landscape ------------------------------------
def test_F_plus_minus_and_sum_and_gradient(thetas):
    for n, th, k in thetas:
        k = k % n
        A, s = a_k(th, k), s_k(th, k)
        fp, fm = shifted_fidelities(th)
        Fp, Fm = f_plus_minus(A, s)
        assert abs(fp[k] - Fp) < 1e-15 and abs(fm[k] - Fm) < 1e-15          # 1, 2
        assert abs(Fp + Fm - A) < 1e-15                                        # 3
        assert abs(exact_gradient(th)[k] - exact_gradient_k(A, s)) < 1e-15     # 4
        assert abs(Fp**2 + Fm**2 - A**2 * (1 + s**2) / 2) < 1e-15


# ---- 5 estimator identity ------------------------------------------------------------------------
def test_estimator_identity_from_counts():
    rng = np.random.default_rng(1)
    for _ in range(50):
        M = int(rng.integers(1, 200)); kp, km = rng.integers(0, M + 1, 2)
        c_plus, c_minus = 1 - kp / M, 1 - km / M
        assert estimator_from_counts(kp, km, M) == pytest.approx(0.5 * (c_plus - c_minus), abs=1e-15)
        assert estimator_from_counts(kp, km, M) == pytest.approx((km - kp) / (2 * M), abs=1e-15)
    assert estimator_from_counts(3, 3, 7) == 0.0


# ---- 6-9 P_zero ---------------------------------------------------------------------------------------
def test_p_zero_matches_bruteforce_small_M():
    rng = np.random.default_rng(2)
    for _ in range(100):
        fp, fm = rng.uniform(0, 1, 2) * rng.choice([1.0, 1e-2, 1e-4]); M = int(rng.integers(1, 40))
        assert abs(float(p_zero_exact(fp, fm, M)) - p_zero_bruteforce(fp, fm, M)) < 1e-12
    assert float(p_zero_exact(0.0, 0.0, 10)) == 1.0 and float(p_zero_exact(1.0, 1.0, 10)) == 1.0
    assert float(p_zero_exact(1.0, 0.0, 10)) == 0.0


def test_p_zero_bounds_and_range():
    rng = np.random.default_rng(3)
    A = 10.0 ** rng.uniform(-12, 0, 300); s = rng.uniform(-1, 1, 300); M = 10.0 ** rng.uniform(0, 12, 300)
    fp, fm = f_plus_minus(A, s)
    pz, pb = p_zero_exact(fp, fm, M), p_both_zero(fp, fm, M)
    assert np.all((pz >= 0) & (pz <= 1)) and np.all(pz >= pb - 1e-12)      # 7, 8
    assert np.all(np.isfinite(pz))


def test_p_zero_to_one_as_A_to_zero_at_fixed_M():
    for M in (1, 100, 10_000):
        vals = [float(p_zero_exact(*f_plus_minus(A, 0.3), M)) for A in (1e-1, 1e-3, 1e-6, 1e-9, 1e-15)]
        assert all(np.diff(vals) >= -1e-12) and vals[-1] > 1 - 1e-9        # 9


def test_p_zero_stable_for_huge_M_and_matches_poisson_limit():
    for A, s, M in ((1e-20, 0.3, 1e20), (1e-30, -0.9, 3e29), (1e-12, 0.0, 1e12), (1e-15, 0.5, 1e17)):
        ex = float(p_zero_exact(*f_plus_minus(A, s), M))
        assert np.isfinite(ex) and abs(ex - float(p_zero_poisson_limit(A, s, M))) < 1e-9


def test_required_shots_exact_is_minimal():
    rng = np.random.default_rng(4)
    A = 10.0 ** rng.uniform(-6, 0, 40); s = rng.uniform(-1, 1, 40)
    fp, fm = f_plus_minus(A, s)
    for q in (0.5, 0.9):
        M, ok = shots_for_nonzero_exact(fp, fm, q)
        assert ok.all() and np.all(M >= 1)
        assert np.all(1 - p_zero_exact(fp, fm, M) >= q - 1e-12)
        below = M > 1
        assert np.all(1 - p_zero_exact(fp[below], fm[below], M[below] - 1) < q + 1e-12)


# ---- 10-12 variance and SNR ----------------------------------------------------------------------------
def test_conditional_variance_simplification(thetas):
    for n, th, k in thetas:
        k = k % n; A, s = a_k(th, k), s_k(th, k); fp, fm = f_plus_minus(A, s)
        for M in (1, 37, 4096):
            assert shot_variance(A, s, M) == pytest.approx(shot_variance_from_F(fp, fm, M), rel=1e-12)   # 10


def test_snr_formula_and_required_shots(thetas):
    for n, th, k in thetas:
        k = k % n; A, s = a_k(th, k), s_k(th, k)
        M = 500
        assert snr_squared(A, s, M) == pytest.approx((A * s / 2) ** 2 / shot_variance(A, s, M), rel=1e-12)   # 11
        for rho in (0.5, 1.0, 2.0, 5.0):
            Ms = shots_for_snr(A, s, rho)
            assert snr_squared(A, s, Ms) == pytest.approx(rho**2, rel=1e-10)                              # 12
    assert np.isinf(shots_for_snr(0.5, 0.0, 1.0)) and snr_squared(0.5, 0.0, 10) == 0.0
    assert np.isinf(shots_for_snr(0.0, 0.5, 1.0)) and snr_squared(0.0, 0.5, 10) == 0.0
    assert shots_for_snr(1.0, 1.0, 1.0) == 0.0  # deterministic estimator (F+ = 0, F- = 1)


# ---- 13-14 typical scaling ------------------------------------------------------------------------------
def test_E_log_cos2_half_and_E_A():
    e = quad(lambda t: np.log(np.cos(t / 2) ** 2), -np.pi, np.pi, limit=200)[0] / (2 * np.pi)
    assert e == pytest.approx(-2 * np.log(2), abs=1e-8) == pytest.approx(E_LOG_COS2_HALF, abs=1e-8)     # 13
    ec2 = quad(lambda t: np.cos(t / 2) ** 2, -np.pi, np.pi)[0] / (2 * np.pi)
    assert ec2 == pytest.approx(0.5, abs=1e-12)
    for n in (2, 5, 11):
        assert arithmetic_mean_a(n) == pytest.approx(0.5 ** (n - 1)) and geometric_typical_a(n) == pytest.approx(4.0 ** (-(n - 1)))  # 14
        assert e_log_a(n) == pytest.approx(-2 * (n - 1) * np.log(2))
        assert geometric_typical_a(n) < arithmetic_mean_a(n) or n == 1
    # Monte Carlo check of both moments at n=6 (SE of mean A ~ 0.03*mean; SE of mean log A ~ 0.02)
    rng = np.random.default_rng(6); th = rng.uniform(-np.pi, np.pi, (40_000, 6))
    A = a_k(th, 0)
    assert abs(np.log(A).mean() - e_log_a(6)) < 0.05 and abs(A.mean() / arithmetic_mean_a(6) - 1) < 0.05


# ---- 15-17 local control -------------------------------------------------------------------------------------
def test_local_gradient_variance_and_snr():
    rng = np.random.default_rng(7)
    for _ in range(30):
        n = int(rng.integers(2, 12)); th = rng.uniform(-np.pi, np.pi, n); k = int(rng.integers(0, n)); M = int(rng.integers(1, 3000))
        s = np.sin(th[k])
        # 15: gradient of C_L = 1 - mean cos^2(theta_j/2)
        h = 1e-6; e = np.zeros(n); e[k] = 1
        cl = lambda t: 1 - np.mean(np.cos(t / 2) ** 2)
        assert local_gradient_k(s, n) == pytest.approx((cl(th + h * e) - cl(th - h * e)) / (2 * h), abs=1e-8)
        pp, pm = local_p_plus_minus(s)
        assert pp + pm == pytest.approx(1.0) and 0 <= pp <= 1
        assert local_estimator_from_counts(3, 5, M, n) == pytest.approx((5 - 3) / (2 * M * n))
        # 16: variance from binomial counts equals cos^2 / (8 M n^2)
        v = (pp * (1 - pp) + pm * (1 - pm)) / (4 * M * n**2)
        assert local_variance(s, M, n) == pytest.approx(v, rel=1e-12) == pytest.approx(np.cos(th[k]) ** 2 / (8 * M * n**2), rel=1e-12)
        # 17: SNR^2 = 2 M tan^2, independent of n
        if abs(np.cos(th[k])) > 1e-6:
            assert local_snr_squared(s, M) == pytest.approx(local_gradient_k(s, n) ** 2 / v, rel=1e-10)
            assert local_snr_squared(s, M) == pytest.approx(2 * M * np.tan(th[k]) ** 2, rel=1e-10)
            assert local_snr_squared(s, local_shots_for_snr(s, 2.0)) == pytest.approx(4.0, rel=1e-10)
    for n1, n2 in ((2, 20), (3, 300)):
        assert local_snr_squared(0.4, 77) == local_snr_squared(0.4, 77)  # trivially n-free: function has no n argument


# ---- 18 binomial Monte Carlo vs P_zero (seeded, tolerant) -----------------------------------------------------
def test_monte_carlo_matches_p_zero_and_variance():
    rng = np.random.default_rng(8)
    for A, s, M in ((1e-2, 0.5, 100), (1e-3, -0.8, 2000), (0.3, 0.1, 8), (1e-4, 0.9, 30_000)):
        fp, fm = f_plus_minus(A, s)
        g = (rng.binomial(M, fm, 40_000) - rng.binomial(M, fp, 40_000)) / (2 * M)
        k = int(np.sum(g == 0)); lo, hi = wilson_ci(k, g.size)
        ex = float(p_zero_exact(fp, fm, M))
        assert lo - 0.01 <= ex <= hi + 0.01
        assert abs(g.var(ddof=1) / shot_variance(A, s, M) - 1) < 0.1
        assert abs(g.mean() - A * s / 2) < 5 * np.sqrt(shot_variance(A, s, M) / g.size)


# ---- 19 PennyLane small-system agreement ------------------------------------------------------------------------
def test_pennylane_projector_shots_match_binomial_model():
    n, M, R = 3, 32, 300
    th = 0.7 * initial_theta(3, n)
    A, s = float(a_k(th, 0)), float(s_k(th, 0)); fp, fm = f_plus_minus(A, s)
    dev = qml.device("default.qubit", wires=n, seed=np.random.default_rng(0))

    @qml.set_shots(shots=M)
    @qml.qnode(dev, diff_method=None)
    def fid(t):
        for j in range(n):
            qml.RX(t[j], wires=j)
        return qml.expval(qml.Projector(np.zeros(n, dtype=int), wires=range(n)))

    e = np.zeros(n); e[0] = np.pi / 2
    g = np.array([0.5 * ((1 - float(fid(th + e))) - (1 - float(fid(th - e)))) for _ in range(R)])
    assert np.allclose(g * 2 * M, np.round(g * 2 * M))  # support {j/(2M)}
    lo, hi = wilson_ci(int(np.sum(g == 0)), R)
    assert lo - 0.02 <= float(p_zero_exact(fp, fm, M)) <= hi + 0.02
    assert abs(g.mean() - A * s / 2) < 5 * np.sqrt(shot_variance(A, s, M) / R)
    assert 0.6 < g.var(ddof=1) / shot_variance(A, s, M) < 1.5
