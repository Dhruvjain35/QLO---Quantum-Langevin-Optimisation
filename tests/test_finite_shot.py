"""Stage 2 tests: finite-shot parameter-shift estimator + Task 10 sanity checks.

Monte-Carlo tests use fixed seeds (deterministic) AND tolerances wide enough
that they would pass for the overwhelming majority of seeds.
"""

import ast
import inspect

import numpy as np
import pennylane as qml
import pytest

import qlo.gradients.exact as exact_mod
import qlo.gradients.finite_shot as fs_mod
from qlo import HardwareEfficientAnsatz, cost_value, global_z_observable, gradient_autograd, make_rng
from qlo.analysis import theoretical_parameter_shift_variance, theoretical_parameter_shift_variance_vector
from qlo.analysis.statistics import summarize_shot_noise
from qlo.gradients.finite_shot import SHIFT, FiniteShotParameterShift
from qlo.utils.indexing import flat_to_multi, multi_to_flat, shifted_params


@pytest.fixture(scope="module")
def system():
    a = HardwareEfficientAnsatz(3, 2, "ring")
    theta = a.init_params(make_rng(7))
    g = gradient_autograd(a, "global", theta).ravel()
    return a, theta, g


# ---- 1. shape -----------------------------------------------------------------------
def test_finite_shot_gradient_shape_matches_exact(system):
    a, theta, g = system
    est = FiniteShotParameterShift(a, "global", shots=16, master_seed=0)
    G = est.gradient(theta)
    assert G.shape == theta.shape == gradient_autograd(a, "global", theta).shape
    assert np.all(np.isfinite(G))
    assert est.n_evaluations == 2 * a.n_params


# ---- 2/3. reproducibility -----------------------------------------------------------
def test_same_master_seed_reproduces_replicates(system):
    a, theta, _ = system
    r1 = FiniteShotParameterShift(a, "global", 32, master_seed=123).replicate_component(theta, 0, 20)
    r2 = FiniteShotParameterShift(a, "global", 32, master_seed=123).replicate_component(theta, 0, 20)
    assert np.array_equal(r1, r2)
    G1 = FiniteShotParameterShift(a, "global", 32, master_seed=[5, 6]).replicate_gradient(theta, 3)
    G2 = FiniteShotParameterShift(a, "global", 32, master_seed=[5, 6]).replicate_gradient(theta, 3)
    assert np.array_equal(G1, G2)


def test_different_master_seed_changes_replicates(system):
    a, theta, _ = system
    r1 = FiniteShotParameterShift(a, "global", 32, master_seed=1).replicate_component(theta, 0, 20)
    r2 = FiniteShotParameterShift(a, "global", 32, master_seed=2).replicate_component(theta, 0, 20)
    assert not np.array_equal(r1, r2)


def test_global_numpy_state_does_not_affect_estimator(system):
    a, theta, _ = system
    np.random.seed(0)
    r1 = FiniteShotParameterShift(a, "global", 32, master_seed=9).replicate_component(theta, 1, 10)
    np.random.seed(12345)
    r2 = FiniteShotParameterShift(a, "global", 32, master_seed=9).replicate_component(theta, 1, 10)
    assert np.array_equal(r1, r2)


# ---- 4. estimator independent of the analytic reference ------------------------------
def test_estimator_module_does_not_import_exact_reference():
    tree = ast.parse(inspect.getsource(fs_mod))
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom):
            assert "exact" not in (node.module or ""), node.module
            assert not any("autograd" in n.name or "parameter_shift" in n.name for n in node.names)
        if isinstance(node, ast.Import):
            assert not any("exact" in n.name for n in node.names)


def test_estimator_runs_when_exact_reference_is_disabled(system, monkeypatch):
    a, theta, _ = system

    def boom(*args, **kwargs):
        raise AssertionError("estimator touched the analytic reference")

    for name in ("gradient_autograd", "gradient_parameter_shift", "gradient_finite_difference", "cost_value", "make_cost"):
        monkeypatch.setattr(exact_mod, name, boom)
    est = FiniteShotParameterShift(a, "global", 8, master_seed=0)
    assert np.isfinite(est.gradient_component(theta, 3))
    assert est.gradient(theta).shape == theta.shape


# ---- 5/6. theoretical variance --------------------------------------------------------
def test_theoretical_variance_nonnegative_and_bounded(system):
    a, theta, _ = system
    for M in (1, 7, 64, 4096):
        v = theoretical_parameter_shift_variance_vector(a, theta, M)
        assert v.shape == (a.n_params,)
        assert np.all(v >= 0.0)
        assert np.all(v <= 0.5 / M + 1e-15)  # max when C_+ = C_- = 0


def test_theoretical_variance_scales_exactly_as_one_over_M(system):
    a, theta, _ = system
    base = theoretical_parameter_shift_variance(a, theta, 4, 1)
    for M in (2, 10, 64, 1000):
        assert theoretical_parameter_shift_variance(a, theta, 4, M) * M == pytest.approx(base, rel=1e-12)


def test_theoretical_variance_matches_explicit_formula(system):
    a, theta, _ = system
    k, M = 5, 100
    cp = cost_value(a, "global", shifted_params(theta, k, +SHIFT))
    cm = cost_value(a, "global", shifted_params(theta, k, -SHIFT))
    assert theoretical_parameter_shift_variance(a, theta, k, M) == pytest.approx((2 - cp**2 - cm**2) / (4 * M), rel=1e-12)


def test_theoretical_variance_rejects_non_global_cost(system):
    a, theta, _ = system
    with pytest.raises(ValueError):
        theoretical_parameter_shift_variance(a, theta, 0, 100, cost="local")


# ---- 7/8. statistical agreement (deterministic seed, wide tolerances) -------------------
def test_empirical_variance_close_to_theory(system):
    # 400 replicates: relative SD of a sample variance ~ sqrt(2/399) ~ 0.07, so
    # [0.65, 1.5] is far outside plausible sampling fluctuation.
    a, theta, g = system
    k, M = 6, 64
    est = FiniteShotParameterShift(a, "global", M, master_seed=[777, k, M])
    g_hat = est.replicate_component(theta, k, 400)
    s = summarize_shot_noise(g_hat, g[k], M, theoretical_parameter_shift_variance(a, theta, k, M))
    assert 0.65 < s["var_ratio_emp_over_theory"] < 1.5, s


def test_mean_shot_noise_compatible_with_zero(system):
    a, theta, g = system
    k, M = 11, 64
    est = FiniteShotParameterShift(a, "global", M, master_seed=[888, k, M])
    g_hat = est.replicate_component(theta, k, 400)
    s = summarize_shot_noise(g_hat, g[k], M, theoretical_parameter_shift_variance(a, theta, k, M))
    assert abs(s["z_mean_xi"]) < 4.0, s  # |z| >= 4 has p ~ 6e-5 under H0
    assert s["ci95_low"] <= 0.0 + 3 * s["se_mean_xi"] and s["ci95_high"] >= 0.0 - 3 * s["se_mean_xi"]


# ---- Task 10 sanity checks --------------------------------------------------------------
def test_A_no_reuse_of_randomness_across_pm_replicates_and_indices(system):
    a, theta, _ = system
    est = FiniteShotParameterShift(a, "global", 16, master_seed=0)
    est.replicate_gradient(theta, 3)  # 3 reps x 12 k x (+,-) = 72 evaluations
    assert est.n_evaluations == 72
    assert len(est.spawn_keys) == 72 and len(set(est.spawn_keys)) == 72
    # and two estimates of the same circuit at the same params differ (fresh streams)
    est2 = FiniteShotParameterShift(a, "global", 64, master_seed=0)
    vals = [est2.cost_estimate(theta) for _ in range(6)]
    assert len(set(vals)) > 1


def test_B_D_finite_shot_actually_uses_M_shots(system):
    # An M-shot mean of ±1 outcomes lies on the lattice {-1, -1+2/M, ..., 1}.
    a, theta, _ = system
    for M in (1, 3, 5, 8):
        est = FiniteShotParameterShift(a, "global", M, master_seed=M)
        for _ in range(5):
            c = est.cost_estimate(theta)
            counts = (c + 1.0) * M / 2.0  # number of +1 outcomes
            assert counts == pytest.approx(round(counts), abs=1e-12) and 0 <= round(counts) <= M
    # the exact value is off-lattice for M=3 (so a shots=None fallback would be caught)
    c_exact = cost_value(a, "global", theta)
    assert (c_exact + 1) * 3 / 2 != pytest.approx(round((c_exact + 1) * 3 / 2), abs=1e-9)


def test_C_analytic_reference_has_no_shot_randomness(system):
    a, theta, _ = system
    before = gradient_autograd(a, "global", theta)
    FiniteShotParameterShift(a, "global", 16, master_seed=3).replicate_gradient(theta, 2)
    np.random.seed(999)
    after = gradient_autograd(a, "global", theta)
    assert np.array_equal(before, after)
    assert "shots" not in inspect.getsource(exact_mod.make_cost).replace("shots=None", "")


def test_E_flatten_indexing_consistent():
    shape = (2, 3, 2)
    for k in range(12):
        m = flat_to_multi(k, shape)
        assert m == tuple(np.unravel_index(k, shape))
        assert multi_to_flat(m, shape) == k
        base = np.arange(12, dtype=float).reshape(shape)
        sh = shifted_params(base, k, 0.5)
        assert sh.ravel()[k] == k + 0.5 and np.sum(sh != base) == 1
    with pytest.raises(IndexError):
        flat_to_multi(12, shape)


def test_E_component_estimator_targets_correct_parameter(system):
    # With many shots, g_hat_k must land near g_exact[k], not some other component.
    a, theta, g = system
    M = 4096
    for k in (0, 5, 11):
        est = FiniteShotParameterShift(a, "global", M, master_seed=[4, k])
        gh = np.mean(est.replicate_component(theta, k, 8))
        sd = np.sqrt(theoretical_parameter_shift_variance(a, theta, k, M) / 8)
        assert abs(gh - g[k]) < 5 * sd + 1e-12


def test_F_global_observable_single_shot_outcomes_are_pm1(system):
    a, theta, _ = system
    dev = qml.device("default.qubit", wires=3, seed=np.random.default_rng(0))

    @qml.set_shots(shots=500)
    @qml.qnode(dev)
    def samples(p):
        a.apply(p)
        return qml.sample(global_z_observable(3))

    v = np.asarray(samples(theta))
    assert v.shape == (500,)
    assert set(np.unique(v)).issubset({-1.0, 1.0})
