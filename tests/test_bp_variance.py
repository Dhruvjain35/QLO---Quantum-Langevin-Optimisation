"""Task 1: Var_theta[dC/dtheta_k] must be a variance ACROSS initializations, per k."""

import numpy as np
import pytest

from qlo import HardwareEfficientAnsatz, gradient_autograd, make_rng
from qlo.analysis import (
    barren_plateau_variance,
    bp_variance_per_parameter,
    bp_variance_summary,
    gradient_samples_across_inits,
)
from qlo.experiments.bp_smoke import SmokeConfig, run_bp_smoke, summarize_gradient


def test_per_parameter_variance_is_across_rows_not_pooled():
    # column 0 constant -> 0; column 1 = [0, 2, 4] -> var 4 (ddof=1); column 2 = [1,1,7] -> 12
    G = np.array([[5.0, 0.0, 1.0], [5.0, 2.0, 1.0], [5.0, 4.0, 7.0]])
    v = bp_variance_per_parameter(G)
    assert v.shape == (3,)
    assert v == pytest.approx([0.0, 4.0, 12.0])
    pooled = np.var(G.ravel(), ddof=1)
    assert not np.isclose(v.mean(), pooled)  # summary across k != pooled variance
    s = bp_variance_summary(v)
    assert s["n_params"] == 3 and s["min_var_k"] == 0.0 and s["max_var_k"] == 12.0


def test_samples_matrix_matches_manual_loop():
    a = HardwareEfficientAnsatz(3, 2)
    G = gradient_samples_across_inits(a, "global", n_init=4, seed=11)
    rng = make_rng(11)
    for i in range(4):
        theta = a.init_params(rng)
        assert np.array_equal(G[i], gradient_autograd(a, "global", theta).ravel())
    assert G.shape == (4, a.n_params)


def test_barren_plateau_variance_single_k_equals_manual():
    a = HardwareEfficientAnsatz(3, 2)
    out = barren_plateau_variance(a, "global", n_init=6, seed=3, k=5)
    rng = make_rng(3)
    manual = [gradient_autograd(a, "global", a.init_params(rng)).ravel()[5] for _ in range(6)]
    assert out["var_k"] == pytest.approx(np.var(manual, ddof=1))
    assert np.array_equal(out["samples_k"], manual)


def test_barren_plateau_variance_all_k_and_determinism():
    a = HardwareEfficientAnsatz(2, 1)
    o1 = barren_plateau_variance(a, "local", n_init=5, seed=0)
    o2 = barren_plateau_variance(a, "local", n_init=5, seed=0)
    assert np.array_equal(o1["per_k"], o2["per_k"])
    assert o1["per_k"].shape == (a.n_params,)
    assert o1["summary"]["mean_var_k"] == pytest.approx(o1["per_k"].mean())


def test_requires_at_least_two_inits():
    with pytest.raises(ValueError):
        gradient_samples_across_inits(HardwareEfficientAnsatz(2, 1), "global", n_init=1, seed=0)


def test_smoke_within_gradient_stat_is_relabelled():
    stats = summarize_gradient(np.array([1.0, -1.0, 3.0]))
    assert "var_grad" not in stats
    assert stats["within_gradient_entry_var"] == pytest.approx(np.var([1.0, -1.0, 3.0]))
    df = run_bp_smoke(SmokeConfig(qubits=(2,), depth=1, n_init=2, seed=0))
    assert "within_gradient_entry_var" in df.columns and "var_grad" not in df.columns
