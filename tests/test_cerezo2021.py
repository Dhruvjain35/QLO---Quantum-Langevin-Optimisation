"""Stage 3 tests: Cerezo-2021 RX-product benchmark. Deterministic identities wherever possible."""

import numpy as np
import pytest
from scipy.integrate import dblquad, quad

from qlo.analysis import bp_variance_per_parameter
from qlo.benchmarks.cerezo2021 import (
    COST_NAMES,
    THEORETICAL_GLOBAL_SEMILOG_SLOPE,
    Cerezo2021Circuit,
    closed_form_global_cost,
    closed_form_global_gradient,
    closed_form_global_gradient_component_logdomain,
    closed_form_gradient_component,
    closed_form_local_cost,
    closed_form_local_gradient,
    factorized_global_second_moment_estimate,
    sample_params,
    theoretical_fourth_moment,
    theoretical_global_variance,
    theoretical_local_variance,
    theoretical_variance,
    theoretical_variance_estimator_rel_se,
)
from qlo.experiments.controlled_bp import ControlledBPConfig, fit_models, gradient_component_stats, run_stage3
from qlo.utils.seeding import make_rng

TIGHT = 1e-12


def _cf_cost(cost, th):
    return closed_form_global_cost(th) if cost == "global" else closed_form_local_cost(th)


def _cf_grad(cost, th):
    return closed_form_global_gradient(th) if cost == "global" else closed_form_local_gradient(th)


# ---- 1/2 sampling -----------------------------------------------------------------------
def test_sampling_reproducible_and_in_range():
    a = sample_params(6, make_rng(5), 100)
    b = sample_params(6, make_rng(5), 100)
    assert np.array_equal(a, b) and a.shape == (100, 6)
    assert np.all(a >= -np.pi) and np.all(a <= np.pi)
    assert not np.array_equal(a, sample_params(6, make_rng(6), 100))
    assert sample_params(3, make_rng(0)).shape == (3,)


# ---- 3-7 PennyLane vs closed form -------------------------------------------------------
@pytest.mark.parametrize("n", [2, 3, 4, 6])
@pytest.mark.parametrize("seed", [0, 1, 2])
@pytest.mark.parametrize("cost", COST_NAMES)
def test_pennylane_cost_and_gradients_match_closed_form(n, seed, cost):
    circ = Cerezo2021Circuit(n)
    th = sample_params(n, make_rng([n, seed]))
    assert abs(circ.cost(cost, th) - _cf_cost(cost, th)) < TIGHT
    g_cf = _cf_grad(cost, th)
    assert np.max(np.abs(circ.gradient(cost, th, "backprop") - g_cf)) < TIGHT
    assert np.max(np.abs(circ.gradient(cost, th, "parameter-shift") - g_cf)) < TIGHT


@pytest.mark.parametrize("cost", COST_NAMES)
def test_batched_parameter_shift_matches_closed_form(cost):
    circ = Cerezo2021Circuit(4)
    TH = sample_params(4, make_rng(9), 64)
    for k in range(4):
        g = circ.gradient_component_parameter_shift_batched(cost, TH, k)
        assert g.shape == (64,)
        assert np.max(np.abs(g - closed_form_gradient_component(cost, TH, k))) < TIGHT


def test_closed_form_component_matches_full_gradient():
    TH = sample_params(5, make_rng(2), 30)
    for cost in COST_NAMES:
        full = _cf_grad(cost, TH)
        for k in range(5):
            assert np.array_equal(full[:, k], closed_form_gradient_component(cost, TH, k))


def test_closed_form_gradient_matches_finite_difference():
    th = sample_params(4, make_rng(3))
    h = 1e-6
    for cost in COST_NAMES:
        fd = np.array([(_cf_cost(cost, th + h * e) - _cf_cost(cost, th - h * e)) / (2 * h) for e in np.eye(4)])
        assert np.max(np.abs(fd - _cf_grad(cost, th))) < 1e-8


# ---- 8/9 theoretical variance by independent quadrature ---------------------------------
def test_theoretical_global_variance_by_quadrature():
    # n=1: E[(sin t / 2)^2] = 1/8 ; n=2: E[(1/2 sin t0 cos^2(t1/2))^2] = 3/64
    v1 = quad(lambda t: (0.5 * np.sin(t)) ** 2, -np.pi, np.pi)[0] / (2 * np.pi)
    v2 = dblquad(lambda t1, t0: (0.5 * np.sin(t0) * np.cos(t1 / 2) ** 2) ** 2, -np.pi, np.pi, -np.pi, np.pi)[0] / (2 * np.pi) ** 2
    assert theoretical_global_variance(1) == pytest.approx(v1, rel=1e-9) == pytest.approx(1 / 8)
    assert theoretical_global_variance(2) == pytest.approx(v2, rel=1e-9) == pytest.approx(3 / 64)
    assert theoretical_global_variance(3) == pytest.approx(9 / 512)
    # mean is exactly zero (odd in theta_0)
    m2 = dblquad(lambda t1, t0: 0.5 * np.sin(t0) * np.cos(t1 / 2) ** 2, -np.pi, np.pi, -np.pi, np.pi)[0]
    assert abs(m2) < 1e-10


def test_theoretical_local_variance_by_quadrature():
    for n in (1, 3, 5):
        v = quad(lambda t: (np.sin(t) / (2 * n)) ** 2, -np.pi, np.pi)[0] / (2 * np.pi)
        assert theoretical_local_variance(n) == pytest.approx(v, rel=1e-9) == pytest.approx(1 / (8 * n**2))
        assert theoretical_variance("local", n) == theoretical_local_variance(n)


# ---- 10/11 exact scaling identities ------------------------------------------------------
def test_global_variance_ratio_is_three_eighths():
    for n in range(1, 30):
        assert theoretical_global_variance(n + 1) / theoretical_global_variance(n) == pytest.approx(3 / 8, rel=1e-13)
    assert THEORETICAL_GLOBAL_SEMILOG_SLOPE == pytest.approx(np.log(3 / 8))


def test_local_variance_scales_as_n_minus_two():
    for n in range(1, 30):
        assert theoretical_local_variance(n) * n**2 == pytest.approx(1 / 8, rel=1e-13)
        assert theoretical_local_variance(2 * n) / theoretical_local_variance(n) == pytest.approx(0.25, rel=1e-13)


def test_fourth_moment_by_quadrature_and_rel_se():
    m4_1 = quad(lambda t: (0.5 * np.sin(t)) ** 4, -np.pi, np.pi)[0] / (2 * np.pi)
    m4_2 = dblquad(lambda t1, t0: (0.5 * np.sin(t0) * np.cos(t1 / 2) ** 2) ** 4, -np.pi, np.pi, -np.pi, np.pi)[0] / (2 * np.pi) ** 2
    assert theoretical_fourth_moment("global", 1) == pytest.approx(m4_1, rel=1e-9)
    assert theoretical_fourth_moment("global", 2) == pytest.approx(m4_2, rel=1e-9)
    assert theoretical_fourth_moment("local", 3) == pytest.approx(quad(lambda t: (np.sin(t) / 6) ** 4, -np.pi, np.pi)[0] / (2 * np.pi), rel=1e-9)
    # local kurtosis ratio is constant 3/2; global grows by 35/18 per qubit
    for n in (1, 5, 20):
        assert theoretical_fourth_moment("local", n) / theoretical_local_variance(n) ** 2 == pytest.approx(1.5)
    r = [theoretical_fourth_moment("global", n) / theoretical_global_variance(n) ** 2 for n in (3, 4)]
    assert r[1] / r[0] == pytest.approx(35 / 18)
    assert theoretical_variance_estimator_rel_se("local", 10, 10_000) == pytest.approx(np.sqrt((1.5 - 9997 / 9999) / 10_000))


# ---- 12 BP statistic is across theta samples ------------------------------------------------
def test_bp_statistic_is_variance_across_theta_samples():
    TH = sample_params(4, make_rng(11), 500)
    G = closed_form_global_gradient(TH)  # (500, 4)
    per_k = bp_variance_per_parameter(G)
    for k in range(4):
        s = gradient_component_stats(G[:, k], theoretical_global_variance(4))
        assert s["empirical_variance"] == pytest.approx(np.var(G[:, k], ddof=1)) == pytest.approx(per_k[k])
        assert s["n_initializations"] == 500
    # the pooled variance over all entries of G is a different statistic and must not be what we report
    assert gradient_component_stats(G[:, 0], 1.0)["empirical_variance"] != np.var(G.ravel(), ddof=1)


# ---- 13 structural-zero audit ------------------------------------------------------------------
@pytest.mark.parametrize("n", [2, 3, 4, 6, 8])
def test_no_parameter_is_structurally_zero(n):
    # Analytic witness: at theta = (pi/2) e_k, dC_G/dtheta_k = 1/2 and dC_L/dtheta_k = 1/(2n) exactly.
    circ = Cerezo2021Circuit(n)
    for k in range(n):
        th = np.zeros(n); th[k] = np.pi / 2
        assert closed_form_global_gradient(th)[k] == pytest.approx(0.5)
        assert closed_form_local_gradient(th)[k] == pytest.approx(1 / (2 * n))
        assert abs(circ.gradient("global", th)[k] - 0.5) < TIGHT
        assert abs(circ.gradient("local", th)[k] - 1 / (2 * n)) < TIGHT
    # Numerical: over random theta every component is nonzero somewhere, for both costs
    TH = sample_params(n, make_rng([13, n]), 200)
    for cost in COST_NAMES:
        assert np.all(np.max(np.abs(_cf_grad(cost, TH)), axis=0) > 1e-3)


# ---- numerical stability + diagnostics ---------------------------------------------------------
def test_logdomain_matches_direct_and_no_underflow_at_n24():
    TH = sample_params(24, make_rng(24), 2000)
    g = closed_form_gradient_component("global", TH, 0)
    assert np.all(np.isfinite(g)) and np.sum(g == 0.0) == 0
    s, la = closed_form_global_gradient_component_logdomain(TH, 0)
    assert np.max(np.abs(s * np.exp(la) - g) / np.abs(g)) < 1e-12


def test_factorized_estimator_is_consistent():
    TH = sample_params(12, make_rng(1), 20_000)
    f = factorized_global_second_moment_estimate(TH, 0)
    # relative error ~ sqrt(0.94*12/20000) ~ 0.024 ; allow 6x
    assert abs(f / theoretical_global_variance(12) - 1) < 0.15
    with pytest.raises(ValueError):
        factorized_global_second_moment_estimate(TH[0], 0)


def test_fit_models_recover_exact_curves():
    n = np.arange(2, 26, 2, dtype=float)
    ge = fit_models(n, (1 / 8) * (3 / 8) ** (n - 1))["exponential"]
    lp = fit_models(n, 1 / (8 * n**2))["power_law"]
    assert ge["slope"] == pytest.approx(np.log(3 / 8), abs=1e-12) and ge["r2"] == pytest.approx(1.0)
    assert lp["slope"] == pytest.approx(-2.0, abs=1e-12) and lp["r2"] == pytest.approx(1.0)


def test_stage3_driver_runs_tiny(tmp_path):
    cfg = ControlledBPConfig(n_values=(2, 3, 4), n_init=300, symmetry_n=(3,), symmetry_n_init=300,
                             pennylane_n=(2, 3), pennylane_n_init_grad=20, pennylane_n_init_batched=100)
    res = run_stage3(cfg, tmp_path, figures=True)
    assert len(res["summary"]) == 6 and set(res["summary"].cost_type) == set(COST_NAMES)
    assert (tmp_path / "controlled_bp_summary.csv").exists() and (tmp_path / "global_fit.json").exists()
    assert (tmp_path / "model_comparison.json").exists() and (tmp_path / "config.json").exists()
    assert len(res["figures"]) == 5 and all(f.stat().st_size > 1000 for f in res["figures"])
    num = res["summary"].select_dtypes(include=[np.number]).drop(columns=["factorized_second_moment", "factorized_ratio"])
    assert np.isfinite(num.to_numpy()).all()
    assert res["crosscheck"].pl_grad_max_abs_diff_vs_closed_form.max() < TIGHT
