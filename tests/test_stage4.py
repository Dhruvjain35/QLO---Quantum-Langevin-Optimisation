"""Stage 4 tests. Statistical tests use fixed seeds and wide, justified tolerances."""

import inspect

import numpy as np
import pytest

from qlo.benchmarks.cerezo2021 import closed_form_global_cost, closed_form_global_gradient
from qlo.stage4 import methods as methods_mod
from qlo.stage4 import shots as shots_mod
from qlo.stage4.landscape import cost, exact_gradient, fidelity, shifted_fidelities, wrap_angles
from qlo.stage4.methods import ConstantLangevin, ExactGD, MatchedGaussian, NoiseOnly, ShotSGD, make_method
from qlo.stage4.seeds import EVALUATION_START_SEEDS, TUNING_START_SEEDS, initial_theta, noise_seed
from qlo.stage4.shots import BinomialProjectorShots, conditional_shot_variance
from qlo.stage4.trajectory import run_trajectory


@pytest.fixture(scope="module")
def theta6():
    return 0.6 * initial_theta(3, 6)  # moderate fidelity so binomial counts are non-degenerate


# ---- landscape ---------------------------------------------------------------------------
def test_landscape_matches_stage3_closed_forms():
    for s in range(5):
        th = initial_theta(s, 7)
        assert abs(cost(th) - closed_form_global_cost(th)) < 1e-15
        assert np.max(np.abs(exact_gradient(th) - closed_form_global_gradient(th))) < 1e-15
        fp, fm = shifted_fidelities(th)
        for k in range(7):
            e = np.zeros(7); e[k] = np.pi / 2
            assert abs(fp[k] - fidelity(th + e)) < 1e-15 and abs(fm[k] - fidelity(th - e)) < 1e-15


def test_wrapping_preserves_cost_and_fidelity():
    th = initial_theta(1, 5) + np.array([0, 2 * np.pi, -4 * np.pi, 6 * np.pi, 10 * np.pi])
    w = wrap_angles(th)
    assert np.all(w >= -np.pi) and np.all(w < np.pi)
    assert abs(fidelity(w) - fidelity(th)) < 1e-12 and abs(cost(w) - cost(th)) < 1e-12
    assert np.max(np.abs(exact_gradient(w) - exact_gradient(th))) < 1e-12


# ---- 1/2 binomial projector estimator mean and variance -------------------------------------
def test_projector_binomial_cost_mean_and_variance(theta6):
    M = 64
    est = BinomialProjectorShots(6, M, np.random.SeedSequence(1))
    F = fidelity(theta6)
    c = np.array([est.cost_estimate(theta6) for _ in range(20_000)])
    assert abs(c.mean() - (1 - F)) < 5 * np.sqrt(F * (1 - F) / M / 20_000)
    assert abs(c.var(ddof=1) / (F * (1 - F) / M) - 1) < 0.1  # rel SD of a variance ratio at 2e4 samples ~ 0.01
    assert np.allclose(c * M, np.round(c * M))  # lattice of K/M


# ---- 3/4 finite-shot gradient unbiased; analytic variance -----------------------------------
def test_shot_gradient_unbiased_and_variance_matches_formula(theta6):
    M = 32
    est = BinomialProjectorShots(6, M, np.random.SeedSequence(2))
    G = np.stack([est.gradient(theta6) for _ in range(20_000)])
    g, v = exact_gradient(theta6), conditional_shot_variance(theta6, M)
    z = (G.mean(0) - g) / np.sqrt(v / 20_000)
    assert np.all(np.abs(z) < 4.5)
    assert np.all(np.abs(G.var(0, ddof=1) / v - 1) < 0.12)
    # formula check against explicit shifted fidelities
    fp, fm = shifted_fidelities(theta6)
    assert np.allclose(v, (fp * (1 - fp) + fm * (1 - fm)) / (4 * M))
    assert np.all(v >= 0) and np.allclose(conditional_shot_variance(theta6, 4 * M) * 4, v)


# ---- 5/6 matched Gaussian ----------------------------------------------------------------------
def test_matched_gaussian_uses_state_dependent_variance_and_matches_covariance(theta6):
    M, eta = 32, 0.3
    m = MatchedGaussian(eta, M, np.random.SeedSequence(3))
    steps = np.stack([(theta6 - m.step(theta6)) / eta for _ in range(20_000)])  # = g + eps
    g, v = exact_gradient(theta6), conditional_shot_variance(theta6, M)
    assert np.all(np.abs((steps.mean(0) - g) / np.sqrt(v / 20_000)) < 4.5)
    C = np.cov(steps, rowvar=False)
    assert np.all(np.abs(np.diag(C) / v - 1) < 0.12)
    off = C[~np.eye(6, dtype=bool)]
    assert np.max(np.abs(off)) < 5 * np.sqrt(np.outer(v, v))[~np.eye(6, dtype=bool)].max() / np.sqrt(20_000) * 3
    # a different theta gives a different variance vector (state dependence)
    assert not np.allclose(conditional_shot_variance(theta6, M), conditional_shot_variance(0.5 * theta6, M))


# ---- 7 independence of shot batches --------------------------------------------------------------
def test_shot_batches_independent_across_shift_parameters_and_iterations():
    th = np.zeros(4)  # at theta=0, F_+ = F_- = 1/2 for every k: maximal binomial variance, identical p everywhere
    M = 16
    est = BinomialProjectorShots(4, M, np.random.SeedSequence(4))
    fp, fm = shifted_fidelities(th)
    assert np.allclose(fp, 0.5) and np.allclose(fm, 0.5)
    G = np.stack([est.gradient(th) for _ in range(4000)])
    # (+,-) independence: Var(g_k) = 2 * (1/4)/(4M) = 1/(8M); a reused batch would give 0
    assert np.all(np.abs(G.var(0, ddof=1) / (1 / (8 * M)) - 1) < 0.1)
    # across parameters: correlation ~ 0 (SE ~ 0.016 at 4000 samples)
    C = np.corrcoef(G, rowvar=False)
    assert np.max(np.abs(C[~np.eye(4, dtype=bool)])) < 0.08
    # across iterations: lag-1 autocorrelation ~ 0
    for k in range(4):
        assert abs(np.corrcoef(G[:-1, k], G[1:, k])[0, 1]) < 0.08
    a, b = est.gradient(th), est.gradient(th)
    assert not np.array_equal(a, b)


def test_estimator_never_imports_exact_gradient():
    src = inspect.getsource(shots_mod)
    assert "exact_gradient" not in src and "closed_form_global_gradient" not in src


# ---- 8/9/10 reproducibility ----------------------------------------------------------------------------
def test_master_seed_reproduces_trajectory_and_noise_seed_changes_it():
    th0 = initial_theta(0, 6)
    for name, hp in (("shot_sgd", {"eta": 0.3, "shots": 16}), ("matched_gaussian", {"eta": 0.3, "shots": 16}),
                     ("constant_langevin", {"eta": 0.3, "sigma_step": 0.1}), ("noise_only", {"eta": 0.3, "shots": 16})):
        r1 = run_trajectory(make_method(name, 6, hp, noise_seed(name, 6, 0, 0)), th0, 50)
        r2 = run_trajectory(make_method(name, 6, hp, noise_seed(name, 6, 0, 0)), th0, 50)
        r3 = run_trajectory(make_method(name, 6, hp, noise_seed(name, 6, 0, 1)), th0, 50)
        assert np.array_equal(r1["theta_final"], r2["theta_final"]) and r1["F_best"] == r2["F_best"]
        assert not np.array_equal(r1["theta_final"], r3["theta_final"])
    # distinct methods with the same (n, seed, replicate) get distinct noise streams
    assert noise_seed("shot_sgd", 6, 0, 0).entropy != noise_seed("matched_gaussian", 6, 0, 0).entropy


def test_exact_gd_is_deterministic_and_monotone_for_small_eta():
    th0 = initial_theta(5, 6)
    r1 = run_trajectory(ExactGD(0.3), th0, 200, record_full=True)
    r2 = run_trajectory(ExactGD(0.3), th0, 200, record_full=True)
    assert np.array_equal(r1["theta_final"], r2["theta_final"])
    assert np.all(np.diff(r1["trajectory_F"]) >= -1e-15)  # gradient descent on C increases F


def test_initial_theta_reproducible_and_in_range():
    a, b = initial_theta(1000, 8), initial_theta(1000, 8)
    assert np.array_equal(a, b) and np.all(a >= -np.pi) and np.all(a <= np.pi)
    assert not np.array_equal(a, initial_theta(1001, 8))


# ---- 12/13 resource accounting --------------------------------------------------------------------------
@pytest.mark.parametrize("n,M", [(3, 7), (6, 64), (10, 1024)])
def test_resource_counter_2n_evaluations_and_2nM_shots(n, M):
    m = ShotSGD(n, 0.1, M, np.random.SeedSequence(0))
    th = initial_theta(0, n)
    for it in range(1, 4):
        th = m.step(th)
        assert m.resources["circuit_evaluations"] == 2 * n * it
        assert m.resources["shots"] == 2 * n * M * it
    r = run_trajectory(ShotSGD(n, 0.1, M, np.random.SeedSequence(0)), initial_theta(0, n), 25)
    assert r["resources_final"] == {"circuit_evaluations": 50 * n, "shots": 50 * n * M, "oracle_gradient_calls": 0}
    for name in ("exact_gd", "matched_gaussian", "constant_langevin", "noise_only"):
        hp = {"eta": 0.1, "shots": M, "sigma_step": 0.1}
        rr = run_trajectory(make_method(name, n, hp, np.random.SeedSequence(1)), initial_theta(0, n), 25)
        assert rr["resources_final"]["shots"] == 0 and rr["resources_final"]["circuit_evaluations"] == 0


# ---- 14 offline evaluator cannot alter optimizer state --------------------------------------------------
def test_offline_evaluator_does_not_alter_trajectory():
    th0 = initial_theta(2, 6)
    for name, hp in (("shot_sgd", {"eta": 0.3, "shots": 16}), ("constant_langevin", {"eta": 0.3, "sigma_step": 0.1}), ("exact_gd", {"eta": 0.3})):
        seed = noise_seed(name, 6, 2, 0) if name != "exact_gd" else None
        with_off = run_trajectory(make_method(name, 6, hp, seed), th0, 300, offline=True)
        without = run_trajectory(make_method(name, 6, hp, seed), th0, 300, offline=False)
        assert np.array_equal(with_off["theta_final"], without["theta_final"])
        assert with_off["resources_final"] == without["resources_final"]
    # the method objects have no attribute that could receive exact fidelity/gradient info from the harness
    for cls in (ShotSGD,):
        assert not any(a in inspect.getsource(cls) for a in ("fidelity(", "exact_gradient("))
    src = inspect.getsource(methods_mod)
    assert "ShotSGD" in src


# ---- 15 seed sets disjoint ---------------------------------------------------------------------------------
def test_tuning_and_evaluation_seed_sets_disjoint():
    assert set(TUNING_START_SEEDS).isdisjoint(EVALUATION_START_SEEDS)
    assert TUNING_START_SEEDS == tuple(range(40)) and EVALUATION_START_SEEDS == tuple(range(1000, 1100))


# ---- trajectory metrics -------------------------------------------------------------------------------------
def test_trajectory_metrics_and_target_eligibility():
    th0 = 0.2 * initial_theta(0, 4)  # F0 above 0.5 -> ineligible for all fixed targets
    r = run_trajectory(ExactGD(0.3), th0, 20)
    assert r["F0"] > 0.5 and not r["eligible_F0.5"] and r["hit_F0.5"] is None
    assert r["F_best"] >= r["F0"] and r["log10_gain_best"] >= 0 and r["C0"] == pytest.approx(1 - r["F0"])
    th1 = initial_theta(0, 4)
    r = run_trajectory(ExactGD(1.0), th1, 2000, record_full=True)
    assert len(r["trajectory_F"]) == 2001 and r["F_best"] == r["trajectory_F"].max()
    if r["hit_F0.1"] is not None:
        assert r["trajectory_F"][r["hit_F0.1"]] >= 0.1 and r["trajectory_F"][r["hit_F0.1"] - 1] < 0.1


def test_langevin_diffusion_coefficient():
    assert ConstantLangevin(0.3, 0.1, 0).diffusion_coefficient == pytest.approx(0.01 / 0.6)
    assert np.isnan(ConstantLangevin(0.0, 0.1, 0).diffusion_coefficient)


def test_noise_only_has_no_drift_on_average(theta6):
    m = NoiseOnly(0.3, 16, np.random.SeedSequence(9))
    d = np.stack([m.step(theta6) - theta6 for _ in range(5000)])
    v = conditional_shot_variance(theta6, 16)
    assert np.all(np.abs(d.mean(0) / (0.3 * np.sqrt(v / 5000))) < 4.5)
