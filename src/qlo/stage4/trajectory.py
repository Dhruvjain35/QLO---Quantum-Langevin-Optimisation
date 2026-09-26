"""Trajectory harness with a strictly OFFLINE exact evaluator.

The method object decides every update using only its own oracle (finite-shot estimates
for ``shot_sgd``). After each update the harness computes the exact ``F``, ``C`` and
``||g||`` from a *copy* of theta for scientific bookkeeping. These values are never passed
to the method, do not change the trajectory (tested: ``offline=False`` gives a
bit-identical final theta), and are not counted as measurements.

Metrics per trajectory: initial / final / best fidelity, cost, gradient norms; log10
fidelity gains; first-hitting iteration for fixed fidelity targets (with the method's
cumulative resources at that moment) and for relative targets 10 F0 and 100 F0. A
trajectory that starts at or above a fixed target is marked ineligible for it.
"""

from __future__ import annotations

import numpy as np

from qlo.stage4.landscape import exact_gradient, fidelity

FIXED_TARGETS = (0.01, 0.1, 0.5)
RELATIVE_TARGETS = (10.0, 100.0)
EPS = 1e-300


def _offline(theta: np.ndarray) -> tuple[float, float]:
    """Exact (F, ||g||) from a copy of theta. Pure function; nothing flows back to the optimizer."""
    t = np.array(theta, dtype=np.float64, copy=True)
    return float(fidelity(t)), float(np.linalg.norm(exact_gradient(t)))


def run_trajectory(method, theta0: np.ndarray, n_iters: int, offline: bool = True, record_full: bool = False) -> dict:
    theta = np.array(theta0, dtype=np.float64, copy=True)
    n = theta.size
    out: dict = {"n_qubits": n, "method": method.name, "n_iters": int(n_iters)}
    if not offline:
        for _ in range(int(n_iters)):
            theta = method.step(theta)
        out["theta_final"] = theta
        out["resources_final"] = method.resources
        return out

    f0, g0 = _offline(theta)
    out.update({"F0": f0, "C0": 1.0 - f0, "gnorm0": g0})
    f_best, it_best, g_max = f0, 0, g0
    hit = {f"hit_F{t}": None for t in FIXED_TARGETS}
    hit_res = {}
    eligible = {t: f0 < t for t in FIXED_TARGETS}
    rel_hit = {f"hit_{int(r)}xF0": None for r in RELATIVE_TARGETS}
    rel_possible = {r: r * f0 <= 1.0 for r in RELATIVE_TARGETS}
    full_F = [f0] if record_full else None
    full_g = [g0] if record_full else None

    for it in range(1, int(n_iters) + 1):
        theta = method.step(theta)
        f, g = _offline(theta)
        if record_full:
            full_F.append(f)
            full_g.append(g)
        if f > f_best:
            f_best, it_best = f, it
        if g > g_max:
            g_max = g
        for t in FIXED_TARGETS:
            key = f"hit_F{t}"
            if hit[key] is None and eligible[t] and f >= t:
                hit[key] = it
                hit_res[key] = dict(method.resources)
        for r in RELATIVE_TARGETS:
            key = f"hit_{int(r)}xF0"
            if rel_hit[key] is None and rel_possible[r] and f >= r * f0:
                rel_hit[key] = it

    out.update({"F_final": f, "C_final": 1.0 - f, "gnorm_final": g,
                "F_best": f_best, "C_best": 1.0 - f_best, "iter_best": it_best, "gnorm_max": g_max,
                "log10_gain_best": float(np.log10((f_best + EPS) / (f0 + EPS))),
                "log10_gain_final": float(np.log10((f + EPS) / (f0 + EPS)))})
    for t in FIXED_TARGETS:
        key = f"hit_F{t}"
        out[f"eligible_F{t}"] = bool(eligible[t])
        out[key] = hit[key]
        res = hit_res.get(key, {})
        out[f"circuit_evals_to_F{t}"] = res.get("circuit_evaluations")
        out[f"shots_to_F{t}"] = res.get("shots")
        out[f"oracle_calls_to_F{t}"] = res.get("oracle_gradient_calls")
    for r in RELATIVE_TARGETS:
        key = f"hit_{int(r)}xF0"
        out[f"possible_{int(r)}xF0"] = bool(rel_possible[r])
        out[key] = rel_hit[key]
    out["theta_final"] = theta
    out["resources_final"] = dict(method.resources)
    if record_full:
        out["trajectory_F"] = np.array(full_F)
        out["trajectory_gnorm"] = np.array(full_g)
    return out
