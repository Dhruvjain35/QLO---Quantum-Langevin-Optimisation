"""Structural-zero audit of every predeclared parameter position (and one positive control).

For each (family, regime, n) the audit draws its own theta sample (stream independent of the master seeds) and, per
position, compares FOUR gradient / fidelity computations:
    production   F_+- from the forward/backward brackets (the path used in the sweeps), g = (F_- - F_+)/2
    re-sim       F(theta +- pi/2 e_k) by full independent forward simulation of the shifted parameter vectors
    analytic     dC/dx from the generator-insertion derivative (no parameter shift)
    PennyLane    -dF/dtheta from autograd on the repo ``HardwareEfficientAnsatz`` (hea_ring) or on an independent
                 PennyLane template of the fallback family (rxry_czbrick)
Agreement is measured as |difference| / S per theta. Any disagreement above tolerance is an implementation error.

Classification (predeclared): STRUCTURAL if every audited theta has |r| <= tol (identically zero derivative);
DEGENERATE if the zero fraction exceeds the VALID limit without being 1, or Var(r) is below its floor; VALID otherwise.
"""

from __future__ import annotations

import numpy as np
import pennylane as qml
from pennylane import numpy as pnp

from qlo.b1 import circuits as C
from qlo.b1 import config as CF
from qlo.b1 import quantities as Q
from qlo.circuits.hardware_efficient import HardwareEfficientAnsatz


def pl_rxry_czbrick(params, n: int, depth: int) -> None:
    """Independent PennyLane template of the fallback family (queued inside a QNode)."""
    for layer in range(depth):
        for q in range(n):
            qml.RX(params[layer, q, 0], wires=q)
            qml.RY(params[layer, q, 1], wires=q)
        for a in range(layer % 2, n - 1, 2):
            qml.CZ(wires=[a, a + 1])


def pl_qnodes(family: str, n: int, depth: int):
    """(fidelity QNode, state QNode) built on PennyLane's default.qubit with backprop."""
    dev = qml.device("default.qubit", wires=n)
    if family == "hea_ring":
        ansatz = HardwareEfficientAnsatz(n, depth, "ring")
        body = ansatz.apply
    elif family == "rxry_czbrick":
        def body(p):
            pl_rxry_czbrick(p, n, depth)
    else:
        raise ValueError(family)

    @qml.qnode(dev, diff_method="backprop", interface="autograd")
    def fid(p):
        body(p)
        return qml.expval(qml.Projector(np.zeros(n, dtype=int), wires=range(n)))

    @qml.qnode(dev)
    def state(p):
        body(p)
        return qml.state()

    return fid, state


def pl_cost_gradients(theta: np.ndarray, family: str, n: int, depth: int) -> np.ndarray:
    """dC/dtheta = -dF/dtheta for every parameter, one autograd call per theta row; shape theta.shape."""
    fid, _ = pl_qnodes(family, n, depth)
    grad = qml.grad(fid)
    return np.stack([-np.asarray(grad(pnp.array(t, requires_grad=True))) for t in theta])


def pl_states(theta: np.ndarray, family: str, n: int, depth: int) -> np.ndarray:
    _, state = pl_qnodes(family, n, depth)
    return np.stack([np.asarray(state(t)) for t in theta])


def classify(zero_fraction: float, var_r: float, cfg: CF.B1Config) -> str:
    if zero_fraction >= 1.0:
        return "STRUCTURAL"
    if zero_fraction > cfg.audit_valid_max_zero_fraction or not var_r > cfg.audit_min_var_r:
        return "DEGENERATE"
    return "VALID"


def fallback_candidates(position, n: int):
    """Predeclared fallback order: the nearest EARLIER rotation of the same gate type (same rotation index j) in
    program order. Program order within a layer is qubit 0, 1, ..., n-1 (each qubit: rotation 0 then rotation 1)."""
    layer, qubit, j = position
    out = [(layer, q, j) for q in range(qubit - 1, -1, -1)]
    for l in range(layer - 1, -1, -1):
        out += [(l, q, j) for q in range(n - 1, -1, -1)]
    return out


def _shifted(theta: np.ndarray, position, delta: float) -> np.ndarray:
    t = theta.copy()
    t[:, position[0], position[1], position[2]] += delta
    return t


def position_audit(theta: np.ndarray, family: str, n: int, depth: int, position, F_plus, F_minus, F_bracket, F_forward,
                   g_pennylane, cfg: CF.B1Config) -> dict:
    """All checks for one position given the production arrays (and PennyLane gradients, or None)."""
    q = Q.per_theta(F_plus, F_minus)
    S = np.maximum(q["S"], 1e-300)
    fp_rs = C.fidelity(_shifted(theta, position, np.pi / 2), family, n, depth)
    fm_rs = C.fidelity(_shifted(theta, position, -np.pi / 2), family, n, depth)
    g_ps_rs = (fm_rs - fp_rs) / 2.0
    g_an = C.analytic_cost_gradient(theta, family, n, depth, position)
    abs_r = np.abs(q["r"])
    zero_fraction = float(np.mean(abs_r <= cfg.audit_zero_tol_r))
    var_r = float(np.var(q["r"]))
    row = {"err_bracket_vs_resim_F": float(np.max(np.maximum(np.abs(F_plus - fp_rs), np.abs(F_minus - fm_rs)) / S)),
           "err_bracket_F_vs_forward": float(np.max(np.abs(F_bracket - F_forward) / S)),
           "err_ps_resim_vs_analytic": float(np.max(np.abs(g_ps_rs - g_an) / S)),
           "err_production_vs_analytic": float(np.max(np.abs(q["g"] - g_an) / S)),
           "err_production_vs_pennylane": float(np.max(np.abs(q["g"] - g_pennylane) / S)) if g_pennylane is not None else float("nan"),
           "pennylane_checked": g_pennylane is not None,
           "zero_fraction": zero_fraction, "var_g": float(np.var(q["g"])), "var_r": var_r,
           "mean_abs_g": float(np.mean(np.abs(q["g"]))), "median_abs_r": float(np.median(abs_r)),
           "max_abs_r": float(np.max(abs_r)), "median_S": float(np.median(q["S"]))}
    errs = [row[k] for k in ("err_bracket_vs_resim_F", "err_bracket_F_vs_forward", "err_ps_resim_vs_analytic", "err_production_vs_analytic")]
    if g_pennylane is not None:
        errs.append(row["err_production_vs_pennylane"])
    row["max_agreement_err"] = float(max(errs))
    row["agreement_ok"] = bool(row["max_agreement_err"] <= cfg.audit_agreement_tol)
    row["classification"] = classify(zero_fraction, var_r, cfg)
    return row


def candidate_classification(theta: np.ndarray, family: str, n: int, depth: int, position, cfg: CF.B1Config) -> dict:
    """Zero fraction, Var(r) and class of an arbitrary rotation (any qubit) by full re-simulation; fallback search."""
    fp = C.fidelity(_shifted(theta, position, np.pi / 2), family, n, depth)
    fm = C.fidelity(_shifted(theta, position, -np.pi / 2), family, n, depth)
    q = Q.per_theta(fp, fm)
    zero_fraction = float(np.mean(np.abs(q["r"]) <= cfg.audit_zero_tol_r))
    var_r = float(np.var(q["r"]))
    return {"position": tuple(position), "zero_fraction": zero_fraction, "var_r": var_r, "classification": classify(zero_fraction, var_r, cfg)}


def fallback(theta: np.ndarray, family: str, n: int, depth: int, position, cfg: CF.B1Config) -> dict:
    """First candidate (nearest earlier rotation of the same gate type) whose audit sample is VALID."""
    for cand in fallback_candidates(position, n):
        c = candidate_classification(theta, family, n, depth, cand, cfg)
        if c["classification"] == "VALID":
            return c
    return {"position": None, "classification": "NO VALID FALLBACK"}


def audit_cell(family: str, regime: str, n: int, cfg: CF.B1Config) -> list[dict]:
    depth = CF.depth_of(regime, n)
    theta = CF.audit_theta(family, regime, n, cfg.audit_samples)
    positions = {**CF.position_map(depth), **CF.control_positions(family, depth)}
    distinct = set(CF.distinct_positions(depth).values()) | set(CF.control_positions(family, depth).values())
    F_forward, sh = C.shifted_fidelities(theta, family, n, depth, distinct)
    G = pl_cost_gradients(theta, family, n, depth) if cfg.audit_pennylane else None
    rows, cache = [], {}
    for label, pos in positions.items():
        if pos not in cache:
            g_pl = G[:, pos[0], pos[1], pos[2]] if G is not None else None
            cache[pos] = position_audit(theta, family, n, depth, pos, sh[pos]["F_plus"], sh[pos]["F_minus"],
                                        sh[pos]["F_bracket"], F_forward, g_pl, cfg)
        gate = C.AXES[family][pos[2]]
        rows.append({"family": family, "regime": regime, "depth": depth, "n": n, "position": label,
                     "layer": pos[0], "qubit": pos[1], "rotation_index": pos[2], "gate": f"R{gate}",
                     "role": "control" if label.startswith("CONTROL") else "predeclared",
                     "same_parameter_as": next((l for l, p in CF.distinct_positions(depth).items() if p == pos and l != label), ""),
                     "audit_samples": int(theta.shape[0]), **cache[pos]})
    return rows
