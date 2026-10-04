"""A2 check: does the Loschmidt vs SWAP shot gap survive without the product form?

Standalone NumPy statevector code (no PennyLane, no repo imports) so it can run anywhere.

For a fidelity objective F(theta) = |<0^n|U(theta)|0^n>|^2 and a gate whose generator is
a Pauli over 2 (RX, RY, RZ), the two term parameter shift rule is exact. With shifted
fidelities F+ and F- and one independent batch of M shots per shift:

  Loschmidt (projector):  Var = [F+(1-F+) + F-(1-F-)] / (4M)
  SWAP (ancilla Z):       Var = [(1-F+^2) + (1-F-^2)] / (4M)
  exact gradient          g   = (F- - F+) / 2   (for C = 1 - F)

so the shots for SNR >= rho are

  M_LE = rho^2 [F+(1-F+) + F-(1-F-)] / (F- - F+)^2
  M_SW = rho^2 [2 - F+^2 - F-^2]     / (F- - F+)^2

These identities need no product structure. The script measures how the typical
(median) values scale with n for three circuit families.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import binom

OUT = Path(__file__).resolve().parent


# ---------------------------------------------------------------- statevector core
def apply_1q(psi, mat, q, n):
    """psi: (B, 2**n). mat: (B, 2, 2). Applies mat on qubit q (qubit 0 = most significant)."""
    B = psi.shape[0]
    psi = psi.reshape(B, 2 ** q, 2, 2 ** (n - q - 1))
    psi = np.einsum("bij,bajc->baic", mat, psi)
    return psi.reshape(B, 2 ** n)


def apply_cnot(psi, c, t, n):
    B = psi.shape[0]
    psi = psi.reshape((B,) + (2,) * n).copy()
    idx = [slice(None)] * (n + 1)
    idx[1 + c] = 1
    sub = psi[tuple(idx)]
    t_axis = 1 + t - (1 if t > c else 0)
    psi[tuple(idx)] = np.flip(sub, axis=t_axis)
    return psi.reshape(B, 2 ** n)


def rx(a):
    c, s = np.cos(a / 2), np.sin(a / 2)
    m = np.empty(a.shape + (2, 2), complex)
    m[..., 0, 0] = c; m[..., 0, 1] = -1j * s; m[..., 1, 0] = -1j * s; m[..., 1, 1] = c
    return m


def ry(a):
    c, s = np.cos(a / 2), np.sin(a / 2)
    m = np.empty(a.shape + (2, 2), complex)
    m[..., 0, 0] = c; m[..., 0, 1] = -s; m[..., 1, 0] = s; m[..., 1, 1] = c
    return m


def rz(a):
    m = np.zeros(a.shape + (2, 2), complex)
    m[..., 0, 0] = np.exp(-0.5j * a); m[..., 1, 1] = np.exp(0.5j * a)
    return m


# ---------------------------------------------------------------- circuit families
def fidelity_rx_product(theta, n):
    # theta: (B, n). Closed form, used as the regression arm.
    return np.prod(np.cos(theta / 2) ** 2, axis=1)


def fidelity_hea(theta, n, depth):
    """theta: (B, depth, n, 2). Layer = RY then RZ on every qubit, then CNOT ring
    (chain for n = 2). Same layout as qlo.circuits.hardware_efficient."""
    B = theta.shape[0]
    psi = np.zeros((B, 2 ** n), complex)
    psi[:, 0] = 1.0
    pairs = [(i, i + 1) for i in range(n - 1)] + ([(n - 1, 0)] if n > 2 else [])
    for l in range(depth):
        for q in range(n):
            psi = apply_1q(psi, ry(theta[:, l, q, 0]), q, n)
            psi = apply_1q(psi, rz(theta[:, l, q, 1]), q, n)
        for c, t in pairs:
            psi = apply_cnot(psi, c, t, n)
    return np.abs(psi[:, 0]) ** 2


def shifted_fidelities(family, n, B, rng, depth=None):
    """Returns F+, F- for the parameter shifted (first rotation on qubit 0)."""
    if family == "rx_product":
        th = rng.uniform(-np.pi, np.pi, size=(B, n))
        tp, tm = th.copy(), th.copy()
        tp[:, 0] += np.pi / 2; tm[:, 0] -= np.pi / 2
        return fidelity_rx_product(tp, n), fidelity_rx_product(tm, n)
    th = rng.uniform(-np.pi, np.pi, size=(B, depth, n, 2))
    tp, tm = th.copy(), th.copy()
    tp[:, 0, 0, 0] += np.pi / 2; tm[:, 0, 0, 0] -= np.pi / 2
    return fidelity_hea(tp, n, depth), fidelity_hea(tm, n, depth)


# ---------------------------------------------------------------- shot formulas
def m_snr(Fp, Fm, rho=1.0):
    d2 = (Fm - Fp) ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        le = rho ** 2 * (Fp * (1 - Fp) + Fm * (1 - Fm)) / d2
        sw = rho ** 2 * (2 - Fp ** 2 - Fm ** 2) / d2
    return le, sw


def conditional_sign_exact(Fp, Fm, M):
    """Exact P(sign(ghat) = sign(g) | ghat != 0) for the Loschmidt estimator,
    ghat = (K- - K+)/(2M), K+- ~ Bin(M, F+-)."""
    r = np.arange(M + 1)
    out = np.empty(len(Fp))
    for i, (a, b) in enumerate(zip(Fp, Fm)):
        pp = binom.pmf(r, M, a); pm = binom.pmf(r, M, b)
        # both tails summed directly; 1 - P(tie) loses all precision when M(F+ + F-) << 1e-8
        p_minus_gt = float(np.dot(pp, binom.sf(r, M, b)))   # K- > K+
        p_plus_gt = float(np.dot(pm, binom.sf(r, M, a)))    # K+ > K-
        nonzero = p_minus_gt + p_plus_gt
        correct = p_minus_gt if b > a else p_plus_gt
        out[i] = correct / nonzero if nonzero > 0 else np.nan
    return out


def fit(ns, ys):
    ns, ys = np.asarray(ns, float), np.asarray(ys, float)
    b, a = np.polyfit(ns, ys, 1)
    pred = a + b * ns
    r2 = 1 - np.sum((ys - pred) ** 2) / np.sum((ys - ys.mean()) ** 2)
    return float(b), float(a), float(r2)


def run(family, ns, B, depth_rule, seed=20261004, chunk=500):
    rows = []
    for n in ns:
        rng = np.random.default_rng([seed, n, len(family)])
        depth = None if depth_rule is None else depth_rule(n)
        Fp_all, Fm_all = [], []
        for start in range(0, B, chunk):
            Fp, Fm = shifted_fidelities(family, n, min(chunk, B - start), rng, depth)
            Fp_all.append(Fp); Fm_all.append(Fm)
        Fp, Fm = np.concatenate(Fp_all), np.concatenate(Fm_all)
        S = Fp + Fm
        rel = np.abs(Fm - Fp) / S                     # relative gradient |r|
        le, sw = m_snr(Fp, Fm)
        ok = np.isfinite(le) & np.isfinite(sw) & (rel > 0)
        ratio = sw[ok] / le[ok]
        # conditional sign law check on a subsample, M = 1024
        sub = np.flatnonzero(ok)[:300]
        cs = conditional_sign_exact(Fp[sub], Fm[sub], 1024)
        law = (1 + rel[sub]) / 2
        good = np.isfinite(cs)
        deep = good & (1024 * S[sub] < 0.05)       # rare event regime, where the law is a limit
        rows.append(dict(
            family=family, n=n, depth=depth, B=int(ok.sum()),
            med_log10_S=float(np.median(np.log10(S[ok]))),
            med_log10_rel=float(np.median(np.log10(rel[ok]))),
            med_log10_M_LE=float(np.median(np.log10(le[ok]))),
            med_log10_M_SW=float(np.median(np.log10(sw[ok]))),
            med_log10_ratio=float(np.median(np.log10(ratio))),
            med_ratio_times_S_over_2=float(np.median(ratio * S[ok] / 2)),
            condsign_deep_count=int(deep.sum()),
            condsign_deep_max_abs_err=(float(np.max(np.abs(cs[deep] - law[deep]))) if deep.any() else None),
            condsign_median_exact=float(np.median(cs[good])) if good.any() else None,
            condsign_median_law=float(np.median(law[good])) if good.any() else None,
        ))
        r = rows[-1]
        print(f"{family:12s} n={n:2d} L={depth} S={r['med_log10_S']:.3f} rel={r['med_log10_rel']:.3f} "
              f"LE={r['med_log10_M_LE']:.3f} SW={r['med_log10_M_SW']:.3f} "
              f"ratio*S/2={r['med_ratio_times_S_over_2']:.4f} deep={r['condsign_deep_count']} "
              f"err={r['condsign_deep_max_abs_err']}", flush=True)
    return rows


def summarize(rows, fit_from):
    ns = [r["n"] for r in rows if r["n"] >= fit_from]
    pick = lambda key: [r[key] for r in rows if r["n"] >= fit_from]
    out = {}
    for key in ("med_log10_S", "med_log10_rel", "med_log10_M_LE", "med_log10_M_SW"):
        b, a, r2 = fit(ns, pick(key))
        out[key] = dict(slope=b, intercept=a, r2=r2)
    out["slope_ratio_SW_over_LE"] = out["med_log10_M_SW"]["slope"] / out["med_log10_M_LE"]["slope"]
    out["predicted_SW_slope_from_LE_plus_S"] = out["med_log10_M_LE"]["slope"] - out["med_log10_S"]["slope"]
    return out


if __name__ == "__main__":
    B = int(sys.argv[1]) if len(sys.argv) > 1 else 4000
    results = {}
    configs = [
        ("rx_product", list(range(2, 21, 2)), None, 4),
        ("hea_shallow", list(range(2, 13)), lambda n: 2, 4),
        ("hea_deep", list(range(2, 13)), lambda n: n, 6),
    ]
    for family, ns, rule, fit_from in configs:
        fam = family.split("_")[0] if family != "rx_product" else family
        rows = run("rx_product" if family == "rx_product" else "hea", ns, B, rule)
        for r in rows:
            r["family"] = family
        results[family] = dict(rows=rows, fit_from_n=fit_from, fits=summarize(rows, fit_from))
        print(json.dumps(results[family]["fits"], indent=1), flush=True)
    (OUT / "general_scaling_results.json").write_text(json.dumps(results, indent=1))
    print("wrote", OUT / "general_scaling_results.json")
