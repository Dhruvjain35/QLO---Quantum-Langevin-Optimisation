"""Seed-level slope fits, A2 prediction residuals, intervals (B3 utilities), concentration and outcome classes.

Slope conventions (B1 spec section 2, used everywhere):
    median_theta log10 S      = a_S  - b_S  n        (b_S  > 0: S decays)
    median_theta log10 |r|    = a_r  - b_r  n        (b_r  > 0: the relative gradient decays)
    median_theta log10 M_LE   = a_LE + b_LE n
    median_theta log10 M_SWAP = a_SW + b_SWAP n
A2 concentrated-regime prediction: b_LE = b_S + 2 b_r,  b_SWAP = 2 b_S + 2 b_r,  b_SWAP - b_LE = b_S.
Residuals per seed: eps_LE = b_LE - (b_S + 2 b_r), eps_SWAP = b_SWAP - (2 b_S + 2 b_r), eps_gap = (b_SWAP - b_LE) - b_S.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd

from qlo.b1 import config as CF
from qlo.hardening import intervals as IV
from qlo.stage6.analysis import linear_fit

SERIES = (("S", "median_log10_S", -1.0), ("r", "median_log10_abs_r", -1.0),
          ("LE", "median_log10_M_LE", 1.0), ("SWAP", "median_log10_M_SWAP", 1.0))
METRICS = ("b_S", "b_r", "b_LE", "b_SWAP", "gap_obs", "eps_LE", "eps_SWAP", "eps_gap", "ratio_obs", "ratio_pred", "ratio_diff")
CELL_KEYS = ["family", "regime", "position", "fit_range"]


def predicted_ratio(b_S, b_r):
    """(2 b_S + 2 b_r) / (b_S + 2 b_r); equals 2 exactly when b_r = 0."""
    b_S, b_r = np.asarray(b_S, dtype=np.float64), np.asarray(b_r, dtype=np.float64)
    return (2.0 * b_S + 2.0 * b_r) / (b_S + 2.0 * b_r)


def prediction_residuals(b_S, b_r, b_LE, b_SWAP) -> dict:
    """Scalars (one seed) or arrays (e.g. theta-bootstrap replicates) of the four slopes."""
    rp = predicted_ratio(b_S, b_r)
    out = {"pred_LE": b_S + 2.0 * b_r, "pred_SWAP": 2.0 * b_S + 2.0 * b_r, "gap_obs": b_SWAP - b_LE, "gap_pred": b_S,
           "eps_LE": b_LE - (b_S + 2.0 * b_r), "eps_SWAP": b_SWAP - (2.0 * b_S + 2.0 * b_r), "eps_gap": (b_SWAP - b_LE) - b_S,
           "ratio_obs": b_SWAP / b_LE, "ratio_pred": rp if np.ndim(rp) else float(rp)}
    out["ratio_diff"] = out["ratio_obs"] - out["ratio_pred"]
    return out


def fit_series(n_values, y) -> dict:
    """Stage 6 ``linear_fit`` (OLS, = np.polyfit) plus the residual vector. Non-finite medians are an error: no point
    is ever dropped from a fit silently."""
    x, y = np.asarray(n_values, dtype=np.float64), np.asarray(y, dtype=np.float64)
    if not np.all(np.isfinite(y)):
        raise ValueError(f"non-finite median in a slope fit: {y}")
    f = linear_fit(x, y)
    f["residuals"] = y - (f["intercept"] + f["slope"] * x)
    return f


def seed_fit(med: pd.DataFrame, n_range) -> dict:
    """``med``: the per-n medians of ONE (cell, seed). Returns slopes in the B1 sign convention, fit diagnostics and
    the prediction residuals."""
    d = med.set_index("n").loc[list(n_range)]
    row = {"n_min": int(min(n_range)), "n_max": int(max(n_range)), "n_points": len(n_range),
           "n_used": " ".join(str(int(n)) for n in n_range), "samples_per_n": int(d.n_rows.min())}
    for key, col, sign in SERIES:
        f = fit_series(d.index.values, d[col].values)
        row[f"b_{key}"] = sign * f["slope"]
        row[f"a_{key}"] = f["intercept"]
        row[f"r2_{key}"] = f["r2"]
        row[f"max_abs_resid_{key}"] = f["max_abs_residual"]
        row[f"resid_{key}"] = " ".join(f"{v:.6g}" for v in f["residuals"])
    row.update(prediction_residuals(row["b_S"], row["b_r"], row["b_LE"], row["b_SWAP"]))
    return row


def seed_fits(medians: pd.DataFrame, ranges_of, group_keys=("family", "regime", "position")) -> pd.DataFrame:
    """One row per (cell, fit range, seed). ``ranges_of(cell_key_tuple)`` -> {range name: n tuple}."""
    rows = []
    for key, g in medians.groupby(list(group_keys), sort=False):
        for rname, n_range in ranges_of(key).items():
            for seed, gs in g.groupby("seed", sort=True):
                rows.append({**dict(zip(group_keys, key)), "fit_range": rname, "seed": int(seed), **seed_fit(gs, n_range)})
    return pd.DataFrame(rows)


def summarize_metrics(fits: pd.DataFrame, cfg: CF.B1Config, metrics=METRICS) -> pd.DataFrame:
    """Seed mean, SD, median, IQR, 95% percentile interval and bootstrap CI of the mean (B3 ``IV.summarize``)."""
    rows = []
    for key, g in fits.groupby(CELL_KEYS, sort=False):
        g = g.sort_values("seed")
        for m in metrics:
            rows.append({**dict(zip(CELL_KEYS, key)), "metric": m, "n_seeds": int(len(g)), "n_min": int(g.n_min.iloc[0]),
                         "n_max": int(g.n_max.iloc[0]), **IV.summarize(g[m].values, cfg.n_boot_seed, cfg.boot_seed)})
    return pd.DataFrame(rows)


def ci_contains(s: dict | pd.Series, value: float = 0.0, kind: str = "boot") -> bool:
    return bool(s[f"{kind}_lo"] <= value <= s[f"{kind}_hi"])


def _halves_slopes(n_values: np.ndarray, m: np.ndarray) -> tuple[float, float]:
    h = math.ceil(len(n_values) / 2)
    lo = -linear_fit(n_values[:h], m[:h])["slope"]
    hi = -linear_fit(n_values[-h:], m[-h:])["slope"]
    return lo, hi


def concentration(medians_cell: pd.DataFrame, fits_cell: pd.DataFrame, bS_summary: pd.Series, n_range, cfg: CF.B1Config) -> dict:
    """Predeclared CLEARLY CONCENTRATED / TRANSITIONAL / NOT CLEARLY CONCENTRATED rule on the seed-mean medians."""
    m = medians_cell.groupby("n").median_log10_S.mean().loc[list(n_range)]
    nv, mv = m.index.values.astype(float), m.values
    S_med = 10.0 ** mv
    b_S = float(bS_summary["mean"])
    lo, hi = _halves_slopes(nv, mv)
    d = {"decay_ci_above_0": bool(bS_summary["boot_lo"] > 0), "monotone_decreasing": bool(np.all(np.diff(mv) < 0)),
         "mean_r2_S": float(fits_cell.r2_S.mean()), "b_S_lower_half": lo, "b_S_upper_half": hi,
         "curvature_rel": float(abs(hi - lo) / b_S) if b_S > 0 else float("inf"),
         "median_S_first_n": float(S_med[0]), "median_S_last_n": float(S_med[-1]), "median_S_max": float(S_med.max())}
    d["r2_ok"] = d["mean_r2_S"] >= cfg.conc_r2_min
    d["curvature_ok"] = d["curvature_rel"] <= cfg.conc_curvature_max
    d["small_S_all_n"] = bool(np.all(S_med <= cfg.conc_median_S_max))
    d["small_S_last_n"] = bool(S_med[-1] <= cfg.conc_median_S_max)
    if d["decay_ci_above_0"] and d["monotone_decreasing"] and d["r2_ok"] and d["curvature_ok"] and d["small_S_all_n"]:
        d["concentration"] = "CLEARLY CONCENTRATED"
    elif d["decay_ci_above_0"] and d["monotone_decreasing"] and d["small_S_last_n"]:
        d["concentration"] = "TRANSITIONAL"
    else:
        d["concentration"] = "NOT CLEARLY CONCENTRATED"
    return d


def cell_outcome(summ: dict, concentration_class: str) -> dict:
    """Predeclared per-cell labels. Rule check: CONSISTENT (eps_gap, eps_LE, eps_SWAP bootstrap CIs all contain 0),
    GAP-ONLY (only eps_gap's CI contains 0), DEVIATES (eps_gap's CI excludes 0). Case mapping: A = consistent in a
    clearly concentrated range; B = deviation outside the clearly concentrated regime (identity holds); D = deviation
    inside a clearly concentrated range; 'consistent (not clearly concentrated)' otherwise."""
    gap_ok, le_ok, sw_ok = (ci_contains(summ[m]) for m in ("eps_gap", "eps_LE", "eps_SWAP"))
    rule = "CONSISTENT" if (gap_ok and le_ok and sw_ok) else ("GAP-ONLY" if gap_ok else "DEVIATES")
    clear = concentration_class == "CLEARLY CONCENTRATED"
    if rule == "CONSISTENT":
        case = "A" if clear else "consistent (not clearly concentrated)"
    else:
        case = "D" if clear else "B"
    return {"gap_ci_contains_0": gap_ok, "le_ci_contains_0": le_ok, "swap_ci_contains_0": sw_ok, "rule_check": rule, "case": case}


def doubling(summ: dict, cfg: CF.B1Config) -> dict:
    """Section 15: test the ratio against 2 only when b_r's CI contains 0 and |b_r| is small relative to b_S."""
    br, bs = summ["b_r"], summ["b_S"]
    eligible = ci_contains(br) and abs(br["mean"]) <= cfg.br_small_fraction * bs["mean"]
    return {"br_ci_contains_0": ci_contains(br), "br_over_bS": float(br["mean"] / bs["mean"]) if bs["mean"] else float("nan"),
            "doubling_test_eligible": bool(eligible),
            "ratio_consistent_with_2": (ci_contains(summ["ratio_obs"], 2.0) if eligible else None),
            "ratio_consistent_with_pred": ci_contains(summ["ratio_diff"])}


def equivalence(s: dict | pd.Series, margin: float) -> bool:
    """Secondary practical-equivalence check: the bootstrap CI of the seed mean lies inside [-margin, +margin]."""
    return bool(-margin <= s["boot_lo"] and s["boot_hi"] <= margin)
