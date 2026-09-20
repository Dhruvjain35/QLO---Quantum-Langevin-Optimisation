"""Statistics for Stage 4 outputs (nonparametric where possible)."""

from __future__ import annotations

import numpy as np
import pandas as pd

from qlo.stage4.trajectory import FIXED_TARGETS, RELATIVE_TARGETS

Z95 = 1.959963984540054


def wilson_ci(k: int, n: int) -> tuple[float, float]:
    """95% Wilson score interval for a binomial proportion."""
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + Z95**2 / n
    c = (p + Z95**2 / (2 * n)) / d
    h = Z95 * np.sqrt(p * (1 - p) / n + Z95**2 / (4 * n**2)) / d
    return (float(c - h), float(c + h))


def summarize_group(df: pd.DataFrame) -> dict:
    """Per-(n, method, config) summary over trajectories."""
    g = df["log10_gain_best"].to_numpy()
    gf = df["log10_gain_final"].to_numpy()
    out = {
        "n_runs": int(len(df)), "n_start_seeds": int(df["start_seed"].nunique()),
        "median_log10_gain_best": float(np.median(g)), "iqr_low_gain_best": float(np.percentile(g, 25)),
        "iqr_high_gain_best": float(np.percentile(g, 75)), "mean_log10_gain_best": float(np.mean(g)),
        "median_log10_gain_final": float(np.median(gf)), "iqr_low_gain_final": float(np.percentile(gf, 25)),
        "iqr_high_gain_final": float(np.percentile(gf, 75)),
        "median_F0": float(df["F0"].median()), "median_F_best": float(df["F_best"].median()),
        "median_F_final": float(df["F_final"].median()),
    }
    for t in FIXED_TARGETS:
        el = df[df[f"eligible_F{t}"]]
        hits = el[f"hit_F{t}"].notna()
        k, m = int(hits.sum()), int(len(el))
        lo, hi = wilson_ci(k, m)
        out.update({f"n_eligible_F{t}": m, f"p_hit_F{t}": k / m if m else float("nan"),
                    f"p_hit_F{t}_ci_low": lo, f"p_hit_F{t}_ci_high": hi,
                    f"median_iter_to_F{t}": float(el.loc[hits, f"hit_F{t}"].median()) if k else float("nan")})
        for col in (f"shots_to_F{t}", f"circuit_evals_to_F{t}"):
            v = el.loc[hits, col].dropna()
            out.update({f"median_{col}": float(v.median()) if len(v) else float("nan"),
                        f"min_{col}": float(v.min()) if len(v) else float("nan"),
                        f"max_{col}": float(v.max()) if len(v) else float("nan")})
    for r in RELATIVE_TARGETS:
        po = df[df[f"possible_{int(r)}xF0"]]
        hits = po[f"hit_{int(r)}xF0"].notna()
        k, m = int(hits.sum()), int(len(po))
        lo, hi = wilson_ci(k, m)
        out.update({f"p_{int(r)}xF0": k / m if m else float("nan"), f"p_{int(r)}xF0_ci_low": lo, f"p_{int(r)}xF0_ci_high": hi,
                    f"median_iter_to_{int(r)}xF0": float(po.loc[hits, f"hit_{int(r)}xF0"].median()) if k else float("nan")})
    return out


def per_seed_mean(df: pd.DataFrame, col: str = "log10_gain_best") -> pd.Series:
    """Average over stochastic replicates -> one value per start seed."""
    return df.groupby("start_seed")[col].mean()


def paired_bootstrap(a: np.ndarray, b: np.ndarray, n_boot: int = 10_000, seed: int = 0) -> dict:
    """Paired differences ``a - b`` (same start seeds). Bootstrap CIs for the median and mean
    difference, fraction of pairs with a > b, and an exact two-sided sign-test p-value."""
    d = np.asarray(a, dtype=np.float64) - np.asarray(b, dtype=np.float64)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, d.size, size=(n_boot, d.size))
    med = np.median(d[idx], axis=1)
    mean = np.mean(d[idx], axis=1)
    from scipy.stats import binomtest

    nz = d[d != 0]
    p_sign = binomtest(int((nz > 0).sum()), int(nz.size), 0.5).pvalue if nz.size else float("nan")
    return {"n_pairs": int(d.size), "median_diff": float(np.median(d)),
            "median_diff_ci95": [float(np.percentile(med, 2.5)), float(np.percentile(med, 97.5))],
            "mean_diff": float(np.mean(d)), "mean_diff_ci95": [float(np.percentile(mean, 2.5)), float(np.percentile(mean, 97.5))],
            "frac_a_greater": float(np.mean(d > 0)), "frac_equal": float(np.mean(d == 0)), "sign_test_p": float(p_sign)}
