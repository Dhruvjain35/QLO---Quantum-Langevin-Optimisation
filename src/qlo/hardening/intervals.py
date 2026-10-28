"""Interval and summary helpers. Every function names its resampling unit; none mixes uncertainty sources."""

from __future__ import annotations

import numpy as np

LEVEL = 0.95


def percentile_interval(x, level: float = LEVEL) -> tuple[float, float]:
    """Central percentile interval of the values themselves (finite values only)."""
    x = np.asarray(x, dtype=np.float64)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return float("nan"), float("nan")
    a = 100.0 * (1.0 - level) / 2.0
    lo, hi = np.percentile(x, [a, 100.0 - a])
    return float(lo), float(hi)


def bootstrap_mean_ci(x, n_boot: int = 10_000, seed: int = 0, level: float = LEVEL) -> dict:
    """Percentile bootstrap CI of the MEAN, resampling the entries of ``x`` (the caller's unit, e.g. seeds)."""
    x = np.asarray(x, dtype=np.float64)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return {"boot_lo": float("nan"), "boot_hi": float("nan"), "boot_se": float("nan")}
    rng = np.random.default_rng(np.random.SeedSequence([3, 0xB3, int(seed)]))
    means = x[rng.integers(0, x.size, (int(n_boot), x.size))].mean(axis=1)
    lo, hi = percentile_interval(means, level)
    return {"boot_lo": lo, "boot_hi": hi, "boot_se": float(means.std(ddof=1))}


def summarize(x, n_boot: int = 10_000, seed: int = 0, level: float = LEVEL) -> dict:
    """mean, sd, median, IQR, min, max, 95% percentile interval, 95% bootstrap CI of the mean."""
    x = np.asarray(x, dtype=np.float64)
    f = x[np.isfinite(x)]
    if f.size == 0:
        keys = ("mean", "sd", "median", "q25", "q75", "min", "max", "pct_lo", "pct_hi", "boot_lo", "boot_hi", "boot_se")
        return {"n_units": 0, **{k: float("nan") for k in keys}}
    lo, hi = percentile_interval(f, level)
    return {"n_units": int(f.size), "mean": float(f.mean()), "sd": float(f.std(ddof=1)) if f.size > 1 else float("nan"),
            "median": float(np.median(f)), "q25": float(np.percentile(f, 25)), "q75": float(np.percentile(f, 75)),
            "min": float(f.min()), "max": float(f.max()), "pct_lo": lo, "pct_hi": hi, **bootstrap_mean_ci(f, n_boot, seed, level)}


def paired_gap(a, b) -> np.ndarray:
    """``b - a`` unit by unit. ``a`` and ``b`` must come from the SAME units in the same order (e.g. LE and SWAP
    slopes of the same seed); the gap is formed before any resampling, so the pairing is never broken."""
    a, b = np.asarray(a, dtype=np.float64), np.asarray(b, dtype=np.float64)
    if a.shape != b.shape:
        raise ValueError("paired gap needs equal-length, unit-aligned inputs")
    return b - a


def ols_slopes(x, Y) -> tuple[np.ndarray, np.ndarray]:
    """OLS slope and intercept of every row of ``Y`` (last axis) against ``x``; identical to ``np.polyfit(x, y, 1)``."""
    x = np.asarray(x, dtype=np.float64)
    Y = np.asarray(Y, dtype=np.float64)
    xc = x - x.mean()
    slope = (Y * xc).sum(axis=-1) / (xc**2).sum()
    return slope, Y.mean(axis=-1) - slope * x.mean()


def ols_slope_se(x, y) -> float:
    """Classical OLS standard error of the slope. DIAGNOSTIC ONLY: it measures scatter about the line, not the
    initialization-sampling uncertainty that the seed and theta intervals measure."""
    x, y = np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)
    if x.size < 3:
        return float("nan")
    b, a = np.polyfit(x, y, 1)
    rss = float(np.sum((y - a - b * x) ** 2))
    return float(np.sqrt(rss / (x.size - 2) / np.sum((x - x.mean()) ** 2)))
