"""Fits and summaries for Stage 5."""

from __future__ import annotations

import numpy as np

LOG10_4 = float(np.log10(4.0))
LOG10_2 = float(np.log10(2.0))


def linear_fit(x, y) -> dict:
    x, y = np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)
    b, a = np.polyfit(x, y, 1)
    pred = a + b * x
    rss = float(np.sum((y - pred) ** 2))
    tss = float(np.sum((y - y.mean()) ** 2))
    return {"slope": float(b), "intercept": float(a), "r2": 1 - rss / tss if tss > 0 else float("nan"),
            "residuals": [float(r) for r in (y - pred)], "n_points": int(x.size)}


def scaling_fit(n, log10M, label: str) -> dict:
    """median log10 M = a + b n; compare b to log10 4 (log-typical A) and log10 2 (arithmetic-mean A)."""
    f = linear_fit(n, log10M)
    return {"statistic": label, **f, "theory_slope_log10_4": LOG10_4, "theory_slope_log10_2": LOG10_2,
            "abs_err_vs_log10_4": abs(f["slope"] - LOG10_4), "abs_err_vs_log10_2": abs(f["slope"] - LOG10_2),
            "closer_to": "log10_4" if abs(f["slope"] - LOG10_4) < abs(f["slope"] - LOG10_2) else "log10_2"}


def wilson_ci(k: int, n: int) -> tuple[float, float]:
    z = 1.959963984540054
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z**2 / n
    c = (p + z**2 / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / d
    return (float(c - h), float(c + h))
