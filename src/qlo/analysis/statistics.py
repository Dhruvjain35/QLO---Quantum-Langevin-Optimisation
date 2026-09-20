"""Small, dependency-free statistics helpers used by the Stage 2 experiment and tests."""

from __future__ import annotations

import numpy as np

Z_95 = 1.959963984540054  # two-sided 95% normal quantile


def summarize_shot_noise(g_hat: np.ndarray, g_exact: float, shots: int, var_theory: float) -> dict[str, float]:
    """Replicate-level summary of ``xi = g_hat - g_exact`` for one ``(k, M)``."""
    g_hat = np.asarray(g_hat, dtype=np.float64)
    n = int(g_hat.size)
    if n < 2:
        raise ValueError("need >= 2 replicates")
    xi = g_hat - float(g_exact)
    mean_xi = float(np.mean(xi))
    var_emp = float(np.var(xi, ddof=1))
    std_emp = float(np.sqrt(var_emp))
    se = std_emp / np.sqrt(n)
    z = mean_xi / se if se > 0 else 0.0
    ratio = var_emp / var_theory if var_theory > 0 else float("nan")
    return {
        "n_replicates": n,
        "shots": int(shots),
        "g_exact": float(g_exact),
        "mean_g_hat": float(np.mean(g_hat)),
        "bias": float(np.mean(g_hat)) - float(g_exact),
        "mean_xi": mean_xi,
        "var_xi_emp": var_emp,
        "std_xi_emp": std_emp,
        "var_xi_theory": float(var_theory),
        "var_ratio_emp_over_theory": ratio,
        "se_mean_xi": float(se),
        "z_mean_xi": float(z),
        "ci95_low": mean_xi - Z_95 * se,
        "ci95_high": mean_xi + Z_95 * se,
        "M_times_var_emp": int(shots) * var_emp,
        "M_times_var_theory": int(shots) * float(var_theory),
    }


def fit_loglog_slope(x: np.ndarray, y: np.ndarray) -> dict[str, float]:
    """Least-squares fit ``log y = a + b log x``; returns slope ``b``, intercept ``a``, ``R^2``."""
    lx, ly = np.log(np.asarray(x, dtype=np.float64)), np.log(np.asarray(y, dtype=np.float64))
    if lx.size < 3:
        raise ValueError("need >= 3 points")
    b, a = np.polyfit(lx, ly, 1)
    pred = a + b * lx
    ss_res = float(np.sum((ly - pred) ** 2))
    ss_tot = float(np.sum((ly - ly.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return {"slope": float(b), "intercept": float(a), "r2": r2, "n_points": int(lx.size)}
