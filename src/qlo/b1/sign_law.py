"""Conditional sign-law spot check (secondary): Loschmidt P(correct sign | g_hat != 0) -> (1 + |r|)/2 for M S << 1.

Exact binomial-difference probabilities from the B3 cancellation-free helper (``direction_accurate``: both tails summed
directly, never 1 - P(>) - P(=)). Loschmidt counts follow the Stage 7 convention: g_hat = (K_- - K_+)/(2M), so
p_pos = F_- and p_neg = F_+. Nothing is fitted: the observed exact probability is compared with (1 + |r|)/2.
"""

from __future__ import annotations

import numpy as np

from qlo.hardening.seed_replicates import direction_accurate


def quantile_rows(S: np.ndarray, quantiles) -> list[tuple[float, int]]:
    """For each predeclared quantile of S, the theta row whose S is closest to that quantile (first on ties)."""
    S = np.asarray(S, dtype=np.float64)
    return [(float(q), int(np.argmin(np.abs(S - np.quantile(S, q))))) for q in quantiles]


def sign_law_rows(F_plus: np.ndarray, F_minus: np.ndarray, quantiles, shots, meta: dict) -> list[dict]:
    rows = []
    for q, i in quantile_rows(F_plus + F_minus, quantiles):
        fp, fm = float(F_plus[i]), float(F_minus[i])
        S = fp + fm
        r = (fm - fp) / S
        M = np.asarray(shots, dtype=np.float64)
        pr = direction_accurate(np.full(M.size, fm), np.full(M.size, fp), M)
        pred = (1.0 + abs(r)) / 2.0
        for j, m in enumerate(M):
            rows.append({**meta, "S_quantile": q, "theta_row": i, "F_plus": fp, "F_minus": fm, "S": S, "r": r, "abs_r": abs(r),
                         "M": int(m), "M_times_S": m * S, "p_nonzero": float(pr["p_nonzero"][j]),
                         "p_correct_given_nonzero": float(pr["p_correct_given_nonzero"][j]), "predicted": pred,
                         "error": float(pr["p_correct_given_nonzero"][j] - pred)})
    return rows
