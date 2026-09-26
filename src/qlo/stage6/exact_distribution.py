"""Exact distribution of a difference of two independent binomial counts.

Both Stage 5 (projector) and Stage 6 (parity) finite-shot parameter-shift estimators have
the form  g_hat proportional to (K_a - K_b)  with  K_a ~ Bin(M, p_a), K_b ~ Bin(M, p_b)
independent. Everything about zeros and directions follows from three probabilities:

    P(K_a = K_b) = sum_r Bin(r;M,p_a) Bin(r;M,p_b)
    P(K_a > K_b) = sum_r Bin(r;M,p_a) * P(K_b <= r-1)
    P(K_a < K_b) = 1 - P(=) - P(>)

Numerics
--------
* log-pmf of K_a on a window around its mean (reusing the Stage 5 windowed log-pmf,
  which is stable to M ~ 1e20 and anchored with ``betaln``);
* the CDF of K_b from the regularized incomplete beta function,
  ``P(K <= r) = I_{1-p}(M-r, r+1)``, which is exact for float M and needs no table.
Truncating K_a's window costs at most its tail mass (<= ~1e-31 at 12 sd), and the CDF
factor is bounded by 1, so the truncation error of P(>) is bounded by that tail mass.
"""

from __future__ import annotations

import numpy as np
from scipy.special import betaincc, logsumexp

from qlo.stage5.theory import _log_binom_pmf_window, _window


def binom_cdf(M, p, r) -> np.ndarray:
    """``P(K <= r)`` for ``K ~ Bin(M, p)`` = ``I_{1-p}(M-r, r+1)`` = ``1 - I_p(r+1, M-r)``; exact for float M.

    Evaluated as the complemented incomplete beta in ``p`` (``betaincc``), never through ``1 - p``:
    for p below ~1e-16 the float ``1 - p`` rounds to exactly 1 and ``I_{1-p}`` would return 1."""
    M, p, r = (np.asarray(x, dtype=np.float64) for x in (M, p, r))
    out = np.where(r < 0, 0.0, np.where(r >= M, 1.0, 0.0))
    mid = (r >= 0) & (r < M)
    if np.any(mid):
        a = np.broadcast_to(M, mid.shape)[mid] - np.broadcast_to(r, mid.shape)[mid]
        b = np.broadcast_to(r, mid.shape)[mid] + 1.0
        x = np.broadcast_to(p, mid.shape)[mid]
        out = np.array(out, dtype=np.float64)
        out[mid] = betaincc(b, a, x)
    return np.clip(out, 0.0, 1.0)


def _chunks(size: int, width: np.ndarray, max_elements: int = 5_000_000, max_rows: int = 4096):
    # ~5M float64 per table keeps each worker near 200 MB even with several live (rows, width) arrays
    order = np.argsort(width)
    i = 0
    while i < order.size:
        j = i + 1
        while j < order.size and j - i < max_rows and (j - i + 1) * int(width[order[j]]) <= max_elements:
            j += 1
        yield order[i:j]
        i = j


def difference_probabilities(p_a, p_b, M, chunk_max_rows: int = 4096) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return ``(P(K_a > K_b), P(K_a = K_b), P(K_a < K_b))``, vectorized and broadcastable."""
    pa, pb, Mv = np.broadcast_arrays(*(np.asarray(x, dtype=np.float64) for x in (p_a, p_b, M)))
    shape = pa.shape
    pa, pb, Mv = pa.ravel().copy(), pb.ravel().copy(), Mv.ravel().copy()
    gt = np.zeros(pa.size)
    eq = np.zeros(pa.size)
    lo, hi = _window(Mv, pa, pa)  # window on K_a only; the K_b factor is a bounded CDF/pmf
    width = np.maximum(hi - lo + 1, 1).astype(np.int64)
    for idx in _chunks(pa.size, width, max_rows=chunk_max_rows):
        w = int(width[idx].max())
        r = lo[idx][:, None] + np.arange(w)[None, :]
        lpa = _log_binom_pmf_window(Mv[idx], pa[idx], lo[idx], w)
        lpa = np.where(r > hi[idx][:, None], -np.inf, lpa)
        pmf_a = np.exp(lpa)
        lpb = _log_binom_pmf_window(Mv[idx], pb[idx], lo[idx], w)
        lpb = np.where(r > hi[idx][:, None], -np.inf, lpb)
        # P(K_b <= r-1) = P(K_b <= lo-1) [one incomplete-beta call per row] + running sum of K_b's pmf
        base = binom_cdf(Mv[idx], pb[idx], lo[idx] - 1.0)
        run = np.cumsum(np.exp(lpb), axis=1)
        cdf_b_lt = base[:, None] + np.concatenate([np.zeros((idx.size, 1)), run[:, :-1]], axis=1)
        gt[idx] = np.sum(pmf_a * np.minimum(cdf_b_lt, 1.0), axis=1)
        eq[idx] = np.exp(logsumexp(lpa + lpb, axis=1))
    gt = np.clip(gt, 0.0, 1.0)
    eq = np.clip(eq, 0.0, 1.0)
    lt = np.clip(1.0 - gt - eq, 0.0, 1.0)
    return gt.reshape(shape), eq.reshape(shape), lt.reshape(shape)


def difference_probabilities_bruteforce(p_a: float, p_b: float, M: int) -> tuple[float, float, float]:
    """Reference for small M: direct double sum with ``scipy.stats.binom``."""
    from scipy.stats import binom

    r = np.arange(int(M) + 1)
    fa, fb = binom.pmf(r, int(M), p_a), binom.pmf(r, int(M), p_b)
    eq = float(np.sum(fa * fb))
    gt = float(np.sum(fa * np.concatenate([[0.0], np.cumsum(fb)[:-1]])))
    return gt, eq, 1.0 - gt - eq


def p_equal_zero_signal(M) -> np.ndarray:
    """Exact ``P(K_a = K_b)`` when ``p_a = p_b = 1/2``: ``C(2M, M)/4^M`` (Vandermonde), ~ 1/sqrt(pi M)."""
    from scipy.special import gammaln

    M = np.asarray(M, dtype=np.float64)
    return np.exp(gammaln(2 * M + 1) - 2 * gammaln(M + 1) - 2 * M * np.log(2.0))


def p_equal_zero_signal_asymptotic(M) -> np.ndarray:
    """Large-M limit ``1/sqrt(pi M)``."""
    return 1.0 / np.sqrt(np.pi * np.asarray(M, dtype=np.float64))
