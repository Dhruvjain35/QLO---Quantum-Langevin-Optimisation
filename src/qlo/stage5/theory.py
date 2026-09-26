"""Exact finite-shot parameter-shift estimator distribution for the projector benchmark.

Benchmark (frozen from Stage 3/4):  F(theta) = prod_j cos^2(theta_j/2),  C = 1 - F.
For component k write  A = A_k = prod_{j!=k} cos^2(theta_j/2),  s = sin(theta_k).  Then

    F_+ = F(theta + pi/2 e_k) = A (1 - s)/2        F_- = A (1 + s)/2          F_+ + F_- = A
    g_k = dC/dtheta_k = A s / 2                    F_+^2 + F_-^2 = A^2 (1 + s^2)/2

With K_+ ~ Bin(M, F_+), K_- ~ Bin(M, F_-) independent and C_hat = 1 - K/M,

    g_hat_k = (K_- - K_+) / (2M)                                              (exact identity)
    P(g_hat_k = 0 | theta, M) = sum_r Bin(r; M, F_+) Bin(r; M, F_-)           (exact)
    P_both0 = (1-F_+)^M (1-F_-)^M  <=  P_zero                                  (lower bound)
    Var(g_hat_k | theta) = [A - A^2 (1+s^2)/2] / (4M)                          (exact)
    SNR_k^2 = M A s^2 / [1 - A (1+s^2)/2]                                      (exact, s != 0)
    M_SNR(rho) = rho^2 [1 - A(1+s^2)/2] / (A s^2)                              (continuous; ceil for integer)

Deep-plateau approximations (A << 1), labelled as such everywhere:
    Poisson limit:  K_± ~ Poisson(M F_±)  =>  P_zero ~= exp(-M A) I_0(M A |cos theta_k|)
    rare-event:     P_zero ~= P_both0 ~= exp(-M A);   M_nonzero(q) ~= -log(1-q) / A

Numerics: binomial log-pmfs are built from the cumulative recurrence
log C(M, r) = sum_{i<=r} log((M-i+1)/i), which is stable for M up to ~1e300 and r << M,
and the sum over r is truncated to the region where both pmfs are non-negligible.
"""

from __future__ import annotations

import numpy as np
from scipy.special import betaln, i0e, logsumexp

LOG2 = float(np.log(2.0))
E_LOG_COS2_HALF = -2.0 * LOG2  # E_theta[log cos^2(theta/2)], theta ~ U[-pi, pi]
E_COS2_HALF = 0.5  # E_theta[cos^2(theta/2)]


# ----------------------------------------------------------------------------------------
# Landscape quantities for component k
# ----------------------------------------------------------------------------------------
def a_k(theta: np.ndarray, k: int) -> np.ndarray:
    """``A_k = prod_{j!=k} cos^2(theta_j/2)``; theta of shape (..., n)."""
    c2 = np.cos(np.asarray(theta, dtype=np.float64) / 2.0) ** 2
    return np.prod(np.delete(c2, k, axis=-1), axis=-1)


def log_a_k(theta: np.ndarray, k: int) -> np.ndarray:
    """``log A_k`` computed as a sum of logs (no underflow for large n)."""
    with np.errstate(divide="ignore"):
        lc2 = np.log(np.cos(np.asarray(theta, dtype=np.float64) / 2.0) ** 2)
    return np.sum(np.delete(lc2, k, axis=-1), axis=-1)


def s_k(theta: np.ndarray, k: int) -> np.ndarray:
    return np.sin(np.asarray(theta, dtype=np.float64)[..., k])


def f_plus_minus(A, s) -> tuple[np.ndarray, np.ndarray]:
    """``F_+ = A(1-s)/2``, ``F_- = A(1+s)/2``."""
    A, s = np.asarray(A, dtype=np.float64), np.asarray(s, dtype=np.float64)
    return A * (1.0 - s) / 2.0, A * (1.0 + s) / 2.0


def exact_gradient_k(A, s) -> np.ndarray:
    """``g_k = A s / 2``."""
    return np.asarray(A, dtype=np.float64) * np.asarray(s, dtype=np.float64) / 2.0


def estimator_from_counts(k_plus, k_minus, M) -> np.ndarray:
    """``g_hat_k = (K_- - K_+) / (2M)``; support is {j/(2M): j = -M..M}."""
    return (np.asarray(k_minus, dtype=np.float64) - np.asarray(k_plus, dtype=np.float64)) / (2.0 * np.asarray(M, dtype=np.float64))


# ----------------------------------------------------------------------------------------
# Conditional variance and SNR (exact)
# ----------------------------------------------------------------------------------------
def shot_variance_from_F(f_plus, f_minus, M) -> np.ndarray:
    """Stage 4 form: ``[F_+(1-F_+) + F_-(1-F_-)] / (4M)``."""
    fp, fm = np.asarray(f_plus, dtype=np.float64), np.asarray(f_minus, dtype=np.float64)
    return (fp * (1 - fp) + fm * (1 - fm)) / (4.0 * np.asarray(M, dtype=np.float64))


def shot_variance(A, s, M) -> np.ndarray:
    """Simplified exact form: ``[A - A^2 (1+s^2)/2] / (4M)``."""
    A, s = np.asarray(A, dtype=np.float64), np.asarray(s, dtype=np.float64)
    return (A - A**2 * (1.0 + s**2) / 2.0) / (4.0 * np.asarray(M, dtype=np.float64))


def snr_squared(A, s, M) -> np.ndarray:
    """``SNR^2 = M A s^2 / [1 - A(1+s^2)/2]``. Returns 0 where s = 0 or A = 0; inf where the
    denominator is 0 (A = 1 and s^2 = 1, i.e. both F_± in {0,1}: a deterministic estimator)."""
    A, s, M = (np.asarray(x, dtype=np.float64) for x in (A, s, M))
    denom = 1.0 - A * (1.0 + s**2) / 2.0
    num = M * A * s**2
    with np.errstate(divide="ignore", invalid="ignore"):
        out = np.where(denom > 0, num / denom, np.where(num > 0, np.inf, 0.0))
    return out


def shots_for_snr(A, s, rho) -> np.ndarray:
    """Continuous ``M_SNR = rho^2 [1 - A(1+s^2)/2] / (A s^2)``; inf where s = 0 or A = 0; 0 where the
    numerator is 0 (deterministic estimator). Integer budget = ceil."""
    A, s, rho = (np.asarray(x, dtype=np.float64) for x in (A, s, rho))
    num = rho**2 * (1.0 - A * (1.0 + s**2) / 2.0)
    den = A * s**2
    with np.errstate(divide="ignore", invalid="ignore"):
        out = np.where(den > 0, num / den, np.where(num > 0, np.inf, 0.0))
    return out


def log10_shots_for_snr(logA, s, rho) -> np.ndarray:
    """Log-space version for tiny A: ``log10 M_SNR = 2 log10 rho + log10(1 - A(1+s^2)/2) - logA/ln10 - log10 s^2``."""
    logA, s, rho = (np.asarray(x, dtype=np.float64) for x in (logA, s, rho))
    A = np.exp(logA)
    with np.errstate(divide="ignore", invalid="ignore"):
        return (2 * np.log10(rho) + np.log10(1.0 - A * (1.0 + s**2) / 2.0) - logA / np.log(10.0) - np.log10(s**2))


# ----------------------------------------------------------------------------------------
# Exact P(g_hat = 0)
# ----------------------------------------------------------------------------------------
def _log_binom_pmf_window(M: np.ndarray, p: np.ndarray, r_lo: np.ndarray, width: int) -> np.ndarray:
    """``log Bin(r; M, p)`` for r = r_lo .. r_lo+width-1 (shape (B, width)). The first column is anchored with
    ``betaln`` (exact 0 when r_lo = 0); subsequent columns use the exact recurrence
    log Bin(r+1) - log Bin(r) = log((M-r)/(r+1)) + log(p/(1-p)). M may be a float (large budgets)."""
    M = np.asarray(M, dtype=np.float64)[:, None]
    p = np.asarray(p, dtype=np.float64)[:, None]
    r = np.asarray(r_lo, dtype=np.float64)[:, None] + np.arange(width, dtype=np.float64)[None, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        logp, log1mp = np.log(p), np.log1p(-p)
        r0 = r[:, :1]
        anchor = np.where(r0 == 0, 0.0, -np.log(M + 1.0) - betaln(M - r0 + 1.0, r0 + 1.0)) + r0 * logp + (M - r0) * log1mp
        inc = np.log(np.maximum(M - r[:, :-1], 0.0)) - np.log(r[:, :-1] + 1.0) + logp - log1mp
        out = np.concatenate([anchor, anchor + np.cumsum(inc, axis=1)], axis=1)
    out = np.where(r > M, -np.inf, out)
    out = np.where(p <= 0.0, np.where(r == 0, 0.0, -np.inf), out)
    out = np.where(p >= 1.0, np.where(r == M, 0.0, -np.inf), out)
    return out


def _window(M, p1, p2, n_sd: float = 12.0, pad: int = 40) -> tuple[np.ndarray, np.ndarray]:
    """Summation window [r_lo, r_hi] = intersection of the two pmfs' (mean ± n_sd sd + pad) windows, clipped to [0, M].
    Outside it the product of the two pmfs is below ~exp(-n_sd^2/2) ~ 1e-31 relative."""
    M = np.asarray(M, dtype=np.float64)
    lo, hi = np.zeros_like(M), M.copy()
    for p in (p1, p2):
        p = np.asarray(p, dtype=np.float64)
        sd = np.sqrt(np.maximum(M * p * (1 - p), 0.0))
        lo = np.maximum(lo, np.floor(M * p - n_sd * sd - pad))
        hi = np.minimum(hi, np.ceil(M * p + n_sd * sd + pad))
    return np.maximum(lo, 0.0), np.minimum(hi, M)


def _p_zero_gaussian_overlap(M, p1, p2) -> np.ndarray:
    """Large-count fallback (both M p >> 1e4 and M > 1e13, where lgamma anchoring loses precision):
    P(K_+ = K_-) ~= (2 pi (s1^2 + s2^2))^{-1/2} exp(-(mu1 - mu2)^2 / (2 (s1^2 + s2^2)))."""
    M, p1, p2 = (np.asarray(x, dtype=np.float64) for x in (M, p1, p2))
    v = M * p1 * (1 - p1) + M * p2 * (1 - p2)
    return np.exp(-(M * (p1 - p2)) ** 2 / (2 * v)) / np.sqrt(2 * np.pi * v)


P_ZERO_GAUSSIAN_FALLBACK_M = 1e13
_MAX_TABLE_ELEMENTS = 20_000_000


def p_zero_exact(f_plus, f_minus, M, chunk: int = 4096) -> np.ndarray:
    """``P(K_+ = K_-) = sum_r Bin(r;M,F_+) Bin(r;M,F_-)``, vectorized over broadcastable inputs.
    Exact (to float rounding) whenever the summation window starts at r = 0 — which is the case for every
    required-shot solve in this study, since there M F_± = O(1..100). Windows starting above 0 are anchored with
    ``betaln`` (log error ~1e-16 * M log M, i.e. < 1e-3 for M <= 1e12); for M > 1e13 with such windows the
    Gaussian-overlap formula is used instead (not needed by any reported quantity)."""
    fp, fm, M = np.broadcast_arrays(*(np.asarray(x, dtype=np.float64) for x in (f_plus, f_minus, M)))
    shape = fp.shape
    fp, fm, M = fp.ravel(), fm.ravel(), M.ravel()
    out = np.zeros(fp.size)
    lo, hi = _window(M, fp, fm)
    width = (hi - lo + 1).astype(np.int64)
    empty = width <= 0
    gauss = (~empty) & (M > P_ZERO_GAUSSIAN_FALLBACK_M) & (lo > 0)
    if gauss.any():
        out[gauss] = _p_zero_gaussian_overlap(M[gauss], fp[gauss], fm[gauss])
    pending = np.flatnonzero(~empty & ~gauss)
    pending = pending[np.argsort(width[pending])]  # group similar widths to keep tables small
    i = 0
    while i < pending.size:
        w = int(width[pending[i]])
        # grow the chunk while the widest row keeps the table within the element budget
        j = i + 1
        while j < pending.size and j - i < chunk and (j - i + 1) * int(width[pending[j]]) <= _MAX_TABLE_ELEMENTS:
            j += 1
        idx = pending[i:j]
        w = int(width[idx].max())
        lp = _log_binom_pmf_window(M[idx], fp[idx], lo[idx], w) + _log_binom_pmf_window(M[idx], fm[idx], lo[idx], w)
        r = lo[idx][:, None] + np.arange(w)[None, :]
        lp = np.where(r > hi[idx][:, None], -np.inf, lp)
        out[idx] = np.exp(logsumexp(lp, axis=1))
        i = j
    return np.clip(out, 0.0, 1.0).reshape(shape)


def p_zero_bruteforce(f_plus: float, f_minus: float, M: int) -> float:
    """Reference: direct sum with scipy.stats.binom (small M only)."""
    from scipy.stats import binom

    r = np.arange(int(M) + 1)
    return float(np.sum(binom.pmf(r, int(M), f_plus) * binom.pmf(r, int(M), f_minus)))


def p_both_zero(f_plus, f_minus, M) -> np.ndarray:
    """``(1-F_+)^M (1-F_-)^M`` — exact probability that both counts are zero; a LOWER bound on P_zero."""
    fp, fm, M = (np.asarray(x, dtype=np.float64) for x in (f_plus, f_minus, M))
    return np.exp(M * (np.log1p(-fp) + np.log1p(-fm)))


def p_zero_poisson_limit(A, s, M) -> np.ndarray:
    """Deep-plateau (A << 1) approximation via K_± ~ Poisson(M F_±):
    ``P_zero ~= exp(-MA) I_0(MA |cos theta_k|)``. Evaluated stably with ``i0e``."""
    A, s, M = (np.asarray(x, dtype=np.float64) for x in (A, s, M))
    x = M * A
    c = np.sqrt(np.clip(1.0 - s**2, 0.0, 1.0))
    return i0e(x * c) * np.exp(-x * (1.0 - c))


def p_zero_rare_event(A, M) -> np.ndarray:
    """Cruder deep-plateau approximation ``exp(-M A)`` (ignores equal non-zero counts)."""
    return np.exp(-np.asarray(M, dtype=np.float64) * np.asarray(A, dtype=np.float64))


# ----------------------------------------------------------------------------------------
# Required shots for P(g_hat != 0) >= q
# ----------------------------------------------------------------------------------------
def shots_for_nonzero_rare_event(A, q) -> np.ndarray:
    """Asymptotic/deep-plateau approximation ``M ~= -log(1-q)/A``."""
    return -np.log1p(-np.asarray(q, dtype=np.float64)) / np.asarray(A, dtype=np.float64)


def log10_shots_for_nonzero_rare_event(logA, q) -> np.ndarray:
    return np.log10(-np.log1p(-np.asarray(q, dtype=np.float64))) - np.asarray(logA, dtype=np.float64) / np.log(10.0)


def shots_for_nonzero_poisson_limit(A, s, q, iters: int = 80) -> np.ndarray:
    """Solve ``exp(-x) I_0(x|c|) = 1 - q`` for x = MA in the Poisson limit (vector bisection), return M = x/A."""
    A, s, q = np.broadcast_arrays(*(np.asarray(x, dtype=np.float64) for x in (A, s, q)))
    c = np.sqrt(np.clip(1.0 - s**2, 0.0, 1.0))
    target = 1.0 - q
    lo = np.zeros_like(A)
    hi = np.full_like(A, 1.0)
    f = lambda x: i0e(x * c) * np.exp(-x * (1.0 - c))
    for _ in range(200):  # grow bracket
        bad = f(hi) > target
        if not bad.any():
            break
        hi = np.where(bad, hi * 2.0, hi)
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        above = f(mid) > target
        lo, hi = np.where(above, mid, lo), np.where(above, hi, mid)
    return hi / A


def shots_for_nonzero_exact(f_plus, f_minus, q, iters: int = 64, chunk: int = 4096) -> tuple[np.ndarray, np.ndarray]:
    """Smallest M (integer where representable) with ``1 - P_zero(M) >= q``, by vector bisection on the exact
    P_zero, assuming P_zero is non-increasing in M (checked numerically in the experiment).
    Returns (M, bracket_ok)."""
    fp, fm, q = np.broadcast_arrays(*(np.asarray(x, dtype=np.float64) for x in (f_plus, f_minus, q)))
    A = fp + fm
    target = 1.0 - q
    # bracket: rare-event guess x 64 and a slow-decay safeguard for s ~ 0 (P_zero ~ 1/sqrt(2 pi M A))
    with np.errstate(divide="ignore"):
        hi = np.maximum(64.0 * shots_for_nonzero_rare_event(A, q), 4.0 / (2 * np.pi * target**2 * A))
    hi = np.ceil(np.maximum(hi, 64.0))
    lo = np.ones_like(hi)
    ok = p_zero_exact(fp, fm, hi, chunk) <= target
    for _ in range(6):
        if ok.all():
            break
        hi = np.where(ok, hi, hi * 8.0)
        ok = p_zero_exact(fp, fm, hi, chunk) <= target
    for _ in range(iters):
        mid = np.floor(0.5 * (lo + hi))
        done = hi - lo <= 1.0
        if done.all():
            break
        pz = p_zero_exact(fp, fm, mid, chunk)
        good = pz <= target
        hi = np.where(~done & good, mid, hi)
        lo = np.where(~done & ~good, np.maximum(mid, lo + 1.0), lo)
    # lo is infeasible (or 1), hi feasible; the answer is hi unless M=1 already works
    one_ok = p_zero_exact(fp, fm, np.ones_like(hi), chunk) <= target
    return np.where(one_ok, 1.0, hi), ok


# ----------------------------------------------------------------------------------------
# Typical-scale constants
# ----------------------------------------------------------------------------------------
def e_log_a(n_qubits: int) -> float:
    """``E[log A_k] = -2 (n-1) log 2``."""
    return E_LOG_COS2_HALF * (n_qubits - 1)


def geometric_typical_a(n_qubits: int) -> float:
    """``exp(E[log A_k]) = 4^{-(n-1)}``."""
    return 4.0 ** (-(n_qubits - 1))


def arithmetic_mean_a(n_qubits: int) -> float:
    """``E[A_k] = (1/2)^{n-1}`` (distinct from the geometric scale)."""
    return 0.5 ** (n_qubits - 1)


# ----------------------------------------------------------------------------------------
# Matched term-wise LOCAL estimator (measurement-efficient, exploits the known decomposition)
# ----------------------------------------------------------------------------------------
def local_gradient_k(s, n_qubits) -> np.ndarray:
    """``g_{L,k} = sin(theta_k) / (2n)``."""
    return np.asarray(s, dtype=np.float64) / (2.0 * n_qubits)


def local_p_plus_minus(s) -> tuple[np.ndarray, np.ndarray]:
    """Shifted single-qubit |0> probabilities ``p_± = (1 ∓ s)/2`` for the term-wise estimator."""
    s = np.asarray(s, dtype=np.float64)
    return (1.0 - s) / 2.0, (1.0 + s) / 2.0


def local_estimator_from_counts(k_plus, k_minus, M, n_qubits) -> np.ndarray:
    """``g_hat_{L,k} = (K_- - K_+) / (2 M n)``."""
    return estimator_from_counts(k_plus, k_minus, M) / n_qubits


def local_variance(s, M, n_qubits) -> np.ndarray:
    """``Var = cos^2(theta_k) / (8 M n^2)``."""
    s = np.asarray(s, dtype=np.float64)
    return (1.0 - s**2) / (8.0 * np.asarray(M, dtype=np.float64) * n_qubits**2)


def local_snr_squared(s, M) -> np.ndarray:
    """``SNR_L^2 = 2 M tan^2(theta_k)`` — no n dependence for this estimator. inf where cos = 0."""
    s, M = np.asarray(s, dtype=np.float64), np.asarray(M, dtype=np.float64)
    with np.errstate(divide="ignore"):
        return 2.0 * M * s**2 / (1.0 - s**2)


def local_shots_for_snr(s, rho) -> np.ndarray:
    """``M = rho^2 (1 - s^2) / (2 s^2)``; inf at s = 0."""
    s, rho = np.asarray(s, dtype=np.float64), np.asarray(rho, dtype=np.float64)
    with np.errstate(divide="ignore"):
        return rho**2 * (1.0 - s**2) / (2.0 * s**2)
