"""Per-sample required shots for both measurement schemes on the common fidelity landscape.

Targets (component k, fixed theta): SNR >= rho (rho = 1, 2), P_correct >= q (q = 0.75, 0.90), P(g_hat != 0) >= 0.90.
All values are integer budgets M >= 1 returned as log10 M, together with a per-sample method code.

Methods
-------
SNR (both)              closed form, log space (``estimators.log10_shots_for_snr``); ceil to an integer.
direction, LOSCHMIDT    exact integer bisection of the exact difference-of-binomials P_correct with a per-sample
                        bracket around max(z_q^2/(A s^2), 1/A) (Stage 6 method); rows with M*A > 1e4 (|s| tiny)
                        use the Skellam/normal limit.
direction, SWAP         exact integer bisection for rows whose normal estimate is <= 10^EXACT_MAX_LOG10 shots;
                        beyond, the continuity-corrected normal inversion in closed form
                            Phi((M d - 1/2)/sqrt(M v)) = q,  d = A|s|/2,  v = [2 - A^2(1+s^2)/2]/4
                            => sqrt(M) = [z sqrt(v) + sqrt(z^2 v + 2d)] / (2d)
                        computed from the exact d (never from q_- - q_+, which is lost next to 1/2 at large n).
non-zero, LOSCHMIDT     exact (Stage 5 ``shots_for_nonzero_exact``).
non-zero, SWAP          exact integer bisection on P_zero (the answer is O(10..100) shots: P_zero -> C(2M,M)/4^M).
"""

from __future__ import annotations

import numpy as np
from scipy.stats import norm

from qlo.stage5.theory import shots_for_nonzero_exact
from qlo.stage6.directional import (
    _bisect_log_shots,
    direction_probabilities,
    shots_for_direction_exact,
    shots_for_direction_skellam,
)
from qlo.stage6.exact_distribution import difference_probabilities
from qlo.stage7.estimators import LN10, log10_shots_for_snr, outcome_probabilities, pair

EXACT_MAX_LOG10 = 4.0   # SWAP: exact integer inversion up to ~1e4 shots
HUGE_MA = 1e4           # Loschmidt: Skellam/normal limit beyond M*A = 1e4
METHOD = {"closed_form": 0, "exact": 1, "normal_cc": 2, "skellam": 3}


def snr_shots(scheme: str, logA, s, rho) -> np.ndarray:
    """log10 of the integer budget ceil(M_SNR), floored at M = 1 (above 1e15, ceil is below float resolution)."""
    l = log10_shots_for_snr(scheme, logA, s, rho)
    small = l < 15.0
    out = np.maximum(l, 0.0)
    with np.errstate(invalid="ignore"):
        out[small] = np.log10(np.maximum(np.ceil(np.power(10.0, l[small])), 1.0))
    return out


def log10_normal_cc_direction(log_d, v, q) -> np.ndarray:
    """Continuity-corrected normal inversion for P_correct = q, from log d (d = |p_pos - p_neg|) and v."""
    z = norm.ppf(q)
    d = np.exp(np.asarray(log_d, dtype=np.float64))
    v = np.asarray(v, dtype=np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        log_x = np.log(z * np.sqrt(v) + np.sqrt(z**2 * v + 2.0 * d)) - np.log(2.0) - np.asarray(log_d)
    return np.maximum(2.0 * log_x / LN10, 0.0)


def swap_direction_shots(logA, s, q) -> tuple[np.ndarray, np.ndarray]:
    logA, s = np.asarray(logA, dtype=np.float64), np.asarray(s, dtype=np.float64)
    A = np.exp(logA)
    with np.errstate(divide="ignore"):
        log_d = logA + np.log(np.abs(s)) - np.log(2.0)
    v = (2.0 - A**2 * (1.0 + s**2) / 2.0) / 4.0
    l_norm = log10_normal_cc_direction(log_d, v, q)
    out, method = l_norm.copy(), np.full(l_norm.shape, METHOD["normal_cc"])
    todo = np.flatnonzero(np.isfinite(l_norm) & (l_norm <= EXACT_MAX_LOG10))
    p_pos, p_neg, _ = pair("swap", A, s)
    for down, up in ((1.0, 1.0), (3.0, 3.0)):
        if todo.size == 0:
            break
        lo, hi = np.maximum(l_norm[todo] - down, 0.0), np.maximum(l_norm[todo] + up, 1.0)
        m_ex, ok = shots_for_direction_exact(p_pos[todo], p_neg[todo], q, lo_log10=lo, hi_log10=hi, n_iter=24)
        at_lo = direction_probabilities(p_pos[todo], p_neg[todo], np.ceil(np.power(10.0, lo)))["p_correct"] >= q
        good = ok & (~at_lo | (lo == 0.0))
        out[todo[good]], method[todo[good]] = np.log10(m_ex[good]), METHOD["exact"]
        todo = todo[~good]
    return out, method


def loschmidt_direction_shots(logA, s, q) -> tuple[np.ndarray, np.ndarray]:
    logA, s = np.asarray(logA, dtype=np.float64), np.asarray(s, dtype=np.float64)
    A = np.exp(logA)
    p_pos, p_neg, _ = pair("loschmidt", A, s)
    z2 = norm.ppf(q) ** 2
    with np.errstate(divide="ignore"):
        l_est = np.maximum(np.log10(z2) - logA / LN10 - np.log10(s**2), -logA / LN10)
    out = np.full(A.shape, np.nan)
    method = np.full(A.shape, -1)
    huge = l_est + logA / LN10 > np.log10(HUGE_MA)
    todo = np.flatnonzero(~huge & np.isfinite(l_est))
    for down, up in ((1.5, 0.7), (4.0, 2.0)):
        if todo.size == 0:
            break
        lo, hi = np.maximum(l_est[todo] - down, 0.0), l_est[todo] + up
        m_ex, ok = shots_for_direction_exact(p_pos[todo], p_neg[todo], q, lo_log10=lo, hi_log10=hi, n_iter=24)
        at_lo = direction_probabilities(p_pos[todo], p_neg[todo], np.ceil(np.power(10.0, lo)))["p_correct"] >= q
        good = ok & (~at_lo | (lo == 0.0))
        out[todo[good]], method[todo[good]] = np.log10(m_ex[good]), METHOD["exact"]
        todo = todo[~good]
    h = np.flatnonzero(huge)
    if h.size:
        m_sk, ok = shots_for_direction_skellam(p_pos[h], p_neg[h], q, lo_log10=np.maximum(l_est[h] - 4.0, 0.0), hi_log10=l_est[h] + 2.0)
        out[h[ok]], method[h[ok]] = np.log10(m_sk[ok]), METHOD["skellam"]
    return out, method


def direction_shots(scheme: str, logA, s, q):
    return loschmidt_direction_shots(logA, s, q) if scheme == "loschmidt" else swap_direction_shots(logA, s, q)


def nonzero_shots(scheme: str, logA, s, q: float = 0.90) -> tuple[np.ndarray, np.ndarray]:
    """Smallest integer M with P(g_hat != 0) >= q (exact for both schemes)."""
    A = np.exp(np.asarray(logA, dtype=np.float64))
    s = np.asarray(s, dtype=np.float64)
    if scheme == "loschmidt":
        fp, fm = outcome_probabilities("loschmidt", A, s)
        m, ok = shots_for_nonzero_exact(fp, fm, q)
        return np.where(ok, np.log10(m), np.nan), np.where(ok, METHOD["exact"], -1)
    p_pos, p_neg, _ = pair("swap", A, s)
    fn = lambda M: 1.0 - difference_probabilities(p_pos, p_neg, M)[1]
    m, ok = _bisect_log_shots(fn, np.full(A.shape, q), 0.0, 4.0, n_iter=20)
    return np.where(ok, np.log10(m), np.nan), np.where(ok, METHOD["exact"], -1)
