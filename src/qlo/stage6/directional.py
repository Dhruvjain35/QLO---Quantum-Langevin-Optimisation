"""Directional correctness of a finite-shot parameter-shift estimator.

Both benchmarks give ``g_hat ∝ (K_pos - K_neg)`` with independent binomial counts, where
``g_true ∝ (p_pos - p_neg)``:

    parity     p_pos = p_+ = (1 + sB)/2,  p_neg = p_- = (1 - sB)/2     (g_hat = (K_+ - K_-)/(2M))
    projector  p_pos = F_- = A(1+s)/2,    p_neg = F_+ = A(1-s)/2       (g_hat = (K_- - K_+)/(2M))

Then, exactly,

    P_zero    = P(K_pos = K_neg)
    P_correct = P(sign(g_hat) = sign(g_true))
    P_wrong   = P(sign(g_hat) = -sign(g_true))       and  P_correct + P_wrong + P_zero = 1.

When g_true = 0 the direction is undefined; by convention P_correct = P_wrong = (1-P_zero)/2
and the row is flagged.

"DIRECTIONAL AMBIGUITY REGION" is an operational term introduced for this project (like the
Stage 5 "dead zone"); it is not established literature terminology. A point is an
(M, q)-directional-ambiguity point for component k when ``P_correct(theta, M, k) <= q``.
"""

from __future__ import annotations

import numpy as np
from scipy.stats import norm

from qlo.stage6.exact_distribution import difference_probabilities

AMBIGUITY_LEVELS = (0.55, 0.60, 0.75)


def direction_probabilities(p_pos, p_neg, M) -> dict[str, np.ndarray]:
    """Exact ``P_correct``, ``P_wrong``, ``P_zero`` (and a ``zero_signal`` flag)."""
    p_pos, p_neg, M = np.broadcast_arrays(*(np.asarray(x, dtype=np.float64) for x in (p_pos, p_neg, M)))
    gt, eq, lt = difference_probabilities(p_pos, p_neg, M)
    pos_signal = p_pos > p_neg
    neg_signal = p_pos < p_neg
    zero_signal = ~(pos_signal | neg_signal)
    correct = np.where(pos_signal, gt, np.where(neg_signal, lt, (1.0 - eq) / 2.0))
    wrong = np.where(pos_signal, lt, np.where(neg_signal, gt, (1.0 - eq) / 2.0))
    return {"p_correct": correct, "p_wrong": wrong, "p_zero": eq, "zero_signal": zero_signal}


def difference_moments(p_pos, p_neg, M) -> tuple[np.ndarray, np.ndarray]:
    """``E[D] = M(p_pos - p_neg)`` and ``Var(D) = M[p_pos(1-p_pos) + p_neg(1-p_neg)]`` for D = K_pos - K_neg."""
    p_pos, p_neg, M = (np.asarray(x, dtype=np.float64) for x in (p_pos, p_neg, M))
    mean = M * (p_pos - p_neg)
    var = M * (p_pos * (1 - p_pos) + p_neg * (1 - p_neg))
    return mean, var


def normal_approximation(p_pos, p_neg, M, continuity: bool = True) -> dict[str, np.ndarray]:
    """Large-M normal approximation. ``P_correct ~= Phi(|E D|/sd)``; with ``continuity`` the
    half-integer correction ``Phi((|E D| - 1/2)/sd)`` is used, which accounts for the tie mass.
    ``|E D|/sd`` equals the exact component SNR of the estimator (the 1/(2M) prefactor cancels)."""
    mean, var = difference_moments(p_pos, p_neg, M)
    sd = np.sqrt(np.maximum(var, 0.0))
    with np.errstate(divide="ignore", invalid="ignore"):
        snr = np.where(sd > 0, np.abs(mean) / sd, np.where(np.abs(mean) > 0, np.inf, 0.0))
        shift = 0.5 if continuity else 0.0
        z = np.where(sd > 0, (np.abs(mean) - shift) / sd, snr)
    return {"snr": snr, "p_correct_normal": norm.cdf(z), "p_zero_normal": np.where(sd > 0, norm.pdf(0.0) / sd, 1.0)}


def shots_for_direction_normal(sB, q, continuity: bool = False) -> np.ndarray:
    """Parity only: invert ``Phi(SNR) = q`` with ``SNR^2 = 2M(sB)^2/(1-(sB)^2)``:
    ``M = z_q^2 (1-(sB)^2)/(2 (sB)^2)``. Ignores the tie correction, so it slightly
    UNDER-estimates the exact requirement; ``shots_for_direction_exact`` quantifies this."""
    sB = np.clip(np.asarray(sB, dtype=np.float64), -1.0, 1.0)
    z = norm.ppf(np.asarray(q, dtype=np.float64))
    num = z**2 * (1.0 - sB**2)
    den = 2.0 * sB**2
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(den > 0, num / den, np.where(num > 0, np.inf, 0.0))


def log10_shots_for_direction_normal(log_abs_sB, q) -> np.ndarray:
    l = np.asarray(log_abs_sB, dtype=np.float64)
    sB2 = np.exp(2 * l)
    z = norm.ppf(q)
    with np.errstate(divide="ignore", invalid="ignore"):
        return 2 * np.log10(abs(z)) + np.log10(np.maximum(1.0 - sB2, 0.0)) - np.log10(2.0) - 2 * l / np.log(10.0)


def skellam_direction_probabilities(p_pos, p_neg, M) -> dict[str, np.ndarray]:
    """Poisson-limit (Skellam) approximation, for the deep-plateau regime where BOTH shifted
    probabilities are tiny (the projector case): K_pos ~ Poisson(M p_pos), K_neg ~ Poisson(M p_neg),
    so D = K_pos - K_neg is Skellam(M p_pos, M p_neg). Exact as p -> 0 at fixed M p; error is O(p).
    Cost is independent of M, which is what makes large-n inversion feasible."""
    import warnings

    from scipy.stats import skellam

    p_pos, p_neg, M = np.broadcast_arrays(*(np.asarray(x, dtype=np.float64) for x in (p_pos, p_neg, M)))
    mu1, mu2 = M * p_pos, M * p_neg
    eq = np.empty(mu1.shape); gt = np.empty(mu1.shape); lt = np.empty(mu1.shape)
    # scipy's skellam returns NaN once the means exceed ~1e18; there Skellam -> Normal(mu1-mu2, mu1+mu2)
    # to machine precision, so use the (continuity-corrected) normal branch for large means.
    big = (mu1 + mu2) > 1e6
    small = ~big
    if np.any(small):
        with warnings.catch_warnings():  # scipy's ncx2 series can fail to converge for very lopsided means
            warnings.simplefilter("ignore")
            e = np.asarray(skellam.pmf(0, mu1[small], mu2[small]), dtype=np.float64)
            c = np.asarray(skellam.cdf(0, mu1[small], mu2[small]), dtype=np.float64)
        bad = ~(np.isfinite(e) & np.isfinite(c))
        if np.any(bad):  # fall back to the normal branch for those entries
            m_, s_ = mu1[small][bad] - mu2[small][bad], np.sqrt(mu1[small][bad] + mu2[small][bad])
            e[bad] = norm.pdf(m_ / s_) / s_
            c[bad] = norm.cdf((0.5 - m_) / s_)
        eq[small] = e
        gt[small] = np.clip(1.0 - c, 0.0, 1.0)
        lt[small] = np.clip(c - e, 0.0, 1.0)
    if np.any(big):
        mean = mu1[big] - mu2[big]
        sd = np.sqrt(mu1[big] + mu2[big])
        eq[big] = norm.pdf(mean / sd) / sd
        gt[big] = norm.sf((0.5 - mean) / sd)
        lt[big] = norm.cdf((-0.5 - mean) / sd)
    gt = np.clip(gt, 0.0, 1.0)
    lt = np.clip(lt, 0.0, 1.0)
    eq = np.clip(eq, 0.0, 1.0)
    pos_signal, neg_signal = p_pos > p_neg, p_pos < p_neg
    correct = np.where(pos_signal, gt, np.where(neg_signal, lt, (1.0 - eq) / 2.0))
    wrong = np.where(pos_signal, lt, np.where(neg_signal, gt, (1.0 - eq) / 2.0))
    return {"p_correct": correct, "p_wrong": wrong, "p_zero": np.clip(eq, 0.0, 1.0)}


def _bisect_log_shots(prob_fn, q, lo_log10, hi_log10, n_iter: int = 34) -> tuple[np.ndarray, np.ndarray]:
    """Smallest integer M with ``prob_fn(M) >= q``, bisecting in log10 M over ``[10^lo, 10^hi]``
    (brackets may be per-row arrays). ``prob_fn`` must be non-decreasing in M (checked
    numerically in the experiment). ``reached`` is False where even ``10^hi`` falls short;
    rows already satisfied at ``10^lo`` return ``10^lo`` (callers choose ``lo`` below the answer)."""
    q = np.asarray(q, dtype=np.float64)
    shape = np.broadcast(q, lo_log10, hi_log10).shape
    lo = np.broadcast_to(np.asarray(lo_log10, dtype=np.float64), shape).copy()
    hi = np.broadcast_to(np.asarray(hi_log10, dtype=np.float64), shape).copy()
    # evaluate ONLY at integer budgets ceil(10^x): shot counts are integers, the float-M binomial formula is
    # just an interpolation between them, and on integers P_correct is non-decreasing (checked), so the
    # predicate is monotone in x and the returned ceil(10^hi) is a budget that was actually tested
    reached = prob_fn(np.ceil(np.power(10.0, hi))) >= q
    for _ in range(n_iter):
        mid = 0.5 * (lo + hi)
        good = prob_fn(np.ceil(np.power(10.0, mid))) >= q
        hi = np.where(good, mid, hi)
        lo = np.where(good, lo, mid)
    return np.ceil(np.power(10.0, hi)), reached


def shots_for_direction_skellam(p_pos, p_neg, q, lo_log10: float = 0.0, hi_log10: float = 60.0) -> tuple[np.ndarray, np.ndarray]:
    """Required shots under the Skellam (Poisson-limit) approximation; cost independent of M."""
    p_pos, p_neg, q = np.broadcast_arrays(*(np.asarray(x, dtype=np.float64) for x in (p_pos, p_neg, q)))
    fn = lambda M: skellam_direction_probabilities(p_pos, p_neg, M)["p_correct"]
    return _bisect_log_shots(fn, q, lo_log10, hi_log10)


def shots_for_direction_exact(p_pos, p_neg, q, lo_log10=0.0, hi_log10=7.0,
                              n_iter: int = 28) -> tuple[np.ndarray, np.ndarray]:
    """Smallest M with exact ``P_correct(M) >= q``, bisecting in log10 M over ``[10^lo, 10^hi]``.
    The exact evaluation cost grows like sqrt(M p) per row, so brackets must keep ``M p`` moderate
    (per-row arrays are accepted); ``reached`` flags rows whose requirement exceeds ``10^hi``."""
    p_pos, p_neg, q = np.broadcast_arrays(*(np.asarray(x, dtype=np.float64) for x in (p_pos, p_neg, q)))
    fn = lambda M: direction_probabilities(p_pos, p_neg, M)["p_correct"]
    return _bisect_log_shots(fn, q, lo_log10, hi_log10, n_iter)
