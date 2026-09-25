"""Required-shot distributions over random initializations, for both benchmarks."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import norm

from qlo.stage5.theory import a_k, f_plus_minus, log_a_k
from qlo.stage5.theory import s_k as proj_s_k
from qlo.stage5.theory import shots_for_snr as proj_shots_for_snr
from qlo.stage6.directional import (
    direction_probabilities,
    log10_shots_for_direction_normal,
    shots_for_direction_exact,
    shots_for_direction_normal,
    shots_for_direction_skellam,
)
from qlo.stage6.finite_shot import log10_shots_for_snr, shots_for_snr
from qlo.stage6.parity_benchmark import b_k, log_abs_b_k, s_k, shifted_probabilities

LN10 = np.log(10.0)


def sample_parity_signal(n_qubits: int, n_samples: int, seed, k: int = 0, chunk: int = 20_000) -> tuple[np.ndarray, np.ndarray]:
    """Returns ``(log|s_k B_k|, sign(s_k B_k))`` over theta ~ U[-pi,pi]^n, streamed in chunks."""
    log_abs = np.empty(n_samples)
    sign = np.empty(n_samples)
    rng = np.random.default_rng(np.random.SeedSequence([6, 100, int(n_qubits), int(seed)]))
    for i in range(0, n_samples, chunk):
        m = min(chunk, n_samples - i)
        th = rng.uniform(-np.pi, np.pi, size=(m, n_qubits))
        s, B = s_k(th, k), b_k(th, k)
        with np.errstate(divide="ignore"):
            log_abs[i:i + m] = np.log(np.abs(s)) + log_abs_b_k(th, k)
        sign[i:i + m] = np.sign(s * B)
    return log_abs, sign


def sample_projector_signal(n_qubits: int, n_samples: int, seed, k: int = 0, chunk: int = 20_000) -> tuple[np.ndarray, np.ndarray]:
    """Returns ``(log A_k, s_k)`` for the frozen Stage 5 projector benchmark."""
    logA, s = np.empty(n_samples), np.empty(n_samples)
    rng = np.random.default_rng(np.random.SeedSequence([6, 101, int(n_qubits), int(seed)]))
    for i in range(0, n_samples, chunk):
        m = min(chunk, n_samples - i)
        th = rng.uniform(-np.pi, np.pi, size=(m, n_qubits))
        logA[i:i + m] = log_a_k(th, k)
        s[i:i + m] = proj_s_k(th, k)
    return logA, s


def _q(x, p):
    return float(np.percentile(x, p))


def _stats(l10: np.ndarray, tag: str) -> dict:
    finite = np.isfinite(l10)
    v = l10[finite]
    return {f"{tag}_finite_frac": float(finite.mean()), f"{tag}_log10M_median": _q(v, 50), f"{tag}_log10M_q25": _q(v, 25),
            f"{tag}_log10M_q75": _q(v, 75), f"{tag}_log10M_q90": _q(v, 90), f"{tag}_log10M_geomean": float(v.mean())}


def parity_required_shots(n_qubits: int, n_samples: int, seed, rho_targets=(1.0, 2.0), q_targets=(0.75, 0.90),
                          table_log10_sB_min: float = -2.5) -> dict:
    """Required shots for SNR and directional-correctness targets, parity benchmark.

    SNR targets are closed form (exact). Directional targets are EXACT (difference-of-binomials
    bisection, via ``parity_direction_table``) for |sB| >= 10^table_log10_sB_min, i.e. up to M ~ 2e4;
    for smaller |sB| the normal inversion ``M = z_q^2 (1-(sB)^2)/(2(sB)^2)`` is used, whose error
    at the table edge is measured and reported (it shrinks further as M grows).
    """
    log_abs, _ = sample_parity_signal(n_qubits, n_samples, seed)
    sB = np.exp(log_abs)
    row = {"n": n_qubits, "benchmark": "parity", "n_samples": int(n_samples),
           "empirical_mean_log_abs_sB": float(log_abs.mean()), "theory_mean_log_abs_sB": -n_qubits * np.log(2.0),
           "empirical_var_log_abs_sB": float(log_abs.var(ddof=1)), "theory_var_log_abs_sB": n_qubits * np.pi**2 / 12.0,
           "empirical_mean_abs_sB": float(sB.mean()), "theory_mean_abs_sB": (2.0 / np.pi) ** n_qubits,
           "empirical_var_g": float(np.mean(0.25 * sB**2)),  # E[g^2]; E[g] = 0 exactly by symmetry
           "theory_var_g": 2.0 ** (-(n_qubits + 2)),
           "log10_abs_g_median": _q((log_abs - np.log(2.0)) / LN10, 50)}
    for rho in rho_targets:
        row.update({**_stats(log10_shots_for_snr(log_abs, rho), f"snr{rho:g}"), f"snr{rho:g}_method": "closed_form_exact"})
    rng = np.random.default_rng(np.random.SeedSequence([6, 102, n_qubits, seed]))
    for q in q_targets:
        l_norm = np.maximum(log10_shots_for_direction_normal(log_abs, q), 0.0)  # physical floor: M >= 1
        grid, grid_l10 = parity_direction_table(q, table_log10_sB_min)
        l_sB = log_abs / LN10
        in_table = l_sB >= grid[0]
        # M_req is non-increasing in |sB|: the grid point at or just below |sB| gives an UPPER bound
        # (conservative), the one just above a lower bound; the gap is reported
        j = np.clip(np.searchsorted(grid, l_sB[in_table], side="right") - 1, 0, grid.size - 1)
        l10 = l_norm.copy()
        l10[in_table] = grid_l10[j]
        upper_nb = grid_l10[np.minimum(j + 1, grid.size - 1)]
        row.update({**_stats(l10, f"dir{q}"),
                    f"dir{q}_method": f"exact_table_for_log10|sB|>={table_log10_sB_min:g}_else_normal_inversion",
                    f"dir{q}_frac_exact": float(in_table.mean()),
                    f"dir{q}_max_table_bracket_log10": float(np.max(grid_l10[j] - upper_nb)) if in_table.any() else float("nan")})
        # approximation audit: exact vs normal on the table itself, overall and in the regime where normal is used
        norm_on_grid = np.maximum(log10_shots_for_direction_normal(grid * LN10, q), 0.0)
        err = grid_l10 - norm_on_grid
        edge = grid <= grid[0] + 0.25  # the table's small-|sB| end: M ~ 1e4, adjacent to the normal regime
        row.update({f"dir{q}_median_log10_ratio_exact_over_normal": float(np.median(err)),
                    f"dir{q}_max_log10_ratio_exact_over_normal": float(np.max(np.abs(err))),
                    f"dir{q}_max_log10_ratio_exact_over_normal_at_table_edge": float(np.max(np.abs(err[edge])))})
    return row


_TABLES: dict = {}


def parity_direction_table(q: float, log10_sB_min: float = -2.5, n_grid: int = 3000) -> tuple[np.ndarray, np.ndarray]:
    """Exact smallest M with P_correct >= q for the parity estimator on a log grid of |sB| in [10^min, 1).

    For parity (p_+- = (1 +- sB)/2) the requirement depends on theta only through |sB|, so one table
    serves every n and sample (cached per process). Solved by exact bisection with a wide bracket.
    """
    key = (float(q), float(log10_sB_min), int(n_grid))
    if key not in _TABLES:
        grid = np.linspace(log10_sB_min, -1e-3, n_grid)
        pp, pm = shifted_probabilities(np.power(10.0, grid), 1.0)
        l_norm = np.maximum(log10_shots_for_direction_normal(grid * LN10, q), 0.0)
        m_ex, ok = shots_for_direction_exact(pp, pm, q, lo_log10=0.0, hi_log10=np.maximum(l_norm + 1.5, 2.0), n_iter=30)
        if not np.all(ok):
            raise RuntimeError("parity direction table: bracket too small")
        _TABLES[key] = (grid, np.log10(m_ex))
    return _TABLES[key]


def projector_required_shots(n_qubits: int, n_samples: int, seed, rho_targets=(1.0, 2.0), q_targets=(0.75, 0.90),
                             exact_subsample: int = 400, exact_max_log10M: float = 6.0) -> dict:
    """Same quantities for the frozen projector benchmark (Stage 5 functions reused unmodified)."""
    logA, s = sample_projector_signal(n_qubits, n_samples, seed)
    A = np.exp(logA)
    row = {"n": n_qubits, "benchmark": "projector", "n_samples": int(n_samples),
           "empirical_mean_logA": float(logA.mean()), "theory_mean_logA": -2 * (n_qubits - 1) * np.log(2.0),
           "log10_abs_g_median": _q(np.log10(np.abs(A * s / 2)), 50)}
    for rho in rho_targets:
        with np.errstate(divide="ignore"):
            l10 = np.log10(proj_shots_for_snr(A, s, rho))
        row.update({**_stats(l10, f"snr{rho:g}"), f"snr{rho:g}_method": "closed_form_exact"})
    # Both shifted probabilities are O(A), so at the requirement M F_+- ~ z_q^2/s^2 = O(1..1e3) even
    # when M itself is ~4^n: the EXACT difference-of-binomials evaluation is cheap if the bisection is
    # bracketed per sample around the normal estimate M* = z_q^2/(A s^2) (never over a global [1, 1e30],
    # whose midpoints put M F ~ 1e14 and summation windows ~1e8 wide). Samples with |s| so small that
    # M* A > HUGE_MA use the Skellam/normal limit (exact as p -> 0; normal to O(1/sqrt(MA)) there);
    # the fraction is recorded. The Skellam limit is also computed on a subsample as a cross-check.
    HUGE_MA = 1e4
    fp_all, fm_all = f_plus_minus(A, s)
    rng = np.random.default_rng(np.random.SeedSequence([6, 103, n_qubits, seed]))
    sub = rng.choice(n_samples, size=min(exact_subsample, n_samples), replace=False)
    for q in q_targets:
        z2 = norm.ppf(q) ** 2
        with np.errstate(divide="ignore"):
            est = np.maximum(z2 / (A * s**2), 1.0 / A)  # >= 1/A: below it the tie mass alone blocks q
        huge = est * A > HUGE_MA
        m_dir = np.full(n_samples, np.nan)
        reached = np.zeros(n_samples, dtype=bool)
        todo = np.flatnonzero(~huge)
        # tight bracket [est/30, 5 est] first; rows whose answer falls outside are re-solved on a wider one
        for down, up in ((1.5, 0.7), (4.0, 2.0)):
            if todo.size == 0:
                break
            hi = np.log10(est[todo]) + up
            lo = np.maximum(np.log10(est[todo]) - down, 0.0)
            m_ex, r_ex = shots_for_direction_exact(fm_all[todo], fp_all[todo], q, lo_log10=lo, hi_log10=hi, n_iter=24)
            at_lo = direction_probabilities(fm_all[todo], fp_all[todo], np.ceil(np.power(10.0, lo)))["p_correct"] >= q
            good = r_ex & (~at_lo | (lo == 0.0))  # answer inside the bracket (or M = 1 already suffices)
            m_dir[todo[good]], reached[todo[good]] = m_ex[good], True
            todo = todo[~good]
        if huge.any():
            l_est = np.log10(est[huge])
            m_sk, r_sk = shots_for_direction_skellam(fm_all[huge], fp_all[huge], q, lo_log10=np.maximum(l_est - 4.0, 0.0),
                                                     hi_log10=l_est + 2.0)
            m_dir[huge], reached[huge] = m_sk, r_sk
        with np.errstate(divide="ignore"):
            l10 = np.where(reached, np.log10(m_dir), np.nan)
        finite = np.isfinite(l10)
        row.update({**_stats(l10[finite], f"dir{q}"), f"dir{q}_method": "exact_bisection_per_sample_bracket",
                    f"dir{q}_exact_reached_frac": float(reached.mean()),
                    f"dir{q}_frac_skellam_normal_limit": float(huge.mean())})
        m_sk, sk_ok = shots_for_direction_skellam(fm_all[sub], fp_all[sub], q)
        both = np.isfinite(l10[sub]) & sk_ok & ~huge[sub]  # compare only against rows solved exactly
        if both.any():
            ratio = np.log10(m_sk[both]) - l10[sub][both]
            row.update({f"dir{q}_skellam_checked": int(both.sum()),
                        f"dir{q}_median_log10_ratio_skellam_over_exact": float(np.median(ratio)),
                        f"dir{q}_max_log10_ratio_skellam_over_exact": float(np.max(np.abs(ratio)))})
    return row
