"""Global-benchmark required-shot distributions, dead-zone fractions and the Stage 4 retrodiction."""

from __future__ import annotations

import numpy as np
import pandas as pd

from qlo.stage5.sampling import sample_log_a_s, theory_log_a_stats
from qlo.stage5.theory import (
    f_plus_minus,
    log10_shots_for_nonzero_rare_event,
    log10_shots_for_snr,
    p_zero_exact,
    shots_for_nonzero_exact,
    shots_for_nonzero_poisson_limit,
)

LN10 = np.log(10.0)


def _q(x, p):
    return float(np.percentile(x, p))


def required_shots_distribution(n_qubits: int, n_samples: int, seed, q_targets=(0.5, 0.9), rho_targets=(1.0, 2.0), k: int = 0) -> tuple[dict, pd.DataFrame]:
    """Per-theta required shots (log10) for nonzero-probability and SNR targets; returns (summary row, per-sample frame)."""
    logA, s = sample_log_a_s(n_qubits, n_samples, seed, k)
    A = np.exp(logA)
    fp, fm = f_plus_minus(A, s)
    per = pd.DataFrame({"logA": logA, "log10A": logA / LN10, "s": s, "log10_abs_g": np.log10(np.abs(A * s / 2)),
                        "log10_var_prefactor": np.log10(A - A**2 * (1 + s**2) / 2)})
    row = {"n": n_qubits, "k": k, "n_samples": int(n_samples),
           "empirical_mean_logA": float(logA.mean()), "theory_mean_logA": theory_log_a_stats(n_qubits)["mean_logA"],
           "empirical_var_logA": float(logA.var(ddof=1)), "theory_var_logA": theory_log_a_stats(n_qubits)["var_logA"],
           "empirical_mean_A": float(A.mean()), "theory_mean_A": theory_log_a_stats(n_qubits)["mean_A"],
           "geometric_A_theory": theory_log_a_stats(n_qubits)["geometric_A"],
           "log10A_median": _q(logA / LN10, 50), "log10A_q25": _q(logA / LN10, 25), "log10A_q75": _q(logA / LN10, 75)}
    for q in q_targets:
        # exact integer solve for every sample; the rare-event and Poisson-limit approximations recorded alongside
        M_exact, ok = shots_for_nonzero_exact(fp, fm, q)
        l10 = np.log10(M_exact)
        l10_rare = log10_shots_for_nonzero_rare_event(logA, q)
        l10_pois = np.log10(shots_for_nonzero_poisson_limit(A, s, q))
        tag = f"nonzero{q}"
        per[f"log10M_{tag}_exact"] = l10
        per[f"log10M_{tag}_rare_event"] = l10_rare
        per[f"log10M_{tag}_poisson"] = l10_pois
        row.update({f"{tag}_method": "exact_bisection_all_samples", f"{tag}_bracket_ok_frac": float(ok.mean()),
                    f"{tag}_log10M_median": _q(l10, 50), f"{tag}_log10M_q25": _q(l10, 25), f"{tag}_log10M_q75": _q(l10, 75),
                    f"{tag}_log10M_q90": _q(l10, 90), f"{tag}_log10M_geomean": float(l10.mean()),
                    f"{tag}_rare_event_log10M_median": _q(l10_rare, 50), f"{tag}_poisson_log10M_median": _q(l10_pois, 50),
                    f"{tag}_median_abs_log10_err_rare_event": float(np.median(np.abs(l10_rare - l10))),
                    f"{tag}_median_abs_log10_err_poisson": float(np.median(np.abs(l10_pois - l10))),
                    f"{tag}_max_abs_log10_err_poisson": float(np.max(np.abs(l10_pois - l10)))})
    for rho in rho_targets:
        l10 = log10_shots_for_snr(logA, s, rho)  # continuous; inf at s = 0 (measure zero)
        finite = np.isfinite(l10)
        tag = f"snr{rho:g}"
        per[f"log10M_{tag}"] = l10
        row.update({f"{tag}_method": "closed_form", f"{tag}_finite_frac": float(finite.mean()),
                    f"{tag}_log10M_median": _q(l10[finite], 50), f"{tag}_log10M_q25": _q(l10[finite], 25),
                    f"{tag}_log10M_q75": _q(l10[finite], 75), f"{tag}_log10M_q90": _q(l10[finite], 90),
                    f"{tag}_log10M_geomean": float(l10[finite].mean())})
    return row, per


def deadzone_fraction(n_qubits: int, n_samples: int, seed, shots, q_levels=(0.5, 0.9, 0.99), k: int = 0,
                      full_gradient_samples: int | None = None) -> list[dict]:
    """Initialization dead-zone fraction  P_theta[ P_zero(theta, M, k) >= q ]  and, optionally, the full-gradient
    version  P_theta[ prod_k P_zero_k >= q ]  (product valid because Stage 4 used independent batches per k)."""
    logA, s = sample_log_a_s(n_qubits, n_samples, seed, k)
    fp, fm = f_plus_minus(np.exp(logA), s)
    rows = []
    full = None
    if full_gradient_samples:
        rng = np.random.default_rng(np.random.SeedSequence([5, 101, int(n_qubits), int(seed)]))
        full = rng.uniform(-np.pi, np.pi, size=(int(full_gradient_samples), n_qubits))
    for M in shots:
        pz = p_zero_exact(fp, fm, float(M))
        row = {"n": n_qubits, "shots": int(M), "k": k, "n_samples": int(n_samples),
               "p_zero_median": _q(pz, 50), "p_zero_q25": _q(pz, 25), "p_zero_q75": _q(pz, 75), "p_zero_mean": float(pz.mean())}
        for q in q_levels:
            row[f"deadzone_frac_q{q}"] = float(np.mean(pz >= q))
        if full is not None:
            pfull = np.ones(full.shape[0])
            for kk in range(n_qubits):
                from qlo.stage5.theory import a_k, s_k

                fpk, fmk = f_plus_minus(a_k(full, kk), s_k(full, kk))
                pfull *= p_zero_exact(fpk, fmk, float(M))
            row.update({"n_samples_full": int(full.shape[0]), "p_full_zero_median": _q(pfull, 50), "p_full_zero_mean": float(pfull.mean())})
            for q in q_levels:
                row[f"full_deadzone_frac_q{q}"] = float(np.mean(pfull >= q))
        rows.append(row)
    return rows


def stage4_prediction(theta0: np.ndarray, shots, n_iters: int = 2000) -> pd.DataFrame:
    """For the actual Stage 4 starting points theta0 (rows), per M: component P_zero (k = 0..n-1), the full-vector
    zero probability prod_k P_zero_k (independent batches per component => product), and the probability the
    trajectory never moves in n_iters iterations = P_full^T (exact: with g_hat = 0 theta is unchanged, so the
    iterations are iid). Returns one row per (theta0 index, M)."""
    from qlo.stage5.theory import a_k, s_k

    n = theta0.shape[1]
    rows = []
    for M in shots:
        pz = np.stack([p_zero_exact(*f_plus_minus(a_k(theta0, kk), s_k(theta0, kk)), float(M)) for kk in range(n)], axis=1)
        pfull = np.prod(pz, axis=1)
        for i in range(theta0.shape[0]):
            rows.append({"theta_index": i, "shots": int(M), "n": n, "F0": float(np.prod(np.cos(theta0[i] / 2) ** 2)),
                         "p_zero_component_min": float(pz[i].min()), "p_zero_component_median": float(np.median(pz[i])),
                         "p_zero_component_max": float(pz[i].max()), "p_full_zero": float(pfull[i]),
                         "p_stuck_all_iters": float(pfull[i] ** n_iters), "n_iters": int(n_iters)})
    return pd.DataFrame(rows)
