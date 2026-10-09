"""Family 0: the frozen Stage 7 RX-product benchmark run through the B1 per-theta pipeline (regression arm).

Same theta distribution (``qlo.stage7.analysis.sample_theta``: iid U[-pi, pi], stream 20 = the Stage 7 required-shots
stream), same objective F = prod_j cos^2(theta_j / 2), target |0^n>, shifted component k = 0 and the same shift
convention (``qlo.stage5.theory``: F_+ = A(1 - s)/2, F_- = A(1 + s)/2). Only the master seeds are new.

B3 published seed intervals for the Loschmidt and SWAP SNR-1 slopes and their gap, but not for b_S or b_r. Those two
are computed here from the B3 theta populations themselves (B3 seeds, B3 sample size, same stream), after checking that
the regenerated populations reproduce B3's published per-seed LE / SWAP slopes.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from qlo.b1 import quantities as Q
from qlo.hardening import seed_replicates as SR
from qlo.stage5.theory import f_plus_minus, log_a_k, s_k
from qlo.stage7 import analysis as AN
from qlo.stage7 import required_shots as R

LN10 = float(np.log(10.0))


def rx_per_theta(n: int, n_samples: int, seed: int, stream: int = 20, k: int = 0, rho: float = 1.0) -> dict:
    th = AN.sample_theta(n, n_samples, seed, stream)
    logA, s = log_a_k(th, k), s_k(th, k)
    fp, fm = f_plus_minus(np.exp(logA), s)
    q = Q.per_theta(fp, fm, rho)
    q["logA"], q["s"] = logA, s
    q["log10_M_LE_ceil"] = Q.ceil_convention(q["log10_M_LE"])
    q["log10_M_SWAP_ceil"] = Q.ceil_convention(q["log10_M_SWAP"])
    return q


def rx_cell_row(n: int, n_samples: int, seed: int, stream: int, k: int, rho: float, zero_tol_r: float) -> tuple[dict, dict]:
    """(sample-summary row with medians for the fits, identity-check row) for one (seed, n)."""
    q = rx_per_theta(n, n_samples, seed, stream, k, rho)
    row = {"family": "rx_product", "regime": "product", "depth": 1, "position": "k=0", "layer": 0, "qubit": k, "gate": "RX", "n": n, "seed": seed,
           **Q.cell_diagnostics(q, zero_tol_r),
           # B3 convention: finite median of the Stage 7 integer-budget values
           "median_log10_M_LE_ceil": SR.finite_median(q["log10_M_LE_ceil"]),
           "median_log10_M_SWAP_ceil": SR.finite_median(q["log10_M_SWAP_ceil"]),
           "n_nonfinite_ceil": int(np.sum(~np.isfinite(q["log10_M_LE_ceil"])) + np.sum(~np.isfinite(q["log10_M_SWAP_ceil"]))),
           "max_abs_r_minus_s": float(np.nanmax(np.abs(q["r"] - q["s"]))),
           "max_abs_log10S_minus_log10A": float(np.max(np.abs(q["log10_S"] - q["logA"] / LN10)))}
    ident = {"family": "rx_product", "regime": "product", "position": "k=0", "n": n, "seed": seed, **Q.identity_check(q)}
    return row, ident


def b3_regenerated_medians(seeds, n_values, n_samples: int, stream: int = SR.REQ_STREAM) -> pd.DataFrame:
    """Per (B3 seed, n): the Stage 7 SNR-1 medians exactly as B3 computed them (``R.snr_shots``), plus the medians of
    log10 S = log10 A and log10 |r| = log10 |s| on the same theta."""
    rows = []
    for seed in seeds:
        for n in n_values:
            th = AN.sample_theta(n, n_samples, seed, stream)
            logA, s = log_a_k(th, 0), s_k(th, 0)
            with np.errstate(divide="ignore"):
                rows.append({"seed": int(seed), "n": int(n), "n_samples": int(n_samples),
                             "median_log10M_loschmidt_snr1": SR.finite_median(R.snr_shots("loschmidt", logA, s, 1.0)),
                             "median_log10M_swap_snr1": SR.finite_median(R.snr_shots("swap", logA, s, 1.0)),
                             "median_log10_S": float(np.median(logA / LN10)),
                             "median_log10_abs_r": float(np.median(np.log10(np.abs(s))))})
    return pd.DataFrame(rows)


def b3_regenerated_slopes(med: pd.DataFrame, ranges: dict) -> pd.DataFrame:
    rows = []
    for seed, g in med.groupby("seed"):
        for rname, n_range in ranges.items():
            d = g.set_index("n").loc[list(n_range)]
            le = SR.fit_series(d.index.values, d.median_log10M_loschmidt_snr1.values)
            sw = SR.fit_series(d.index.values, d.median_log10M_swap_snr1.values)
            fs = SR.fit_series(d.index.values, d.median_log10_S.values)
            fr = SR.fit_series(d.index.values, d.median_log10_abs_r.values)
            rows.append({"seed": int(seed), "fit_range": rname, "slope_loschmidt_snr1": le["slope"], "slope_swap_snr1": sw["slope"],
                         "gap_swap_minus_le_snr1": sw["slope"] - le["slope"], "b_S": -fs["slope"], "b_r": -fr["slope"]})
    return pd.DataFrame(rows)
