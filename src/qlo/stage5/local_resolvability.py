"""Matched term-wise LOCAL estimator control.

Local cost C_L = 1 - (1/n) sum_j cos^2(theta_j/2). Only the j = k term depends on theta_k, so a
measurement-efficient estimator measures just that qubit's |0> projector at the two shifts:
p_± = (1 ∓ s)/2, K_± ~ Bin(M, p_±), g_hat_{L,k} = (K_- - K_+)/(2 M n). This EXPLOITS the known
decomposition of this particular cost; it is not a generic local-cost measurement strategy.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from qlo.stage5.theory import local_p_plus_minus, local_shots_for_snr, p_zero_exact, shots_for_nonzero_exact

LN10 = np.log(10.0)


def _q(x, p):
    return float(np.percentile(x, p))


def local_required_shots_distribution(n_qubits: int, n_samples: int, seed, q_targets=(0.5, 0.9), rho_targets=(1.0, 2.0)) -> dict:
    rng = np.random.default_rng(np.random.SeedSequence([5, 102, int(n_qubits), int(seed)]))
    s = np.sin(rng.uniform(-np.pi, np.pi, size=int(n_samples)))
    pp, pm = local_p_plus_minus(s)
    row = {"n": n_qubits, "n_samples": int(n_samples), "estimator": "term_wise_local_projector"}
    for q in q_targets:
        M, ok = shots_for_nonzero_exact(pp, pm, q)
        l10 = np.log10(M)
        tag = f"nonzero{q}"
        row.update({f"{tag}_log10M_median": _q(l10, 50), f"{tag}_log10M_q25": _q(l10, 25), f"{tag}_log10M_q75": _q(l10, 75),
                    f"{tag}_log10M_q90": _q(l10, 90), f"{tag}_bracket_ok_frac": float(ok.mean())})
    for rho in rho_targets:
        l10 = np.log10(local_shots_for_snr(s, rho))
        f = np.isfinite(l10)
        tag = f"snr{rho:g}"
        row.update({f"{tag}_log10M_median": _q(l10[f], 50), f"{tag}_log10M_q25": _q(l10[f], 25), f"{tag}_log10M_q75": _q(l10[f], 75),
                    f"{tag}_log10M_q90": _q(l10[f], 90)})
    for M in (16, 64, 256, 1024):
        row[f"p_zero_median_M{M}"] = _q(p_zero_exact(pp, pm, float(M)), 50)
    return row
