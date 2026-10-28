"""A2 check of the general conditional sign law for the Loschmidt estimator.

Claim: for any circuit, as M (F+ + F-) -> 0,
  P(sign ghat = sign g | ghat != 0) -> (1 + |r|) / 2,   r = (F- - F+) / (F+ + F-).
On the RX product circuit r = sin(theta_k), which is the Stage 7 law.
Exact binomial sums, compared point by point.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from general_scaling import conditional_sign_exact, shifted_fidelities

OUT = Path(__file__).resolve().parent
rng = np.random.default_rng(20261004)
rows = []
for family, n, depth in [("rx_product", 12, None), ("rx_product", 20, None),
                         ("hea", 10, 2), ("hea", 12, 2), ("hea", 10, 10), ("hea", 12, 12)]:
    Fp, Fm = shifted_fidelities(family, n, 400, rng, depth)
    S = Fp + Fm
    rel = np.abs(Fm - Fp) / S
    law = (1 + rel) / 2
    for M in (16, 1024):
        cs = conditional_sign_exact(Fp, Fm, M)
        ok = np.isfinite(cs)
        x = M * S
        bands = {}
        for lo, hi in [(0, 1e-3), (1e-3, 1e-2), (1e-2, 1e-1), (1e-1, 1)]:
            m = ok & (x >= lo) & (x < hi)
            bands[f"MS in [{lo:g},{hi:g})"] = (
                dict(count=int(m.sum()), max_abs_err=float(np.max(np.abs(cs[m] - law[m]))))
                if m.any() else dict(count=0, max_abs_err=None))
        row = dict(family=family, n=n, depth=depth, M=M,
                   median_exact=float(np.median(cs[ok])), median_law=float(np.median(law[ok])),
                   median_M_times_S=float(np.median(x)), bands=bands)
        rows.append(row)
        print(json.dumps(row), flush=True)
(OUT / "condsign_results.json").write_text(json.dumps(rows, indent=1))
