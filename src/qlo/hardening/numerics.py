"""Acceptance checks for required-shot points beyond the Stage 7 range (n = 22, 24).

Predeclared rule: an extension n is ACCEPTED into the (separately labelled) extended fit only if every check below
passes; otherwise it is reported as numerically unresolved and the fit keeps the Stage 7 range. The primary fit is
always n = 2..20 regardless. No method is changed to make a point pass.

Per-point checks (``validate_required_point``):
  all_finite                 every per-theta log10 M of all 10 series is finite (no unsolved rows)
  no_underflow               log A finite and A = exp(log A) > 1e-280 for every theta (no float underflow)
  swap_exact_rows_no_tie     no SWAP row solved by EXACT bisection has float q_- == q_+ (P_correct would be
                             silently pinned to (1 - P_zero)/2 ~ 1/2)
  le_direction_inversion     on a subsample of LE exact-bisection rows: P_correct(M) >= q and P_correct(M_b) < q,
                             M_b = min(M - 1, floor(M (1 - INV_REL_TOL)))  (the Stage 6/7 bisection is in log10 M with a
                             finite number of halvings, so the returned M is minimal only to ~1e-7 relative)
  le_nonzero_inversion       the same for P(g != 0 | M) >= 0.9
                             (P(M) >= q is tested with tolerance PROB_TOL = 1e-9: above 2^53 shots M is not an exact integer
                             in float64 and the bisection stops on the threshold to ~1e-12)
  Rows returned as M = 2 whose true answer is M = 1 are a known property of the Stage 6/7 bisection (it never
  evaluates exactly 10^lo = 1); they are counted separately as ``floor_artifact`` and do not fail the point.
Method check (``swap_normal_handover_check``, n-independent): the continuity-corrected normal inversion used for
SWAP direction targets beyond 1e4 shots agrees with exact integer inversion to <= HANDOVER_TOL decades at every
tested budget decade, and its error does not grow with M.
"""

from __future__ import annotations

import numpy as np
from scipy.stats import norm

from qlo.stage6.directional import direction_probabilities, shots_for_direction_exact
from qlo.stage7 import required_shots as R
from qlo.stage7.estimators import pair
from qlo.hardening.seed_replicates import SERIES

MIN_A = 1e-280
HANDOVER_TOL = 1e-3   # decades; predeclared (Stage 7 measured <= 1.3e-4 at its handover)
N_VERIFY = 300
INV_REL_TOL = 1e-5
PROB_TOL = 1e-9      # at M > 2^53 the bisection lands on the threshold to ~1e-12; integers are not representable there


def validate_required_point(res: dict, n_verify: int = N_VERIFY, seed: int = 0) -> dict:
    vals, meth, logA, s = res["vals"], res["meth"], res["logA"], res["s"]
    A = np.exp(logA)
    chk = {"n": res["n"], "seed": res["seed"]}
    chk["all_finite"] = bool(np.all(np.isfinite(vals)))
    chk["min_log10_A"] = float(np.min(logA) / np.log(10.0)) if np.all(np.isfinite(logA)) else float("-inf")
    chk["no_underflow"] = bool(np.all(np.isfinite(logA)) and np.all(A > MIN_A))
    rng = np.random.default_rng(np.random.SeedSequence([0xB3, 9, int(res["n"]), int(res["seed"])]))
    ties = 0
    for i, (sc, tag) in enumerate(SERIES):
        if sc == "swap" and tag.startswith("dir"):
            ex = meth[i] == R.METHOD["exact"]
            p_pos, p_neg, _ = pair("swap", A[ex], s[ex])
            ties += int(np.sum(p_pos == p_neg))
    chk["swap_exact_rows_float_ties"] = ties
    chk["swap_exact_rows_no_tie"] = ties == 0
    worst_ok = True
    for i, (sc, tag) in enumerate(SERIES):
        if sc != "loschmidt" or tag == "snr1" or tag == "snr2":
            continue
        rows = np.flatnonzero((meth[i] == R.METHOD["exact"]) & np.isfinite(vals[i]))
        if rows.size == 0:
            chk[f"le_{tag}_verified"] = 0
            continue
        rows = rng.choice(rows, min(n_verify, rows.size), replace=False)
        M = np.round(np.power(10.0, vals[i][rows]))
        Mb = np.maximum(np.minimum(M - 1, np.floor(M * (1 - INV_REL_TOL))), 1)
        p_pos, p_neg, _ = pair("loschmidt", A[rows], s[rows])
        if tag.startswith("dir"):
            q = 0.75 if tag == "dir0.75" else 0.90
            prob = lambda m: direction_probabilities(p_pos, p_neg, m)["p_correct"]
        else:
            q = 0.90
            prob = lambda m: 1.0 - direction_probabilities(p_pos, p_neg, m)["p_zero"]
        at, below = prob(M), prob(Mb)
        floor = (M == 2) & (below >= q)                     # true answer is M = 1 (Stage 6/7 bisection floor)
        ok = (at >= q - PROB_TOL) & ((M == 1) | (below < q) | floor)
        chk[f"le_{tag}_verified"] = int(rows.size)
        chk[f"le_{tag}_floor_artifact"] = int(np.sum(floor))
        chk[f"le_{tag}_inversion_failures"] = int(np.sum(~ok))
        worst_ok &= bool(ok.all())
    chk["le_inversions_ok"] = worst_ok
    chk["accepted"] = bool(chk["all_finite"] and chk["no_underflow"] and chk["swap_exact_rows_no_tie"] and worst_ok)
    return chk


def swap_normal_handover_check(decades=(3.5, 4.5, 5.5, 6.5), n_rows: int = 150, seed: int = 0) -> dict:
    """Normal-cc vs exact integer inversion of SWAP P_correct at increasing budgets. For each target budget
    decade D, A and s are chosen so that the normal answer is ~10^D; the exact inversion brackets it by +-0.3 dec."""
    rng = np.random.default_rng(np.random.SeedSequence([0xB3, 10, int(seed)]))
    out = []
    for q in (0.75, 0.90):
        for D in decades:
            s = rng.uniform(0.2, 1.0, n_rows)
            # d = A|s|/2 with M ~ z^2 v / d^2, v ~ 1/2  =>  A ~ 2 z sqrt(1/2) / (|s| sqrt(10^D))
            A = np.clip(2 * norm.ppf(q) * np.sqrt(0.5) / (s * np.sqrt(10.0**D)), 1e-300, 0.5)
            v = (2 - A**2 * (1 + s**2) / 2) / 4
            ln = R.log10_normal_cc_direction(np.log(A * s / 2), v, q)
            pp, pn, _ = pair("swap", A, s)
            me, ok = shots_for_direction_exact(pp, pn, q, lo_log10=ln - 0.3, hi_log10=ln + 0.3, n_iter=26)
            err = np.abs(np.log10(me) - ln)
            out.append({"q": q, "target_log10M": D, "rows": n_rows, "all_bracketed": bool(ok.all()),
                        "max_abs_log10_err": float(err[ok].max()) if ok.any() else float("nan"),
                        "median_abs_log10_err": float(np.median(err[ok])) if ok.any() else float("nan")})
    ok = all(r["all_bracketed"] and r["max_abs_log10_err"] <= HANDOVER_TOL for r in out)
    for q in (0.75, 0.90):
        e = [r["max_abs_log10_err"] for r in out if r["q"] == q]
        ok &= all(e[i + 1] <= e[i] + 2e-5 for i in range(len(e) - 1))   # 2e-5: integer-ceil granularity at 10^3.5
    return {"rows": out, "tolerance_decades": HANDOVER_TOL, "passed": bool(ok)}
