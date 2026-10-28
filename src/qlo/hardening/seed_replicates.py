"""Per-theta Stage 7 quantities for one (seed, n), and seed-level slope fits.

Every scientific number is produced by the frozen Stage 7 functions (``qlo.stage7.required_shots``,
``qlo.stage6.directional``); this module only keeps the PER-THETA arrays that Stage 7 summarized and threw away,
so that they can be replicated over seeds and resampled.

Pairing: one theta array per (seed, n), drawn by ``qlo.stage7.analysis.sample_theta`` with the Stage 7 stream ids,
is used for BOTH measurement schemes. Nothing is ever drawn per scheme.
"""

from __future__ import annotations

import numpy as np

from qlo.stage5.theory import log_a_k, s_k
from qlo.stage6.analysis import linear_fit
from qlo.stage6.directional import direction_probabilities
from qlo.stage6.exact_distribution import difference_probabilities
from qlo.stage7 import analysis as AN
from qlo.stage7 import required_shots as R
from qlo.stage7.estimators import SCHEMES, pair
from qlo.hardening.intervals import ols_slope_se

STAGE7_SEED = 0                                  # the Stage 7 master seed (Stage7Config.seed)
B3_SEEDS = tuple(range(10_000, 10_020))          # 20 independent master seeds, disjoint from Stage 7
REQ_STREAM, FIXED_STREAM = 20, 30                # Stage 7 stream ids (required shots, fixed-M probabilities)
TARGETS = AN.TARGETS                             # (tag, kind, x): snr1, snr2, dir0.75, dir0.9, nz0.9
SERIES = tuple((sc, tag) for sc in SCHEMES for tag, _, _ in TARGETS)   # 10 series, fixed order
PRIMARY_N = (2, 4, 6, 8, 10, 12, 14, 16, 18, 20)  # the Stage 7 fit range; the primary fit never changes
EXTENSION_N = (22, 24)                            # used only in a separately labelled extended fit, if accepted


def check_seeds(seeds, forbidden=(STAGE7_SEED,)) -> tuple[int, ...]:
    seeds = tuple(int(s) for s in seeds)
    if len(set(seeds)) != len(seeds):
        raise ValueError("duplicate master seeds")
    bad = set(seeds) & set(int(f) for f in forbidden)
    if bad:
        raise ValueError(f"master seed(s) {sorted(bad)} reuse the Stage 7 seed")
    return seeds


def series_name(sc: str, tag: str) -> str:
    return f"{sc}_{tag}"


# ---- required shots, per theta ------------------------------------------------------------------------------
def required_shots_per_theta(n: int, n_samples: int, seed: int) -> dict:
    """Per-theta log10 M for all 10 (scheme, target) series on ONE theta array (exactly paired).

    Same calls, same order and same theta stream as ``qlo.stage7.analysis.required_shots_for_n``; its medians
    are therefore reproduced exactly (tested)."""
    th = AN.sample_theta(n, n_samples, seed, REQ_STREAM)
    logA, s = log_a_k(th, 0), s_k(th, 0)
    vals = np.empty((len(SERIES), n_samples))
    meth = np.empty((len(SERIES), n_samples), dtype=np.int8)
    for i, (sc, tag) in enumerate(SERIES):
        kind, x = {t: (k, v) for t, k, v in TARGETS}[tag]
        if kind == "snr":
            v, m = R.snr_shots(sc, logA, s, x), np.full(n_samples, R.METHOD["closed_form"])
        elif kind == "dir":
            v, m = R.direction_shots(sc, logA, s, x)
        else:
            v, m = R.nonzero_shots(sc, logA, s, x)
        vals[i], meth[i] = v, m
    return {"n": int(n), "seed": int(seed), "n_samples": int(n_samples), "logA": logA, "s": s, "vals": vals, "meth": meth}


def finite_median(x) -> float:
    """Median over finite values; same definition as Stage 7 (``analysis._q(x, 50)``)."""
    x = np.asarray(x, dtype=np.float64)
    x = x[np.isfinite(x)]
    return float(np.percentile(x, 50)) if x.size else float("nan")


def medians_and_methods(res: dict) -> list[dict]:
    rows = []
    for i, (sc, tag) in enumerate(SERIES):
        v, m = res["vals"][i], res["meth"][i]
        row = {"seed": res["seed"], "n": res["n"], "scheme": sc, "target": tag, "n_samples": res["n_samples"],
               "median_log10M": finite_median(v), "unsolved": int(np.sum(~np.isfinite(v)))}
        for name, code in R.METHOD.items():
            row[f"frac_{name}"] = float(np.mean(m == code))
        rows.append(row)
    return rows


def fit_series(n_values, medians) -> dict:
    """Stage 7 fit (``qlo.stage6.analysis.linear_fit``: median log10 M = a + b n) plus a diagnostic OLS slope SE."""
    n_values, medians = np.asarray(n_values, dtype=np.float64), np.asarray(medians, dtype=np.float64)
    ok = np.isfinite(medians)
    f = linear_fit(n_values[ok], medians[ok])
    f["ols_slope_se_diagnostic"] = ols_slope_se(n_values[ok], medians[ok])
    f["n_min"], f["n_max"] = int(n_values[ok].min()), int(n_values[ok].max())
    f["n_used"] = " ".join(str(int(x)) for x in n_values[ok])
    return f


# ---- fixed-M probabilities, per theta -----------------------------------------------------------------------------
FIXED_KEYS = ("p_zero", "p_correct", "p_correct_given_nonzero")


def fixed_shot_per_theta(n: int, shots, n_samples: int, seed: int) -> dict:
    """Exact per-theta P_zero, P_correct, P(correct | nonzero) for both schemes at every M, on ONE theta array.
    Same theta stream and same calls as ``qlo.stage7.analysis.paired_fixed_shots`` (medians reproduced, tested)."""
    th = AN.sample_theta(n, n_samples, seed, FIXED_STREAM)
    logA, s = log_a_k(th, 0), s_k(th, 0)
    A = np.exp(logA)
    out = {"n": int(n), "seed": int(seed), "n_samples": int(n_samples), "shots": tuple(int(M) for M in shots),
           "logA": logA, "s": s}
    for sc in SCHEMES:
        p_pos, p_neg, _ = pair(sc, A, s)
        arr = np.empty((len(shots), len(FIXED_KEYS), n_samples))
        tie = np.empty((len(shots), n_samples), dtype=bool)
        for j, M in enumerate(shots):
            pr = direction_probabilities(p_pos, p_neg, float(M))
            nz = 1.0 - pr["p_zero"]
            with np.errstate(invalid="ignore", divide="ignore"):
                cond = np.where(nz > 0, pr["p_correct"] / nz, np.nan)
            arr[j, 0], arr[j, 1], arr[j, 2] = pr["p_zero"], pr["p_correct"], cond
            tie[j] = pr["zero_signal"] & (s != 0)
        out[sc], out[f"{sc}_float_tie"] = arr, tie
    # cancellation-free LE P(correct | nonzero), for the sign-law check at large n (see direction_accurate)
    p_pos, p_neg, _ = pair("loschmidt", A, s)
    out["loschmidt_cond_accurate"] = np.stack([direction_accurate(p_pos, p_neg, float(M))["p_correct_given_nonzero"] for M in shots])
    return out


def direction_accurate(p_pos, p_neg, M) -> dict:
    """P_correct and P(nonzero) with BOTH tails summed directly: P(K_pos > K_neg) and P(K_neg > K_pos) each come from
    ``difference_probabilities`` (the second with the arguments swapped), never as ``1 - P(>) - P(=)``.

    Why: in Stage 6/7 ``P(<)`` (and hence P_correct whenever the signal is negative) and ``1 - P_zero`` are formed by
    subtraction from ~1. When P(nonzero) ~ M A is below ~1e-10 that subtraction keeps only a few significant digits,
    which matters for P(correct | nonzero) deep in the Loschmidt dead zone (n >= 20, small M). Stage 7's own
    definition is kept unchanged elsewhere; this is used only to measure and avoid that rounding."""
    gt, eq, _ = difference_probabilities(p_pos, p_neg, M)
    lt, _, _ = difference_probabilities(p_neg, p_pos, M)
    p_pos, p_neg = np.asarray(p_pos), np.asarray(p_neg)
    correct = np.where(p_pos > p_neg, gt, np.where(p_pos < p_neg, lt, 0.5 * (gt + lt)))
    nz = gt + lt
    with np.errstate(invalid="ignore", divide="ignore"):
        cond = np.where(nz > 0, correct / nz, np.nan)
    return {"p_correct": correct, "p_nonzero": nz, "p_zero": eq, "p_correct_given_nonzero": cond}
