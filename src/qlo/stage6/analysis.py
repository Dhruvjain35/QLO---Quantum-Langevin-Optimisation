"""Matched-signal comparison, failure fractions, information distances, vector reliability."""

from __future__ import annotations

import numpy as np
import pandas as pd

from qlo.stage5.theory import a_k, f_plus_minus
from qlo.stage5.theory import s_k as proj_s_k
from qlo.stage6.directional import direction_probabilities
from qlo.stage6.parity_benchmark import b_k, s_k, shifted_probabilities

Z95 = 1.959963984540054


def wilson_ci(k: int, n: int) -> tuple[float, float]:
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + Z95**2 / n
    c = (p + Z95**2 / (2 * n)) / d
    h = Z95 * np.sqrt(p * (1 - p) / n + Z95**2 / (4 * n**2)) / d
    # the Wilson bounds are exactly 0 at k = 0 and exactly 1 at k = n; c -+ h only rounds to them
    return (0.0 if k == 0 else float(c - h), 1.0 if k == n else float(c + h))


def linear_fit(x, y) -> dict:
    x, y = np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)
    b, a = np.polyfit(x, y, 1)
    pred = a + b * x
    rss = float(np.sum((y - pred) ** 2))
    tss = float(np.sum((y - y.mean()) ** 2))
    return {"slope": float(b), "intercept": float(a), "r2": 1 - rss / tss if tss > 0 else float("nan"),
            "max_abs_residual": float(np.max(np.abs(y - pred))), "n_points": int(x.size)}


def scaling_fit(n, log10M, label: str) -> dict:
    f = linear_fit(n, log10M)
    refs = {"log10_4": float(np.log10(4)), "log10_2": float(np.log10(2))}
    closer = min(refs, key=lambda k: abs(f["slope"] - refs[k]))
    return {"statistic": label, **f, **{f"theory_slope_{k}": v for k, v in refs.items()},
            "abs_err_vs_log10_4": abs(f["slope"] - refs["log10_4"]), "abs_err_vs_log10_2": abs(f["slope"] - refs["log10_2"]),
            "closer_to": closer}


# ---- benchmark-agnostic per-sample probability pairs -------------------------------------------
def parity_pairs(theta: np.ndarray, k: int = 0) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(p_pos, p_neg, g_true)`` for the parity estimator (g_hat ∝ K_+ - K_-)."""
    s, B = s_k(theta, k), b_k(theta, k)
    pp, pm = shifted_probabilities(s, B)
    return pp, pm, 0.5 * s * B


def projector_pairs(theta: np.ndarray, k: int = 0) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(p_pos, p_neg, g_true)`` for the projector estimator (g_hat ∝ K_- - K_+)."""
    A, s = a_k(theta, k), proj_s_k(theta, k)
    fp, fm = f_plus_minus(A, s)
    return fm, fp, 0.5 * A * s


PAIRS = {"parity": parity_pairs, "projector": projector_pairs}
BENCH_ID = {"parity": 1, "projector": 2}  # deterministic seed ids (never hash(): Python string hashing is salted)


# ---- matched-signal analysis ---------------------------------------------------------------------
def matched_signal_table(n_values, n_samples: int, shots, seed: int = 0, k: int = 0,
                         bin_edges=None, min_count: int = 200) -> pd.DataFrame:
    """Pool theta from both benchmarks over ``n_values``, bin by ``log10 |g_true|``, and compare
    P_zero / P_correct / conditional SD at the SAME |g| bin and the SAME M.

    Bins are shared between benchmarks; a (bin, benchmark, M) cell is reported only when the bin
    holds at least ``min_count`` samples for BOTH benchmarks, so no cell mixes incompatible bins.
    """
    if bin_edges is None:
        bin_edges = np.arange(-12.0, 0.51, 0.5)
    rows = []
    store = {}
    for bench, fn in PAIRS.items():
        lg, pp, pn = [], [], []
        for n in n_values:
            rng = np.random.default_rng(np.random.SeedSequence([6, 200, int(n), int(seed), BENCH_ID[bench]]))
            th = rng.uniform(-np.pi, np.pi, size=(int(n_samples), int(n)))
            a, b, g = fn(th, k)
            with np.errstate(divide="ignore"):
                lg.append(np.log10(np.abs(g)))
            pp.append(a)
            pn.append(b)
        store[bench] = (np.concatenate(lg), np.concatenate(pp), np.concatenate(pn))
    idx = {bench: np.digitize(store[bench][0], bin_edges) for bench in PAIRS}
    for bi in range(1, len(bin_edges)):
        counts = {bench: int(np.sum(idx[bench] == bi)) for bench in PAIRS}
        if min(counts.values()) < min_count:
            continue
        for bench in PAIRS:
            lg, pp, pn = store[bench]
            m = idx[bench] == bi
            # cap the per-cell work; a random subsample of the bin is unbiased for the bin mean
            sel = np.flatnonzero(m)
            if sel.size > 4000:
                sel = np.random.default_rng(bi).choice(sel, 4000, replace=False)
            for M in shots:
                d = direction_probabilities(pp[sel], pn[sel], float(M))
                sd = np.sqrt((pp[sel] * (1 - pp[sel]) + pn[sel] * (1 - pn[sel])) / (4.0 * M))
                rows.append({"bin_low": float(bin_edges[bi - 1]), "bin_high": float(bin_edges[bi]),
                             "log10_abs_g_mid": float(0.5 * (bin_edges[bi - 1] + bin_edges[bi])), "benchmark": bench,
                             "shots": int(M), "n_in_bin": counts[bench], "n_used": int(sel.size),
                             "mean_p_zero": float(d["p_zero"].mean()), "median_p_zero": float(np.median(d["p_zero"])),
                             "mean_p_correct": float(d["p_correct"].mean()), "median_p_correct": float(np.median(d["p_correct"])),
                             "mean_p_wrong": float(d["p_wrong"].mean()),
                             "mean_estimator_sd": float(sd.mean()), "mean_abs_g": float(np.mean(10.0 ** lg[sel])),
                             "mean_p_pos": float(pp[sel].mean()), "mean_p_neg": float(pn[sel].mean())})
    return pd.DataFrame(rows)


# ---- fixed-shot failure fractions -------------------------------------------------------------------
def failure_fractions(n_values, shots, n_samples: int, seed: int = 0, k: int = 0,
                      p_zero_level: float = 0.9, p_correct_level: float = 0.60) -> pd.DataFrame:
    rows = []
    for bench, fn in PAIRS.items():
        for n in n_values:
            rng = np.random.default_rng(np.random.SeedSequence([6, 300, int(n), int(seed), BENCH_ID[bench]]))
            th = rng.uniform(-np.pi, np.pi, size=(int(n_samples), int(n)))
            pp, pn, g = fn(th, k)
            for M in shots:
                d = direction_probabilities(pp, pn, float(M))
                rows.append({"benchmark": bench, "n": n, "shots": int(M), "n_samples": int(n_samples),
                             "frac_p_zero_ge": float(np.mean(d["p_zero"] >= p_zero_level)),
                             "frac_p_correct_le": float(np.mean(d["p_correct"] <= p_correct_level)),
                             "frac_both": float(np.mean((d["p_zero"] >= p_zero_level) & (d["p_correct"] <= p_correct_level))),
                             "frac_resolved": float(np.mean((d["p_zero"] < p_zero_level) & (d["p_correct"] > p_correct_level))),
                             "median_p_zero": float(np.median(d["p_zero"])), "median_p_correct": float(np.median(d["p_correct"])),
                             "p_zero_level": p_zero_level, "p_correct_level": p_correct_level})
    return pd.DataFrame(rows)


# ---- information distances between the two shifted Bernoulli distributions ---------------------------
def information_distances(p_a, p_b) -> dict[str, np.ndarray]:
    """Per-shot distances between ``Bernoulli(p_a)`` and ``Bernoulli(p_b)``."""
    pa, pb = np.asarray(p_a, dtype=np.float64), np.asarray(p_b, dtype=np.float64)
    tv = np.abs(pa - pb)
    hell2 = 1.0 - np.sqrt(pa * pb) - np.sqrt((1 - pa) * (1 - pb))
    with np.errstate(divide="ignore", invalid="ignore"):
        kl = np.where((pa > 0) & (pa < 1) & (pb > 0) & (pb < 1),
                      pa * np.log(pa / pb) + (1 - pa) * np.log((1 - pa) / (1 - pb)), np.inf)
        # symmetric chi^2-type / local Fisher scale: (p_a - p_b)^2 / [p(1-p)] at the midpoint
        mid = 0.5 * (pa + pb)
        chi2_mid = np.where((mid > 0) & (mid < 1), (pa - pb) ** 2 / (mid * (1 - mid)), np.inf)
    return {"tv": tv, "hellinger2": np.clip(hell2, 0.0, 1.0), "kl": kl, "chi2_mid": chi2_mid}


def information_distance_table(n_values, n_samples: int, seed: int = 0, k: int = 0) -> pd.DataFrame:
    rows = []
    for bench, fn in PAIRS.items():
        for n in n_values:
            rng = np.random.default_rng(np.random.SeedSequence([6, 400, int(n), int(seed), BENCH_ID[bench]]))
            th = rng.uniform(-np.pi, np.pi, size=(int(n_samples), int(n)))
            pp, pn, g = fn(th, k)
            d = information_distances(pp, pn)
            row = {"benchmark": bench, "n": n, "n_samples": int(n_samples),
                   "median_abs_g": float(np.median(np.abs(g))), "median_p_pos": float(np.median(pp))}
            for name, v in d.items():
                finite = np.isfinite(v)
                row[f"median_{name}"] = float(np.median(v[finite]))
                row[f"mean_{name}"] = float(np.mean(v[finite]))
            row["median_shots_1_over_hellinger2"] = float(np.median(1.0 / np.maximum(d["hellinger2"], 1e-300)))
            row["median_shots_1_over_tv2"] = float(np.median(1.0 / np.maximum(d["tv"] ** 2, 1e-300)))
            row["median_shots_4_over_chi2"] = float(np.median(4.0 / np.maximum(d["chi2_mid"], 1e-300)))
            rows.append(row)
    return pd.DataFrame(rows)


# ---- full-vector reliability -----------------------------------------------------------------------
def vector_reliability(n_values, shots, n_theta: int, replicates: int, seed: int = 0) -> pd.DataFrame:
    """Full-gradient finite-shot estimates with one independent batch per (k, shift).

    Reports cosine similarity with the exact gradient, P(dot > 0) (descent direction),
    and P(g_hat vector = 0). Components share theta, so their errors are NOT independent
    across k in general; no independence is assumed here — everything is measured directly.
    """
    rows = []
    for bench, fn in PAIRS.items():
        for n in n_values:
            rng = np.random.default_rng(np.random.SeedSequence([6, 500, int(n), int(seed), BENCH_ID[bench]]))
            th = rng.uniform(-np.pi, np.pi, size=(int(n_theta), int(n)))
            P_pos = np.empty((n_theta, n))
            P_neg = np.empty((n_theta, n))
            G = np.empty((n_theta, n))
            for k in range(n):
                P_pos[:, k], P_neg[:, k], G[:, k] = fn(th, k)
            gnorm = np.linalg.norm(G, axis=1)
            for M in shots:
                cos_all, dot_all, zero_all = [], [], []
                for _ in range(replicates):
                    kp = rng.binomial(int(M), P_pos)
                    km = rng.binomial(int(M), P_neg)
                    ghat = (kp - km) / (2.0 * M)  # sign convention already folded into (p_pos, p_neg)
                    hn = np.linalg.norm(ghat, axis=1)
                    dot = np.sum(ghat * G, axis=1)
                    with np.errstate(invalid="ignore", divide="ignore"):
                        # a zero estimate means "no step": cos is undefined; stored as NaN, and the all-draws
                        # summaries below count it as cos = 0 (neither aligned nor anti-aligned)
                        cos = np.where((hn > 0) & (gnorm > 0), dot / (hn * gnorm), np.nan)
                    cos_all.append(cos)
                    dot_all.append(dot)
                    zero_all.append(hn == 0)
                cos = np.concatenate(cos_all)
                dot = np.concatenate(dot_all)
                zero = np.concatenate(zero_all)
                nz = ~zero
                cos0 = np.where(np.isnan(cos), 0.0, cos)  # zero step -> cos 0 (all-draws summaries)
                rows.append({"benchmark": bench, "n": n, "shots": int(M), "n_theta": int(n_theta), "replicates": int(replicates),
                             "p_vector_zero": float(zero.mean()), "p_dot_gt_0": float(np.mean(dot > 0)),
                             "p_dot_le_0": float(np.mean(dot <= 0)), "p_dot_lt_0": float(np.mean(dot < 0)),
                             "median_cos": float(np.median(cos0)), "mean_cos": float(np.mean(cos0)),
                             "cos_q25": float(np.percentile(cos0, 25)), "cos_q75": float(np.percentile(cos0, 75)),
                             "p_cos_gt_0p5": float(np.mean(cos0 > 0.5)), "p_abs_cos_lt_0p1_given_nonzero":
                                 float(np.mean(np.abs(cos[nz]) < 0.1)) if nz.any() else float("nan"),
                             "p_dot_gt_0_given_nonzero": float(np.mean(dot[nz] > 0)) if nz.any() else float("nan"),
                             "median_cos_given_nonzero": float(np.nanmedian(cos[nz])) if nz.any() else float("nan")})
    return pd.DataFrame(rows)
