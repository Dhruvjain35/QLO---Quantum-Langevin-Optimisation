"""Stage 7 analyses. Every comparison uses the SAME theta for both measurement schemes (paired by construction)."""

from __future__ import annotations

import numpy as np
import pandas as pd

from qlo.stage5.theory import log_a_k, s_k
from qlo.stage6.analysis import wilson_ci
from qlo.stage6.directional import direction_probabilities
from qlo.stage7 import required_shots as R
from qlo.stage7.estimators import SCHEMES, pair, pair_delta, snr_squared

LN10 = float(np.log(10.0))
SCHEME_ID = {"loschmidt": 1, "swap": 2, "swap_null": 3}


def sample_theta(n: int, n_samples: int, seed: int, stream: int) -> np.ndarray:
    return np.random.default_rng(np.random.SeedSequence([7, stream, int(n), int(seed)])).uniform(-np.pi, np.pi, (int(n_samples), int(n)))


def _q(x, p):
    x = x[np.isfinite(x)]
    return float(np.percentile(x, p)) if x.size else float("nan")


def _stats(x: np.ndarray, tag: str) -> dict:
    f = np.isfinite(x)
    v = x[f]
    return {f"{tag}_finite_frac": float(f.mean()), f"{tag}_median": _q(v, 50), f"{tag}_q25": _q(v, 25), f"{tag}_q75": _q(v, 75),
            f"{tag}_q90": _q(v, 90), f"{tag}_geomean": float(v.mean()) if v.size else float("nan")}


# ---- common landscape -------------------------------------------------------------------------------
def common_landscape_check(n: int, n_samples: int, seed: int) -> dict:
    """The exact gradient implied by each scheme's outcome probabilities, c (p_pos - p_neg), must equal A s / 2 for
    every theta, i.e. the two schemes estimate the same quantity on the same landscape."""
    th = sample_theta(n, n_samples, seed, 10)
    A, s = np.exp(log_a_k(th, 0)), s_k(th, 0)
    g = A * s / 2.0
    out = {"n": n, "n_samples": n_samples}
    for sc in SCHEMES:
        p_pos, p_neg, c = pair(sc, A, s)
        out[f"{sc}_max_abs_implied_g_minus_g"] = float(np.max(np.abs(c * (p_pos - p_neg) - g)))
        out[f"{sc}_max_abs_c_delta_minus_g"] = float(np.max(np.abs(c * pair_delta(sc, A, s) - g)))
    return out


# ---- required shots (Tasks 9-11) --------------------------------------------------------------------------
TARGETS = (("snr1", "snr", 1.0), ("snr2", "snr", 2.0), ("dir0.75", "dir", 0.75), ("dir0.9", "dir", 0.90), ("nz0.9", "nz", 0.90))


def required_shots_for_n(n: int, n_samples: int, seed: int, keep: int = 2000) -> dict:
    """All targets for both schemes on the same theta. Returns scheme summaries, paired-difference summaries,
    and a small per-sample subsample (for figures)."""
    th = sample_theta(n, n_samples, seed, 20)
    logA, s = log_a_k(th, 0), s_k(th, 0)
    vals, meth = {}, {}
    for sc in SCHEMES:
        for tag, kind, x in TARGETS:
            if kind == "snr":
                v, m = R.snr_shots(sc, logA, s, x), np.full(n_samples, R.METHOD["closed_form"])
            elif kind == "dir":
                v, m = R.direction_shots(sc, logA, s, x)
            else:
                v, m = R.nonzero_shots(sc, logA, s, x)
            vals[(sc, tag)], meth[(sc, tag)] = v, m
    rows, paired = [], {"n": n, "n_samples": n_samples}
    for sc in SCHEMES:
        row = {"n": n, "scheme": sc, "n_samples": n_samples,
               "empirical_mean_logA": float(logA.mean()), "theory_mean_logA": -2.0 * (n - 1) * np.log(2.0),
               "log10_abs_g_median": _q(np.log10(np.exp(logA) * np.abs(s) / 2.0), 50)}
        for tag, _, _ in TARGETS:
            row.update(_stats(vals[(sc, tag)], f"{tag}_log10M"))
            m = meth[(sc, tag)]
            for name, code in R.METHOD.items():
                row[f"{tag}_frac_{name}"] = float(np.mean(m == code))
            row[f"{tag}_unsolved"] = int(np.sum(~np.isfinite(vals[(sc, tag)])))
        rows.append(row)
    for tag, _, _ in TARGETS:
        d = vals[("swap", tag)] - vals[("loschmidt", tag)]
        paired.update(_stats(d, f"{tag}_log10_swap_over_le"))
        paired[f"{tag}_frac_swap_needs_more"] = float(np.mean(d[np.isfinite(d)] > 1e-12))
        paired[f"{tag}_frac_swap_needs_fewer"] = float(np.mean(d[np.isfinite(d)] < -1e-12))
    rng = np.random.default_rng(np.random.SeedSequence([7, 21, n, seed]))
    sub = rng.choice(n_samples, size=min(keep, n_samples), replace=False)
    samples = pd.DataFrame({"n": n, "log10_A": logA[sub] / LN10, "s": s[sub],
                            **{f"{sc}_{tag}": vals[(sc, tag)][sub] for sc in SCHEMES for tag, _, _ in TARGETS}})
    return {"rows": rows, "paired": paired, "samples": samples}


# ---- fixed-M paired probabilities (Tasks 6, 7, 11) --------------------------------------------------------
def paired_fixed_shots(n_values, shots, n_samples: int, seed: int) -> pd.DataFrame:
    rows = []
    for n in n_values:
        th = sample_theta(n, n_samples, seed, 30)
        logA, s = log_a_k(th, 0), s_k(th, 0)
        A = np.exp(logA)
        for M in shots:
            d = {}
            for sc in SCHEMES:
                p_pos, p_neg, _ = pair(sc, A, s)
                pr = direction_probabilities(p_pos, p_neg, float(M))
                nz = 1.0 - pr["p_zero"]
                with np.errstate(invalid="ignore", divide="ignore"):
                    cond = np.where(nz > 0, pr["p_correct"] / nz, np.nan)
                d[sc] = {"p_zero": pr["p_zero"], "p_correct": pr["p_correct"], "p_wrong": pr["p_wrong"], "p_correct_given_nonzero": cond,
                         "snr": np.sqrt(snr_squared(sc, A, s, float(M))), "float_tie": pr["zero_signal"] & (s != 0)}
            row = {"n": n, "shots": int(M), "n_samples": int(n_samples)}
            for sc in SCHEMES:
                for key in ("p_zero", "p_correct", "p_wrong", "p_correct_given_nonzero", "snr"):
                    row[f"{sc}_median_{key}"] = _q(d[sc][key], 50)
                    row[f"{sc}_mean_{key}"] = float(np.nanmean(d[sc][key]))
                row[f"{sc}_frac_p_zero_ge_0.9"] = float(np.mean(d[sc]["p_zero"] >= 0.9))
                row[f"{sc}_frac_p_correct_le_0.6"] = float(np.mean(d[sc]["p_correct"] <= 0.6))
                row[f"{sc}_frac_float_tie"] = float(np.mean(d[sc]["float_tie"]))
            for key in ("p_zero", "p_correct", "p_correct_given_nonzero"):
                diff = d["swap"][key] - d["loschmidt"][key]
                row[f"paired_median_diff_{key}"] = _q(diff, 50)
                row[f"paired_frac_swap_higher_{key}"] = float(np.nanmean(diff > 0))
            with np.errstate(divide="ignore", invalid="ignore"):
                r = np.log10(d["swap"]["snr"] / d["loschmidt"]["snr"])
            row["paired_median_log10_snr_swap_over_le"] = _q(r, 50)
            row["paired_frac_swap_snr_higher"] = float(np.nanmean(r > 0))
            rows.append(row)
    return pd.DataFrame(rows)


# ---- information distances (Task 14) ------------------------------------------------------------------------
def bernoulli_distances(p_a, p_b, delta) -> dict:
    """Per-shot distances between Bernoulli(p_a) and Bernoulli(p_b) with the EXACT difference delta = p_a - p_b
    supplied (cancellation-free for p_a ~ p_b ~ 1/2)."""
    pa, pb, dl = (np.asarray(x, dtype=np.float64) for x in (p_a, p_b, delta))
    with np.errstate(divide="ignore", invalid="ignore"):
        r1 = dl / (np.sqrt(pa) + np.sqrt(pb))
        r2 = -dl / (np.sqrt(1 - pa) + np.sqrt(1 - pb))
        h2 = 0.5 * (r1**2 + r2**2)
        bhat = -np.log1p(-np.minimum(h2, 1.0 - 1e-16))
        pbar = 0.5 * (pa + pb)
        chi2 = dl**2 / (pbar * (1 - pbar))
        small = np.abs(dl) < 1e-6 * np.minimum(pb, 1 - pb)
        kl_exact = pa * np.log1p(dl / pb) + (1 - pa) * np.log1p(-dl / (1 - pb))
        kl = np.where(small, dl**2 / (2 * pb * (1 - pb)), kl_exact)
        kl = np.where((pb > 0) & (pb < 1), kl, np.inf)
    return {"tv": np.abs(dl), "hellinger2": h2, "bhattacharyya": bhat, "kl": kl, "chi2_mid": chi2}


def information_distance_table(n_values, n_samples: int, seed: int) -> pd.DataFrame:
    rows = []
    for n in n_values:
        th = sample_theta(n, n_samples, seed, 40)
        logA, s = log_a_k(th, 0), s_k(th, 0)
        A = np.exp(logA)
        for sc in SCHEMES:
            p_pos, p_neg, _ = pair(sc, A, s)
            dist = bernoulli_distances(p_pos, p_neg, pair_delta(sc, A, s))
            row = {"n": n, "scheme": sc, "n_samples": n_samples}
            for k, v in dist.items():
                f = np.isfinite(v) & (v > 0)
                row[f"median_{k}"] = _q(v[f], 50)
                row[f"median_log10_{k}"] = _q(np.log10(v[f]), 50)
            row["median_log10_shots_1_over_hellinger2"] = -row["median_log10_hellinger2"]
            row["median_log10_shots_1_over_tv2"] = -2 * row["median_log10_tv"]
            rows.append(row)
    return pd.DataFrame(rows)


# ---- prior-work replication (Task 18) ----------------------------------------------------------------------
def prior_work_replication(n_values, shots, n_samples: int, seed: int, exact_tv_samples: int = 400) -> pd.DataFrame:
    """REPLICATIONS of known behaviour (Thanasilp et al. 2024; Aghaei Saem et al. 2026), at the level of the
    FIDELITY estimate: (1) the Loschmidt estimate is mostly exactly 0; (2) the SWAP estimate's distribution
    approaches the parameter-independent F = 0 distribution (ancilla outcomes ~ 50/50)."""
    from scipy.stats import binom

    rows = []
    for n in n_values:
        th = sample_theta(n, n_samples, seed, 50)
        logF = np.sum(np.log(np.cos(th / 2.0) ** 2), axis=1)
        F = np.exp(logF)
        sub = np.random.default_rng(np.random.SeedSequence([7, 51, n, seed])).choice(n_samples, min(exact_tv_samples, n_samples), replace=False)
        for M in shots:
            p0 = np.exp(M * np.log1p(-F))  # P(F_hat_LE = 0) = (1 - F)^M
            snr_swap = F / np.sqrt((1 - F**2) / M)  # |E F_hat| / sd(F_hat) for SWAP
            snr_le = np.where(F < 1, F / np.sqrt(np.maximum(F * (1 - F), 1e-300) / M), np.inf)
            r = np.arange(M + 1)
            null = binom.pmf(r, M, 0.5)
            tv = np.array([0.5 * np.abs(binom.pmf(r, M, (1 + f) / 2) - null).sum() for f in F[sub]])
            # normal-limit TV between Bin(M,(1+F)/2) and Bin(M,1/2): 2 Phi(F sqrt(M)/2) - 1 (valid for small F)
            from scipy.stats import norm

            tv_norm = 2 * norm.cdf(F[sub] * np.sqrt(M) / 2) - 1
            rows.append({"n": n, "shots": int(M), "n_samples": int(n_samples), "median_F": _q(F, 50),
                         "le_median_p_fhat_zero": _q(p0, 50), "le_frac_p_fhat_zero_ge_0.9": float(np.mean(p0 >= 0.9)),
                         "le_median_mean_over_sd": _q(snr_le, 50),
                         "swap_median_mean_over_sd": _q(snr_swap, 50), "swap_frac_mean_over_sd_lt_0.1": float(np.mean(snr_swap < 0.1)),
                         "swap_median_tv_to_null_exact": _q(tv, 50), "swap_q90_tv_to_null_exact": _q(tv, 90),
                         "swap_median_tv_to_null_normal": _q(tv_norm, 50), "tv_samples": int(sub.size)})
    return pd.DataFrame(rows)


# ---- full-gradient vector reliability (Task 13) ---------------------------------------------------------------
def _component_quantities(th: np.ndarray):
    """A_k, s_k for every k (prefix/suffix products, no division)."""
    c2 = np.cos(th / 2.0) ** 2
    n = th.shape[-1]
    pre = np.ones(th.shape[:-1] + (n + 1,))
    suf = np.ones(th.shape[:-1] + (n + 1,))
    pre[..., 1:] = np.cumprod(c2, axis=-1)
    suf[..., :-1] = np.cumprod(c2[..., ::-1], axis=-1)[..., ::-1]
    return pre[..., :n] * suf[..., 1:], np.sin(th)


def sample_gradient(scheme: str, th: np.ndarray, M: int, rng) -> tuple[np.ndarray, np.ndarray]:
    """One full finite-shot gradient per row of theta (independent batch per component and shift) and the exact gradient."""
    A, s = _component_quantities(th)
    if scheme == "swap_null":  # random-walk control: SWAP counts with the signal removed (q_+ = q_- = 1/2)
        ghat = (rng.binomial(int(M), 0.5, A.shape) - rng.binomial(int(M), 0.5, A.shape)) / float(M)
        return ghat, A * s / 2.0
    p_pos, p_neg, c = pair(scheme, A, s)
    ghat = c * (rng.binomial(int(M), p_pos) - rng.binomial(int(M), p_neg)) / float(M)
    return ghat, A * s / 2.0


def vector_reliability(n_values, shots, n_theta: int, replicates: int, seed: int) -> pd.DataFrame:
    rows = []
    for n in n_values:
        th = sample_theta(n, n_theta, seed, 60)
        for M in shots:
            for sc in SCHEMES:
                rng = np.random.default_rng(np.random.SeedSequence([7, 61, n, M, seed, SCHEME_ID[sc]]))
                cos, dot, zero, ratio, czero, csign = [], [], [], [], [], []
                for _ in range(replicates):
                    gh, g = sample_gradient(sc, th, M, rng)
                    hn, gn = np.linalg.norm(gh, axis=1), np.linalg.norm(g, axis=1)
                    dt = np.sum(gh * g, axis=1)
                    with np.errstate(divide="ignore", invalid="ignore"):
                        cos.append(np.where((hn > 0) & (gn > 0), dt / (hn * gn), np.nan))
                        ratio.append(hn / gn)
                    dot.append(dt); zero.append(hn == 0)
                    czero.append(np.mean(gh == 0, axis=1))
                    csign.append(np.mean(np.sign(gh) == np.sign(g), axis=1))
                cos, dot, zero = np.concatenate(cos), np.concatenate(dot), np.concatenate(zero)
                ratio, czero, csign = np.concatenate(ratio), np.concatenate(czero), np.concatenate(csign)
                cos0 = np.where(np.isnan(cos), 0.0, cos)
                nz = ~zero
                rows.append({"n": n, "shots": int(M), "scheme": sc, "n_theta": n_theta, "replicates": replicates,
                             "p_vector_zero": float(zero.mean()), "p_dot_gt_0": float(np.mean(dot > 0)), "p_dot_le_0": float(np.mean(dot <= 0)),
                             "p_dot_lt_0": float(np.mean(dot < 0)),
                             "median_cos": float(np.median(cos0)), "mean_cos": float(np.mean(cos0)),
                             "cos_q10": float(np.percentile(cos0, 10)), "cos_q25": float(np.percentile(cos0, 25)),
                             "cos_q75": float(np.percentile(cos0, 75)), "cos_q90": float(np.percentile(cos0, 90)),
                             "median_cos_given_nonzero": float(np.nanmedian(cos[nz])) if nz.any() else float("nan"),
                             "p_dot_gt_0_given_nonzero": float(np.mean(dot[nz] > 0)) if nz.any() else float("nan"),
                             "median_log10_norm_ratio": float(np.median(np.log10(ratio[np.isfinite(ratio) & (ratio > 0)]))) if np.any(ratio > 0) else float("nan"),
                             "mean_frac_components_zero": float(czero.mean()), "mean_frac_component_signs_correct": float(csign.mean())})
    return pd.DataFrame(rows)


# ---- same-landscape optimization diagnostic (Task 17) ------------------------------------------------------------
def optimization_diagnostic(n_values, variants, n_seeds: int, iters: int, eta: float, seed: int, curve_every: int = 5):
    """Fixed-eta GD, theta <- theta - eta g_hat, identical start theta per (n, seed) for EVERY variant (exact and both
    schemes), no tuning. ``variants``: list of (label, scheme or None for exact, M). Returns (per-run DataFrame,
    per-iteration median-curve DataFrame)."""
    runs, curves = [], []
    for n in n_values:
        th0 = sample_theta(n, n_seeds, seed, 70)
        F0 = np.prod(np.cos(th0 / 2.0) ** 2, axis=1)
        for label, sc, M in variants:
            rng = np.random.default_rng(np.random.SeedSequence([7, 71, n, int(M), seed, SCHEME_ID.get(sc, 0)]))
            th = th0.copy()
            Fbest = F0.copy()
            zero_it = np.zeros(n_seeds); wrong_it = np.zeros(n_seeds); cos_sum = np.zeros(n_seeds); cos_cnt = np.zeros(n_seeds)
            for t in range(iters):
                if sc is None:
                    A, s = _component_quantities(th)
                    g = A * s / 2.0
                    gh = g
                else:
                    gh, g = sample_gradient(sc, th, M, rng)
                hn, gn = np.linalg.norm(gh, axis=1), np.linalg.norm(g, axis=1)
                dt = np.sum(gh * g, axis=1)
                zero_it += hn == 0
                wrong_it += (hn > 0) & (dt <= 0)
                ok = (hn > 0) & (gn > 0)
                cos_sum[ok] += dt[ok] / (hn[ok] * gn[ok]); cos_cnt[ok] += 1
                th = (th - eta * gh + np.pi) % (2 * np.pi) - np.pi
                F = np.prod(np.cos(th / 2.0) ** 2, axis=1)
                Fbest = np.maximum(Fbest, F)
                if t % curve_every == 0 or t == iters - 1:
                    curves.append({"n": n, "variant": label, "iter": t + 1, "median_log10_F": float(np.median(np.log10(np.maximum(F, 1e-300)))),
                                   "q25_log10_F": float(np.percentile(np.log10(np.maximum(F, 1e-300)), 25)),
                                   "q75_log10_F": float(np.percentile(np.log10(np.maximum(F, 1e-300)), 75))})
            for i in range(n_seeds):
                runs.append({"n": n, "variant": label, "scheme": sc or "exact", "shots": int(M), "seed": i, "iters": iters, "eta": eta,
                             "total_measurement_shots": int(2 * n * M * iters) if sc else 0,
                             "total_state_copies": int(2 * n * M * iters * (2 if sc in ("swap", "swap_null") else 1)) if sc else 0,
                             "F0": float(F0[i]), "F_final": float(F[i]), "F_best": float(Fbest[i]),
                             "log10_gain_best": float(np.log10(max(Fbest[i], 1e-300)) - np.log10(max(F0[i], 1e-300))),
                             "frac_zero_iters": float(zero_it[i] / iters), "frac_wrong_iters": float(wrong_it[i] / iters),
                             "mean_cos_nonzero": float(cos_sum[i] / cos_cnt[i]) if cos_cnt[i] else float("nan")})
    return pd.DataFrame(runs), pd.DataFrame(curves)


def summarize_optimization(runs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (n, label), g in runs.groupby(["n", "variant"], sort=False):
        ex = runs[(runs.n == n) & (runs.variant == "exact")].set_index("seed")
        gi = g.set_index("seed")
        k = int((gi.F_best >= 0.5).sum())
        lo, hi = wilson_ci(k, len(gi))
        rows.append({"n": n, "variant": label, "scheme": g.scheme.iloc[0], "shots": int(g.shots.iloc[0]), "seeds": len(gi),
                     "median_F0": float(gi.F0.median()), "median_F_final": float(gi.F_final.median()), "median_F_best": float(gi.F_best.median()),
                     "median_log10_gain": float(gi.log10_gain_best.median()), "p_F_best_ge_0.5": k / len(gi), "p_ci_low": lo, "p_ci_high": hi,
                     "mean_frac_zero_iters": float(gi.frac_zero_iters.mean()), "mean_frac_wrong_iters": float(gi.frac_wrong_iters.mean()),
                     "mean_cos_nonzero": float(gi.mean_cos_nonzero.mean()),
                     "paired_frac_better_than_exact": float(np.mean(gi.F_best > ex.F_best + 1e-9)),
                     "paired_frac_worse_than_exact": float(np.mean(gi.F_best < ex.F_best - 1e-9)),
                     # F_best is a running maximum, which ANY random walk inflates; F_final is the fair endpoint
                     "paired_frac_final_better_than_exact": float(np.mean(gi.F_final > ex.F_final + 1e-9)),
                     "paired_frac_final_worse_than_exact": float(np.mean(gi.F_final < ex.F_final - 1e-9)),
                     "p_F_final_ge_0.5": float(np.mean(gi.F_final >= 0.5)),
                     "total_measurement_shots_per_run": int(gi.total_measurement_shots.iloc[0]),
                     "total_state_copies_per_run": int(gi.total_state_copies.iloc[0])})
    return pd.DataFrame(rows)
