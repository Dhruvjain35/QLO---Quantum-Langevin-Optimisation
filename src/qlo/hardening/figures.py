"""B3 figures: uncertainty-aware versions of the key Stage 7 figures. Plain matplotlib; every seed is drawn
(no selected examples). SEED and THETA intervals are always drawn as separate, labelled marks."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from qlo.hardening.seed_replicates import SERIES, TARGETS

COLOR = {"loschmidt": "C0", "swap": "C3"}
LABEL = {"loschmidt": "Loschmidt / projector", "swap": "SWAP test"}
SLOPE_TAGS = ("snr1", "snr2", "dir0.75", "dir0.9")


def _plt():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def _save(fig, path: Path, plt) -> Path:
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def slope_seed_distribution(fits: pd.DataFrame, path: Path) -> Path:
    plt = _plt()
    f = fits[fits.fit_range == "primary_2_20"]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, sc, ref in ((axes[0], "loschmidt", np.log10(4)), (axes[1], "swap", np.log10(16))):
        tags = [t for t, _, _ in TARGETS if not (sc == "swap" and t.startswith("nz"))]
        for i, t in enumerate(tags):
            v = f[(f.scheme == sc) & (f.target == t)].slope.values
            ax.plot(np.full(v.size, i) + np.linspace(-0.15, 0.15, v.size), v, "o", ms=3.5, color=COLOR[sc], alpha=0.75)
            ax.plot([i - 0.25, i + 0.25], [np.median(v)] * 2, color="black", lw=1.2)
        ax.axhline(ref, color="grey", ls=":", lw=1, label=f"reference {ref:.4f}")
        ax.set_xticks(range(len(tags)), tags)
        ax.set_ylabel("fitted slope (median log10 M per qubit), n = 2..20")
        ax.set_title(f"{LABEL[sc]}: one point per independent seed (line = median)", fontsize=9)
        ax.grid(True, alpha=0.3); ax.legend(fontsize=7)
    return _save(fig, path, plt)


def slope_confidence_intervals(chk: pd.DataFrame, path: Path) -> Path:
    plt = _plt()
    c = chk[(chk.fit_range == "primary_2_20") & chk.metric.str.startswith("slope_") & ~chk.metric.str.endswith("nz0.9")
            | (chk.metric == "slope_loschmidt_nz0.9") & (chk.fit_range == "primary_2_20")].copy()
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))
    for ax, sc in zip(axes, ("loschmidt", "swap")):
        d = c[c.metric.str.startswith(f"slope_{sc}_")].reset_index(drop=True)
        for i, r in d.iterrows():
            ax.plot([i - 0.12, i - 0.12], [r.seed_pct_lo, r.seed_pct_hi], color=COLOR[sc], lw=3, label="SEED 95% percentile interval" if i == 0 else None)
            ax.plot([i + 0.12, i + 0.12], [r.theta_lo, r.theta_hi], color="black", lw=3, label="THETA-bootstrap 95% CI" if i == 0 else None)
            ax.plot(i, float(r.stage7_estimate), "D", color="gold", mec="black", ms=6, label="Stage 7 estimate" if i == 0 else None)
        ax.set_xticks(range(len(d)), [m.split("_", 2)[2] for m in d.metric])
        ax.set_title(f"{LABEL[sc]}: Stage 7 slope vs hardened intervals (n = 2..20)", fontsize=9)
        ax.set_ylabel("slope"); ax.grid(True, alpha=0.3); ax.legend(fontsize=7)
    return _save(fig, path, plt)


def exponent_gap(fits: pd.DataFrame, ssum: pd.DataFrame, tsum: pd.DataFrame, st7: dict, path: Path) -> Path:
    plt = _plt()
    f = fits[fits.fit_range == "primary_2_20"]
    fig, ax = plt.subplots(figsize=(7.5, 4.4))
    for i, t in enumerate(SLOPE_TAGS):
        a = f[(f.scheme == "loschmidt") & (f.target == t)].sort_values("seed")
        b = f[(f.scheme == "swap") & (f.target == t)].sort_values("seed")
        gap = b.slope.values - a.slope.values
        ax.plot(np.full(gap.size, i) + np.linspace(-0.2, 0.0, gap.size), gap, "o", ms=3.5, color="C2", alpha=0.75,
                label="per-seed paired gap" if i == 0 else None)
        s = ssum[(ssum.fit_range == "primary_2_20") & (ssum.metric == f"gap_swap_minus_le_{t}")].iloc[0]
        th = tsum[(tsum.fit_range == "primary_2_20") & (tsum.metric == f"gap_swap_minus_le_{t}")].iloc[0]
        ax.plot([i + 0.1] * 2, [s.boot_lo, s.boot_hi], color="C2", lw=3, label="SEED: 95% CI of mean gap" if i == 0 else None)
        ax.plot([i + 0.2] * 2, [th.pct_lo, th.pct_hi], color="black", lw=3, label="THETA-bootstrap 95% CI" if i == 0 else None)
        v = st7.get(("primary_2_20", f"gap_swap_minus_le_{t}"))
        if v is not None:
            ax.plot(i + 0.15, v, "D", color="gold", mec="black", ms=6, label="Stage 7" if i == 0 else None)
    ax.axhline(np.log10(4), color="grey", ls=":", lw=1, label="log10 4 = 0.602")
    ax.set_xticks(range(len(SLOPE_TAGS)), SLOPE_TAGS); ax.set_ylabel("γ_SWAP − γ_LE (same seed, same θ)")
    ax.set_title("Paired exponent gap, n = 2..20", fontsize=9); ax.grid(True, alpha=0.3); ax.legend(fontsize=7)
    return _save(fig, path, plt)


def required_shots_with_intervals(med: pd.DataFrame, seeds, ref_seed: int, boot_meds: np.ndarray, fit_n, path: Path) -> Path:
    plt = _plt()
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4), sharey=True)
    for ax, tag in zip(axes, ("snr1", "dir0.9")):
        for sc in ("loschmidt", "swap"):
            d = med[(med.scheme == sc) & (med.target == tag)]
            s = d[d.seed.isin(seeds)].groupby("n").median_log10M
            n = np.array(sorted(s.groups))
            lo, hi = s.quantile(0.025).loc[n].values, s.quantile(0.975).loc[n].values
            ax.fill_between(n, lo, hi, color=COLOR[sc], alpha=0.25, label=f"{LABEL[sc]}: SEED 2.5-97.5%")
            r = d[d.seed == ref_seed].sort_values("n")
            ax.plot(r.n, r.median_log10M, "o-", color=COLOR[sc], ms=3, lw=1, label=f"{LABEL[sc]}: reference (Stage 7 θ)")
            i = SERIES.index((sc, tag))
            b = boot_meds[:, i, :]
            ax.errorbar(np.asarray(fit_n) + 0.25, np.median(b, axis=0), yerr=[np.median(b, 0) - np.percentile(b, 2.5, 0), np.percentile(b, 97.5, 0) - np.median(b, 0)],
                        fmt="none", ecolor="black", elinewidth=1, capsize=2, label="THETA-bootstrap 95%" if sc == "swap" else None)
        ax.set_xlabel("n"); ax.set_title(f"median required shots, target {tag}", fontsize=9); ax.grid(True, alpha=0.3); ax.legend(fontsize=6)
    axes[0].set_ylabel("median log10 M")
    return _save(fig, path, plt)


def _prob_panel(fsum: pd.DataFrame, quantity: str, ax, shots=(64, 1024, 16384, 65536)):
    for M, ls in zip(shots, ("-", "--", ":", "-.")):
        for sc in ("loschmidt", "swap"):
            d = fsum[(fsum.scheme == sc) & (fsum.quantity == quantity) & (fsum.shots == M)].sort_values("n")
            if d.empty:
                continue
            ax.plot(d.n, d.reference_median, ls=ls, color=COLOR[sc], lw=1, label=f"{LABEL[sc]}, M={M}")
            ax.fill_between(d.n, d.SEED_pct_lo, d.SEED_pct_hi, color=COLOR[sc], alpha=0.15)
            ax.errorbar(d.n, d.reference_median, yerr=[d.reference_median - d.THETA_lo, d.THETA_hi - d.reference_median], fmt="none", ecolor="black", elinewidth=0.8, capsize=1.5)


def zero_probability_with_intervals(fsum: pd.DataFrame, path: Path) -> Path:
    plt = _plt()
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    _prob_panel(fsum, "p_zero", ax)
    ax.set_xlabel("n"); ax.set_ylabel("median P(ĝ_k = 0 | θ, M)"); ax.set_ylim(-0.02, 1.02); ax.grid(True, alpha=0.3); ax.legend(fontsize=6)
    ax.set_title("Exact-zero probability: band = SEED 95% percentile, bars = THETA-bootstrap 95% CI", fontsize=8)
    return _save(fig, path, plt)


def directional_correctness_with_intervals(fsum: pd.DataFrame, path: Path) -> Path:
    plt = _plt()
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6), sharey=True)
    _prob_panel(fsum, "p_correct", axes[0])
    _prob_panel(fsum, "p_correct_given_nonzero", axes[1])
    axes[0].set_ylabel("median probability")
    axes[0].set_title("P_correct", fontsize=9); axes[1].set_title("P(correct | ĝ ≠ 0)  (Stage 7 subtraction form)", fontsize=9)
    for ax in axes:
        ax.set_xlabel("n"); ax.set_ylim(-0.02, 1.02); ax.grid(True, alpha=0.3)
    axes[1].legend(fontsize=6)
    return _save(fig, path, plt)


def vector_alignment_with_intervals(vsum: pd.DataFrame, path: Path) -> Path:
    plt = _plt()
    metrics = (("p_vector_zero", "P(full ĝ = 0)"), ("median_cos", "median cos(ĝ, g) (zero vector = 0)"),
               ("median_cos_given_nonzero", "median cos | ĝ ≠ 0"), ("p_dot_gt_0", "P(ĝ·g > 0)"),
               ("mean_frac_component_signs_correct", "component sign accuracy"), ("median_log10_norm_ratio", "median log10 ‖ĝ‖/‖g‖"))
    fig, axes = plt.subplots(2, 3, figsize=(13, 7))
    for ax, (m, lab) in zip(axes.ravel(), metrics):
        for sc in ("loschmidt", "swap"):
            for M, ls in ((64, "--"), (1024, "-")):
                d = vsum[(vsum.scheme == sc) & (vsum.metric == m) & (vsum.shots == M)].sort_values("n")
                ax.plot(d.n, d["median"], ls=ls, marker="o", ms=3, color=COLOR[sc], label=f"{LABEL[sc]}, M={M}")
                ax.fill_between(d.n, d.pct_lo, d.pct_hi, color=COLOR[sc], alpha=0.15)
                ax.plot(d.n, d.stage7_value, "x", color="black", ms=4)
        ax.set_title(lab, fontsize=9); ax.set_xlabel("n"); ax.grid(True, alpha=0.3)
    axes[0, 0].legend(fontsize=6)
    fig.suptitle("Full-gradient reliability over 20 independent seeds (band: SEED 95% percentile; ×: Stage 7)", fontsize=9)
    return _save(fig, path, plt)


def conditional_sign_validation(sl: pd.DataFrame, path: Path) -> Path:
    plt = _plt()
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))
    ax = axes[0]
    for M, c in zip(sorted(sl.shots.unique()), plt.cm.viridis(np.linspace(0, 0.9, sl.shots.nunique()))):
        d = sl[sl.shots == M].sort_values("n")
        ax.errorbar(d.n + 0.05 * np.log2(M / 16), d.reference_observed_median_accurate,
                    yerr=[d.reference_observed_median_accurate - d.THETA_lo, d.THETA_hi - d.reference_observed_median_accurate],
                    fmt="o", ms=3, color=c, capsize=2, label=f"M={M} (THETA CI)")
        ax.plot(d.n + 0.05 * np.log2(M / 16) + 0.03, d.SEED_mean, "s", ms=2.5, color=c, alpha=0.6)
        ax.vlines(d.n + 0.05 * np.log2(M / 16) + 0.03, d.SEED_pct_lo, d.SEED_pct_hi, color=c, alpha=0.6, lw=1)
    ax.axhline(sl.prediction_fixed.iloc[0], color="black", ls=":", label="fixed prediction (1+1/√2)/2")
    ax.set_xlabel("n"); ax.set_ylabel("median P(correct | ĝ_LE ≠ 0)"); ax.grid(True, alpha=0.3); ax.legend(fontsize=6)
    ax.set_title("Loschmidt conditional sign law (squares/lines: SEED mean, 95% percentile)", fontsize=8)
    ax = axes[1]
    deep = sl[sl.deep_regime]
    for M, c in zip(sorted(sl.shots.unique()), plt.cm.viridis(np.linspace(0, 0.9, sl.shots.nunique()))):
        d = deep[deep.shots == M].sort_values("n")
        if d.empty:
            continue
        ax.errorbar(d.n + 0.05 * np.log2(M / 16), d.per_theta_disc_median, yerr=[d.per_theta_disc_median - d.per_theta_disc_q05, d.per_theta_disc_q95 - d.per_theta_disc_median],
                    fmt="o", ms=3, color=c, capsize=2, label=f"M={M}")
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xlabel("n"); ax.set_ylabel("per-θ exact − (1+|s|)/2  (median, 5-95%)"); ax.set_yscale("symlog", linthresh=1e-8)
    ax.set_title("Per-θ discrepancy, deep-regime cells only (reference median P_zero ≥ 0.99)", fontsize=8); ax.grid(True, alpha=0.3); ax.legend(fontsize=6)
    return _save(fig, path, plt)
