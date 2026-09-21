"""Stage 5 figures. Plain matplotlib."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from qlo.stage5.analysis import LOG10_4


def _plt():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def fig_zero_probability_vs_shots(rep: pd.DataFrame, path: Path, mc: pd.DataFrame | None = None) -> Path:
    plt = _plt()
    keys = rep[["n", "label"]].drop_duplicates().values.tolist()
    fig, axes = plt.subplots(1, len(keys), figsize=(3.4 * len(keys), 3.6), squeeze=False, sharey=True)
    for ax, (n, lab) in zip(axes[0], keys):
        d = rep[(rep.n == n) & (rep.label == lab)].sort_values("shots")
        ax.semilogx(d.shots, d.p_zero_exact, "-", label="exact P(ĝ_k=0)")
        ax.semilogx(d.shots, d.p_both_zero, "--", label="both-counts-zero bound")
        ax.semilogx(d.shots, d.p_zero_rare_event, ":", label="rare-event exp(−MA)")
        ax.semilogx(d.shots, d.p_zero_poisson, "-.", label="Poisson limit e^{−MA} I₀(MA|cos θ_k|)")
        if mc is not None:
            e = mc[(mc.n == n) & (mc.label == lab)]
            ax.errorbar(e.shots, e.p_zero_empirical, yerr=np.vstack([e.p_zero_empirical - e.p_zero_ci_low, e.p_zero_ci_high - e.p_zero_empirical]), fmt="o", ms=4, capsize=2, color="black", label=f"binomial MC ({int(e.replicates.iloc[0])} reps, Wilson 95%)" if len(e) else None)
        ax.set_title(f"n={n}, {lab}\nlog10A={d.log10A.iloc[0]:.2f}, s={d.s.iloc[0]:.2f}", fontsize=8); ax.set_xlabel("M"); ax.grid(True, which="both", alpha=0.3)
    axes[0][0].set_ylabel("P(ĝ_k = 0 | θ, M)"); axes[0][0].legend(fontsize=6)
    fig.tight_layout(); fig.savefig(path, dpi=140); plt.close(fig)
    return path


def _band(ax, n, med, lo, hi, label, marker="o"):
    ax.plot(n, med, marker=marker, label=label)
    ax.fill_between(n, lo, hi, alpha=0.15)


def fig_shots_vs_n(summary: pd.DataFrame, path: Path, tags: list[tuple[str, str]], title: str, fits: dict | None = None) -> Path:
    plt = _plt()
    fig, ax = plt.subplots(figsize=(6, 4.5))
    s = summary.sort_values("n")
    for tag, lab in tags:
        _band(ax, s.n, s[f"{tag}_log10M_median"], s[f"{tag}_log10M_q25"], s[f"{tag}_log10M_q75"], lab)
        if fits and tag in fits:
            f = fits[tag]; ax.plot(s.n, f["intercept"] + f["slope"] * s.n, "--", color="gray", alpha=0.7, label=f"fit {lab}: slope {f['slope']:.3f}" if tag == tags[0][0] else None)
    ax.plot(s.n, s[f"{tags[0][0]}_log10M_median"].iloc[0] + LOG10_4 * (s.n - s.n.iloc[0]), ":", color="black", label="slope log10 4 (log-typical A) reference")
    ax.set_xlabel("n"); ax.set_ylabel("log10 M required (median, IQR band; 100k θ/n)")
    ax.set_title(title, fontsize=8); ax.grid(True, alpha=0.3); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)
    return path


def fig_global_vs_local(g: pd.DataFrame, l: pd.DataFrame, path: Path) -> Path:
    plt = _plt()
    fig, ax = plt.subplots(figsize=(6, 4.5))
    g, l = g.sort_values("n"), l.sort_values("n")
    for tag, lab, m in (("nonzero0.9", "P(ĝ≠0) ≥ 0.9", "o"), ("snr2", "SNR ≥ 2", "s")):
        _band(ax, g.n, g[f"{tag}_log10M_median"], g[f"{tag}_log10M_q25"], g[f"{tag}_log10M_q75"], f"global projector, {lab}", m)
        _band(ax, l.n, l[f"{tag}_log10M_median"], l[f"{tag}_log10M_q25"], l[f"{tag}_log10M_q75"], f"local term-wise, {lab}", m)
    ax.set_xlabel("n"); ax.set_ylabel("log10 M required (median, IQR)")
    ax.set_title("Shot requirement: global projector vs matched term-wise local estimator\nθ ~ U[-π,π]^n, k=0", fontsize=8)
    ax.grid(True, alpha=0.3); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)
    return path


def fig_deadzone_fraction(dz: pd.DataFrame, path: Path) -> Path:
    plt = _plt()
    qs = [0.5, 0.9, 0.99]
    fig, axes = plt.subplots(1, len(qs), figsize=(4 * len(qs), 4), sharey=True, squeeze=False)
    for ax, q in zip(axes[0], qs):
        for M in sorted(dz.shots.unique()):
            d = dz[dz.shots == M].sort_values("n")
            ax.plot(d.n, d[f"deadzone_frac_q{q}"], marker="o", label=f"M={M}")
            if f"full_deadzone_frac_q{q}" in d and d[f"full_deadzone_frac_q{q}"].notna().any():
                ax.plot(d.n, d[f"full_deadzone_frac_q{q}"], marker="x", linestyle="--", alpha=0.6, color=ax.lines[-1].get_color())
        ax.set_title(f"q = {q}   (solid: component k=0; dashed x: full gradient)", fontsize=8); ax.set_xlabel("n"); ax.set_ylim(-0.02, 1.02); ax.grid(True, alpha=0.3)
    axes[0][0].set_ylabel("initialization dead-zone fraction  P_θ[P(ĝ=0|θ,M) ≥ q]"); axes[0][0].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)
    return path


def fig_stage4_prediction(pred: pd.DataFrame, emp: pd.DataFrame, dist: pd.DataFrame, path: Path) -> Path:
    plt = _plt()
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    ax = axes[0]
    p = pred.groupby("shots").agg(pred_stuck=("p_stuck_all_iters", "mean"), pred_full_zero_iter1=("p_full_zero", "mean")).reset_index()
    ax.semilogx(p.shots, p.pred_stuck, "o-", label="predicted P(never moves in 2000 it.) = E_θ0[P_full^2000]")
    ax.semilogx(p.shots, p.pred_full_zero_iter1, "s--", label="predicted P(full ĝ = 0 at one iteration) = E_θ0[P_full]")
    e = emp.sort_values("shots")
    ax.errorbar(e.shots, e.empirical_frac_never_moved, yerr=np.vstack([e.empirical_frac_never_moved - e.ci_low, e.ci_high - e.empirical_frac_never_moved]), fmt="^", capsize=3, color="black", label="Stage 4 observed (500 runs/M, Wilson 95%)")
    ax.set_xlabel("M"); ax.set_ylabel("probability"); ax.set_ylim(-0.02, 1.02); ax.set_title("n=12, Stage 4 start seeds 1000–1099: predicted vs observed zero-gradient trajectories", fontsize=8)
    ax.grid(True, which="both", alpha=0.3); ax.legend(fontsize=7)
    ax = axes[1]
    for M in sorted(dist.shots.unique()):
        d = dist[dist.shots == M]
        ax.hist(d.p_zero, bins=40, range=(0, 1), histtype="step", label=f"M={M}")
    ax.set_xlabel("P(ĝ_0 = 0 | θ, M)"); ax.set_ylabel("count (random θ, n=12)"); ax.set_title("Distribution of component zero probability over random initializations, n=12", fontsize=8)
    ax.grid(True, alpha=0.3); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)
    return path


def fig_a_distribution(summary: pd.DataFrame, path: Path) -> Path:
    plt = _plt()
    s = summary.sort_values("n")
    fig, ax = plt.subplots(figsize=(6, 4.5))
    _band(ax, s.n, s.log10A_median, s.log10A_q25, s.log10A_q75, "empirical median log10 A (IQR band)")
    ax.plot(s.n, np.log10(s.geometric_A_theory), "--", color="gray", label="geometric scale 4^{-(n-1)} = exp E[log A]")
    ax.plot(s.n, np.log10(s.theory_mean_A), ":", color="black", label="arithmetic mean E[A] = 2^{-(n-1)}")
    ax.plot(s.n, np.log10(s.empirical_mean_A), "x", color="black", label="empirical mean A")
    ax.set_xlabel("n"); ax.set_ylabel("log10 A_k"); ax.set_title("A_k = ∏_{j≠k} cos²(θ_j/2) over θ ~ U[-π,π]^n: typical vs mean scale", fontsize=8)
    ax.grid(True, alpha=0.3); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)
    return path
