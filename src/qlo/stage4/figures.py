"""Stage 4 figures. Plain matplotlib, no decorative styling."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from qlo.stage4.trajectory import FIXED_TARGETS

LABELS = {"exact_gd": "A exact GD", "shot_sgd": "B finite-shot SGD", "matched_gaussian": "C matched Gaussian",
          "constant_langevin": "D constant Langevin", "noise_only": "E noise-only (diagnostic)"}
ORDER = ["exact_gd", "shot_sgd", "matched_gaussian", "constant_langevin", "noise_only"]


def _plt():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def fig_snr(snr: pd.DataFrame, path: Path) -> Path:
    plt = _plt()
    fig, ax = plt.subplots(figsize=(6, 4.5))
    for M in sorted(snr.shots.unique()):
        s = snr[snr.shots == M].sort_values("n")
        ax.semilogy(s.n, s.nsr_median, "o-", label=f"M={M} (median)")
        ax.fill_between(s.n, s.nsr_q25, s.nsr_q75, alpha=0.15)
    ax.axhline(1, color="black", linewidth=1, label="noise = signal")
    ax.set_xlabel("n"); ax.set_ylabel("sqrt(trace Σ_shot) / ||g_exact||")
    ax.set_title("Initial shot-noise-to-gradient ratio vs n (median, IQR band)\n5000 θ ~ U[-π,π]^n per n; analytic conditional shot variance", fontsize=8)
    ax.grid(True, which="both", alpha=0.3); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)
    return path


def fig_trajectories(traj: pd.DataFrame, path: Path, primary_config: dict) -> Path:
    plt = _plt()
    ns = sorted(traj.n.unique()); seeds = sorted(traj.start_seed.unique())
    fig, axes = plt.subplots(len(ns), len(seeds), figsize=(3.6 * len(seeds), 2.6 * len(ns)), squeeze=False, sharex=True)
    for i, n in enumerate(ns):
        for j, s in enumerate(seeds):
            ax = axes[i][j]
            d = traj[(traj.n == n) & (traj.start_seed == s)]
            for m in ORDER:
                dm = d[d.method == m].sort_values("iter")
                if len(dm):
                    ax.semilogy(dm["iter"], dm["F"], linewidth=1, label=LABELS[m])
            ax.set_title(f"n={n}, start seed {s}, replicate 0", fontsize=8); ax.grid(True, which="both", alpha=0.3)
            if i == len(ns) - 1: ax.set_xlabel("iteration")
            if j == 0: ax.set_ylabel("exact fidelity F (offline)")
    axes[0][0].legend(fontsize=6)
    fig.suptitle(f"Representative trajectories (predetermined seeds), selected hyperparameters {primary_config}", fontsize=8)
    fig.tight_layout(); fig.savefig(path, dpi=130); plt.close(fig)
    return path


def fig_gain(runs: pd.DataFrame, path: Path) -> Path:
    plt = _plt()
    ns = sorted(runs.n.unique())
    fig, axes = plt.subplots(1, len(ns), figsize=(3.2 * len(ns), 4), sharey=True, squeeze=False)
    for ax, n in zip(axes[0], ns):
        data = [runs[(runs.n == n) & (runs.method == m)]["log10_gain_best"].to_numpy() for m in ORDER]
        ax.boxplot(data, tick_labels=["A", "B", "C", "D", "E"], showfliers=True, flierprops={"markersize": 2})
        ax.axhline(0, color="black", linewidth=1); ax.set_title(f"n={n}", fontsize=9); ax.grid(True, alpha=0.3)
    axes[0][0].set_ylabel("log10(F_best / F0)")
    fig.suptitle("Best-fidelity gain by method and n (held-out seeds 1000-1099; B/C/D/E: 5 replicates each)\n" + ", ".join(f"{k}={v}" for k, v in LABELS.items()), fontsize=7)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)
    return path


def fig_target_success(summary: pd.DataFrame, path: Path) -> Path:
    plt = _plt()
    fig, axes = plt.subplots(1, len(FIXED_TARGETS), figsize=(4 * len(FIXED_TARGETS), 4), sharey=True, squeeze=False)
    for ax, t in zip(axes[0], FIXED_TARGETS):
        for m in ORDER:
            s = summary[summary.method == m].sort_values("n")
            if len(s):
                yerr = np.clip(np.vstack([s[f"p_hit_F{t}"] - s[f"p_hit_F{t}_ci_low"], s[f"p_hit_F{t}_ci_high"] - s[f"p_hit_F{t}"]]), 0, None)
                ax.errorbar(s.n, s[f"p_hit_F{t}"], yerr=yerr, marker="o", capsize=3, label=LABELS[m])  # Wilson interval can exclude p at p in {0,1}
        ax.set_title(f"target F ≥ {t}", fontsize=9); ax.set_xlabel("n"); ax.set_ylim(-0.02, 1.02); ax.grid(True, alpha=0.3)
    axes[0][0].set_ylabel("P(reach target within 2000 iterations) ± Wilson 95% CI"); axes[0][0].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)
    return path


def fig_shot_cost(runs: pd.DataFrame, path: Path) -> Path:
    plt = _plt()
    b = runs[runs.method == "shot_sgd"]
    fig, axes = plt.subplots(1, len(FIXED_TARGETS), figsize=(4 * len(FIXED_TARGETS), 4), squeeze=False)
    for ax, t in zip(axes[0], FIXED_TARGETS):
        for n in sorted(b.n.unique()):
            v = b[(b.n == n) & b[f"eligible_F{t}"]][f"shots_to_F{t}"].dropna().to_numpy()
            if len(v):
                ax.hist(np.log10(v), bins=20, alpha=0.5, label=f"n={n} ({len(v)} hits)")
        ax.set_title(f"finite-shot SGD: shots to reach F ≥ {t}", fontsize=9); ax.set_xlabel("log10(cumulative shots)"); ax.grid(True, alpha=0.3); ax.legend(fontsize=7)
    axes[0][0].set_ylabel("count of successful runs")
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)
    return path


def fig_shot_budget(runs: pd.DataFrame, path: Path, n_iters: int) -> Path:
    plt = _plt()
    b = runs[runs.method == "shot_sgd"]
    fig, ax = plt.subplots(figsize=(6, 4.5))
    for n in sorted(b.n.unique()):
        s = b[b.n == n].groupby("shots").agg(gain=("log10_gain_best", "median"), q25=("log10_gain_best", lambda x: np.percentile(x, 25)),
                                             q75=("log10_gain_best", lambda x: np.percentile(x, 75))).reset_index()
        total = 2 * n * s.shots * n_iters
        ax.errorbar(total, s.gain, yerr=np.vstack([s.gain - s.q25, s.q75 - s.gain]), marker="o", capsize=3, label=f"n={n}")
        for M, x, y in zip(s.shots, total, s.gain):
            ax.annotate(f"M={M}", (x, y), fontsize=6, xytext=(3, 3), textcoords="offset points")
    ax.set_xscale("log"); ax.set_xlabel(f"total shots consumed over {n_iters} iterations (2 n M per iteration)")
    ax.set_ylabel("median log10(F_best/F0) (IQR bars)"); ax.axhline(0, color="black", linewidth=1)
    ax.set_title("Finite-shot SGD: improvement vs total shot consumption (held-out seeds, selected η, all M)", fontsize=8)
    ax.grid(True, which="both", alpha=0.3); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)
    return path
