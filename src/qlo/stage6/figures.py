"""Stage 6 figures. Plain matplotlib, no decorative styling."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

COLOR = {"projector": "C0", "parity": "C1"}
LABEL = {"projector": "global projector (Stages 3-5)", "parity": "global parity (Stage 6)"}


def _plt():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def fig_two_failure_modes(ff: pd.DataFrame, path: Path, shots: int) -> Path:
    plt = _plt()
    d = ff[ff.shots == shots]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), sharex=True)
    for b in ("projector", "parity"):
        s = d[d.benchmark == b].sort_values("n")
        axes[0].plot(s.n, s.median_p_zero, "o-", color=COLOR[b], label=LABEL[b])
        axes[1].plot(s.n, s.median_p_correct, "o-", color=COLOR[b], label=LABEL[b])
    axes[0].axhline(1.0, color="black", lw=1, ls=":")
    axes[0].set_ylabel("median P(ĝ_k = 0 | θ, M)"); axes[0].set_title("count starvation: exact-zero probability", fontsize=9)
    axes[1].axhline(0.5, color="black", lw=1, ls=":", label="random direction")
    axes[1].set_ylabel("median P(correct sign | θ, M)"); axes[1].set_title("directional ambiguity: sign reliability", fontsize=9)
    for ax in axes:
        ax.set_xlabel("n"); ax.set_ylim(-0.02, 1.02); ax.grid(True, alpha=0.3); ax.legend(fontsize=7)
    fig.suptitle(f"Two finite-shot failure modes at fixed M = {shots} (θ ~ U[-π,π]^n, k=0)", fontsize=9)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)
    return path


def fig_zero_probability_comparison(ff: pd.DataFrame, path: Path) -> Path:
    plt = _plt()
    shots = sorted(ff.shots.unique())
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), sharey=True)
    for ax, b in zip(axes, ("projector", "parity")):
        for M in shots:
            s = ff[(ff.benchmark == b) & (ff.shots == M)].sort_values("n")
            ax.plot(s.n, s.median_p_zero, "o-", label=f"M={M}")
        ax.set_title(LABEL[b], fontsize=9); ax.set_xlabel("n"); ax.set_ylim(-0.02, 1.02); ax.grid(True, alpha=0.3)
    axes[0].set_ylabel("median P(ĝ_k = 0 | θ, M)"); axes[0].legend(fontsize=7)
    fig.suptitle("Exact-zero estimator probability vs n (median over θ ~ U[-π,π]^n, k=0); reference: zero-signal parity limit C(2M,M)/4^M ≈ 1/√(πM)", fontsize=9)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)
    return path


def fig_directional_correctness(ff: pd.DataFrame, path: Path) -> Path:
    plt = _plt()
    shots = sorted(ff.shots.unique())
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), sharey=True)
    for ax, b in zip(axes, ("projector", "parity")):
        for M in shots:
            s = ff[(ff.benchmark == b) & (ff.shots == M)].sort_values("n")
            ax.plot(s.n, s.median_p_correct, "o-", label=f"M={M}")
        ax.axhline(0.5, color="black", lw=1, ls=":")
        ax.set_title(LABEL[b], fontsize=9); ax.set_xlabel("n"); ax.grid(True, alpha=0.3)
    axes[0].set_ylabel("median P(correct sign | θ, M)"); axes[0].legend(fontsize=7)
    fig.suptitle("Directional correctness vs n (median over θ, k=0). Dotted 0.5 = coin-flip sign.\nP_correct + P_wrong + P_zero = 1: a low P_correct can come from zeros OR wrong signs (see P_zero figure).", fontsize=9)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)
    return path


def _band(ax, s, tag, label, color=None, marker="o"):
    ax.plot(s.n, s[f"{tag}_log10M_median"], marker=marker, color=color, label=label)
    ax.fill_between(s.n, s[f"{tag}_log10M_q25"], s[f"{tag}_log10M_q75"], alpha=0.15, color=color)


def fig_shots_for_direction(summary: pd.DataFrame, path: Path, fits: dict) -> Path:
    plt = _plt()
    fig, ax = plt.subplots(figsize=(6.5, 4.6))
    for b, ls in (("parity", "-"), ("projector", "--")):
        s = summary[summary.benchmark == b].sort_values("n")
        for tag, m in (("dir0.75", "o"), ("dir0.9", "s")):
            if f"{tag}_log10M_median" in s and s[f"{tag}_log10M_median"].notna().any():
                d = s.dropna(subset=[f"{tag}_log10M_median"])
                ax.plot(d.n, d[f"{tag}_log10M_median"], marker=m, linestyle=ls, color=COLOR[b],
                        label=f"{b}, P_correct ≥ {tag[3:]}")
                ax.fill_between(d.n, d[f"{tag}_log10M_q25"], d[f"{tag}_log10M_q75"], alpha=0.12, color=COLOR[b])
    n0 = summary.n.min()
    ax.plot(sorted(summary.n.unique()), [1.0 + np.log10(4) * (n - n0) for n in sorted(summary.n.unique())], ":", color="black", label="slope log₁₀4 reference")
    ax.set_xlabel("n"); ax.set_ylabel("log₁₀ M for directional correctness (median, IQR)")
    ax.set_title("Shots for a reliable gradient SIGN (θ ~ U[-π,π]^n, k=0)\nparity: exact for |sB| ≥ 10^-2.5, normal beyond\nprojector: exact bisection (Skellam limit for |s| ≲ 0.02)", fontsize=8)
    ax.grid(True, alpha=0.3); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)
    return path


def fig_snr_shot_scaling(summary: pd.DataFrame, path: Path, fits: dict) -> Path:
    plt = _plt()
    fig, ax = plt.subplots(figsize=(6.5, 4.6))
    for b in ("projector", "parity"):
        s = summary[summary.benchmark == b].sort_values("n")
        for tag, m in (("snr1", "o"), ("snr2", "s")):
            _band(ax, s, tag, f"{b}, SNR ≥ {tag[3:]}", COLOR[b], m)
            f = fits.get(f"{b}_{tag}")
            if f:
                ax.plot(s.n, f["intercept"] + f["slope"] * s.n, ":", color=COLOR[b], alpha=0.6)
    ax.set_xlabel("n"); ax.set_ylabel("log₁₀ M for component SNR (median, IQR)")
    ax.set_title("Shots for component SNR: both benchmarks grow exponentially\n(dotted: fitted median slope)", fontsize=8)
    ax.grid(True, alpha=0.3); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)
    return path


def fig_matched_gradient(ms: pd.DataFrame, path: Path) -> Path:
    plt = _plt()
    shots = sorted(ms.shots.unique())
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
    M0 = shots[len(shots) // 2]
    d = ms[ms.shots == M0]
    for b in ("projector", "parity"):
        s = d[d.benchmark == b].sort_values("log10_abs_g_mid")
        axes[0].plot(s.log10_abs_g_mid, s.mean_p_zero, "o-", color=COLOR[b], label=LABEL[b])
        axes[1].plot(s.log10_abs_g_mid, s.mean_p_correct, "o-", color=COLOR[b], label=LABEL[b])
        axes[2].semilogy(s.log10_abs_g_mid, s.mean_p_pos, "o-", color=COLOR[b], label=f"{b}: p_pos")
        axes[2].semilogy(s.log10_abs_g_mid, s.mean_p_neg, "s--", color=COLOR[b], alpha=0.6, label=f"{b}: p_neg")
    axes[0].set_ylabel(f"mean P(ĝ=0) at M={M0}"); axes[0].set_title("exact-zero probability at matched |g|", fontsize=9)
    axes[1].axhline(0.5, color="black", lw=1, ls=":"); axes[1].set_ylabel(f"mean P(correct sign) at M={M0}")
    axes[1].set_title("directional correctness at matched |g|", fontsize=9)
    axes[2].set_ylabel("shifted Bernoulli probabilities"); axes[2].set_title("why: where the two Bernoullis sit", fontsize=9)
    for ax in axes:
        ax.set_xlabel("log₁₀ |g_true| (shared bins)"); ax.grid(True, alpha=0.3); ax.legend(fontsize=6)
    fig.suptitle("Matched-signal comparison: bins of equal |g_true| pooled over n, identical M", fontsize=9)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)
    return path


def fig_vector_alignment(vr: pd.DataFrame, path: Path) -> Path:
    plt = _plt()
    shots = sorted(vr.shots.unique())
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), sharex=True)
    M0 = shots[len(shots) // 2]
    d = vr[vr.shots == M0]
    for b in ("projector", "parity"):
        s = d[d.benchmark == b].sort_values("n")
        axes[0].plot(s.n, s.p_vector_zero, "o-", color=COLOR[b], label=LABEL[b])
        axes[1].plot(s.n, s.median_cos, "o-", color=COLOR[b], label=f"{b} (all draws; zero step = 0)")
        axes[1].plot(s.n, s.median_cos_given_nonzero, "s--", color=COLOR[b], alpha=0.6, label=f"{b} (nonzero ĝ only)")
        axes[2].plot(s.n, s.p_dot_gt_0, "o-", color=COLOR[b], label=f"{b} (all draws; zero step = 0)")
        axes[2].plot(s.n, s.p_dot_gt_0_given_nonzero, "s--", color=COLOR[b], alpha=0.6, label=f"{b} (nonzero ĝ only)")
    axes[0].set_ylabel("P(ĝ vector = 0)"); axes[0].set_title("no vector at all", fontsize=9)
    axes[1].axhline(0.0, color="black", lw=1, ls=":"); axes[1].set_ylabel("median cos(ĝ, g)"); axes[1].set_title("alignment with the exact gradient", fontsize=9)
    axes[2].axhline(0.5, color="black", lw=1, ls=":"); axes[2].set_ylabel("P(ĝ·g > 0)"); axes[2].set_title("is it a descent direction?", fontsize=9)
    for ax in axes:
        ax.set_xlabel("n"); ax.grid(True, alpha=0.3); ax.legend(fontsize=6)
    fig.suptitle(f"Full-vector reliability at M = {M0} (one independent batch per component and shift)", fontsize=9)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)
    return path


def fig_mechanism_map(ff: pd.DataFrame, path: Path) -> Path:
    plt = _plt()
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), sharey=True)
    ns = sorted(ff.n.unique())
    shots = sorted(ff.shots.unique())
    for ax, b in zip(axes, ("projector", "parity")):
        grid = np.zeros((len(ns), len(shots)))
        for i, n in enumerate(ns):
            for j, M in enumerate(shots):
                r = ff[(ff.benchmark == b) & (ff.n == n) & (ff.shots == M)].iloc[0]
                # 0 = resolved, 1 = directional ambiguity only, 2 = count starvation (dead zone)
                grid[i, j] = 2 if r.frac_p_zero_ge >= 0.5 else (1 if r.frac_p_correct_le >= 0.5 else 0)
        im = ax.imshow(grid, origin="lower", aspect="auto", vmin=0, vmax=2, cmap="viridis",
                       extent=[-0.5, len(shots) - 0.5, -0.5, len(ns) - 0.5])
        ax.set_xticks(range(len(shots))); ax.set_xticklabels(shots, fontsize=7)
        ax.set_yticks(range(len(ns))); ax.set_yticklabels(ns, fontsize=7)
        ax.set_xlabel("M"); ax.set_title(LABEL[b], fontsize=9)
        for i in range(len(ns)):
            for j in range(len(shots)):
                ax.text(j, i, f"{int(grid[i, j])}", ha="center", va="center", color="white", fontsize=7)
    axes[0].set_ylabel("n")
    cb = fig.colorbar(im, ax=axes, ticks=[0, 1, 2], fraction=0.03)
    cb.ax.set_yticklabels(["0 resolved", "1 directional\nambiguity", "2 count\nstarvation"], fontsize=7)
    fig.suptitle("Dominant failure region by (n, M): majority of random initializations\n"
                 "2 if ≥50% have P(ĝ=0) ≥ 0.9; else 1 if ≥50% have P_correct ≤ 0.60; else 0", fontsize=9, y=1.04)
    fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
    return path
