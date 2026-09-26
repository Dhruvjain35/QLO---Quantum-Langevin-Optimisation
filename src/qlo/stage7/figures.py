"""Stage 7 figures. Plain matplotlib, no decorative styling."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binom

from qlo.stage5.theory import a_k, s_k
from qlo.stage7 import estimators as E
from qlo.stage7.analysis import sample_theta

COLOR = {"loschmidt": "C0", "swap": "C3"}
LABEL = {"loschmidt": "Loschmidt / projector", "swap": "SWAP test"}


def _plt():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def _estimator_pmf(scheme: str, A: float, s: float, M: int):
    """Exact pmf of g_hat = c (K_pos - K_neg)/M by convolving the two binomial pmfs."""
    p_pos, p_neg, c = E.pair(scheme, A, s)
    r = np.arange(M + 1)
    pmf = np.convolve(binom.pmf(r, M, p_pos), binom.pmf(r, M, p_neg)[::-1])  # index i <-> D = i - M
    D = np.arange(-M, M + 1)
    return c * D / M, pmf


def fig_same_landscape(path: Path, M: int = 1024) -> Path:
    plt = _plt()
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    for ax, n in zip(axes, (4, 8, 12)):
        pool = sample_theta(n, 400, 0, 90)
        g_all = a_k(pool, 0) * s_k(pool, 0) / 2
        i = int(np.argmin(np.abs(np.log(np.abs(g_all)) - np.median(np.log(np.abs(g_all))))))
        A, s = float(a_k(pool[i], 0)), float(s_k(pool[i], 0))
        g = A * s / 2
        for sc in E.SCHEMES:
            x, pmf = _estimator_pmf(sc, A, s, M)
            keep = pmf > 1e-6 * pmf.max()
            ax.vlines(x[keep], 0, pmf[keep], color=COLOR[sc], alpha=0.7, lw=1.2, label=f"{LABEL[sc]} (P(ĝ=0)={pmf[x == 0].sum():.3f})")
            if sc == "loschmidt":  # its mass is concentrated on a few atoms (mostly exactly 0): mark them
                ax.plot(x[keep], pmf[keep], "o", color=COLOR[sc], ms=5, zorder=5)
        ax.axvline(g, color="black", ls=":", lw=1.2, label=f"exact g = {g:.2e}")
        ax.set_title(f"n = {n}, median-|g| θ, M = {M}", fontsize=9)
        ax.set_xlabel("ĝ_k (same θ, same exact gradient)"); ax.set_yscale("log"); ax.set_ylim(1e-6, 1.2); ax.grid(True, alpha=0.3)
        ax.legend(fontsize=6)
    axes[0].set_ylabel("exact probability")
    fig.suptitle("One landscape, two measurement schemes: exact finite-shot parameter-shift gradient distributions", fontsize=9)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)
    return path


def fig_zero_probability(pf: pd.DataFrame, path: Path) -> Path:
    plt = _plt()
    fig, ax = plt.subplots(figsize=(6.8, 4.4))
    for M, ls in zip(sorted(pf.shots.unique()), ("-", "--", ":")):
        d = pf[pf.shots == M].sort_values("n")
        for sc in E.SCHEMES:
            ax.plot(d.n, d[f"{sc}_median_p_zero"], marker="o", ls=ls, color=COLOR[sc], label=f"{LABEL[sc]}, M={M}")
        ax.axhline(1 / np.sqrt(np.pi * M), color="grey", lw=0.8, ls=ls)
    ax.set_xlabel("n"); ax.set_ylabel("median P(ĝ_k = 0 | θ, M)"); ax.set_ylim(-0.02, 1.02); ax.grid(True, alpha=0.3); ax.legend(fontsize=6)
    ax.set_title("Exact-zero probability on the SAME landscape (grey: 1/√(πM))", fontsize=9)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)
    return path


def fig_directional(pf: pd.DataFrame, path: Path) -> Path:
    plt = _plt()
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4), sharey=True)
    for M, ls in zip(sorted(pf.shots.unique()), ("-", "--", ":")):
        d = pf[pf.shots == M].sort_values("n")
        for sc in E.SCHEMES:
            axes[0].plot(d.n, d[f"{sc}_median_p_correct"], marker="o", ls=ls, color=COLOR[sc], label=f"{LABEL[sc]}, M={M}")
            axes[1].plot(d.n, d[f"{sc}_median_p_correct_given_nonzero"], marker="o", ls=ls, color=COLOR[sc], label=f"{LABEL[sc]}, M={M}")
    for ax, t in zip(axes, ("unconditional P(correct sign)", "P(correct sign | ĝ ≠ 0)")):
        ax.axhline(0.5, color="black", lw=1, ls=":"); ax.set_xlabel("n"); ax.set_title(t, fontsize=9); ax.grid(True, alpha=0.3)
    axes[0].set_ylabel("median over θ"); axes[1].legend(fontsize=6)
    fig.suptitle("Directional correctness, same θ for both schemes", fontsize=9)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)
    return path


def fig_required(summary: pd.DataFrame, tags, path: Path, what: str) -> Path:
    plt = _plt()
    fig, ax = plt.subplots(figsize=(7, 4.8))
    markers = ("o", "s", "^")
    for tag, m in zip(tags, markers):
        for sc in E.SCHEMES:
            d = summary[summary.scheme == sc].sort_values("n")
            ax.plot(d.n, d[f"{tag}_log10M_median"], marker=m, color=COLOR[sc], label=f"{LABEL[sc]}: {tag}")
            ax.fill_between(d.n, d[f"{tag}_log10M_q25"], d[f"{tag}_log10M_q75"], color=COLOR[sc], alpha=0.10)
    ax.set_xlabel("n"); ax.set_ylabel("log₁₀ M required (median, IQR)"); ax.grid(True, alpha=0.3); ax.legend(fontsize=7)
    ax.set_title(f"Required shots per component for {what} (θ ~ U[-π,π]^n, k = 0)", fontsize=9)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)
    return path


def fig_scaling(summary: pd.DataFrame, fits: dict, path: Path) -> Path:
    plt = _plt()
    fig, ax = plt.subplots(figsize=(7, 4.8))
    ns = np.array(sorted(summary.n.unique()))
    for sc in E.SCHEMES:
        d = summary[summary.scheme == sc].sort_values("n")
        for tag, m in (("snr1", "o"), ("dir0.9", "s")):
            ax.plot(d.n, d[f"{tag}_log10M_median"], m, color=COLOR[sc], label=f"{LABEL[sc]} {tag} (slope {fits[f'{sc}_{tag}']['slope']:.3f})")
            f = fits[f"{sc}_{tag}"]
            ax.plot(ns, f["intercept"] + f["slope"] * ns, "-", color=COLOR[sc], alpha=0.5)
    for ref, val, ls in (("4ⁿ", np.log10(4), ":"), ("16ⁿ", np.log10(16), "--")):
        ax.plot(ns, val * (ns - ns[0]), ls, color="black", lw=1, label=f"reference slope {ref} ({val:.3f})")
    ax.set_xlabel("n"); ax.set_ylabel("median log₁₀ M"); ax.grid(True, alpha=0.3); ax.legend(fontsize=6)
    ax.set_title("Same landscape, fitted median shot scaling vs predeclared references", fontsize=9)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)
    return path


def fig_vector(vr: pd.DataFrame, path: Path, M: int = 1024) -> Path:
    plt = _plt()
    d = vr[vr.shots == M]
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    for sc in E.SCHEMES:
        s = d[d.scheme == sc].sort_values("n")
        axes[0].plot(s.n, s.p_vector_zero, "o-", color=COLOR[sc], label=LABEL[sc])
        axes[1].plot(s.n, s.median_cos, "o-", color=COLOR[sc], label=f"{LABEL[sc]} (all; zero step = 0)")
        axes[1].fill_between(s.n, s.cos_q10, s.cos_q90, color=COLOR[sc], alpha=0.12)
        axes[1].plot(s.n, s.median_cos_given_nonzero, "s--", color=COLOR[sc], alpha=0.6, label=f"{LABEL[sc]} (nonzero only)")
        axes[2].plot(s.n, s.p_dot_gt_0, "o-", color=COLOR[sc], label=LABEL[sc])
    axes[0].set_ylabel("P(ĝ vector = 0)"); axes[1].set_ylabel("cos(ĝ, g): median, 10–90% band"); axes[2].set_ylabel("P(ĝ·g > 0)")
    axes[1].axhline(0, color="black", lw=1, ls=":"); axes[2].axhline(0.5, color="black", lw=1, ls=":")
    for ax in axes:
        ax.set_xlabel("n"); ax.grid(True, alpha=0.3); ax.legend(fontsize=6)
    fig.suptitle(f"Full-gradient reliability at M = {M}, same θ for both schemes", fontsize=9)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)
    return path


def fig_trajectories(curves: pd.DataFrame, path: Path) -> Path:
    plt = _plt()
    ns = sorted(curves.n.unique())
    fig, axes = plt.subplots(1, len(ns), figsize=(4.2 * len(ns), 4.2), sharey=False)
    style = {"exact": ("black", "-")}
    for ax, n in zip(np.atleast_1d(axes), ns):
        for label, c in curves[curves.n == n].groupby("variant", sort=False):
            col = "black" if label == "exact" else (COLOR["loschmidt"] if label.startswith("LE") else ("grey" if "random walk" in label else COLOR["swap"]))
            ls = "-" if label == "exact" else ("--" if "M=64" in label or "M=32" in label else "-.")
            if "matched" in label:
                ls = ":"
            ax.plot(c.iter, c.median_log10_F, ls, color=col, label=label)
        ax.set_title(f"n = {n}", fontsize=9); ax.set_xlabel("iteration"); ax.grid(True, alpha=0.3)
    np.atleast_1d(axes)[0].set_ylabel("median log₁₀ F (exact fidelity)")
    np.atleast_1d(axes)[-1].legend(fontsize=6)
    fig.suptitle("Same starts, same η, same landscape: exact GD vs Loschmidt vs SWAP finite-shot GD (diagnostic only)", fontsize=9)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)
    return path


def fig_information(info: pd.DataFrame, path: Path) -> Path:
    plt = _plt()
    fig, ax = plt.subplots(figsize=(7, 4.6))
    for sc in E.SCHEMES:
        d = info[info.scheme == sc].sort_values("n")
        ax.plot(d.n, d.median_log10_hellinger2, "o-", color=COLOR[sc], label=f"{LABEL[sc]}: squared Hellinger")
        ax.plot(d.n, 2 * d.median_log10_tv, "s--", color=COLOR[sc], alpha=0.6, label=f"{LABEL[sc]}: TV²")
    ax.set_xlabel("n"); ax.set_ylabel("median log₁₀ (per-shot distance)"); ax.grid(True, alpha=0.3); ax.legend(fontsize=7)
    ax.set_title("Distinguishability of the + and − shifted single-shot distributions", fontsize=9)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)
    return path
