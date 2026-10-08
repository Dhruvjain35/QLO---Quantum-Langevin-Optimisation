"""B2 figures (plain matplotlib, read from the result CSVs)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

LABEL = {"loschmidt": "Loschmidt", "swap_ancilla": "ancilla SWAP", "swap_destructive": "destructive SWAP", "hadamard_u": "Hadamard U-stat"}
STYLE = {"loschmidt": ("C0", "o", "-"), "swap_ancilla": ("C3", "s", "-"), "swap_destructive": ("k", "x", ":"), "hadamard_u": ("C2", "^", "-")}


def _plt():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def _save(fig, path, plt):
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def all_figures(out: Path, cfg) -> None:
    plt = _plt()
    figs = out / "figures"
    med = pd.read_csv(out / "required_shots_rx.csv")
    m0 = med[(med.seed == cfg.stage7_seed) & (med.rho == cfg.rho_primary)]
    summ = pd.read_csv(out / "slope_summary.csv")
    fb = pd.read_csv(out / "fixed_budget_failure_modes.csv")
    est = list(cfg.estimators)

    # 1 required shots
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    for e in est:
        g = m0[m0.estimator == e].sort_values("n")
        c, mk, ls = STYLE[e]
        ax.plot(g.n, g.median_log10_total_executions, ls, marker=mk, color=c, label=LABEL[e])
    ax.set_xlabel("n")
    ax.set_ylabel("median log10 total circuit executions (SNR 1)")
    ax.set_title("RX-product benchmark: executions per gradient component", fontsize=10)
    ax.legend(fontsize=8)
    _save(fig, figs / "four_estimator_required_shots.png", plt)

    # 2 slopes with replicate CIs
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for ax, fr in zip(axes, ("primary_2_20", "n_ge_8")):
        for j, q in enumerate(("median_log10_total_executions", "median_log10_M_internal")):
            s = summ[(summ.fit_range == fr) & (summ.quantity == q)].set_index("estimator").loc[est]
            x = np.arange(len(est)) + (j - 0.5) * 0.25
            ax.errorbar(x, s.replicate_mean, yerr=[s.replicate_mean - s.replicate_pct_lo, s.replicate_pct_hi - s.replicate_mean],
                        fmt="o", capsize=3, label=("total executions" if j == 0 else "internal M") + " (10-seed mean, 95% pct)")
            ax.plot(x, s.stage7_population_slope, "kx")
        ax.axhline(np.log10(4), color="0.5", ls=":", lw=0.8)
        ax.axhline(np.log10(16), color="0.5", ls="--", lw=0.8)
        ax.set_xticks(np.arange(len(est)))
        ax.set_xticklabels([LABEL[e] for e in est], fontsize=7)
        ax.set_title(f"slope of median log10 required ({fr}); x = Stage 7 population", fontsize=9)
    axes[0].legend(fontsize=7)
    _save(fig, figs / "four_estimator_slopes.png", plt)

    # 3 ancilla vs destructive SWAP
    sub = pd.read_csv(out / "required_shots_rx_theta_subsample.csv.gz")
    ctrl = pd.read_csv(out / "swap_destructive_vs_ancilla_rx.csv")
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].plot(sub.log10_total_swap_ancilla, sub.log10_total_swap_destructive, ".", ms=1, alpha=0.3)
    lo, hi = sub.log10_total_swap_ancilla.min(), sub.log10_total_swap_ancilla.max()
    axes[0].plot([lo, hi], [lo, hi], "k-", lw=0.8)
    axes[0].set_xlabel("ancilla SWAP log10 total executions (per theta)")
    axes[0].set_ylabel("destructive SWAP")
    axes[1].semilogy(ctrl.n, np.maximum(ctrl.median_abs_dlog10, 1e-17), "o-", label="median |d log10|")
    axes[1].semilogy(ctrl.n, np.maximum(ctrl.q999_abs_dlog10, 1e-17), "s-", label="99.9% quantile")
    axes[1].semilogy(ctrl.n, np.maximum(ctrl.max_abs_dlog10, 1e-17), "^-", label="max (tail: one pair overlap ~1e-19)")
    axes[1].set_xlabel("n")
    axes[1].legend(fontsize=7)
    axes[1].set_title("per-theta difference (float precision only)", fontsize=9)
    _save(fig, figs / "ancilla_vs_destructive_swap.png", plt)

    # 4 Hadamard variance validation
    hv = pd.read_csv(out / "hadamard_variance_validation.csv").dropna(subset=["var_formula"])
    hv = hv[hv.gradient_case.isna()] if "gradient_case" in hv else hv
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    for (x, y), g in hv.groupby(["x", "y"]):
        line, = ax.loglog(g.M, g.var_formula, "-", lw=0.8, label=f"x={x}, y={y}")
        ax.loglog(g.M, g.mc_var, "o", color=line.get_color(), mfc="none")
        ex = g.dropna(subset=["exact_var"])
        ax.loglog(ex.M, ex.exact_var, "x", color=line.get_color())
    ax.loglog([2, 256], [4 / (2 * 1), 4 / (256 * 255)], "k:", lw=0.8, label="4/[M(M-1)] (x=y=0)")
    ax.set_xlabel("M (shots per quadrature)")
    ax.set_ylabel("Var F_hat_HT")
    ax.set_title("line = formula, x = exact enumeration, o = Monte Carlo", fontsize=9)
    ax.legend(fontsize=6)
    _save(fig, figs / "hadamard_variance_validation.png", plt)

    # 5 Hadamard required shots vs Loschmidt
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    h, l = m0[m0.estimator == "hadamard_u"].sort_values("n"), m0[m0.estimator == "loschmidt"].sort_values("n")
    axes[0].plot(h.n, h.median_log10_M_internal, "^-", color="C2", label="Hadamard internal M (per quadrature per shift)")
    axes[0].plot(h.n, h.median_log10_total_executions, "^--", color="C2", label="Hadamard total = 4M")
    axes[0].plot(l.n, l.median_log10_total_executions, "o-", color="C0", label="Loschmidt total = 2M")
    axes[0].set_xlabel("n")
    axes[0].set_ylabel("median log10")
    axes[0].legend(fontsize=7)
    axes[1].plot(h.n, 10 ** (h.median_log10_total_executions.values - l.median_log10_total_executions.values), "o-")
    axes[1].axhspan(8, 2 * 2 * (1 + np.sqrt(3)), color="0.9", label="deep-limit prediction 8 to 10.9")
    axes[1].set_xlabel("n")
    axes[1].set_ylabel("Hadamard / Loschmidt median total executions")
    axes[1].legend(fontsize=7)
    _save(fig, figs / "hadamard_required_shots.png", plt)

    # 6-8 fixed budgets
    for fname, col, ylab in (("fixed_budget_sign_accuracy.png", "median_P_correct", "median P(correct sign)"),
                             ("fixed_budget_snr.png", "median_SNR", "median SNR")):
        fig, axes = plt.subplots(1, len(cfg.fixed_n), figsize=(3.2 * len(cfg.fixed_n), 3.6), sharey=True)
        for ax, n in zip(axes, cfg.fixed_n):
            for e in est:
                g = fb[(fb.n == n) & (fb.estimator == e)]
                c, mk, ls = STYLE[e]
                ax.semilogx(g.total_executions, g[col], ls, marker=mk, color=c, label=LABEL[e])
            if col == "median_SNR":
                ax.set_yscale("log")
            ax.set_title(f"n = {n}", fontsize=9)
            ax.set_xlabel("total executions")
        axes[0].set_ylabel(ylab)
        axes[0].legend(fontsize=6)
        _save(fig, figs / fname, plt)
    fig, axes = plt.subplots(2, len(cfg.fixed_n), figsize=(3.2 * len(cfg.fixed_n), 6), sharey="row")
    for j, n in enumerate(cfg.fixed_n):
        for e in est:
            g = fb[(fb.n == n) & (fb.estimator == e)]
            c, mk, ls = STYLE[e]
            axes[0, j].semilogx(g.total_executions, g.median_P_zero, ls, marker=mk, color=c, label=LABEL[e])
            axes[1, j].semilogx(g.total_executions, g.median_P_correct_given_nonzero, ls, marker=mk, color=c)
        axes[0, j].set_title(f"n = {n}", fontsize=9)
        axes[1, j].set_xlabel("total executions")
    axes[0, 0].set_ylabel("median P(g_hat = 0)")
    axes[1, 0].set_ylabel("median P(correct | nonzero)")
    axes[0, 0].legend(fontsize=6)
    _save(fig, figs / "failure_mode_comparison.png", plt)

    # 9 resource-normalized
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    for e in est:
        g = m0[m0.estimator == e].sort_values("n")
        c, mk, ls = STYLE[e]
        ax.plot(g.n, g.median_log10_total_executions.values - l.median_log10_total_executions.values, ls, marker=mk, color=c, label=LABEL[e])
    ax.set_xlabel("n")
    ax.set_ylabel("log10 (median total executions / Loschmidt)")
    ax.set_title("Executions relative to Loschmidt (circuit width, copies and control cost NOT included)", fontsize=8)
    ax.legend(fontsize=7)
    _save(fig, figs / "resource_normalized_comparison.png", plt)

    # 10 Hadamard distribution examples
    de = pd.read_csv(out / "hadamard_distribution_examples.csv")
    fig, axes = plt.subplots(1, de.total_executions.nunique(), figsize=(12, 3.6))
    for ax, (T, g) in zip(axes, de.groupby("total_executions")):
        ax.bar(g.g_hat, g.prob, width=(g.g_hat.max() - g.g_hat.min()) / 150 + 1e-12)
        ax.axvline(g.g_true.iloc[0], color="r", lw=1, label="true g")
        ax.axvline(0, color="k", lw=0.5)
        ax.set_title(f"T = {T} executions (M = {T // 4}), n = 10 median-S theta", fontsize=8)
        ax.set_xlabel("g_hat (exact pmf)")
    axes[0].legend(fontsize=7)
    _save(fig, figs / "hadamard_distribution_examples.png", plt)
