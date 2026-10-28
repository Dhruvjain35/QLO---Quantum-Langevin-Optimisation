"""B1 figures. Plain matplotlib, primary fit range, every seed drawn as a small point next to the seed mean and its
bootstrap 95% CI (B3 utilities). Reference lines are guides, never fits."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

FAM_SHORT = {"rx_product": "RX", "hea_ring": "HEA", "rxry_czbrick": "BRICK"}
POS_COLOR = {"EARLY": "C0", "MIDDLE": "C2", "LATE": "C3", "k=0": "k"}
FAM_MARK = {"rx_product": "s", "hea_ring": "o", "rxry_czbrick": "^"}


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


def cell_label(r) -> str:
    if r["family"] == "rx_product":
        return "RX product"
    return f"{FAM_SHORT[r['family']]} {r['regime']} {r['position']}"


def _primary(main: pd.DataFrame, fits: pd.DataFrame):
    m = main[main.fit_range == "primary"].reset_index(drop=True)
    f = fits[fits.fit_range == "primary"]
    return m, f


def _seed_values(f: pd.DataFrame, r, col: str) -> np.ndarray:
    g = f[(f.family == r["family"]) & (f.regime == r["regime"]) & (f.position == r["position"])]
    return g[col].values


def _cell_errorbar(ax, m: pd.DataFrame, f: pd.DataFrame, col: str, refs=(), title: str = "", xlabel: str = "") -> None:
    y = np.arange(len(m))[::-1]
    for yi, (_, r) in zip(y, m.iterrows()):
        c = POS_COLOR[r["position"]]
        ax.plot(_seed_values(f, r, col), np.full(len(_seed_values(f, r, col)), yi), ".", color=c, alpha=0.35, ms=4)
        ax.errorbar(r[col], yi, xerr=[[r[col] - r[f"{col}_ci_lo"]], [r[f"{col}_ci_hi"] - r[col]]], fmt=FAM_MARK[r["family"]], color=c,
                    capsize=3, ms=5)
    for v, lab, ls in refs:
        ax.axvline(v, color="0.5", ls=ls, lw=0.8, label=lab)
    ax.set_yticks(y)
    ax.set_yticklabels([cell_label(r) for _, r in m.iterrows()], fontsize=7)
    ax.set_title(title, fontsize=10)
    ax.set_xlabel(xlabel)
    if refs:
        ax.legend(fontsize=7, loc="best")


def bS_by_family(main, fits, path: Path) -> Path:
    plt = _plt()
    m, f = _primary(main, fits)
    fig, ax = plt.subplots(figsize=(7, 0.28 * len(m) + 1.5))
    _cell_errorbar(ax, m, f, "b_S", ((np.log10(4), "log10 4", ":"), (np.log10(2), "log10 2", "--")),
                   "b_S: median log10 S = a_S - b_S n (primary range; dots = seeds, bar = bootstrap 95% CI)", "b_S (decades / qubit)")
    return _save(fig, path, plt)


def br_by_family_and_position(main, fits, path: Path) -> Path:
    plt = _plt()
    m, f = _primary(main, fits)
    fig, ax = plt.subplots(figsize=(7, 0.28 * len(m) + 1.5))
    _cell_errorbar(ax, m, f, "b_r", ((0.0, "b_r = 0", "-"),),
                   "b_r: median log10 |r| = a_r - b_r n (primary range)", "b_r (decades / qubit)")
    return _save(fig, path, plt)


def le_swap_slopes(main, fits, path: Path) -> Path:
    plt = _plt()
    m, f = _primary(main, fits)
    fig, axes = plt.subplots(1, 2, figsize=(11, 0.28 * len(m) + 1.5), sharey=True)
    _cell_errorbar(axes[0], m, f, "b_LE", (), "Loschmidt: median log10 M_LE = a + b_LE n", "b_LE (decades / qubit)")
    _cell_errorbar(axes[1], m, f, "b_SWAP", (), "SWAP: median log10 M_SWAP = a + b_SWAP n", "b_SWAP (decades / qubit)")
    return _save(fig, path, plt)


def _obs_vs_pred(main, fits, xcol_fn, ycol: str, xlabel: str, ylabel: str, title: str, path: Path) -> Path:
    plt = _plt()
    m, f = _primary(main, fits)
    fig, ax = plt.subplots(figsize=(5.6, 5.2))
    xs = []
    for _, r in m.iterrows():
        g = f[(f.family == r["family"]) & (f.regime == r["regime"]) & (f.position == r["position"])]
        x = xcol_fn(g)
        ax.plot(x, g[ycol], ".", color=POS_COLOR[r["position"]], alpha=0.25, ms=3)
        ax.errorbar(x.mean(), r[ycol], yerr=[[r[ycol] - r[f"{ycol}_ci_lo"]], [r[f"{ycol}_ci_hi"] - r[ycol]]], fmt=FAM_MARK[r["family"]],
                    color=POS_COLOR[r["position"]], ms=5, capsize=2)
        xs += list(x) + list(g[ycol])
    lo, hi = min(xs), max(xs)
    pad = 0.05 * (hi - lo)
    ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], "k-", lw=0.8, label="identity (prediction)")
    for fam, mk in FAM_MARK.items():
        ax.plot([], [], mk, color="0.3", label=FAM_SHORT[fam])
    for p, c in POS_COLOR.items():
        if p != "k=0":
            ax.plot([], [], "s", color=c, label=p)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=10)
    ax.legend(fontsize=7)
    return _save(fig, path, plt)


def observed_vs_predicted_le(main, fits, path: Path) -> Path:
    return _obs_vs_pred(main, fits, lambda g: (g.b_S + 2 * g.b_r).values, "b_LE", "b_S + 2 b_r (per seed; marker = seed mean)", "b_LE",
                        "Loschmidt slope vs A2 prediction (primary range)", path)


def observed_vs_predicted_swap(main, fits, path: Path) -> Path:
    return _obs_vs_pred(main, fits, lambda g: (2 * g.b_S + 2 * g.b_r).values, "b_SWAP", "2 b_S + 2 b_r (per seed; marker = seed mean)",
                        "b_SWAP", "SWAP slope vs A2 prediction (primary range)", path)


def gap_test(main, fits, path: Path) -> Path:
    plt = _plt()
    m, f = _primary(main, fits)
    fig, axes = plt.subplots(1, 2, figsize=(12, max(5.0, 0.28 * len(m) + 1.5)), gridspec_kw={"width_ratios": [1, 1.1]})
    ax = axes[0]
    xs = []
    for _, r in m.iterrows():
        g = f[(f.family == r["family"]) & (f.regime == r["regime"]) & (f.position == r["position"])]
        ax.plot(g.b_S, g.gap_obs, ".", color=POS_COLOR[r["position"]], alpha=0.25, ms=3)
        ax.plot(r["b_S"], r["gap_obs"], FAM_MARK[r["family"]], color=POS_COLOR[r["position"]], ms=5)
        xs += list(g.b_S) + list(g.gap_obs)
    lo, hi = min(xs), max(xs)
    ax.plot([lo, hi], [lo, hi], "k-", lw=0.8, label="b_SWAP - b_LE = b_S")
    for fam, mk in FAM_MARK.items():
        ax.plot([], [], mk, color="0.3", label=FAM_SHORT[fam])
    for p, c in POS_COLOR.items():
        if p != "k=0":
            ax.plot([], [], "s", color=c, label=p)
    ax.set_xlabel("b_S")
    ax.set_ylabel("b_SWAP - b_LE")
    ax.set_title("Gap rule (primary range; dots = seeds)", fontsize=10)
    ax.legend(fontsize=7)
    _cell_errorbar(axes[1], m, f, "eps_gap", ((0.0, "eps_gap = 0", "-"),), "eps_gap = (b_SWAP - b_LE) - b_S, bootstrap 95% CI",
                   "eps_gap (decades / qubit)")
    return _save(fig, path, plt)


def slope_ratio_vs_br(main, fits, path: Path) -> Path:
    plt = _plt()
    m, f = _primary(main, fits)
    fig, ax = plt.subplots(figsize=(6, 4.6))
    t = []
    for _, r in m.iterrows():
        g = f[(f.family == r["family"]) & (f.regime == r["regime"]) & (f.position == r["position"])]
        x = (g.b_r / g.b_S).values
        ax.plot(x, g.ratio_obs, ".", color=POS_COLOR[r["position"]], alpha=0.25, ms=3)
        ax.errorbar(x.mean(), r["ratio_obs"], yerr=[[r["ratio_obs"] - r["ratio_obs_ci_lo"]], [r["ratio_obs_ci_hi"] - r["ratio_obs"]]],
                    fmt=FAM_MARK[r["family"]], color=POS_COLOR[r["position"]], ms=5, capsize=2)
        t += list(x)
    tt = np.linspace(min(min(t), 0.0) - 0.02, max(t) + 0.02, 200)
    ax.plot(tt, (2 + 2 * tt) / (1 + 2 * tt), "k-", lw=0.8, label="(2 b_S + 2 b_r) / (b_S + 2 b_r)")
    ax.axhline(2.0, color="0.5", ls=":", lw=0.8, label="2 (b_r = 0)")
    ax.set_xlabel("b_r / b_S")
    ax.set_ylabel("b_SWAP / b_LE")
    ax.set_title("Exponent ratio vs relative-gradient slope (primary range)", fontsize=10)
    ax.legend(fontsize=7)
    return _save(fig, path, plt)


def parameter_position_effect(pos: pd.DataFrame, path: Path) -> Path:
    plt = _plt()
    p = pos[pos.fit_range == "primary"]
    groups = list(p.groupby(["family", "regime"], sort=False))
    fig, axes = plt.subplots(2, len(groups), figsize=(2.6 * len(groups), 5.6), sharey="row", squeeze=False)
    for j, ((fam, reg), g) in enumerate(groups):
        x = np.arange(len(g))
        for i, (col, lo, hi, ylab, c) in enumerate((("b_r_mean", "b_r_ci_lo", "b_r_ci_hi", "b_r", "C0"),
                                                     ("eps_gap_mean", "eps_gap_ci_lo", "eps_gap_ci_hi", "eps_gap", "C3"))):
            ax = axes[i, j]
            ax.errorbar(x, g[col], yerr=[g[col] - g[lo], g[hi] - g[col]], fmt="o", color=c, capsize=3)
            ax.axhline(0.0, color="0.5", lw=0.8)
            ax.set_xticks(x)
            ax.set_xticklabels(g.position, fontsize=7)
            if j == 0:
                ax.set_ylabel(ylab)
            if i == 0:
                ax.set_title(f"{FAM_SHORT[fam]} {reg}", fontsize=9)
    fig.suptitle("Parameter position: b_r and eps_gap, seed mean with bootstrap 95% CI (primary range)", fontsize=10)
    return _save(fig, path, plt)


def ratio_identity_error(iv: pd.DataFrame, path: Path) -> Path:
    plt = _plt()
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    for (fam, reg), g in iv.groupby(["family", "regime"], sort=False):
        h = g.groupby("n").max_rel_err.max()
        ax.semilogy(h.index, np.maximum(h.values, 1e-18), FAM_MARK[fam] + "-", ms=4, label=f"{FAM_SHORT[fam]} {reg}")
    ax.set_xlabel("n")
    ax.set_ylabel("max |M_SWAP/M_LE - (2-Q)/(S-Q)| / [(2-Q)/(S-Q)]")
    ax.set_title("Exact shot-ratio identity, every theta (max over seeds and positions)", fontsize=10)
    ax.legend(fontsize=7, ncol=2)
    return _save(fig, path, plt)


def conditional_sign_error(df: pd.DataFrame, path: Path) -> Path:
    plt = _plt()
    fig, ax = plt.subplots(figsize=(6.6, 4.6))
    for (fam, reg), g in df.groupby(["family", "regime"], sort=False):
        ax.semilogx(g.M_times_S, g.error, FAM_MARK[fam], ms=3, alpha=0.6, label=f"{FAM_SHORT[fam]} {reg}")
    ax.axhline(0.0, color="k", lw=0.8)
    ax.set_xlabel("M S")
    ax.set_ylabel("P(correct | g_hat != 0) - (1 + |r|)/2")
    ax.set_title("Loschmidt conditional sign law, exact binomial differences (n, seed, M, S-quantiles predeclared)", fontsize=9)
    ax.legend(fontsize=7, ncol=2)
    return _save(fig, path, plt)


def median_S_vs_n(med: pd.DataFrame, path: Path) -> Path:
    plt = _plt()
    fig, ax = plt.subplots(figsize=(7, 4.8))
    for (fam, reg, pos), g in med.groupby(["family", "regime", "position"], sort=False):
        h = g.groupby("n").median_log10_S.mean()
        ax.plot(h.index, h.values, FAM_MARK[fam] + "-", color=POS_COLOR[pos], ms=3, lw=0.8, alpha=0.8)
    for fam, mk in FAM_MARK.items():
        ax.plot([], [], mk + "-", color="0.3", label=FAM_SHORT[fam])
    for p, c in POS_COLOR.items():
        ax.plot([], [], "-", color=c, label=p)
    ax.axhline(-1.0, color="0.5", ls=":", lw=0.8, label="S = 0.1")
    ax.set_xlabel("n")
    ax.set_ylabel("seed mean of median log10 S")
    ax.set_title("Typical (median) S vs n, every cell", fontsize=10)
    ax.legend(fontsize=7, ncol=2)
    return _save(fig, path, plt)


def all_figures(med, fits, main, pos, iv, figs: Path) -> list[Path]:
    return [bS_by_family(main, fits, figs / "bS_by_family.png"),
            br_by_family_and_position(main, fits, figs / "br_by_family_and_position.png"),
            le_swap_slopes(main, fits, figs / "LE_SWAP_slopes.png"),
            observed_vs_predicted_le(main, fits, figs / "observed_vs_predicted_LE.png"),
            observed_vs_predicted_swap(main, fits, figs / "observed_vs_predicted_SWAP.png"),
            gap_test(main, fits, figs / "gap_test.png"),
            slope_ratio_vs_br(main, fits, figs / "slope_ratio_vs_br.png"),
            parameter_position_effect(pos, figs / "parameter_position_effect.png"),
            ratio_identity_error(iv, figs / "ratio_identity_error.png"),
            median_S_vs_n(med, figs / "median_S_vs_n.png")]
