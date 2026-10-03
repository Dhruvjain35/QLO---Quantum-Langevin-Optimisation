"""Paired theta bootstrap (SECONDARY uncertainty) and the optional two-level bootstrap.

A replicate draws ONE vector of resampling counts per (population, n) and applies it to every series of that
population (both schemes, all targets), so the LE / SWAP pairing on identical theta is preserved inside every
replicate. Different n use different theta populations and are resampled independently.

Medians of a resample are computed from the resampling counts on pre-sorted values (``CountMedian``); this gives
exactly ``np.percentile(v[idx], 50)`` over the finite values (tested) at O(N) cost per series.
"""

from __future__ import annotations

import numpy as np

from qlo.hardening.intervals import ols_slopes


class CountMedian:
    """Median of a with-replacement resample of each row of ``vals`` (S, N), given per-index counts (N,)."""

    def __init__(self, vals: np.ndarray):
        vals = np.atleast_2d(np.asarray(vals, dtype=np.float64))
        self.rows = []
        for v in vals:
            fin = np.flatnonzero(np.isfinite(v))
            order = fin[np.argsort(v[fin], kind="stable")]
            self.rows.append((order, v[order]))

    def medians(self, counts: np.ndarray) -> np.ndarray:
        out = np.empty(len(self.rows))
        for i, (order, sv) in enumerate(self.rows):
            cum = np.cumsum(counts[order])
            tot = cum[-1] if cum.size else 0
            if tot == 0:
                out[i] = np.nan
                continue
            pos = 0.5 * (tot - 1)
            lo, hi = int(np.floor(pos)), int(np.ceil(pos))
            a = sv[np.searchsorted(cum, lo + 1)]
            b = sv[np.searchsorted(cum, hi + 1)]
            out[i] = a + (pos - lo) * (b - a)
        return out


def resample_counts(rng, n_items: int) -> np.ndarray:
    """Counts of a size-N with-replacement resample of N items."""
    return np.bincount(rng.integers(0, n_items, n_items), minlength=n_items)


def theta_bootstrap(pop: dict, n_values, n_boot: int, seed: int) -> np.ndarray:
    """``pop[n]`` = (S, N) per-theta values of one population. Returns medians (n_boot, S, len(n_values)); for every
    replicate and n the SAME counts are applied to all S series (pairing kept)."""
    rng = np.random.default_rng(np.random.SeedSequence([0xB3, 1, int(seed)]))
    cms = [CountMedian(pop[n]) for n in n_values]
    sizes = [pop[n].shape[1] for n in n_values]
    out = np.empty((int(n_boot), cms[0].rows.__len__(), len(n_values)))
    for b in range(int(n_boot)):
        for j, (cm, N) in enumerate(zip(cms, sizes)):
            out[b, :, j] = cm.medians(resample_counts(rng, N))
    return out


def hierarchical_bootstrap(pops: dict, seeds, n_values, n_boot: int, seed: int) -> np.ndarray:
    """TWO-LEVEL bootstrap. ``pops[seed][n]`` = (S, N). Each replicate resamples master seeds with replacement,
    then, for every selected seed and n, resamples theta (same counts for all series). Returns the replicate MEAN
    over selected seeds of each seed's fitted slope: array (n_boot, S)."""
    rng = np.random.default_rng(np.random.SeedSequence([0xB3, 2, int(seed)]))
    seeds = list(seeds)
    cms = {sd: [CountMedian(pops[sd][n]) for n in n_values] for sd in seeds}
    S = len(cms[seeds[0]][0].rows)
    out = np.empty((int(n_boot), S))
    for b in range(int(n_boot)):
        pick = rng.integers(0, len(seeds), len(seeds))
        slopes = np.empty((len(seeds), S))
        for i, k in enumerate(pick):
            sd = seeds[k]
            med = np.empty((S, len(n_values)))
            for j, n in enumerate(n_values):
                med[:, j] = cms[sd][j].medians(resample_counts(rng, pops[sd][n].shape[1]))
            slopes[i] = ols_slopes(n_values, med)[0]
        out[b] = slopes.mean(axis=0)
    return out


def bootstrap_fixed_medians(arrs: dict, n_boot: int, seed: int) -> dict:
    """``arrs[key]`` = (S, N) per-theta probabilities for ONE (n, M) cell; keys share the theta array. Returns
    {key: (n_boot, S) medians} with one count vector per replicate shared by every key (paired)."""
    rng = np.random.default_rng(np.random.SeedSequence([0xB3, 3, int(seed)]))
    cms = {k: CountMedian(v) for k, v in arrs.items()}
    N = next(iter(arrs.values())).shape[1]
    out = {k: np.empty((int(n_boot), len(cm.rows))) for k, cm in cms.items()}
    for b in range(int(n_boot)):
        c = resample_counts(rng, N)
        for k, cm in cms.items():
            out[k][b] = cm.medians(c)
    return out
