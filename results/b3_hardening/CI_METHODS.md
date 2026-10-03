# B3 confidence-interval methods

Every interval in `results/b3_hardening/` is labelled with its uncertainty source. "95%" never appears
without one of the labels below. The sources are reported **separately**, never pooled into one number,
except in the clearly labelled HIERARCHICAL rows.

## SEED: independent-seed replicate interval (PRIMARY)

**What varies.** Everything that depends on the master seed. Each of the 20 master seeds `10000..10019`
draws a completely fresh θ population for every n (25,000 θ per n for required shots, 10,000 θ per (n, M) for
fixed-shot probabilities, 400 θ × 25 shot replicates for full-gradient vectors, 50 paired starts per n for the
optimization diagnostic), and the whole Stage 7 pipeline runs on it. The Stage 7 seed (0) is excluded
(`check_seeds`, tested). LE and SWAP use the same θ within a seed.

**What question it answers.** Would the headline numbers change materially on a fresh random draw?

**What is reported** (`seed_summary.csv`, the `SEED_*` columns elsewhere), over the 20 seed-level values:
- mean, SD, median, IQR, min, max;
- **95% percentile interval**: the 2.5–97.5% percentiles of the 20 seed values, i.e. the replicate-to-replicate
  spread. With 20 values this is close to [min, max] (linear interpolation between the two smallest and the two largest);
- **95% bootstrap CI of the mean**: percentile bootstrap, 10,000 resamples, **seed is the resampling unit**.
  This is the uncertainty of the *average* slope, so it is narrower than the percentile interval.

**Exponent gap.** `γ_SWAP − γ_LE` is formed **inside each seed first** (`intervals.paired_gap`) and then summarized
over seeds, so the between-seed variation that LE and SWAP share (same θ) cancels. LE and SWAP seed results are
never treated as independent.

## THETA: paired θ-bootstrap interval (SECONDARY)

**What varies.** Only which θ points are drawn, inside one fixed reference population. The reference is the
regenerated Stage 7 population: master seed 0 with the Stage 7 θ streams, 100,000 θ per n for required shots
(20,000 per n for fixed-shot probabilities). Stage 7 did not save per-θ values, so they were recomputed
deterministically; the regenerated medians equal the Stage 7 medians (see `run_meta.json`, `stage7_reproduction`).

**What question it answers.** How much of the estimate is finite-initialization sampling noise in a population
of this size?

**Procedure.** 2,000 replicates. In each replicate and for each n, one vector of with-replacement resampling
counts is drawn and applied to **every series of that population** (both schemes, all targets), so LE and SWAP stay
paired on identical θ. Medians are recomputed at every n, the slope is refitted across n, and the gap is formed
within the replicate. Reported: bootstrap median, 2.5–97.5% percentile CI, bootstrap SE.

**Caveat.** Because the reference is the Stage 7 population, its plug-in estimate *is* the Stage 7 estimate. So
"Stage 7 inside the THETA CI" is close to automatic and is **not** evidence of robustness. The useful information is
the width of the CI.

## HIERARCHICAL: two-level bootstrap (optional, separately labelled)

**What varies.** Both sources. 2,000 replicates: resample the 20 seeds with replacement, then for every selected
seed and n resample its θ (one count vector shared by all series), refit each selected seed's slopes, and average
over the selected seeds. **Estimand:** the mean slope over seeds. Reported only for the primary range n = 2..20.
It does not replace the SEED or THETA rows.

## Finite-shot Monte Carlo (separate again)

`mc_validation.csv` checks the exact P_zero / P_correct formulas against binomial simulation (20,000 replicates per
cell) at the extended M and n. Its intervals are Wilson intervals over MC replicates and say nothing about θ or seed
variation.

## Not used as primary

Ordinary least-squares slope standard errors (`ols_slope_se_diagnostic`) are reported only as a diagnostic. They
measure scatter of the 10 medians around a line, not initialization-sampling uncertainty. An R² near 0.9999 is not by
itself evidence of robustness.

## Predeclared rule for the Stage 7 point-estimate check

A Stage 7 slope is one replicate computed on 100k θ. "Inside the seed interval" is therefore judged against the
SEED **95% percentile interval**. Inside / outside the stricter SEED bootstrap CI of the mean is reported alongside, plus
a z diagnostic in which the seed SD is rescaled from 25k to 100k θ (assuming variance ∝ 1/N). A miss in any column is
flagged and investigated. Neither intervals nor methods were changed after the results were seen.
