
## run 2026-10-07 17:18:06 CDT  phases=audit  git HEAD=735a8dfba93c8cb500a1205a70969395959ef5de  config sha256=2b7306ae18f9f6caf833670faf6ab55d48d11dcbb394d13930ab70bc5d77e497
audit: 40 (family, regime, n) cells, 140 position rows, 704s
  max agreement error (|diff|/S) over all rows: 3.17e-13; PennyLane-checked rows: 140/140
  positive control CONTROL_FINAL_RZ: {'STRUCTURAL': 20}
  predeclared positions: {'VALID': 120}

## run 2026-10-07 17:29:54 CDT  phases=rx  git HEAD=735a8dfba93c8cb500a1205a70969395959ef5de  config sha256=2b7306ae18f9f6caf833670faf6ab55d48d11dcbb394d13930ab70bc5d77e497
rx: B3 populations regenerated (20 seeds x 10 n x 25000 theta); max |slope diff| vs B3 seed_slopes.csv = 2.22e-16 over 80 values
  primary_2_20: regenerated LE summary mean 0.6053152901 vs published 0.6053152901; pct [0.6033809, 0.6069591] vs [0.6033809, 0.6069591]
  n_ge_8: regenerated LE summary mean 0.6016619291 vs published 0.6016619291; pct [0.5973119, 0.6050237] vs [0.5973119, 0.6050237]
rx: B1 arm 10 seeds x 10 n x 25000 theta (10s)
  [PRIMARY] primary           b_LE_ceil   B1 mean +0.604932 (boot +0.604370..+0.605592)  B3 pct +0.603381..+0.606959  inside=True overlap=True
  [PRIMARY] primary           b_SWAP_ceil B1 mean +1.200465 (boot +1.199505..+1.201510)  B3 pct +1.198529..+1.204970  inside=True overlap=True
  [PRIMARY] primary           gap_ceil    B1 mean +0.595533 (boot +0.594936..+0.596104)  B3 pct +0.594463..+0.598182  inside=True overlap=True
  [PRIMARY] primary           b_S         B1 mean +0.601301 (boot +0.600810..+0.601803)  B3 pct +0.600251..+0.603712  inside=True overlap=True
  [PRIMARY] primary           b_r         B1 mean -0.000023 (boot -0.000100..+0.000048)  B3 pct -0.000145..+0.000149  inside=True overlap=True
  [secondary] secondary_n_ge_8  b_LE_ceil   B1 mean +0.601549 (boot +0.600614..+0.602587)  B3 pct +0.597312..+0.605024  inside=True overlap=True
  [secondary] secondary_n_ge_8  b_SWAP_ceil B1 mean +1.202949 (boot +1.201297..+1.204772)  B3 pct +1.195894..+1.209470  inside=True overlap=True
  [secondary] secondary_n_ge_8  gap_ceil    B1 mean +0.601400 (boot +0.600410..+0.602351)  B3 pct +0.597628..+0.604885  inside=True overlap=True
  [secondary] secondary_n_ge_8  b_S         B1 mean +0.602035 (boot +0.601052..+0.603031)  B3 pct +0.598404..+0.605414  inside=True overlap=True
  [secondary] secondary_n_ge_8  b_r         B1 mean -0.000123 (boot -0.000227..-0.000024)  B3 pct -0.000336..+0.000322  inside=True overlap=True
  RX identity check: max rel err 1.75e-11, excluded rows 0

## run 2026-10-07 17:30:08 CDT  phases=sweep  git HEAD=735a8dfba93c8cb500a1205a70969395959ef5de  config sha256=2b7306ae18f9f6caf833670faf6ab55d48d11dcbb394d13930ab70bc5d77e497

## note 2026-10-07 17:33:05 CDT (sweep running; no entangling result and no A2 number read yet)
A2 comparison rule, fixed before any A2 number is read (applied only in the final a2 phase):
- corresponding cells = same family, same depth rule, same parameter. A2 shifted the first rotation on qubit 0 of
  the repo HEA, which is B1 hea_ring EARLY (RY, layer 0, qubit 0). Only the A2 depths that match a B1 regime
  (2, 4 or n) are compared; others are listed as "no corresponding B1 cell".
- if A2's n range differs from the B1 primary range, B1 slopes are ALSO recomputed on A2's n range from B1's own
  per-seed medians (post-hoc, labelled); the B1 primary-range values are shown alongside.
- agreement = the A2 value lies inside B1's 10-seed 95% percentile interval for that slope. A2 is one replicate
  (one seed, 4,000 theta per n), the same per-seed size as B1, so the replicate spread is the right yardstick
  (the same rule B3 used to check Stage 7). The B1 bootstrap CI of the seed mean is reported alongside.
Post-freeze code changes so far (no configuration change): pooled-over-seeds per-cell diagnostics
(cell_diagnostics_pooled.csv), eps_gap CIs in parameter_position_summary.csv and its figure, a legend in gap_test.png,
markdown table renderer (qlo.b1.report_tables), tests.

## note 2026-10-07 17:44:26 CDT sweep restarted
The first sweep attempt (6 workers, 64 MiB state batches) was stopped after 13 min with no task finished: the
machine was under heavy external memory pressure (swap 9.4/10 GB, internal disk 0.68 GB free, load ~125) and the
workers were thrashing. Restarted with 16 MiB state batches (qlo.b1.circuits.MAX_STATE_ELEMENTS = 2^20) and 4
workers. Batch size is an implementation detail (rows are simulated independently); no configuration change,
no output had been produced or inspected.

## run 2026-10-07 17:44:34 CDT  phases=sweep  git HEAD=735a8dfba93c8cb500a1205a70969395959ef5de  config sha256=2b7306ae18f9f6caf833670faf6ab55d48d11dcbb394d13930ab70bc5d77e497

## note 2026-10-07 17:56:11 CDT simulator kernel change before any sweep output
The relaunched sweep was still ~9x slower than the preflight (machine shared with other jobs; the M4's 4 performance
cores busy, workers on efficiency cores, memory bandwidth contended) and had finished no task after ~20 min. It was
stopped. qlo.b1.circuits now applies the rotations of up to 4 neighbouring qubits as one 16x16 Kronecker matrix per
row through BLAS matmul (FUSE_QUBITS = 4) instead of one 2x2 pass per qubit; the audited per-qubit path is kept as
fuse=1. Validation (fused_kernel_validation.csv): on the audit theta of all 40 (family, regime, n) cells, fused vs
the audited reference path agree to max |dF_pm|/S = 4.5e-14 and max |dF|/F = 7.5e-14; B1 unit tests (PennyLane
states and autograd gradients) pass on the fused path. The configuration is unchanged; no sweep output existed.

## run 2026-10-07 17:56:21 CDT  phases=sweep  git HEAD=735a8dfba93c8cb500a1205a70969395959ef5de  config sha256=2b7306ae18f9f6caf833670faf6ab55d48d11dcbb394d13930ab70bc5d77e497

## note 2026-10-07 19:44:07 CDT sweep resumed after the controlling session ended (87/400 tasks cached; cached tasks reused unchanged)

## run 2026-10-07 19:44:18 CDT  phases=sweep  git HEAD=735a8dfba93c8cb500a1205a70969395959ef5de  config sha256=2b7306ae18f9f6caf833670faf6ab55d48d11dcbb394d13930ab70bc5d77e497
sweep: 400 tasks (313 computed, 87 cached); compute 5.03 CPU-h; wall 418s; max |F(theta) from brackets - forward| / S = 1.95e-12

## run 2026-10-07 19:51:19 CDT  phases=analyse  git HEAD=735a8dfba93c8cb500a1205a70969395959ef5de  config sha256=2b7306ae18f9f6caf833670faf6ab55d48d11dcbb394d13930ab70bc5d77e497
analyse: identity check over 6740000 theta rows: max abs err 2.96e+02, max rel err 1.75e-11, excluded 0
  degenerate-row counts over all cells: {'n_S_zero': 0, 'n_delta_zero': 0, 'n_var_le_zero': 0, 'n_underflow': 0, 'n_r_numerically_zero': 0, 'n_nonfinite_log10_M_LE': 0, 'n_nonfinite_log10_M_SWAP': 0}; NaN rows in fitted medians: 0
analyse: 340 seed fits, 34 cells x ranges (38s)

## run 2026-10-07 19:52:17 CDT  phases=sign  git HEAD=735a8dfba93c8cb500a1205a70969395959ef5de  config sha256=2b7306ae18f9f6caf833670faf6ab55d48d11dcbb394d13930ab70bc5d77e497
sign: 720 (cell, theta, M) rows; max |error| where M S <= 1e-2: 4.84e-04 (265 rows); max |error| overall 3.28e-01

## run 2026-10-07 19:53:27 CDT  phases=a2  git HEAD=735a8dfba93c8cb500a1205a70969395959ef5de  config sha256=2b7306ae18f9f6caf833670faf6ab55d48d11dcbb394d13930ab70bc5d77e497
a2: 18 comparisons over 3 corresponding cells; A2 inside B1 seed 95% interval: 13/18; max |A2 - B1 mean| = 0.0127

## run 2026-10-07 19:54:06 CDT  phases=analyse  git HEAD=735a8dfba93c8cb500a1205a70969395959ef5de  config sha256=2b7306ae18f9f6caf833670faf6ab55d48d11dcbb394d13930ab70bc5d77e497
analyse: identity check over 6740000 theta rows: max abs err 2.96e+02, max rel err 1.75e-11, excluded 0
  degenerate-row counts over all cells: {'n_S_zero': 0, 'n_delta_zero': 0, 'n_var_le_zero': 0, 'n_underflow': 0, 'n_r_numerically_zero': 0, 'n_nonfinite_log10_M_LE': 0, 'n_nonfinite_log10_M_SWAP': 0}; NaN rows in fitted medians: 0
analyse: 340 seed fits, 34 cells x ranges (70s)
