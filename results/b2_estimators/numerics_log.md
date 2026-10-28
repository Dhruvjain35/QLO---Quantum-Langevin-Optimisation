
## run 2026-10-07 20:26:24 CDT phase=validate config sha256=9dc965c2ede0a00145b75827382fd47333de4b6b11d17f5749cf3e4e0c58cf4a
validate: destructive/ancilla SWAP circuits max |P(+1) - (1+F)/2| = 4.44e-16; finite-shot max |z| = 2.42 (36 cells); var x M rel err max 0.080
validate: Hadamard circuits (numpy + PennyLane qml.ctrl, complex overlaps) max error = 6.66e-16
validate: Hadamard U-statistic exact (M <= 256) max |bias|,|var err|,|pmf-1| = 1.78e-15; gradient exact bias/var err = 8.88e-16; MC max |mean z| = 2.19, max |var rel err| = 0.0128 (5s)

## run 2026-10-07 20:26:40 CDT phase=rx config sha256=9dc965c2ede0a00145b75827382fd47333de4b6b11d17f5749cf3e4e0c58cf4a
rx: Stage 7 LE/SWAP median regression max |diff| = 1.78e-15 over 40 cells
rx: slopes of median log10 total executions (Stage 7 population, n = 2..20): loschmidt 0.6059, swap_ancilla 1.2025, swap_destructive 1.2025, hadamard_u 0.6066 (34s)
rx: ancilla vs destructive SWAP slope difference max 9.55e-15

## run 2026-10-07 20:27:27 CDT phase=fixed config sha256=9dc965c2ede0a00145b75827382fd47333de4b6b11d17f5749cf3e4e0c58cf4a

## run 2026-10-07 20:30:52 CDT phase=spot config sha256=9dc965c2ede0a00145b75827382fd47333de4b6b11d17f5749cf3e4e0c58cf4a
spot: SECONDARY ENTANGLING SPOT CHECK slopes (HEA depth n, EARLY, n 6..12): loschmidt 0.3023, swap_ancilla 0.6042, swap_destructive 0.6042, hadamard_u 0.3017 (47s)

## run 2026-10-07 20:31:41 CDT phase=figures config sha256=9dc965c2ede0a00145b75827382fd47333de4b6b11d17f5749cf3e4e0c58cf4a

## run 2026-10-07 20:32:09 CDT phase=figures config sha256=9dc965c2ede0a00145b75827382fd47333de4b6b11d17f5749cf3e4e0c58cf4a
note 20:32 hadamard_distribution_examples budgets changed to (32,128,512): exact joint pmf at T=2048 exceeded memory (illustration only; not a config item)

## run 2026-10-07 20:32:30 CDT phase=figures config sha256=9dc965c2ede0a00145b75827382fd47333de4b6b11d17f5749cf3e4e0c58cf4a
fixed: 96 (n, budget, estimator) cells written (the run then stopped in the distribution-example step; rerun as above)
note: post-freeze implementation change (before the sweep): hadamard_theory.t_pmf made sparse over reachable S^2 values (same exact pmf; dense convolution was too slow at M = 512). Re-validated: exact bias/var errors <= 1.8e-15.
