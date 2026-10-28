"""B3: statistical hardening of the Stage 7 results (Satyabrat work allocation, delegated track).

Stage 7 (commit ba51866) is frozen. Nothing here changes a Stage 7 definition; every scientific quantity is
computed by calling the Stage 5-7 functions. B3 only adds uncertainty, kept as SEPARATE sources:

  seed_replicates  PRIMARY. Independent master seeds (10000..10019): a complete fresh theta population and
                   pipeline per seed; LE and SWAP evaluated on the same theta within a seed.
  bootstrap        SECONDARY. Paired theta bootstrap inside one fixed reference population (the regenerated
                   Stage 7 100k-theta set); the same resampled indices are used for LE and SWAP.
                   Also an optional, separately labelled two-level (seed, then theta) bootstrap.
  intervals        interval / summary helpers and a vectorized OLS slope.
  numerics         acceptance checks for extending n beyond 20 (rejects, never silently accepts).
  figures          uncertainty-aware versions of the key Stage 7 figures.
"""
