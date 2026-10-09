"""Stage 7: same landscape, two measurement schemes.

The Stage 3-5 fidelity landscape (|psi(theta)> = prod_j RX(theta_j)|0^n>, F = prod_j cos^2(theta_j/2),
C = 1 - F) is estimated at finite shots in two physically different ways:

  loschmidt  projector / Loschmidt-echo readout: one Bernoulli(F) outcome per shot  (Stage 5 estimator)
  swap       SWAP test against |0^n>: ancilla Z = +1 with probability (1 + F)/2

Both estimate the SAME C, the SAME parameter-shift gradient g_k = A_k s_k / 2, at the SAME theta, so any
difference in estimator behaviour comes from the measurement statistics alone (no landscape confound).

estimators       per-scheme shifted outcome probabilities, gradient algebra, variance, SNR, resource accounting
required_shots   per-sample required shots (SNR, sign reliability, non-zero estimate), exact where feasible
analysis         paired same-theta tables, vector reliability, information distances, prior-work replication,
                 same-landscape optimization diagnostic
pennylane_check  explicit projector and ancilla SWAP-test circuits vs the binomial model
figures          plots
"""
