"""Stage 5: exact characterization of the finite-shot gradient "dead zone" on the Stage 3/4
global projector benchmark, and the matched term-wise local control. No optimizer here.

theory              exact distribution of the finite-shot parameter-shift estimator, P(g_hat=0),
                    conditional variance / SNR, required shots, deep-plateau approximations
sampling            vectorized log-space sampling of A_k, s_k over theta ~ U[-pi,pi]^n
global_deadzone     required-shot distributions, dead-zone fractions, Stage 4 retrodiction
local_resolvability matched term-wise local estimator control
analysis / figures  fits, summaries, plots
"""
