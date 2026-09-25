"""Stage 6: do all finite-shot barren plateaus fail the same way?

Compares two analytically controlled global barren plateaus on the same circuit
(|psi> = prod_j RX(theta_j)|0^n>, theta ~ U[-pi,pi]^n):

  A. GLOBAL PROJECTOR  C_G = 1 - prod_j cos^2(theta_j/2)      (Stages 3-5, frozen)
  B. GLOBAL PARITY     C_Z = [1 - prod_j cos(theta_j)] / 2    (new here)

Both have exponentially vanishing gradient variance. Their finite-shot parameter-shift
estimators are both differences of two binomial counts, but with shifted success
probabilities in different places: near 0 for the projector (both ~ A_k) and near 1/2
for parity. Stage 6 asks whether that changes the *mode* of failure.

parity_benchmark     exact parity cost / gradient / B_k / theory moments
finite_shot          binomial parity estimator, conditional variance, SNR, required shots
exact_distribution   difference-of-binomials: P(=), P(>), P(<) for either benchmark
directional          P_correct / P_wrong / P_zero, ambiguity regions, normal approximation
scaling              required-shot distributions and fits over random initializations
analysis             matched-signal binning, failure fractions, information distances
figures              plots
"""
