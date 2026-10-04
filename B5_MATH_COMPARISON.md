# B5 mathematical comparison: Stage 7 formulas vs the two prior papers

Purpose: place every Stage 7 formula beside what Thanasilp et al. 2024 (paper **A**) and Aghaei Saem et al. 2026
(paper **B**, reviewed as arXiv:2507.22054v2) explicitly contain.

**This document does not decide novelty.** For each Stage 7 expression it records one of:
- **exact prior-work match located**;
- **related formula located**;
- **no explicit formula located**.

Algebra that connects a prior formula to a Stage 7 formula is collected separately in §6, *Derived implications of
prior work*. That algebra is ours: neither source performs it. It carries no novelty label of any kind.

Notation (Stage 3–7, frozen): |ψ(θ)⟩ = ⊗_j RX(θ_j)|0ⁿ⟩, θ_j ~ U[−π, π], F = ∏_j cos²(θ_j/2), C = 1 − F. For component k:
- A = ∏_{j≠k} cos²(θ_j/2), s = sin θ_k;
- F± = F(θ ± (π/2)e_k) = A(1 ∓ s)/2;
- g = ∂C/∂θ_k = As/2;
- M shots per shifted circuit, one independent batch per (k, ±).

Evidence ids (A-xx, B-xx) refer to [`B5_EVIDENCE_LEDGER.md`](B5_EVIDENCE_LEDGER.md). † = present only in arXiv v2 of
paper B (absent from v1; published-version status unverified).

## 1. Landscape quantities

| Stage 7 expression | Status in A | Status in B |
|---|---|---|
| F = ∏ cos²(θ_j/2) concentrates over θ ~ U[−π, π]ⁿ | **Exact match (distributional)**: Prop. 3 / SI Eqs. (207)–(216). Fidelity of product single-qubit rotations with uniform angles, Var ≤ E[κ²] = (3/8)ⁿ, E[κ] = 2^(−n) (A-03). Stated for a kernel κ(x, x′) with y-axis rotations; the per-factor law equals Stage 7's x-rotation factor (our observation; it would not hold for z-rotations). | Related: same circuit used with global-Z costs (B-08c). Fidelity concentration stated only for a 2-design (B-09a†). |
| Log-typical scale A_typ = 4^(−(n−1)) (E log cos²(θ/2) = −2 log 2, Stage 5 Prop. 5) | No explicit formula located (A uses the mean 2^(−n) and second moment only). | No explicit formula located. |
| F± = A(1 ∓ s)/2, g = As/2 (shifted fidelities, exact gradient) | No explicit formula located (no parameter shift in A). | No explicit formula located (generic shift rule Eq. (6) only; B-01). |

## 2. Loschmidt (projector) readout

| # | Stage 7 expression | Status in A | Status in B |
|---|---|---|---|
| 2.1 | One shot ~ Bernoulli(F); F̂ = K/M | **Exact match located (fidelity level)**: main Eq. (12) + LE paragraph p. 4; SI Eq. (17) (A-05, A-08b). | **Exact match located (fidelity level)**: POVM {\|0⟩⟨0\|^⊗n, 𝟙 − \|0⟩⟨0\|^⊗n}, Fig. 5a (B-09a†). |
| 2.2 | P(F̂ = 0) = (1 − F)^M | **Exact match located**: (1 − s)^N inside SI Eqs. (26)–(27); (1 − μ)^N ≈ 1 − Nμ on p. 4 (A-05, A-06a). | Related: fixed distribution (0, 1) ⇒ estimate zero w.h.p. (B-09a†). No formula. |
| 2.3 | Var(F̂_LE) = F(1 − F)/M | Related: follows in one line from the stated outcome model; not written (A-08b). | **Exact match located (per shot)**: Var^(LE) = F(1 − F), §IV p. 10 (B-10a†). |
| 2.4 | ĝ_LE = (Ĉ₊ − Ĉ₋)/2 = (K₋ − K₊)/(2M) | No explicit formula located. | Related: generic finite-shot parameter-shift update, Eq. (6) (B-01). The Loschmidt gradient estimator is not written. |
| 2.5 | Var(ĝ_LE) = [F₊(1 − F₊) + F₋(1 − F₋)]/(4M) = [A − A²(1+s²)/2]/(4M) | No explicit formula located. | **Fidelity-estimator variance located, gradient-estimator variance not located** (B-10a†). Algebra: §6.2. |
| 2.6 | SNR²_LE = MAs²/[1 − A(1+s²)/2] | No explicit formula located. | No explicit formula located. |
| 2.7 | P(ĝ_LE = 0) = Σ_r Bin(r; M, F₊) Bin(r; M, F₋) (≥ (1 − F₊)^M(1 − F₋)^M) | Related: single-estimate (1 − s)^N only (A-06b). | Related: fidelity-level zero w.h.p. only (B-09b†). |
| 2.8 | P_correct, P_wrong (exact difference of binomials); P_correct → 0 through zeros | No explicit formula located. | No explicit formula located. Indistinguishability results are related (B-06c), but they concern the Pauli/±1 case. |
| 2.9 | P(correct \| ĝ_LE ≠ 0) → (1 + \|s\|)/2 as M·A → 0; population median (1 + 1/√2)/2 | No explicit formula located (A-NL13). | No explicit formula located (B-NL13). No prior formula from which it follows was located. Stage 7 §9 derives it from the single-count argument F₋/(F₊ + F₋) = (1 + s)/2. |
| 2.10 | Required shots: M_SNR ≈ ρ²/(As²), median ∝ 4ⁿ; fitted slopes 0.600–0.607 decades/qubit | Related: kernel-level numerics, "N ∈ Ω(2ⁿ)" for a fixed non-zero fraction (SI Fig. 4a; A-10b). | Related, implied only: Eq. (12) + Var^(LE) (B-10c†). Algebra: §6.3. |

## 3. SWAP-test readout

| # | Stage 7 expression | Status in A | Status in B |
|---|---|---|---|
| 3.1 | One shot: ancilla +1 w.p. q = (1 + F)/2; F̂ = 2K/M − 1 | **Exact match located (fidelity level)**: p₊ = 1/2 + κ/2, main p. 4; SI Eq. (46) (A-08). | **Exact match located (fidelity level)**: ancilla POVM {\|0⟩⟨0\|, \|1⟩⟨1\|}, Fig. 5b (B-09a†). |
| 3.2 | q± = (1 + F±)/2 | No explicit formula located (no shifts). | No explicit formula located. |
| 3.3 | Var(F̂_SWAP) = (1 − F²)/M | Related: follows from the stated outcome model; not written (A-08b). | **Exact match located (per shot)**: Var^(SWAP) = 1 − F², §IV p. 10 (B-10a†). |
| 3.4 | ĝ_SWAP = (Ĉ₊ − Ĉ₋)/2 = (K₋ − K₊)/M | No explicit formula located. | Related: Eq. (B18) writes the Pauli (±1-outcome) parameter-shift update as a difference of two ±1 means. Same structure with ℓ = F (B-06b). |
| 3.5 | Var(ĝ_SWAP) = [2 − F₊² − F₋²]/(4M) = [2 − A²(1+s²)/2]/(4M) | No explicit formula located. | **Fidelity-estimator variance located, gradient-estimator variance not located** (B-10a†). Algebra: §6.2. |
| 3.6 | SNR²_SWAP = MA²s²/[2 − A²(1+s²)/2] | No explicit formula located. | No explicit formula located. |
| 3.7 | P(ĝ_SWAP = 0) → C(2M, M)/4^M ≈ 1/√(πM) (n-independent) | Related: null estimate κ̂^(rand) = mean of ±1 equiprobable outcomes (Eq. 14; A-08). Tie probability not computed. | Related: null outcomes z = ±1 equiprobable (B15, B19, B26; B-06b). Tie probability not computed. Algebra: §6.5. |
| 3.8 | P_correct → (1 − P_zero)/2 (coin-flip sign) | No explicit formula located. | Related: indistinguishable from a sign-symmetric null (B-06c). No formula. |
| 3.9 | Required shots: M_SNR ≈ 2ρ²/(A²s²), median ∝ 16ⁿ; fitted slopes 1.197–1.203 decades/qubit | Related: kernel-level numerics, "at least exponentially" (SI Fig. 4b; A-10b). | Related, implied only (B-10c†). Algebra: §6.3. |

## 4. Paired (same-θ) readout comparison

| Stage 7 expression | Status in A | Status in B |
|---|---|---|
| SNR²_SWAP / SNR²_LE = A(1 − A(1+s²)/2)/(2 − A²(1+s²)/2) → A/2 | No explicit formula located. | No explicit formula located (implied by B-10a† after §6.2). |
| Exponent gap ≈ log₁₀4 per qubit (0.600–0.607 vs 1.197–1.203) | No explicit formula located (A-NL17). | No explicit formula located (B-NL17). |
| Per-shot TV: A\|s\| (LE) vs A\|s\|/2 (SWAP). Squared Hellinger ∝ A (LE) vs ∝ A² (SWAP) | Related: 1-norm (TV) hypothesis-testing bounds only (A-11). | Related: 1-norm bounds only (B-16). |
| Zero-vs-random failure split survives on an identical landscape | Related: the split is shown for kernels (same kernel, two tests; A-09a, A-10a). | Related: the split is stated for one fidelity objective under a 2-design (B-09a†). |

## 5. Full-vector and trajectory quantities

| Stage 7 expression | Status in A | Status in B |
|---|---|---|
| P(full ĝ = 0) = ∏_k P₀(θ, M, k) (independent batches) | Related: Pr[K̂ = 𝟙] = ∏_{i<j} Pr[κ̂_ij = 0] for an estimated Gram matrix (SI Eqs. 33–36; A-07). | Related: fidelity-level zero only (B-09b†). |
| Median cos(ĝ, g); P(ĝ·g > 0); component sign accuracy | No explicit formula located. | Related: update vector indistinguishable from a parameter-independent random vector (B-06d). No metric computed. |
| ‖ĝ‖/‖g‖ (SWAP ≈ 10^4.5 at n = 12, M = 1024) | No explicit formula located. | Related: null update magnitude is parameter-independent (B15; B-06e). Algebra: §6.6. |
| SWAP-driven GD vs signal-free random walk from identical starts; F_final as endpoint | No explicit formula located (A-NL24); random-matrix model control only (A-09b). | Related/explicit: random-walk reference on displacement statistics, Fig. 3b–c (B-07b); start matching not stated (B-07c). |
| Loschmidt-driven GD "frozen" (58–94% exactly-zero steps at n ≥ 10) | No explicit formula located. | No explicit statement located. Follows from applying Theorem 2 / Corollary 3 to the (0, 1) fixed distribution (§6.7). |

---

## 6. Derived implications of prior work

The algebra below is **ours**. Each step starts from a formula located in a source and states every assumption added.
None of it is stated in either source. It is recorded so that the principal audit can see exactly how far the
located prior formulas reach. It is not a judgment of novelty.

### 6.1 From A's single-estimate zero probability to the r = 0 term of the gradient zero probability

Located (A-06a): P(κ̂ = 0 | κ = s) = (1 − s)^N for one Loschmidt estimate from N shots.

Added assumptions: two shifted fidelities F± estimated with M shots each, independent batches; the parameter-shift
estimate is zero iff the two counts are equal.

Then P(ĝ_LE = 0) = P(K₊ = K₋) = Σ_{r=0}^{M} Bin(r; M, F₊) Bin(r; M, F₋). Its r = 0 term,
(1 − F₊)^M (1 − F₋)^M, is the product of two of A's factors and is a lower bound on P(ĝ = 0). The r ≥ 1 terms
(equal nonzero counts) have no counterpart in A. Stage 5 quantified them: the cruder "both counts zero" form misses
≈ 0.18 of probability at typical |s|, and the Poisson-limit form e^(−MA) I₀(MA|cos θ_k|) includes them.

### 6.2 From B's per-shot variances to Stage 7's gradient-estimator variances

Located (B-10a†): Var^(LE) = F(1 − F), Var^(SWAP) = 1 − F² (one shot, fidelity / infidelity estimate).

Added assumptions: Ĉ± = 1 − F̂± from M shots each; ĝ = (Ĉ₊ − Ĉ₋)/2; independent batches (Stage 7 design); shift
identities F₊ + F₋ = A and F₊² + F₋² = A²(1+s²)/2 (Stage 5 Prop. 1, not in either source).

- Var(ĝ_LE) = [F₊(1 − F₊) + F₋(1 − F₋)]/(4M) = [A − A²(1+s²)/2]/(4M)
- Var(ĝ_SWAP) = [(1 − F₊²) + (1 − F₋²)]/(4M) = [2 − A²(1+s²)/2]/(4M)
- With g = As/2: SNR²_LE = MAs²/[1 − A(1+s²)/2] and SNR²_SWAP = MA²s²/[2 − A²(1+s²)/2]. Their ratio → A/2 as A → 0.

These are Stage 7's §5–§6 expressions. The step from B's formulas is one line plus the landscape identities, which
are not in B.

### 6.3 From B's resolution ratio ε_N to shot exponents — the result depends on the criterion and the ensemble

Located (B-10b†, B-10c†): ε_N(α) = Var_ρ(α)[ℓ̂(α)]/(N Var_α[ℓ(α)]). Resolution needs ε_N ≲ 1, together with the two
per-shot variances. The source computes no exponent.

**(a) The source's own setting (2-design, fidelity level).** For a Haar-random state and fixed target in dimension
d = 2ⁿ:
- F ~ Beta(1, d − 1), so E F = 1/d, Var F = (d − 1)/(d²(d + 1)) ≈ d^(−2), and median F ≈ (ln 2)/d (checked by
  sampling at n = 6, 10).
- Loschmidt: N ≳ F(1 − F)/Var F ≈ d = 2ⁿ at typical F.
- SWAP: N ≳ (1 − F²)/Var F ≈ d² = 4ⁿ.

The SWAP requirement is approximately the square of the Loschmidt one (ratio ≈ 1/F). This is a readout-dependent
exponent **at the fidelity level under a 2-design**, implied by the source's formulas but not stated.

**(b) The same ε_N criterion on Stage 7's product landscape (fidelity level).** Var_θ F = (3/8)ⁿ − 4^(−n) ≈ (3/8)ⁿ
(checked by sampling), log-typical F = 4^(−n).
- Loschmidt: N ≳ F/Var_θ F, which is (2/3)ⁿ at the log-typical F (no exponential budget) or (4/3)ⁿ at the mean F = 2^(−n).
- SWAP: N ≳ 1/Var_θ F = (8/3)ⁿ.

These differ from Stage 7's exponents and depend on which scale of F is inserted. Var_θ F is dominated by rare near-target θ (heavy tail: (3/8)ⁿ ≫ F_typ² =
16^(−n)), so ε_N measures resolution against the landscape-wide spread, not against the local signal. The source calls
ε_N a rule of thumb.

**(c) A per-point signal-to-noise criterion (the criterion Stage 7 uses).**
- Fidelity level: N ≳ Var_1/F², i.e. ≈ 1/F (LE) vs ≈ 1/F² (SWAP).
- Gradient level (§6.2): M ≈ ρ²/(As²) (LE) vs ≈ 2ρ²/(A²s²) (SWAP).
- With the log-typical A = 4^(−(n−1)) (Stage 5 Prop. 5): median log₁₀ M grows by log₁₀4 = 0.602 (LE) and
  log₁₀16 = 1.204 (SWAP) decades per qubit. These are the references for Stage 7's fitted 0.600–0.607 and
  1.197–1.203.

**What is needed beyond the located formulas to reach Stage 7's gradient-level 4ⁿ vs 16ⁿ:**
1. B's per-shot variances (B-10a†).
2. The parameter-shift difference with independent batches.
3. The product-landscape shift identities.
4. A per-θ criterion (SNR ≥ ρ or P_correct ≥ q) instead of B's ε_N.
5. The log-typical value of A (not in A or B).
6. For sign targets, exact difference-of-binomial inversion (Stage 6/7; not located).

Neither source performs steps 2–6, and neither compares the readouts on a same-θ grid.

### 6.4 Why the located 1-norm (total-variation) bounds do not separate the readouts

Located: A's SI Eq. (14) (success ≤ 1/2 + N|ε|/2 for binary distributions; A-11) and B's Eq. (A19) (success ≤ 1/2 +
N‖P₀ − P₀′‖₁/4; B-16). Both tests in the sources compare an outcome distribution with a fixed one.

Our application to the two shifted per-shot distributions:
- TV_LE = |F₋ − F₊| = A|s|
- TV_SWAP = |q₋ − q₊| = A|s|/2

Both are linear in A, so a 1-norm bound of this form gives the same 1/(A|s|) threshold for both readouts. The squared
Hellinger distance, which governs the sample complexity of telling the two apart, scales differently:
- H²_LE ≈ A(√(1+s) − √(1−s))²/4 ∝ A
- H²_SWAP ≈ (As/2)²/2 ∝ A²

This was checked numerically at A = 10^(−4), 10^(−6) and s = 0.5, and it matches Stage 7 §15. The readout dependence
is therefore invisible to the located 1-norm machinery and appears only with a variance- or Hellinger-type criterion.

### 6.5 SWAP null tie probability from the located null distributions

Located: A's κ̂^(rand) = (1/N) Σ λ̃_m and B's z = ±1 means (B19, B26). Each is 2K/N − 1 with K ~ Bin(N, 1/2).

Added assumption: two independent such means, one per shift.

Then P(equal) = Σ_r C(N, r)² 4^(−N) = C(2N, N)/4^N ≈ 1/√(πN) (Vandermonde), independent of n. Checked: 0.017839 at
N = 1000 vs 1/√(1000π) = 0.017841. Stage 6/7's SWAP P_zero plateau. Not computed in either source.

### 6.6 Magnitude of B's null update

Located: B15, [Δ_N]_k = −(η/2) Σ_i c_i (Z̄_ik − Z̄′_ik), with Z̄ means of N equiprobable ±1 values.

Then E[Δ_k] = 0 and E‖Δ_N‖² = (η²/4) · N_p · (2/N) · Σ_i c_i². This is parameter-independent, while the exact step
η‖∇L‖ is exponentially small on a barren plateau, so the ratio of estimated to exact step norms grows as ‖∇L‖ shrinks.
Stage 7's measured SWAP per-component noise sd ≈ 1/√(2M) and median log₁₀ ‖ĝ‖/‖g‖ ≈ 4.5 at n = 12, M = 1024 are the
corresponding Stage 7 quantities. B states no norm comparison.

### 6.7 Loschmidt "frozen" dynamics from B's framework (not stated in B)

Located: B's Theorem 2 / Corollary 3, and the Loschmidt fixed distribution P_fixed = (0, 1) (B-09a†).

Our application: under P_fixed every outcome is "not all-zero", so every fixed-distribution fidelity estimate is
exactly 0 and the corresponding parameter-shift update is exactly zero. The fixed-distribution dynamics is therefore
a stationary point ("frozen"), not the ±1 random walk of B15, which is written for Pauli/±1 outcomes. Stage 7 §17
observed Loschmidt GD frozen (58–94% exactly-zero steps at n ≥ 10). Neither source states this contrast.

### 6.8 A's shot count and its constant-fluctuation assumption

Located (A-12a): A's sufficient shot count N ≥ 2‖O‖²_∞ log(2/p)/(ε̃² Var_α[X]) (SI Supp. Prop. 5) is Hoeffding-based, so
it uses only the outcome range. A's Gram-matrix corollary (SI Eq. 150) additionally assumes that "statistical
fluctuations associated with individual measurement outcomes stay constant". The source frames the count as a
requirement.

Our observation: replacing the range by the actual per-shot variance σ² (CLT/Bernstein-type reasoning) gives
N ≈ σ²/(ε̃² Var_α[X]).
- SWAP (±1 outcomes): σ² = 1 − F² ≈ 1, so the constant-fluctuation assumption holds.
- Loschmidt (0/1 outcomes): σ² = F(1 − F) ≈ F is exponentially small, so the assumption fails and the variance-aware
  count is smaller by the factor F.

This is the same mechanism that separates Stage 7's per-θ requirements (∝ 1/A for Loschmidt vs ∝ 1/A² for SWAP).
A states the assumption but does not connect it to the choice of readout.

## 7. One-line summary per candidate

- **C1 (exact zero/sign probabilities).**
  - Located: fidelity-level zero probability (A, B†) and the ±1 difference structure of the Pauli update (B).
  - Not located: the gradient-level difference-of-binomials pmf, P(ĝ = 0) with equal-count terms, P_correct and
    P_wrong.
- **C2 (conditional Loschmidt sign law).** No formula located in either source, and no located formula implies it.
- **C3 (readout-dependent gradient exponent on one landscape).**
  - Located: per-shot fidelity variances and a resolution ratio (B†), and fidelity-level same-kernel numerics (A).
  - Gradient-level exponents follow only via §6.2–§6.3 with added assumptions; the 4ⁿ vs 16ⁿ pair is not stated.
- **C4 (vector / trajectory consequences).**
  - Located: random-walk indistinguishability and displacement-statistics comparison (B), and a product-form joint
    zero probability for Gram matrices (A).
  - Not located: cosine, norm-ratio, dot-sign and component-sign metrics, and matched-start comparisons with a
    zero-signal control on fidelity endpoints.

This document is evidence for the principal novelty audit. It does not make the novelty decision.
