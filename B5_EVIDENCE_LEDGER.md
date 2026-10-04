# B5 evidence ledger

This ledger is the evidence packet for the principal novelty audit (A1). It records what two prior papers
explicitly contain, where, and under which assumptions. **It makes no novelty determination.**

Papers:

- **A** — S. Thanasilp, S. Wang, M. Cerezo, Z. Holmes, *Exponential concentration in quantum kernel methods*,
  Nature Communications **15**, 5200 (2024), doi:10.1038/s41467-024-49287-w. Reviewed: the published article
  (13 pp., online 2024-06-18) and its Supplementary Information PDF (`41467_2024_49287_MOESM1_ESM.pdf`, 51 pp.).
  Main-text locators are the published equation/figure numbers. Locators prefixed **SI** use the Supplementary
  Information's own numbering.
- **B** — R. Aghaei Saem, B. Tafreshi, Z. Holmes, S. Thanasilp, *Pitfalls when tackling the exponential concentration
  of parameterized quantum models*, Quantum Science and Technology **11**, 015049 (2026),
  doi:10.1088/2058-9565/ae2202 (online 2026-01-30). **The published (IOP) text could not be retrieved:** every IOP URL
  returned a bot-protection page, and no repository copy of the published version was found. **Every B locator below
  therefore refers to arXiv:2507.22054v2 (2026-06-04)**, read in full from its PDF and LaTeX source. Where the
  item is absent from arXiv:2507.22054v1 (2025-07-29), this is stated. Published-version page, equation and
  section numbers have **not** been verified (see "Open items" at the end).

Allowed classifications:

- `EXPLICIT` — the source states it.
- `PARTIAL / RELATED` — the source contains a related but non-identical statement.
- `PARTIAL / IMPLIED` — requested for candidate 3: the source's formulas imply the item only after extra algebra
  that the source does not perform.
- `NOT LOCATED` — not found in the reviewed version after the listed searches. It means only that, **not**
  that the item is absent from the literature.

One entry is one (source location, claim) pair. Where a single location supports several claims at different
strengths, it appears as separate entries (e.g. A-06a and A-06b). Matrix rows refer to `B5_PRIOR_WORK_MATRIX.md`.
Candidate ids: C1–C4 are the four Stage 7 candidate contributions; F-A…F-H are the "already expected prior" facts
of the B5 brief; CTX is context.

Stage 7 notation used in "Relation to Stage 7":
- F = ∏ cos²(θ_j/2), A = ∏_{j≠k} cos²(θ_j/2), s = sin θ_k, F± = A(1∓s)/2, g = As/2;
- Loschmidt (LE) counts K± ~ Bin(M, F±);
- SWAP probabilities q± = (1+F±)/2.

---

## Paper A — Thanasilp et al., Nat. Commun. 15, 5200 (2024)

<a id="A-01"></a>
### A-01 · Definition 1 (exponential concentration), Eqs. (9)–(11)
- **Paper:** A
- **Candidate:** F-A
- **Matrix rows:** 1
- **Claim category:** Exponential concentration of measured quantities (definition)
- **Classification:** EXPLICIT
- **Section:** Results
- **Subsection:** Why exponential concentration is problematic
- **Equation:** (9), (10), (11)
- **Figure:** —
- **Appendix:** SI Eqs. (18)–(19) restate it
- **Page:** 4
- **Source version:** Published version (Nat. Commun. 15, 5200; online 2024-06-18)
- **Source URL:** https://doi.org/10.1038/s41467-024-49287-w
- **Short quote:** —
- **Source statement (paraphrase):** A quantity X(α) measured as an expectation value is deterministically exponentially concentrated if |X(α) − μ| ≤ β ∈ O(1/bⁿ) for all α, and probabilistically concentrated if Pr_α[|X(α) − μ| ≥ δ] ≤ β/δ² with β ∈ O(1/bⁿ), b > 1. By Chebyshev, an exponentially small Var_α[X(α)] suffices. For kernels, α = {x, x′} and X = κ(x, x′).
- **Mathematical expression:** Pr_α[|X(α) − μ| ≥ δ] ≤ β/δ², β ∈ O(1/bⁿ); Var_α[X(α)] ∈ O(1/bⁿ)
- **Assumptions:** X(α) is the expectation of an observable; α drawn from a stated distribution.
- **Scope:** General definition; applied to fidelity and projected quantum kernels.
- **Relation to Stage 7:** Stage 7's fidelity F and gradient g over θ ~ U[−π, π]ⁿ concentrate in this sense (Stage 3 verified Var_θ[∂_k C] = (1/8)(3/8)^(n−1)).
- **Does NOT establish:** Anything about finite-shot parameter-shift gradient estimators, readout-dependent shot exponents, or zero/sign statistics.

<a id="A-02"></a>
### A-02 · Theorem 1 (expressivity-induced concentration), Eqs. (20)–(22)
- **Paper:** A
- **Candidate:** F-A
- **Matrix rows:** 1
- **Claim category:** Exponential concentration of the fidelity kernel (expressive / near-2-design embeddings)
- **Classification:** EXPLICIT
- **Section:** Results
- **Subsection:** Sources of exponential concentration — 1. Expressivity-induced concentration
- **Equation:** (17)–(22)
- **Figure:** Fig. 6
- **Appendix:** SI Note IV (proof)
- **Page:** 6–7
- **Source version:** Published version (Nat. Commun. 15, 5200; online 2024-06-18)
- **Source URL:** https://doi.org/10.1038/s41467-024-49287-w
- **Short quote:** —
- **Source statement (paraphrase):** Kernel values concentrate around their mean with a bound G_n(ε_U)/δ² set by an expressivity measure ε_U. For the fidelity kernel, G_n includes β_Haar = 1/(2^(n−1)(2ⁿ+1)). Close to a 2-design the kernel exponentially concentrates (fidelity-kernel mean 1/2ⁿ) and exponentially many shots are needed.
- **Mathematical expression:** Pr[|κ − E κ| ≥ δ] ≤ G_n(ε_U)/δ²; G_n = β_Haar + ε_U(ε_U + 2√β_Haar); β_Haar = 1/(2^(n−1)(2ⁿ+1))
- **Assumptions:** x, x′ drawn from the same distribution (relaxed in SI Note IV A); pure input states.
- **Scope:** Data-embedding ensembles; Haar / 2-design limit.
- **Relation to Stage 7:** Background. The Stage 7 landscape is a product (non-2-design) landscape (see A-03).
- **Does NOT establish:** Product-landscape typical (log-mean) values, finite-shot gradient statistics, readout comparison.

<a id="A-03"></a>
### A-03 · Proposition 3 (global-measurement-induced concentration), Eq. (26); SI Note VI Eqs. (207)–(216)
- **Paper:** A
- **Candidate:** F-A; CTX
- **Matrix rows:** 1, 27
- **Claim category:** Concentration of the fidelity between tensor-product single-qubit-rotation states
- **Classification:** EXPLICIT
- **Section:** Results
- **Subsection:** Sources of exponential concentration — 3. Global-measurement-induced concentration
- **Equation:** (26); SI (207)–(216), esp. SI (213), (214), (216)
- **Figure:** Fig. 7
- **Appendix:** SI Note VI (proof)
- **Page:** 8–9 (main; Fig. 7 on p. 9); SI 35
- **Source version:** Published version (Nat. Commun. 15, 5200) + Supplementary Information (MOESM1)
- **Source URL:** https://doi.org/10.1038/s41467-024-49287-w
- **Short quote:** —
- **Source statement (paraphrase):** Setting: fidelity kernel, product embedding U(x) = ⊗_k U_k(x_k) of single-qubit y-axis rotations, components independent and uniform on [−π, π], input |0⟩^⊗n. Result: Pr[|κ − 1/2ⁿ| ≥ δ] ≤ (3/8)ⁿ/δ². The proof shows E[κ²] = (3/8)ⁿ (per-qubit factor 3/8, SI Eq. 214) and E[κ] = 1/2ⁿ (SI Eq. 216). Fig. 7 confirms the exponential variance decay numerically for single layers of R_x, R_y and H·R_z.
- **Mathematical expression:** Var[κ^FQ] ≤ E[(κ^FQ)²] = (3/8)ⁿ; E[κ^FQ] = 2^(−n)
- **Assumptions:** U_k(x_k) = e^(−i x_k Y) (SI); x_k iid uniform on [−π, π]; product input state.
- **Scope:** The fidelity value between two encoded product states (a kernel), not a cost gradient; y-axis rotations only (U_k = e^(−i x_k Y)). The same single-qubit overlap law (cos² of a uniform angle) also holds for x-axis rotations acting on |0⟩ (our observation, not stated in the source) but not for z-axis rotations (overlap ≡ 1). The third case in Fig. 7 is Hadamard followed by R_z.
- **Relation to Stage 7:** The per-qubit factor is cos² of a uniformly distributed angle. For the R_x circuit of Stage 3–7 the same law holds (our observation): cos²(θ_j/2), θ_j ~ U[−π, π], has mean 1/2 and second moment 3/8. The Stage 7 fidelity landscape therefore belongs to this concentrating family. Stage 5/7 additionally used E[log cos²(θ/2)] = −2 log 2 (log-typical A = 4^(−(n−1))), which this proposition does not state.
- **Does NOT establish:** The gradient variance (1/8)(3/8)^(n−1) (Stage 3; Cerezo et al. 2021), log-typical values, shifted fidelities F±, or any finite-shot quantity.

<a id="A-04"></a>
### A-04 · SI Supplemental Proposition 6 (generalised global-measurement concentration), SI Eqs. (217)–(218)
- **Paper:** A
- **Candidate:** F-A
- **Matrix rows:** 1
- **Claim category:** Fidelity-kernel concentration for arbitrary local (product) unitaries
- **Classification:** EXPLICIT
- **Section:** SI Note VI
- **Subsection:** A. Extension to arbitrary local unitaries
- **Equation:** SI (217)–(228)
- **Figure:** —
- **Appendix:** SI Note VI A
- **Page:** SI 36–37
- **Source version:** Supplementary Information (MOESM1) to the published version
- **Source URL:** https://doi.org/10.1038/s41467-024-49287-w
- **Short quote:** —
- **Source statement (paraphrase):** For a general product embedding of independent local unitaries, the fidelity kernel concentrates with bound ∏_k G₁^(k)(ε_k)/δ², where G₁^(k) = 1/3 + ε_k(ε_k + √(4/3)). For fully random single-qubit unitaries this bound is 1/3ⁿ.
- **Mathematical expression:** Pr[|κ − μ| ≥ δ] ≤ ∏_k G₁^(k)/δ²; G₁^(k) = 1/3 + ε_k(ε_k + √(4/3))
- **Assumptions:** Independent components; product initial state; x and x′ drawn from the same distribution (used in the proof, SI p. 37).
- **Scope:** Fidelity kernel; product embeddings.
- **Relation to Stage 7:** Background generalisation of A-03.
- **Does NOT establish:** Any finite-shot or gradient-level statement.

<a id="A-05"></a>
### A-05 · Loschmidt Echo estimate and Proposition 1, Eq. (13)
- **Paper:** A
- **Candidate:** F-B
- **Matrix rows:** 2
- **Claim category:** Loschmidt (overlap-test) fidelity estimates collapse to zero
- **Classification:** EXPLICIT
- **Section:** Results
- **Subsection:** Why exponential concentration is problematic
- **Equation:** (12), (13)
- **Figure:** Fig. 1 (caption), Fig. 2
- **Appendix:** SI Supplemental Proposition 2 (full version)
- **Page:** 4 (Fig. 1 on p. 2; Fig. 2 on p. 5)
- **Source version:** Published version (Nat. Commun. 15, 5200; online 2024-06-18)
- **Source URL:** https://doi.org/10.1038/s41467-024-49287-w
- **Short quote:** —
- **Source statement (paraphrase):** In the Loschmidt Echo test, the kernel is the probability of the all-zero bitstring, with outcome 1 for all-zero and 0 otherwise. If the kernel concentrates to an exponentially small μ, the chance of never seeing all-zero in N shots is (1 − μ)^N ≈ 1 − Nμ, so it is likely that the statistical estimate is zero. Proposition 1: with N ∈ O(poly(n)) shots and a training set of size N_s, the estimated Gram matrix equals the identity with probability ≥ 1 − δ′, δ′ ∈ O(c^(−n)).
- **Mathematical expression:** P(no all-zero outcome in N shots) = (1 − μ)^N ≈ 1 − Nμ; Pr[K̂ = 𝟙] ≥ 1 − δ′, δ′ ∈ O(c^(−n))
- **Assumptions:** Kernel exponentially concentrates to an exponentially small μ (Def. 1); N ∈ O(poly(n)); independent shots; N_s ∈ O(poly(n)) for the Gram-matrix statement (SI pp. 7–8).
- **Scope:** Fidelity-kernel estimates (one estimate per data pair).
- **Relation to Stage 7:** Same outcome model as Stage 7's Loschmidt readout (Bernoulli(F) per shot). Stage 7 §18 lists the fidelity-level zero estimate as a replication of prior work.
- **Does NOT establish:** The zero probability of a parameter-shift gradient (a difference of two shifted estimates), the equal-nonzero-count (tie) contribution, or sign statistics.

<a id="A-06a"></a>
### A-06a · SI Supplemental Proposition 2, SI Eqs. (20)–(32)
- **Paper:** A
- **Candidate:** F-B
- **Matrix rows:** 2
- **Claim category:** Exact single-estimate Loschmidt zero probability (fidelity level)
- **Classification:** EXPLICIT
- **Section:** SI Note III
- **Subsection:** A. Fidelity quantum kernel — 1. Loschmidt Echo test
- **Equation:** SI (20), (21); proof SI (26)–(32)
- **Figure:** —
- **Appendix:** SI Note III A 1
- **Page:** SI 7–8
- **Source version:** Supplementary Information (MOESM1) to the published version
- **Source URL:** https://doi.org/10.1038/s41467-024-49287-w
- **Short quote:** —
- **Source statement (paraphrase):** With polynomially many shots, the Loschmidt estimate of the fidelity kernel is zero with probability ≥ 1 − δ, δ ∈ O(c^(−n)). The proof writes P(κ̂ = 0) as an integral of the exact conditional probability (1 − s)^N over the kernel-value distribution, then lower-bounds it by (1 − N(μ + β^(1/4)))(1 − √β).
- **Mathematical expression:** Pr[κ̂ = 0] = ∫₀¹ (1 − s)^N Pr[κ = s] ds ≥ (1 − N(μ + β^(1/4)))(1 − √β)
- **Assumptions:** Probabilistic exponential concentration with μ ∈ O(1/b′ⁿ), β ∈ O(1/bⁿ); N ∈ O(poly(n)).
- **Scope:** A single fidelity (kernel) estimate.
- **Relation to Stage 7:** The conditional factor (1 − s)^N is the exact zero probability of one Loschmidt fidelity estimate. Stage 7 §18 reports P(F̂ = 0) = (1 − F)^M at the fidelity level.
- **Does NOT establish:** Gradient-level P(ĝ = 0) (see A-06b).

<a id="A-06b"></a>
### A-06b · SI Supplemental Proposition 2 proof, SI Eqs. (26)–(27), read against the gradient-level question
- **Paper:** A
- **Candidate:** C1
- **Matrix rows:** 9
- **Claim category:** Gradient-level exact zero probability P(ĝ = 0)
- **Classification:** PARTIAL / RELATED
- **Section:** SI Note III
- **Subsection:** A. Fidelity quantum kernel — 1. Loschmidt Echo test
- **Equation:** SI (26)–(27)
- **Figure:** —
- **Appendix:** SI Note III A 1
- **Page:** SI 7
- **Source version:** Supplementary Information (MOESM1) to the published version
- **Source URL:** https://doi.org/10.1038/s41467-024-49287-w
- **Short quote:** —
- **Source statement (paraphrase):** Only the zero probability of a single fidelity estimate is derived, P(κ̂ = 0 | κ = s) = (1 − s)^N. No parameter-shift difference of two estimates is considered anywhere in the paper.
- **Mathematical expression:** P(κ̂ = 0 | κ = s) = (1 − s)^N
- **Assumptions:** As A-06a.
- **Scope:** Fidelity level, not gradient level.
- **Relation to Stage 7:** Stage 5/7 gradient-level law: P(ĝ = 0) = Σ_r Bin(r; M, F₊) Bin(r; M, F₋). Its r = 0 term, (1 − F₊)^M (1 − F₋)^M, is a product of two factors of the source's form. The r ≥ 1 equal-count terms have no counterpart in the source. Algebra in B5_MATH_COMPARISON §6.1.
- **Does NOT establish:** The difference-of-binomials law, the tie (r ≥ 1) mass, the Poisson-limit e^(−MA) I₀(MA|cos θ_k|) form, or any statement about gradients.

<a id="A-07"></a>
### A-07 · SI Eqs. (33)–(36): joint all-zero probability of the estimated Gram matrix
- **Paper:** A
- **Candidate:** C4
- **Matrix rows:** 19
- **Claim category:** Joint probability that a whole collection of Loschmidt estimates is exactly zero
- **Classification:** PARTIAL / RELATED
- **Section:** SI Note III
- **Subsection:** A. Fidelity quantum kernel — 1. Loschmidt Echo test (second half of the proof of Supp. Prop. 2)
- **Equation:** SI (33)–(36)
- **Figure:** —
- **Appendix:** SI Note III A 1
- **Page:** SI 8
- **Source version:** Supplementary Information (MOESM1) to the published version
- **Source URL:** https://doi.org/10.1038/s41467-024-49287-w
- **Short quote:** —
- **Source statement (paraphrase):** Kernel estimates for different data pairs are independent. So the probability that every off-diagonal Gram entry is estimated as exactly zero (K̂ = 𝟙) is the product of per-entry zero probabilities, bounded below by (1 − δ)^(N_s(N_s − 1)/2) ≥ 1 − N_s(N_s − 1)δ/2.
- **Mathematical expression:** Pr[K̂ = 𝟙] = ∏_{i<j} Pr[κ̂_ij = 0] ≥ (1 − δ)^(N_s(N_s−1)/2)
- **Assumptions:** Independent estimates per pair; per-entry concentration.
- **Scope:** Gram matrix of kernel values, not a gradient vector.
- **Relation to Stage 7:** Structurally the same as Stage 5 §10 and Stage 7 §14 (full-gradient zero probability = ∏_k P₀ under independent batches). The source applies it to kernel matrices, not gradient vectors.
- **Does NOT establish:** The full-gradient-vector zero probability, its values at the Stage 7 settings, or the conditional alignment of nonzero vectors.

<a id="A-08"></a>
### A-08 · SWAP-test estimate and Proposition 2, Eq. (14); SI Eqs. (46)–(49), (54)
- **Paper:** A
- **Candidate:** F-C
- **Matrix rows:** 3
- **Claim category:** SWAP-test fidelity estimates become data-independent random variables
- **Classification:** EXPLICIT
- **Section:** Results
- **Subsection:** Why exponential concentration is problematic
- **Equation:** (14); SI (46)–(49), SI (54)
- **Figure:** Fig. 1 (caption), Fig. 2; SI Fig. 3
- **Appendix:** SI Note III A 2 (Supplemental Lemma 3, Supplemental Corollary 2)
- **Page:** 4 (main; Fig. 1 on p. 2, Fig. 2 on p. 5); SI 9–13
- **Source version:** Published version (Nat. Commun. 15, 5200) + Supplementary Information (MOESM1)
- **Source URL:** https://doi.org/10.1038/s41467-024-49287-w
- **Short quote:** —
- **Source statement (paraphrase):** SWAP-test outcomes are +1 with probability p₊ = 1/2 + κ/2 and −1 otherwise, so the kernel is a perturbation of the uniform distribution. With polynomially many shots, each estimate is, with probability exponentially close to 1 over input pairs, statistically indistinguishable from κ̂^(rand) = (1/N) Σ λ̃_m with λ̃_m = ±1 equiprobable. For a polynomial-size training set the estimated Gram matrix is indistinguishable from a data-independent random matrix.
- **Mathematical expression:** P_κ = {(1+κ)/2, (1−κ)/2}; P₀ = {1/2, 1/2}; κ̂_N^(rand) = (1/N) Σ_m λ̃_m
- **Assumptions:** Exponential concentration to an exponentially small μ; N ∈ O(poly(n)); statements hold with probability ≥ 1 − δ_κ over input pairs; N_s ∈ O(poly(n)) for the Gram-matrix and model statements (SI Supp. Cors. 1–2, union bound SI Eqs. 59–61); statistical indistinguishability defined at success probability ≤ 0.51 (SI Def. 2).
- **Scope:** Fidelity-kernel estimates.
- **Relation to Stage 7:** Same outcome model as Stage 7's SWAP readout (q = (1+F)/2). Stage 7 §18 lists the data-independent SWAP fidelity estimate as a replication of prior work.
- **Does NOT establish:** A gradient-level SWAP estimator, its zero/tie probability (central binomial), its sign statistics, or its shot exponent.

<a id="A-08b"></a>
### A-08b · Outcome models of the two tests, Eqs. (12), (14); SI Eqs. (17), (46) — read against estimator variances
- **Paper:** A
- **Candidate:** F-D
- **Matrix rows:** 4
- **Claim category:** Different estimator variances for Loschmidt vs SWAP
- **Classification:** PARTIAL / RELATED
- **Section:** Results; SI Note III A
- **Subsection:** Why exponential concentration is problematic; SI III A 1–2
- **Equation:** (12), (14); SI (17), (46), (47)
- **Figure:** SI Fig. 1
- **Appendix:** SI Note III A
- **Page:** 4; SI 5–6, 9
- **Source version:** Published version + Supplementary Information (MOESM1)
- **Source URL:** https://doi.org/10.1038/s41467-024-49287-w
- **Short quote:** —
- **Source statement (paraphrase):** Both tests are written as empirical means of single-shot outcomes. Loschmidt: outcome 1/0 with P(1) = κ. SWAP: ±1 with P(+1) = (1+κ)/2. The two outcome distributions are fully specified, but the per-shot estimator variances (κ(1 − κ) and 1 − κ²) are not written down or compared.
- **Mathematical expression:** LE: λ ∈ {1, 0}, P(1) = κ; SWAP: λ ∈ {+1, −1}, P(+1) = (1+κ)/2
- **Assumptions:** —
- **Scope:** Fidelity-kernel level.
- **Relation to Stage 7:** The variances follow from these outcome models in one line (B5_MATH_COMPARISON §2.3, §3.3). Paper B states them explicitly (B-10a).
- **Does NOT establish:** The variance expressions themselves, or any consequence of their difference for shot exponents.

<a id="A-09a"></a>
### A-09a · Fig. 1 caption, Corollary 1 Eqs. (15)–(16), Figs. 2–3: same-kernel Loschmidt vs SWAP behaviour
- **Paper:** A
- **Candidate:** F-B; F-C; C3
- **Matrix rows:** 14a
- **Claim category:** Same-objective Loschmidt vs SWAP comparison at the fidelity/kernel-estimate level
- **Classification:** EXPLICIT
- **Section:** Introduction; Results
- **Subsection:** Why exponential concentration is problematic
- **Equation:** (15), (16)
- **Figure:** Fig. 1, Fig. 2, Fig. 3
- **Appendix:** SI Supplemental Corollaries 1–2
- **Page:** 2 (Fig. 1); 5 (Corollary 1, Figs. 2–3)
- **Source version:** Published version (Nat. Commun. 15, 5200; online 2024-06-18)
- **Source URL:** https://doi.org/10.1038/s41467-024-49287-w
- **Short quote:** —
- **Source statement (paraphrase):** For the same kernels and data, the two readouts behave differently. Loschmidt gives an identity estimated Gram matrix and zero predictions on unseen data. SWAP gives a random matrix and predictions that fluctuate around zero. Fig. 3 shows this numerically for an engineered 40-qubit dataset (tensor-product encoding, N = 1000 shots per kernel value, 10 repetitions), comparing training with exact kernels, LE estimates, SWAP estimates and a random matrix.
- **Mathematical expression:** a₀(y, λ) = y/(1 − λ) (LE); a_rand(y, λ) = (K̂_N^(rand) − λ𝟙)^(−1) y (SWAP)
- **Assumptions:** Kernel concentration; polynomial shots; kernel ridge regression; statements hold with probability exponentially close to 1; N_s ∈ O(poly(n)).
- **Scope:** Kernel estimates and downstream kernel-ridge-regression models; no parameter-shift gradients.
- **Relation to Stage 7:** A same-objective, two-readout comparison exists here at the kernel (fidelity) level. Stage 7's comparison is at the parameter-shift-gradient level on a paired θ grid.
- **Does NOT establish:** Gradient-level comparison, shot-exponent comparison, zero/sign probabilities.

<a id="A-09b"></a>
### A-09b · Fig. 3 "random" curve: a data-independent random-matrix baseline
- **Paper:** A
- **Candidate:** C4
- **Matrix rows:** 25
- **Claim category:** Signal-free control
- **Classification:** PARTIAL / RELATED
- **Section:** Results
- **Subsection:** Why exponential concentration is problematic
- **Equation:** —
- **Figure:** Fig. 3 (legend entry "random"); SI Fig. 6 (projected kernel analogue)
- **Appendix:** SI Note III B 2
- **Page:** 5–6; SI 19
- **Source version:** Published version + Supplementary Information (MOESM1)
- **Source URL:** https://doi.org/10.1038/s41467-024-49287-w
- **Short quote:** —
- **Source statement (paraphrase):** Fig. 3 includes a model trained on a random matrix whose off-diagonal entries are data-independent random variables. With SWAP estimates, the relative test loss follows that random-matrix control.
- **Mathematical expression:** —
- **Assumptions:** As A-09a.
- **Scope:** Kernel-model generalisation; not an optimisation trajectory.
- **Relation to Stage 7:** A signal-free control exists, applied to a trained kernel model rather than to a gradient-descent trajectory with matched starts.
- **Does NOT establish:** A random-walk optimisation control, matched-start trajectories, or alignment metrics.

<a id="A-10a"></a>
### A-10a · SI Supplementary Fig. 4 and SI Note III A 3: numerical Loschmidt vs SWAP on product-R_y kernels
- **Paper:** A
- **Candidate:** C3; F-E
- **Matrix rows:** 5, 14a
- **Claim category:** Numerical same-kernel comparison of the two readouts and the shot burden
- **Classification:** EXPLICIT
- **Section:** SI Note III
- **Subsection:** A. Fidelity quantum kernel — 3. Numerical simulation
- **Equation:** —
- **Figure:** SI Supplementary Fig. 4 (a) Loschmidt Echo, (b) SWAP
- **Appendix:** SI Note III A 3
- **Page:** SI 13–14
- **Source version:** Supplementary Information (MOESM1) to the published version
- **Source URL:** https://doi.org/10.1038/s41467-024-49287-w
- **Short quote:** "exponentially many measurement shots are required i.e., N ∈ Ω(2ⁿ)"
- **Source statement (paraphrase):** Setup: N_s = 25 inputs with components uniform on [0, 2π], tensor product of single-qubit R_y encodings, n = 5–40, shots per kernel value up to ~2×10⁶. Loschmidt: the fraction of zero estimates falls with N, and a fixed non-zero fraction (~0.75) needs N ∈ Ω(2ⁿ). At 30 and 40 qubits every estimate is zero even at 2×10⁶ shots. SWAP: the fraction of estimates passing a binomial test against the uniform distribution (p < 0.01) needs N to grow at least exponentially. Vertical lines mark 2ⁿ.
- **Mathematical expression:** Empirical: zero ratio (LE) and binomial-test success ratio (SWAP) vs N for n ∈ {5, 7, 10, 15, 20, 30, 40}
- **Assumptions:** Product-R_y encoding with uniform data; independent shots; binomial test at p < 0.01.
- **Scope:** Kernel (fidelity) estimates; numerical; no fitted exponents.
- **Relation to Stage 7:** Closest located analogue of Stage 5–7's per-θ shot requirements, on the same product-rotation fidelity family, but at the fidelity (not gradient) level and with different success metrics.
- **Does NOT establish:** Gradient-level shot requirements, a fitted exponent for either readout, or a comparison of the two exponents.

<a id="A-10b"></a>
### A-10b · SI Supplementary Fig. 4, read against the gradient-exponent questions
- **Paper:** A
- **Candidate:** C3
- **Matrix rows:** 15, 16
- **Claim category:** Shot-complexity exponent for each readout
- **Classification:** PARTIAL / RELATED
- **Section:** SI Note III
- **Subsection:** A. Fidelity quantum kernel — 3. Numerical simulation
- **Equation:** —
- **Figure:** SI Supplementary Fig. 4
- **Appendix:** SI Note III A 3
- **Page:** SI 14
- **Source version:** Supplementary Information (MOESM1) to the published version
- **Source URL:** https://doi.org/10.1038/s41467-024-49287-w
- **Short quote:** —
- **Source statement (paraphrase):** For Loschmidt, the source states N ∈ Ω(2ⁿ) is needed for a fixed non-zero fraction, a lower-bound statement at the kernel level based on numerics. For SWAP it states "at least exponentially". Neither exponent is fitted and the two are not compared.
- **Mathematical expression:** N_LE ∈ Ω(2ⁿ) (numerical, kernel level); N_SWAP: exponential (unspecified base)
- **Assumptions:** As A-10a.
- **Scope:** Kernel (fidelity) level; numerical.
- **Relation to Stage 7:** Stage 5/7 report per-θ median requirements at the gradient level: ≈4ⁿ for LE (log-typical A = 4^(−(n−1))) and ≈16ⁿ for SWAP. The SI Fig. 4 caption ties the 2ⁿ vertical lines to the Hilbert-space dimension. Elsewhere the source gives the mean fidelity kernel 2^(−n) for this family (Prop. 3). Reading 2ⁿ as an arithmetic-mean scale, as opposed to Stage 5's log-typical 4^(−(n−1)), is ours. Ω(2ⁿ) is a lower bound and does not contradict ≈4ⁿ.
- **Does NOT establish:** Exponent values, the 4ⁿ vs 16ⁿ pair, gradient-level results.

<a id="A-11"></a>
### A-11 · SI Note II: Supplemental Lemmas 1–2, Supplemental Proposition 1 (SI Eq. 14); SI Definitions 2–3
- **Paper:** A
- **Candidate:** F-F
- **Matrix rows:** 6
- **Claim category:** Outcome-distribution hypothesis testing / statistical indistinguishability as the finite-shot object
- **Classification:** EXPLICIT
- **Section:** SI Note II; SI Note III A 2
- **Subsection:** A. One sample; B. Many samples; SWAP test
- **Equation:** SI (1), (9), (14), (15), (48)
- **Figure:** SI Fig. 2
- **Appendix:** SI Notes II–III
- **Page:** SI 2–4, 10, 12
- **Source version:** Supplementary Information (MOESM1) to the published version
- **Source URL:** https://doi.org/10.1038/s41467-024-49287-w
- **Short quote:** —
- **Source statement (paraphrase):** One-sample success probability ≤ 1/2 + ‖P − Q‖₁/4. N-sample product bound ‖P^⊗N − Q^⊗N‖₁ ≤ N‖P − Q‖₁. For binary P₀ = (p₀, 1 − p₀) vs P_ε = (p₀ + ε, 1 − p₀ − ε), success ≤ 1/2 + N|ε|/2. Statistical indistinguishability is defined as success ≤ 0.51 (distributions, Def. 2) and extended to any map of the samples (outputs, Def. 3).
- **Mathematical expression:** Pr[right decision] ≤ 1/2 + N|ε|/2 (binary); ‖P^⊗N − Q^⊗N‖₁ ≤ N‖P − Q‖₁
- **Assumptions:** Independent samples; equal priors.
- **Scope:** Binary outcome distributions (kernel tests); 1-norm (total-variation) bounds.
- **Relation to Stage 7:** Stage 7 §15 notes that per-shot total variation is ∝ A for both readouts (SWAP's is exactly half of Loschmidt's), whereas squared Hellinger scales as A (LE) vs A² (SWAP). A 1-norm-based bound of this form gives the same 1/A threshold for both readouts (B5_MATH_COMPARISON §6.4).
- **Does NOT establish:** Hellinger/variance-based sample complexity or a readout-dependent exponent.

<a id="A-12a"></a>
### A-12a · SI Supplemental Proposition 5 (SI Eqs. 142–149) and Supplemental Corollary 5 (SI Eq. 150)
- **Paper:** A
- **Candidate:** F-E
- **Matrix rows:** 5
- **Claim category:** Exponential shot count sufficient for resolving concentrated quantities
- **Classification:** EXPLICIT
- **Section:** SI Note III
- **Subsection:** D. Sufficient condition to resolve kernel values
- **Equation:** SI (142)–(150)
- **Figure:** —
- **Appendix:** SI Note III D
- **Page:** SI 28–29
- **Source version:** Supplementary Information (MOESM1) to the published version
- **Source URL:** https://doi.org/10.1038/s41467-024-49287-w
- **Short quote:** "exponential scaling in measurement shots is indeed required"
- **Source statement (paraphrase):** With relative error ε̃ = ε/√Var_α[X], Hoeffding's inequality makes N ≥ 2‖O‖²_∞ log(2/p)/(ε̃² Var_α[X]) shots sufficient. Under exponential concentration this scales as Ω(b^(2n)/ε̃²), assuming ‖O‖_∞ ∈ O(1). For a full Gram matrix it is Ω(N_s² b^(2n)/ε̃²). The source presents these counts as requirements (SI p. 28; Cor. 5 "the number of measurement shots N required"), although the argument is a Hoeffding sufficiency bound.
- **Mathematical expression:** N ≥ 2‖O‖²_∞ log(2/p) / (ε̃² Var_α[X(α)]) ∈ Ω(b^(2n)/ε̃²)
- **Assumptions:** Bounded observable (‖O‖_∞ ∈ O(1)); Hoeffding (range-based) concentration; relative-error criterion ε̃ ≲ 1; Supp. Cor. 5 additionally assumes "statistical fluctuations associated with individual measurement outcomes stay constant" (SI p. 29).
- **Scope:** Any expectation-value estimate; the bound depends on the outcome range ‖O‖_∞, not on the measurement scheme's variance. The constant-fluctuation assumption of Supp. Cor. 5 holds for ±1 (SWAP-type) outcomes but not for the Loschmidt 0/1 outcome at small F (per-shot variance F(1 − F)). The source does not discuss this.
- **Relation to Stage 7:** This is the source's shot-count statement. It is range-based, so it cannot separate Loschmidt from SWAP exponents. The per-outcome fluctuation that Cor. 5 assumes constant is exactly what differs between the readouts in Stage 7 (per-shot variance F(1 − F) vs 1 − F²), giving ∝ 1/A vs ∝ 1/A² requirements. The source does not draw this connection.
- **Does NOT establish:** A measurement-dependent exponent, gradient-level requirements, or a proven necessary (lower-bound) shot count (the count is framed as required but derived from a sufficiency bound).

<a id="A-13"></a>
### A-13 · Main-text statements of the exponential shot burden
- **Paper:** A
- **Candidate:** F-E
- **Matrix rows:** 5
- **Claim category:** Exponential measurement burden under concentration
- **Classification:** EXPLICIT
- **Section:** Introduction; Results; Discussion
- **Subsection:** Introduction (p. 2); after Theorem 1 (p. 6); Discussion (p. 11)
- **Equation:** —
- **Figure:** —
- **Appendix:** —
- **Page:** 2, 6, 11
- **Source version:** Published version (Nat. Commun. 15, 5200; online 2024-06-18)
- **Source URL:** https://doi.org/10.1038/s41467-024-49287-w
- **Short quote:** —
- **Source statement (paraphrase):** When kernels exponentially concentrate, resolving them needs exponentially many shots. With polynomially many shots the trained model is independent of the input data.
- **Mathematical expression:** —
- **Assumptions:** Exponential concentration.
- **Scope:** Kernel estimation.
- **Relation to Stage 7:** The general fact Stage 5–7 take as given. Stage 7 adds readout-specific gradient-level exponents (see A-10b and B-10c for the closest prior statements).
- **Does NOT establish:** Readout-dependent exponents; gradient-level statements.

<a id="A-14"></a>
### A-14 · Corollary 1 (Eqs. 15–16); SI Definition 3, Supplemental Corollaries 1–2 (SI Eqs. 37–38, 54–57)
- **Paper:** A
- **Candidate:** F-H
- **Matrix rows:** 26
- **Claim category:** Post-processing (training/prediction) cannot recover information lost to concentration
- **Classification:** EXPLICIT
- **Section:** Results; SI Note III A
- **Subsection:** Why exponential concentration is problematic; SI III A 1–2
- **Equation:** (15), (16); SI (37)–(45), (54)–(61)
- **Figure:** SI Fig. 3
- **Appendix:** SI Note III A
- **Page:** 5; SI 8–13
- **Source version:** Published version + Supplementary Information (MOESM1)
- **Source URL:** https://doi.org/10.1038/s41467-024-49287-w
- **Short quote:** —
- **Source statement (paraphrase):** Any map applied to samples from indistinguishable distributions gives indistinguishable outputs (SI Def. 3). Optimal kernel-ridge parameters and predictions trained on concentrated estimates are, with probability exponentially close to 1, equal to (Loschmidt) or indistinguishable from (SWAP) data-independent quantities.
- **Mathematical expression:** Pr[a_opt = a₀(y, λ)] ≥ 1 − δ (LE); outputs of Φ on indistinguishable samples are indistinguishable (Def. 3)
- **Assumptions:** Kernel concentration; polynomial shots; polynomial training set.
- **Scope:** Kernel methods.
- **Relation to Stage 7:** The no-post-processing fact that Stage 7 does not claim.
- **Does NOT establish:** Gradient-level statements.

<a id="A-15"></a>
### A-15 · SI Note III C: Supplemental Lemma 6, Proposition 4, Corollary 4 (SI Eqs. 134–141)
- **Paper:** A
- **Candidate:** F-H
- **Matrix rows:** 26
- **Claim category:** Multi-copy coherent processing does not distinguish exponentially close states
- **Classification:** EXPLICIT
- **Section:** SI Note III
- **Subsection:** C. Indistinguishability of concentrated quantum states
- **Equation:** SI (134)–(141)
- **Figure:** —
- **Appendix:** SI Note III C
- **Page:** SI 27–28
- **Source version:** Supplementary Information (MOESM1) to the published version
- **Source URL:** https://doi.org/10.1038/s41467-024-49287-w
- **Short quote:** —
- **Source statement (paraphrase):** Given m copies of ρ or σ, the optimal success probability is ≤ 1/2 + m‖ρ − σ‖₁/4. For exponentially close states and polynomial m, they remain indistinguishable even with coherent processing. Hence multi-copy error-mitigation schemes cannot remove concentration-induced data independence.
- **Mathematical expression:** Pr[right decision] ≤ 1/2 + m‖ρ − σ‖₁/4
- **Assumptions:** ‖ρ − σ‖₁ ∈ O(1/bⁿ); m ∈ O(poly(n)).
- **Scope:** State-level distinguishability.
- **Relation to Stage 7:** Background; Stage 7 makes no state-discrimination claim.
- **Does NOT establish:** Anything about specific readouts' gradient exponents.

<a id="A-16"></a>
### A-16 · SI Note VIII (error mitigation), Supplemental Theorem 1 (SI Eqs. 274–278); main text p. 10
- **Paper:** A
- **Candidate:** F-H
- **Matrix rows:** 26
- **Claim category:** Error mitigation cannot remove (noise-induced) kernel concentration
- **Classification:** EXPLICIT
- **Section:** SI Note VIII; Results
- **Subsection:** SI VIII; main text "4. Noise-induced concentration"
- **Equation:** SI (274)–(278)
- **Figure:** —
- **Appendix:** SI Note VIII
- **Page:** SI 42–43; main 10
- **Source version:** Published version + Supplementary Information (MOESM1)
- **Source URL:** https://doi.org/10.1038/s41467-024-49287-w
- **Short quote:** —
- **Source statement (paraphrase):** Quoting a result of its Ref. [18] (Wang et al.), common error-mitigation strategies within a unified framework concentrate onto a state-independent fixed point at linear depth. So they cannot mitigate noise-induced exponential concentration of kernel values, and can even impair resolvability.
- **Mathematical expression:** |C_m − F₀| ∈ O(2^(−bn) a_max |T_EM| M_max)
- **Assumptions:** Local depolarising noise; circuit depth Ω(n); the cited theorem's assumptions.
- **Scope:** Noisy kernels; error mitigation.
- **Relation to Stage 7:** Background (Stage 7 is noiseless; B4 is Satyabrat's).
- **Does NOT establish:** Anything about noiseless finite-shot gradient readouts.

<a id="A-17"></a>
### A-17 · Proposition 4 (Eq. 37), Fig. 9; SI Notes IX–X: trainable embeddings
- **Paper:** A
- **Candidate:** CTX
- **Matrix rows:** —
- **Claim category:** Exponentially flat training landscape of kernel target alignment (trainability)
- **Classification:** EXPLICIT
- **Section:** Results
- **Subsection:** Training parameterized quantum kernels
- **Equation:** (35)–(37); SI (280)–(342)
- **Figure:** Fig. 9
- **Appendix:** SI Notes IX–X
- **Page:** 10–11; SI 43–50
- **Source version:** Published version + Supplementary Information (MOESM1)
- **Source URL:** https://doi.org/10.1038/s41467-024-49287-w
- **Short quote:** —
- **Source statement (paraphrase):** The deviation probability of the kernel target alignment over variational parameters is "approximately bounded" by the variances of the parameterised kernels. If those vanish exponentially, the source concludes that the alignment landscape is exponentially flat and untrainable with polynomially many shots. Fig. 9 shows the exponential decay of Var_θ[TA] (exact variances over 500 initialisations of a single R_y layer + HEE; no finite-shot simulation).
- **Mathematical expression:** Pr_θ[|TA(θ) − E TA| ≥ δ] ≲ M Σ_{ij} Var_θ[κ_θ(x_i, x_j)]/δ² ("approximately bounded", Eq. (37))
- **Assumptions:** As stated in Proposition 4.
- **Scope:** Trainability of parameterised kernels; no finite-shot gradient-estimator analysis, no trajectories.
- **Relation to Stage 7:** Context only; listed so that the absence of a random-walk statement in paper A (A-NL08) is not mistaken for absence of any trainability content.
- **Does NOT establish:** Random-walk behaviour, parameter-shift finite-shot statistics, or trajectory comparisons.

<a id="A-18"></a>
### A-18 · Source's own survey of prior work: Introduction (p. 2) and SI Note I
- **Paper:** A
- **Candidate:** CTX
- **Matrix rows:** —
- **Claim category:** Source's positioning of earlier work (for follow-up only)
- **Classification:** EXPLICIT
- **Section:** Introduction; SI Note I
- **Subsection:** SI I.A–I.D
- **Equation:** —
- **Figure:** —
- **Appendix:** SI Note I
- **Page:** 2; SI 1–2
- **Source version:** Published version + Supplementary Information (MOESM1)
- **Source URL:** https://doi.org/10.1038/s41467-024-49287-w
- **Short quote:** —
- **Source statement (paraphrase):** The source says its refs. 8 and 16 rigorously study the shots needed to train the fidelity kernel without addressing concentration. SI I.D says its SI refs. [9, 10] study shot noise in kernel methods without concentration. SI I.B says its SI ref. [3] suggests, without proving, that exponentially many shots are needed for certain embeddings.
- **Mathematical expression:** —
- **Assumptions:** —
- **Scope:** Bibliographic.
- **Relation to Stage 7:** Feeds B5_FOLLOWUP_SOURCES.md (shot-complexity studies of fidelity-kernel estimation).
- **Does NOT establish:** Anything about those works' contents (not reviewed in B5).

### Paper A — not-located search log

All entries below: whole paper reviewed (main text pp. 1–13 and SI pp. 1–51), by reading every page plus text
searches of the extracted text. Common search terms for all entries: gradient, parameter shift, parameter-shift,
derivative, sign, direction, conditional, tie, binomial, Skellam, cosine, angle, inner product, norm, random walk,
trajectory, exponent, 4^n, 16^n, sample complexity, shots, SNR, signal-to-noise.

<a id="A-NL07"></a>
### A-NL07 · Parameter-shift finite-shot analysis
- **Paper:** A
- **Candidate:** C1; C3
- **Matrix rows:** 7
- **Claim category:** Parameter-shift finite-shot analysis
- **Classification:** NOT LOCATED
- **Section:** Whole paper (main text + Supplementary Information)
- **Subsection:** —
- **Equation:** —
- **Figure:** —
- **Appendix:** —
- **Page:** main 1–13; SI 1–51
- **Source version:** Published version + Supplementary Information (MOESM1)
- **Source URL:** https://doi.org/10.1038/s41467-024-49287-w
- **Short quote:** —
- **Source statement (paraphrase):** Not located in the reviewed version after searches for: parameter shift, parameter-shift, gradient estimator, derivative, finite-difference. "Gradient" appears only in barren-plateau background and reference titles. Kernel target alignment (A-17) is treated through variances, not finite-shot gradient estimators.
- **Mathematical expression:** —
- **Assumptions:** —
- **Scope:** Search covered main text, all SI notes and figure captions.
- **Relation to Stage 7:** —
- **Does NOT establish:** Absence from the wider literature.

<a id="A-NL08"></a>
### A-NL08 · Parameter-shift random-walk behaviour
- **Paper:** A
- **Candidate:** C4; F-G
- **Matrix rows:** 8
- **Claim category:** Random-walk behaviour of finite-shot gradient descent
- **Classification:** NOT LOCATED
- **Section:** Whole paper (main text + Supplementary Information)
- **Subsection:** —
- **Equation:** —
- **Figure:** —
- **Appendix:** —
- **Page:** main 1–13; SI 1–51
- **Source version:** Published version + Supplementary Information (MOESM1)
- **Source URL:** https://doi.org/10.1038/s41467-024-49287-w
- **Short quote:** —
- **Source statement (paraphrase):** Not located in the reviewed version after searches for: random walk, trajectory, gradient descent, optimization step, update rule. The closest content is the random-matrix Gram estimate (A-08, A-09b) and the kernel-target-alignment trainability statement (A-17), neither of which concerns gradient-descent trajectories.
- **Mathematical expression:** —
- **Assumptions:** —
- **Scope:** As above.
- **Relation to Stage 7:** —
- **Does NOT establish:** Absence from the wider literature.

<a id="A-NL10"></a>
### A-NL10 · Exact gradient P_correct
- **Paper:** A
- **Candidate:** C1
- **Matrix rows:** 10
- **Claim category:** Exact finite-shot gradient correct-sign probability
- **Classification:** NOT LOCATED
- **Section:** Whole paper (main text + Supplementary Information)
- **Subsection:** —
- **Equation:** —
- **Figure:** —
- **Appendix:** —
- **Page:** main 1–13; SI 1–51
- **Source version:** Published version + Supplementary Information (MOESM1)
- **Source URL:** https://doi.org/10.1038/s41467-024-49287-w
- **Short quote:** —
- **Source statement (paraphrase):** Not located in the reviewed version after searches for: sign, correct direction, direction, P_correct, probability of correct sign, binomial difference. No sign or direction statistics are given for any estimate.
- **Mathematical expression:** —
- **Assumptions:** —
- **Scope:** As above.
- **Relation to Stage 7:** —
- **Does NOT establish:** Absence from the wider literature.

<a id="A-NL11"></a>
### A-NL11 · Exact gradient P_wrong
- **Paper:** A
- **Candidate:** C1
- **Matrix rows:** 11
- **Claim category:** Exact finite-shot gradient wrong-sign probability
- **Classification:** NOT LOCATED
- **Section:** Whole paper (main text + Supplementary Information)
- **Subsection:** —
- **Equation:** —
- **Figure:** —
- **Appendix:** —
- **Page:** main 1–13; SI 1–51
- **Source version:** Published version + Supplementary Information (MOESM1)
- **Source URL:** https://doi.org/10.1038/s41467-024-49287-w
- **Short quote:** —
- **Source statement (paraphrase):** Not located in the reviewed version after searches for: wrong sign, sign, direction, P_wrong. No sign statistics appear.
- **Mathematical expression:** —
- **Assumptions:** —
- **Scope:** As above.
- **Relation to Stage 7:** —
- **Does NOT establish:** Absence from the wider literature.

<a id="A-NL12"></a>
### A-NL12 · Exact difference-of-binomials gradient law
- **Paper:** A
- **Candidate:** C1
- **Matrix rows:** 12
- **Claim category:** Distribution of a difference of two binomial counts (parameter-shift estimator)
- **Classification:** NOT LOCATED
- **Section:** Whole paper (main text + Supplementary Information)
- **Subsection:** —
- **Equation:** —
- **Figure:** —
- **Appendix:** —
- **Page:** main 1–13; SI 1–51
- **Source version:** Published version + Supplementary Information (MOESM1)
- **Source URL:** https://doi.org/10.1038/s41467-024-49287-w
- **Short quote:** —
- **Source statement (paraphrase):** Not located in the reviewed version after searches for: binomial, difference, Skellam, tie, equal counts. "Binomial" occurs only for the binomial hypothesis test of SI Figs. 4–5. All estimates are single empirical means; no difference of two estimates is analysed.
- **Mathematical expression:** —
- **Assumptions:** —
- **Scope:** As above.
- **Relation to Stage 7:** —
- **Does NOT establish:** Absence from the wider literature.

<a id="A-NL13"></a>
### A-NL13 · Conditional Loschmidt sign law
- **Paper:** A
- **Candidate:** C2
- **Matrix rows:** 13
- **Claim category:** P(correct sign | nonzero estimate) for the Loschmidt readout
- **Classification:** NOT LOCATED
- **Section:** Whole paper (main text + Supplementary Information)
- **Subsection:** —
- **Equation:** —
- **Figure:** —
- **Appendix:** —
- **Page:** main 1–13; SI 1–51
- **Source version:** Published version + Supplementary Information (MOESM1)
- **Source URL:** https://doi.org/10.1038/s41467-024-49287-w
- **Short quote:** —
- **Source statement (paraphrase):** Not located in the reviewed version after searches for: conditional, given nonzero, sign, direction, success probability, single count, (1+|sin|)/2. Conditional probabilities appear only as the zero-estimate integrand P(κ̂ = 0 | κ = s) (A-06a) and in hypothesis-test success probabilities (SI Eqs. (2)–(3), (50)). None is a sign law.
- **Mathematical expression:** —
- **Assumptions:** —
- **Scope:** As above.
- **Relation to Stage 7:** —
- **Does NOT establish:** Absence from the wider literature.

<a id="A-NL14b"></a>
### A-NL14b · Same-objective Loschmidt vs SWAP comparison at the parameter-shift-gradient level
- **Paper:** A
- **Candidate:** C3
- **Matrix rows:** 14b
- **Claim category:** Same-θ, same-gradient comparison of two readouts
- **Classification:** NOT LOCATED
- **Section:** Whole paper (main text + Supplementary Information)
- **Subsection:** —
- **Equation:** —
- **Figure:** —
- **Appendix:** —
- **Page:** main 1–13; SI 1–51
- **Source version:** Published version + Supplementary Information (MOESM1)
- **Source URL:** https://doi.org/10.1038/s41467-024-49287-w
- **Short quote:** —
- **Source statement (paraphrase):** Not located in the reviewed version after searches for: gradient, parameter shift, same landscape, same objective. The two-readout comparison exists only for kernel estimates (A-09a, A-10a).
- **Mathematical expression:** —
- **Assumptions:** —
- **Scope:** As above.
- **Relation to Stage 7:** —
- **Does NOT establish:** Absence from the wider literature.

<a id="A-NL17"></a>
### A-NL17 · Explicit 4ⁿ vs 16ⁿ comparison
- **Paper:** A
- **Candidate:** C3
- **Matrix rows:** 17
- **Claim category:** Explicit pair of shot exponents 4ⁿ (LE) and 16ⁿ (SWAP)
- **Classification:** NOT LOCATED
- **Section:** Whole paper (main text + Supplementary Information)
- **Subsection:** —
- **Equation:** —
- **Figure:** —
- **Appendix:** —
- **Page:** main 1–13; SI 1–51
- **Source version:** Published version + Supplementary Information (MOESM1)
- **Source URL:** https://doi.org/10.1038/s41467-024-49287-w
- **Short quote:** —
- **Source statement (paraphrase):** Not located in the reviewed version after searches for: 4^n, 16^n, 4n, 16n, exponent, base of exponential, decades per qubit, slope. Exponential statements are Ω(2ⁿ) (numerical, LE, kernel level; A-10b), Ω(b^(2n)) (scheme-independent sufficient bound; A-12a), and "at least exponentially" (SWAP; A-10b).
- **Mathematical expression:** —
- **Assumptions:** —
- **Scope:** As above.
- **Relation to Stage 7:** —
- **Does NOT establish:** Absence from the wider literature.

<a id="A-NL18"></a>
### A-NL18 · Measurement scheme changing the gradient-resolution exponent
- **Paper:** A
- **Candidate:** C3
- **Matrix rows:** 18
- **Claim category:** Readout-dependent exponent of the resolution shot count
- **Classification:** NOT LOCATED
- **Section:** Whole paper (main text + Supplementary Information)
- **Subsection:** Nearest content: SI Note III D (A-12a)
- **Equation:** Nearest content: SI (143)–(144)
- **Figure:** —
- **Appendix:** —
- **Page:** main 1–13; SI 1–51
- **Source version:** Published version + Supplementary Information (MOESM1)
- **Source URL:** https://doi.org/10.1038/s41467-024-49287-w
- **Short quote:** —
- **Source statement (paraphrase):** Not located in the reviewed version after searches for: exponent, measurement strategy, readout, POVM, shot scaling, resolution. The source's resolution bound (SI Supp. Prop. 5) is range-based and identical for both tests, and no statement ties the exponent to the measurement scheme.
- **Mathematical expression:** —
- **Assumptions:** —
- **Scope:** As above.
- **Relation to Stage 7:** —
- **Does NOT establish:** Absence from the wider literature.

<a id="A-NL20"></a>
### A-NL20 · Estimated/exact gradient cosine similarity
- **Paper:** A
- **Candidate:** C4
- **Matrix rows:** 20
- **Claim category:** Cosine similarity of estimated vs exact gradient vectors
- **Classification:** NOT LOCATED
- **Section:** Whole paper (main text + Supplementary Information)
- **Subsection:** —
- **Equation:** —
- **Figure:** —
- **Appendix:** —
- **Page:** main 1–13; SI 1–51
- **Source version:** Published version + Supplementary Information (MOESM1)
- **Source URL:** https://doi.org/10.1038/s41467-024-49287-w
- **Short quote:** —
- **Source statement (paraphrase):** Not located in the reviewed version after searches for: cosine, angle, alignment, inner product, gradient vector. "Angle" occurs only for rotation angles of the encodings.
- **Mathematical expression:** —
- **Assumptions:** —
- **Scope:** As above.
- **Relation to Stage 7:** —
- **Does NOT establish:** Absence from the wider literature.

<a id="A-NL21"></a>
### A-NL21 · Gradient norm inflation
- **Paper:** A
- **Candidate:** C4
- **Matrix rows:** 21
- **Claim category:** Norm of estimated vs exact gradient
- **Classification:** NOT LOCATED
- **Section:** Whole paper (main text + Supplementary Information)
- **Subsection:** —
- **Equation:** —
- **Figure:** —
- **Appendix:** —
- **Page:** main 1–13; SI 1–51
- **Source version:** Published version + Supplementary Information (MOESM1)
- **Source URL:** https://doi.org/10.1038/s41467-024-49287-w
- **Short quote:** —
- **Source statement (paraphrase):** Not located in the reviewed version after searches for: norm (hits are Schatten/1-norms of states or distributions, the feature-space norm ‖a‖_H of the kernel-ridge regulariser (main p. 3) and operator norms ‖O‖_∞ (SI Supp. Prop. 5)), gradient norm, magnitude.
- **Mathematical expression:** —
- **Assumptions:** —
- **Scope:** As above.
- **Relation to Stage 7:** —
- **Does NOT establish:** Absence from the wider literature.

<a id="A-NL22"></a>
### A-NL22 · Probability that ĝ·g > 0
- **Paper:** A
- **Candidate:** C4
- **Matrix rows:** 22
- **Claim category:** Probability that the estimated gradient is a descent direction
- **Classification:** NOT LOCATED
- **Section:** Whole paper (main text + Supplementary Information)
- **Subsection:** —
- **Equation:** —
- **Figure:** —
- **Appendix:** —
- **Page:** main 1–13; SI 1–51
- **Source version:** Published version + Supplementary Information (MOESM1)
- **Source URL:** https://doi.org/10.1038/s41467-024-49287-w
- **Short quote:** —
- **Source statement (paraphrase):** Not located in the reviewed version after searches for: inner product, dot product, descent direction, sign.
- **Mathematical expression:** —
- **Assumptions:** —
- **Scope:** As above.
- **Relation to Stage 7:** —
- **Does NOT establish:** Absence from the wider literature.

<a id="A-NL23"></a>
### A-NL23 · Component sign accuracy
- **Paper:** A
- **Candidate:** C4
- **Matrix rows:** 23
- **Claim category:** Fraction of gradient components with correct sign
- **Classification:** NOT LOCATED
- **Section:** Whole paper (main text + Supplementary Information)
- **Subsection:** —
- **Equation:** —
- **Figure:** —
- **Appendix:** —
- **Page:** main 1–13; SI 1–51
- **Source version:** Published version + Supplementary Information (MOESM1)
- **Source URL:** https://doi.org/10.1038/s41467-024-49287-w
- **Short quote:** —
- **Source statement (paraphrase):** Not located in the reviewed version after searches for: sign, component, sign accuracy.
- **Mathematical expression:** —
- **Assumptions:** —
- **Scope:** As above.
- **Relation to Stage 7:** —
- **Does NOT establish:** Absence from the wider literature.

<a id="A-NL24"></a>
### A-NL24 · Matched-start finite-shot optimisation trajectories
- **Paper:** A
- **Candidate:** C4
- **Matrix rows:** 24
- **Claim category:** Trajectories from identical starts under different readouts/shot budgets
- **Classification:** NOT LOCATED
- **Section:** Whole paper (main text + Supplementary Information)
- **Subsection:** —
- **Equation:** —
- **Figure:** —
- **Appendix:** —
- **Page:** main 1–13; SI 1–51
- **Source version:** Published version + Supplementary Information (MOESM1)
- **Source URL:** https://doi.org/10.1038/s41467-024-49287-w
- **Short quote:** —
- **Source statement (paraphrase):** Not located in the reviewed version after searches for: trajectory, iteration, initialization, paired, matched. Kernel ridge regression is solved in closed form. Fig. 9 reports variances over 500 random initialisations, not trajectories.
- **Mathematical expression:** —
- **Assumptions:** —
- **Scope:** As above.
- **Relation to Stage 7:** —
- **Does NOT establish:** Absence from the wider literature.

<a id="A-NL28"></a>
### A-NL28 · Global-Z (parity) parameter-shift training numerics on a single rotation layer
- **Paper:** A
- **Candidate:** CTX
- **Matrix rows:** 28
- **Claim category:** Finite-shot training numerics with a global Pauli-Z cost on a single-qubit-rotation layer
- **Classification:** NOT LOCATED
- **Section:** Whole paper (main text + Supplementary Information)
- **Subsection:** —
- **Equation:** —
- **Figure:** —
- **Appendix:** —
- **Page:** main 1–13; SI 1–51
- **Source version:** Published version + Supplementary Information (MOESM1)
- **Source URL:** https://doi.org/10.1038/s41467-024-49287-w
- **Short quote:** —
- **Source statement (paraphrase):** Not located in the reviewed version after searches for: Pauli-Z, parity, global observable, training curve. Single rotation layers appear only as kernel embeddings (Prop. 3, Fig. 7, main Fig. 3, SI Fig. 4) and as the trainable R_y layer of Fig. 9 (Var_θ[TA] only, no training curves).
- **Mathematical expression:** —
- **Assumptions:** —
- **Scope:** As above.
- **Relation to Stage 7:** —
- **Does NOT establish:** Absence from the wider literature.

---

## Paper B — Aghaei Saem, Tafreshi, Holmes, Thanasilp, QST 11, 015049 (2026)

All B locators: **arXiv:2507.22054v2** (dated 2026-06-04). The published (IOP) version was not inspected.
"v1" = arXiv:2507.22054v1 (2025-07-29).

<a id="B-01"></a>
### B-01 · §II Framework: procedure P, Eqs. (1)–(4), (6)
- **Paper:** B
- **Candidate:** C1; C3
- **Matrix rows:** 7
- **Claim category:** Parameter-shift gradient descent with finite-shot loss estimates, inside the general procedure
- **Classification:** EXPLICIT
- **Section:** II. Framework
- **Subsection:** Gradient-based and non-gradient based training
- **Equation:** (1)–(4), (6)
- **Figure:** Fig. 2
- **Appendix:** App. B Eq. (B17) (formal re-statement)
- **Page:** 3–4 (arXiv v2)
- **Source version:** arXiv:2507.22054v2 (2026-06-04); same equation numbers in v1
- **Source URL:** https://arxiv.org/abs/2507.22054v2
- **Short quote:** —
- **Source statement (paraphrase):** Each quantity ℓ_i is estimated from N POVM outcomes (outcome probabilities p_k = Tr[ρ M_k]). For vanilla gradient descent with the parameter-shift rule, the k-th update uses the difference of estimated losses at θ ± (π/2)ê_k. For a single Pauli observable this needs N_ℓ = 2N_p estimated quantities.
- **Mathematical expression:** [Φ_P]_k = θ_k − (η/2)[L̂(θ + (π/2)ê_k) − L̂(θ − (π/2)ê_k)]
- **Assumptions:** Circuits obeying the parameter-shift rule; polynomially many quantities.
- **Scope:** General loss functions; the finite-shot statistics are analysed only through the indistinguishability results (B-04–B-06).
- **Relation to Stage 7:** Stage 7's estimator ĝ = [Ĉ(θ + π/2 e_k) − Ĉ(θ − π/2 e_k)]/2 with independent batches is this update's gradient term with an infidelity loss.
- **Does NOT establish:** The distribution, variance, zero or sign probabilities of the parameter-shift estimate for any specific readout.

<a id="B-02"></a>
### B-02 · §II "Polynomial POVMs in disguise": global Pauli-Z as a two-outcome parity POVM
- **Paper:** B
- **Candidate:** CTX
- **Matrix rows:** —
- **Claim category:** Global Pauli-Z estimation uses the two-element parity POVM {Π₊, Π₋}
- **Classification:** EXPLICIT
- **Section:** II. Framework
- **Subsection:** Polynomial POVMs in disguise
- **Equation:** —
- **Figure:** —
- **Appendix:** —
- **Page:** 5 (arXiv v2)
- **Source version:** arXiv:2507.22054v2 (2026-06-04); the global-Z paragraph is absent from v1
- **Source URL:** https://arxiv.org/abs/2507.22054v2
- **Short quote:** —
- **Source statement (paraphrase):** Estimating ⟨⊗Z_i⟩ requires only the eigenvalue sector of each bitstring, so the relevant POVM is {Π₊, Π₋}: projectors onto even- and odd-parity subspaces.
- **Mathematical expression:** M = {Π₊, Π₋}
- **Assumptions:** —
- **Scope:** POVM bookkeeping.
- **Relation to Stage 7:** Matches the two-outcome parity estimator of the Stage 6 parity benchmark.
- **Does NOT establish:** Finite-shot gradient statistics for that POVM.

<a id="B-03"></a>
### B-03 · §III Definition 1 (outcome probability concentration), Eq. (8)
- **Paper:** B
- **Candidate:** F-A; F-F
- **Matrix rows:** 1, 6
- **Claim category:** Concentration at the level of POVM outcome probabilities
- **Classification:** EXPLICIT
- **Section:** III. Exponential concentration
- **Subsection:** Definition 1
- **Equation:** (8)
- **Figure:** Fig. 1
- **Appendix:** App. B (used throughout)
- **Page:** 5–7; Fig. 1 on p. 2 (arXiv v2)
- **Source version:** arXiv:2507.22054v2 (2026-06-04); v1 has Definition 1/Eq. (8) without the explicit "drawn from D" qualifier and the following remark
- **Source URL:** https://arxiv.org/abs/2507.22054v2
- **Short quote:** —
- **Source statement (paraphrase):** A POVM's outcome probabilities concentrate with respect to α ~ D if, for every element, Pr(|p_k(α) − μ_k| ≥ δ) ≤ β/δ² with β ∈ O(exp(−n)) and α-independent μ_k. The mechanisms are those of expectation-value concentration, since p_k is the expectation of a POVM element.
- **Mathematical expression:** Pr_{α~D}(|p_k(α) − μ_k| ≥ δ) ≤ β/δ², β ∈ O(exp(−n))
- **Assumptions:** Distribution D over the variables.
- **Scope:** Any POVM.
- **Relation to Stage 7:** Stage 5–7 analyse exactly these outcome probabilities (F±, q±) at fixed θ and their distribution over θ.
- **Does NOT establish:** Finite-shot gradient laws or exponents.

<a id="B-04"></a>
### B-04 · Theorem 1 (informal) and Theorem 2 (formal), Eqs. (B1)–(B13)
- **Paper:** B
- **Candidate:** F-F; F-E
- **Matrix rows:** 5, 6
- **Claim category:** Concentrated polynomial-size POVM ⇒ samples indistinguishable from a fixed distribution under polynomial shots
- **Classification:** EXPLICIT
- **Section:** III. Exponential concentration; App. B
- **Subsection:** Theorem 1; Theorem 2
- **Equation:** (B1)–(B13)
- **Figure:** Fig. 1
- **Appendix:** App. B
- **Page:** 7; 19–20; Fig. 1 on p. 2 (arXiv v2)
- **Source version:** arXiv:2507.22054v2 (2026-06-04); same statement and (B)-numbering in v1
- **Source URL:** https://arxiv.org/abs/2507.22054v2
- **Short quote:** —
- **Source statement (paraphrase):** With |M| ∈ O(poly(n)) and N ∈ O(poly(n)), with probability ≥ 1 − δ over α the samples from P_α cannot be told apart from samples of P_fixed with success above 1/2 + ε. Here δ = |M|√β and ε = N|M|β^(1/4)/4 are exponentially small.
- **Mathematical expression:** Pr[right decision] ≤ 1/2 + N|M|β^(1/4)/4 with prob. ≥ 1 − |M|√β
- **Assumptions:** Definition 1 for all outcomes; polynomial POVM and shots; 1-norm hypothesis-testing bound (App. A).
- **Scope:** Any procedure with polynomial POVMs.
- **Relation to Stage 7:** The general indistinguishability fact; Stage 7 does not claim it.
- **Does NOT establish:** Readout-specific exponents (the bound is 1-norm based), gradient zero/sign laws.

<a id="B-05"></a>
### B-05 · Corollary 1 (informal, Eq. 9) and Corollary 3 (formal)
- **Paper:** B
- **Candidate:** F-H
- **Matrix rows:** 26
- **Claim category:** No classical post-processing removes indistinguishability
- **Classification:** EXPLICIT
- **Section:** III. Exponential concentration; App. B
- **Subsection:** Corollary 1; Corollary 3
- **Equation:** (9)
- **Figure:** Fig. 1
- **Appendix:** App. B
- **Page:** 7; 20; Fig. 1 on p. 2 (arXiv v2)
- **Source version:** arXiv:2507.22054v2 (2026-06-04); same in v1
- **Source URL:** https://arxiv.org/abs/2507.22054v2
- **Short quote:** —
- **Source statement (paraphrase):** Applying any map Φ′ to concentrated measurement outcomes yields, with high probability, an estimate statistically indistinguishable from the α-independent random variable ℓ̂_fixed = Φ′(S_N,fixed).
- **Mathematical expression:** ℓ̂_fixed = Φ′(S_N,fixed)
- **Assumptions:** As Theorem 2.
- **Scope:** Arbitrary post-processing.
- **Relation to Stage 7:** A fact Stage 7 does not claim. Stage 7's gradient estimators are one such map.
- **Does NOT establish:** The specific form of the fixed-distribution gradient laws.

<a id="B-06a"></a>
### B-06a · Corollary 2 (Eq. 10), Corollary 4 (Eqs. B14–B15) and proof (B16)–(B30): random walk
- **Paper:** B
- **Candidate:** C4; F-G
- **Matrix rows:** 8
- **Claim category:** Polynomial-shot parameter-shift gradient descent on a concentrated loss is statistically a random walk
- **Classification:** EXPLICIT
- **Section:** III. Exponential concentration; App. B
- **Subsection:** Corollary 2; Corollary 4 and proof
- **Equation:** (10); (B14)–(B30)
- **Figure:** Fig. 3
- **Appendix:** App. B
- **Page:** 7; 20–23; Fig. 3 on p. 6 (arXiv v2)
- **Source version:** arXiv:2507.22054v2 (2026-06-04); v1 has the same corollaries with wording "results in a random walk" (see Version differences V-05)
- **Source URL:** https://arxiv.org/abs/2507.22054v2
- **Short quote:** —
- **Source statement (paraphrase):** Setting: loss Σ_i c_i Tr[ρ(θ)O_i] with Pauli O_i, parameter-shift rule, random initialisation, polynomially many iterations and shots. With probability ≥ 1 − c, c ∈ O(exp(−n)), each update θ_new = θ_current + Δ_N is indistinguishable from one with parameter-independent components [Δ_N]_k = −(η/2) Σ_i c_i [Σ_j z_ijk/N − Σ_j z′_ijk/N], with z, z′ = ±1 equiprobable. A union bound over steps extends this to the whole trajectory.
- **Mathematical expression:** [Δ_N]_k = −(η/2) Σ_i c_i [Σ_{j=1}^N z_ijk/N − Σ_{j=1}^N z′_ijk/N], z, z′ ∈ {±1} equiprobable
- **Assumptions:** Pauli observables (two-outcome POVMs, ±1 outcomes); every POVM exponentially concentrated; N, N_ℓ, N_step ∈ O(poly(n)); union bound.
- **Scope:** Concentration point 1/2 for both outcomes (Pauli/SWAP-like). The Loschmidt projector POVM, whose fixed distribution is (0, 1) (B-09a), is not the case written in (B15).
- **Relation to Stage 7:** Stage 7 §17–18 lists SWAP-driven GD at n ≥ 10 matching a signal-free random walk as a replication of this corollary.
- **Does NOT establish:** Zero/sign probabilities, alignment, norm or trajectory metrics (B-06b–B-06e), or the Loschmidt "frozen" behaviour.

<a id="B-06b"></a>
### B-06b · Eqs. (B18), (B26): update written as a difference of two empirical means of ±1 outcomes
- **Paper:** B
- **Candidate:** C1
- **Matrix rows:** 12
- **Claim category:** Difference-of-binomials form of the finite-shot parameter-shift estimate
- **Classification:** PARTIAL / RELATED
- **Section:** App. B
- **Subsection:** Proof of Corollary 4, step (1)
- **Equation:** (B18), (B19), (B26)
- **Figure:** —
- **Appendix:** App. B
- **Page:** 22–23 (arXiv v2)
- **Source version:** arXiv:2507.22054v2 (2026-06-04); same numbering in v1
- **Source URL:** https://arxiv.org/abs/2507.22054v2
- **Short quote:** —
- **Source statement (paraphrase):** Each update component is −(η/2) Σ_i c_i [Σ_q λ_i1kq/N − Σ_q λ_i2kq/N], with λ = +1 w.p. (1 + ℓ_i(α))/2 and −1 otherwise. Each sum is an affine function of a binomial count (our rewriting; the word "binomial" does not appear in the source, which uses only ±1-outcome sums). In the concentrated limit both are replaced by z = ±1 equiprobable (B26). The distribution of the difference is not derived.
- **Mathematical expression:** (our rewriting) Σ_q λ_q/N = 2K/N − 1, K ~ Bin(N, (1 + ℓ)/2)
- **Assumptions:** As B-06a.
- **Scope:** Pauli ±1 outcomes.
- **Relation to Stage 7:** The SWAP gradient estimator ĝ_SWAP = (K₋ − K₊)/M has exactly this structure with ℓ = F (ancilla ⟨Z⟩ = F). Stage 5–7 derive its exact pmf-based P_zero/P_correct.
- **Does NOT establish:** The pmf of the difference, P(ĝ = 0), P_correct/P_wrong, tie probability, or the Loschmidt (0/1-outcome) version.

<a id="B-06c"></a>
### B-06c · Corollaries 2/4 read against gradient sign probabilities
- **Paper:** B
- **Candidate:** C1
- **Matrix rows:** 10, 11
- **Claim category:** Correct/wrong-sign probabilities of the finite-shot gradient
- **Classification:** PARTIAL / RELATED
- **Section:** III. Exponential concentration; App. B
- **Subsection:** Corollary 2; Corollary 4
- **Equation:** (10), (B15), (B26)
- **Figure:** —
- **Appendix:** App. B
- **Page:** 7; 21–23 (arXiv v2)
- **Source version:** arXiv:2507.22054v2 (2026-06-04)
- **Source URL:** https://arxiv.org/abs/2507.22054v2
- **Short quote:** —
- **Source statement (paraphrase):** The update is indistinguishable from a parameter-independent variable that is symmetric under sign flip. This implies the sign carries no usable information about the landscape with high probability, but no probability of the correct or wrong sign is stated or computed.
- **Mathematical expression:** —
- **Assumptions:** As B-06a.
- **Scope:** Pauli/±1 outcomes; asymptotic indistinguishability statement.
- **Relation to Stage 7:** Stage 6–7 give exact finite-M P_correct, P_wrong and P_zero (e.g. SWAP P_correct → (1 − P_zero)/2; Loschmidt P_correct → 0 through zeros). Neither the values nor the Loschmidt asymmetry appear in the source.
- **Does NOT establish:** Exact or asymptotic P_correct/P_wrong values; their M- or n-dependence.

<a id="B-06d"></a>
### B-06d · Corollaries 2/4 read against full-vector direction metrics
- **Paper:** B
- **Candidate:** C4
- **Matrix rows:** 20, 22, 23
- **Claim category:** Direction of the estimated gradient vector (cosine, ĝ·g > 0, component sign accuracy)
- **Classification:** PARTIAL / RELATED
- **Section:** III. Exponential concentration; App. B
- **Subsection:** Corollary 2; Corollary 4
- **Equation:** (10), (B14)–(B15)
- **Figure:** Fig. 3
- **Appendix:** App. B
- **Page:** 7; 21; Fig. 3 on p. 6 (arXiv v2)
- **Source version:** arXiv:2507.22054v2 (2026-06-04)
- **Source URL:** https://arxiv.org/abs/2507.22054v2
- **Short quote:** —
- **Source statement (paraphrase):** The whole update vector is indistinguishable from a parameter-independent random vector. The source states this at the level of trajectories (random walk); it does not compute any alignment measure (cosine, ĝ·g sign, fraction of correct component signs).
- **Mathematical expression:** —
- **Assumptions:** As B-06a.
- **Scope:** Trajectory-level statement.
- **Relation to Stage 7:** Stage 7 §14 reports these metrics directly: SWAP median cosine ≈ 0.009, P(ĝ·g > 0) ≈ 0.51; Loschmidt conditional cosine ≈ 0.74 at n = 12, M = 1024. "Optimizer behaves like a random walk" and "cosine ≈ 0" are related but not identical statements.
- **Does NOT establish:** Any alignment value, the Loschmidt "zero-or-aligned" pattern, or finite-n behaviour.

<a id="B-06e"></a>
### B-06e · Corollary 4 read against gradient-norm inflation
- **Paper:** B
- **Candidate:** C4
- **Matrix rows:** 21
- **Claim category:** Norm of the estimated gradient relative to the exact gradient
- **Classification:** PARTIAL / RELATED
- **Section:** App. B
- **Subsection:** Corollary 4
- **Equation:** (B15)
- **Figure:** —
- **Appendix:** App. B
- **Page:** 21 (arXiv v2)
- **Source version:** arXiv:2507.22054v2 (2026-06-04)
- **Source URL:** https://arxiv.org/abs/2507.22054v2
- **Short quote:** —
- **Source statement (paraphrase):** The null update (B15) has a parameter-independent magnitude set by η, N and the coefficients c_i. The source makes no statement comparing the estimated-gradient norm with the exact-gradient norm.
- **Mathematical expression:** (implied, not stated) E‖Δ_N‖² = (η²/4) Σ_k Σ_i c_i² (2/N)
- **Assumptions:** As B-06a; the expectation shown is our algebra on (B15) (B5_MATH_COMPARISON §6.6).
- **Scope:** Null (concentrated) update only.
- **Relation to Stage 7:** Stage 7 §14 measures SWAP median log₁₀‖ĝ‖/‖g‖ ≈ 4.5 at n = 12, M = 1024.
- **Does NOT establish:** Any norm-ratio statement or value.

<a id="B-07a"></a>
### B-07a · Fig. 3 and §III text: 15-qubit training with 150 / 2¹⁵ / infinite shots
- **Paper:** B
- **Candidate:** C4; F-G
- **Matrix rows:** 8, 28
- **Claim category:** Numerical random-walk behaviour of finite-shot training on a barren plateau
- **Classification:** EXPLICIT
- **Section:** III. Exponential concentration
- **Subsection:** Text following Corollary 2
- **Equation:** —
- **Figure:** Fig. 3 (a)–(c)
- **Appendix:** App. C (setup)
- **Page:** 6–7 (arXiv v2)
- **Source version:** arXiv:2507.22054v2 (2026-06-04); v1 has the figure, but its caption lacks the definition of the plotted quantity
- **Source URL:** https://arxiv.org/abs/2507.22054v2
- **Short quote:** —
- **Source statement (paraphrase):** For n = 15, a single layer of X rotations and the global Pauli-Z observable: trajectories with 2¹⁵ and infinite shots converge, while trajectories with 150 shots wander randomly. The mean and variance of (1/N_p)‖θ^(t) − θ^(0)‖₁ over initialisations under 150 shots closely resemble a random walk.
- **Mathematical expression:** Plotted statistic: (1/N_p)‖θ^(t) − θ^(0)‖₁ (mean and variance over initialisations)
- **Assumptions:** Global-Z cost; "random initialization" (Cor. 2; the Fig. 3 caption says "different parameter initializations"); learning rate η (value and schedule not stated in the reviewed text).
- **Scope:** One system size, one cost (parity), displacement statistics.
- **Relation to Stage 7:** Same circuit family as Stage 3–7 but the parity cost (Stage 6's parity benchmark). Stage 7's random-walk diagnostic uses the fidelity cost with the SWAP readout.
- **Does NOT establish:** Fidelity/SWAP or Loschmidt trajectories, endpoint-fidelity comparisons, alignment metrics, or start matching (B-07c).

<a id="B-07b"></a>
### B-07b · Fig. 3 (b)–(c): "Random Walk" reference curve
- **Paper:** B
- **Candidate:** C4
- **Matrix rows:** 25
- **Claim category:** Signal-free random-walk control
- **Classification:** EXPLICIT
- **Section:** III. Exponential concentration
- **Subsection:** Text following Corollary 2
- **Equation:** —
- **Figure:** Fig. 3 (b), (c) — legend "Random Walk"
- **Appendix:** —
- **Page:** 6–7 (arXiv v2)
- **Source version:** arXiv:2507.22054v2 (2026-06-04)
- **Source URL:** https://arxiv.org/abs/2507.22054v2
- **Short quote:** —
- **Source statement (paraphrase):** A random-walk reference is plotted alongside the 150-shot and 2¹⁵-shot runs for the mean and variance of parameter displacement. How the reference curve is generated is not specified in the reviewed text.
- **Mathematical expression:** —
- **Assumptions:** As B-07a.
- **Scope:** Displacement statistics only; parity cost; n = 15.
- **Relation to Stage 7:** Stage 7's control replaces SWAP counts with zero-signal counts (q₊ = q₋ = 1/2) under the identical update rule and identical starts, and compares fidelity endpoints and alignment. The source's comparison is on displacement statistics.
- **Does NOT establish:** Endpoint (loss/fidelity) comparison with a control, matched starts, or the construction of the reference.

<a id="B-07c"></a>
### B-07c · Fig. 3 (a) and Fig. 4: trajectories at different shot budgets — start matching not stated
- **Paper:** B
- **Candidate:** C4
- **Matrix rows:** 24
- **Claim category:** Matched-start finite-shot optimisation trajectories
- **Classification:** PARTIAL / RELATED
- **Section:** III. Exponential concentration; IV. Practical step-by-step guidelines
- **Subsection:** Text following Corollary 2; discussion of Fig. 4
- **Equation:** —
- **Figure:** Fig. 3 (a); Fig. 4
- **Appendix:** App. C
- **Page:** 6–8 (arXiv v2)
- **Source version:** arXiv:2507.22054v2 (2026-06-04)
- **Source URL:** https://arxiv.org/abs/2507.22054v2
- **Short quote:** —
- **Source statement (paraphrase):** Trajectories under different shot budgets are shown together (Fig. 3a, a PCA cut) and training curves are compared across budgets (Fig. 4). Neither the text nor the captions say whether initial parameters are identical across budgets. Fig. 3 (b)–(c) average over different initialisations.
- **Mathematical expression:** —
- **Assumptions:** —
- **Scope:** Parity / global-Z costs.
- **Relation to Stage 7:** Stage 7 pairs every arm (exact, LE, SWAP, random walk) from identical θ₀ by construction (Stage 6 found an unpaired comparison gave a false positive).
- **Does NOT establish:** That starts were matched; any paired statistic.

<a id="B-08a"></a>
### B-08a · Fig. 4, §IV text and App. C: training curves at 10n vs 2ⁿ vs infinite shots
- **Paper:** B
- **Candidate:** F-E; F-G
- **Matrix rows:** 5
- **Claim category:** Exponential shot budget needed for training on a barren plateau (numerics)
- **Classification:** EXPLICIT
- **Section:** IV. Practical step-by-step guidelines
- **Subsection:** Discussion of Fig. 4
- **Equation:** (C1)–(C11)
- **Figure:** Fig. 4 (a)–(d)
- **Appendix:** App. C
- **Page:** 6, 8; 23–25 (arXiv v2)
- **Source version:** arXiv:2507.22054v2 (2026-06-04); also in v1
- **Source URL:** https://arxiv.org/abs/2507.22054v2
- **Short quote:** —
- **Source statement (paraphrase):** For quantum natural gradient, CVaR, classical-NN initialisation and the rescaled parameter-shift rule (n = 9–17), training succeeds with 2ⁿ shots but does not move meaningfully with 10 × n shots.
- **Mathematical expression:** —
- **Assumptions:** Single X-rotation layer (C1); global-Z observables (C8) or a sum of four global-Z terms (C9); uniform random initialisation for QNG, CVaR and the rescaled parameter-shift rule; classical-NN-output initialisation for the NN method.
- **Scope:** Global costs; numerical.
- **Relation to Stage 7:** General exponential-burden fact; not readout-specific.
- **Does NOT establish:** Readout-dependent exponents or gradient-level laws.

<a id="B-08b"></a>
### B-08b · App. C, Eqs. (C1), (C8), (C9): single X-rotation layer with global Pauli-Z costs
- **Paper:** B
- **Candidate:** CTX
- **Matrix rows:** 28
- **Claim category:** Numerical setup on the single-qubit-rotation product circuit
- **Classification:** EXPLICIT
- **Section:** App. C
- **Subsection:** Further details of numerical simulation
- **Equation:** (C1), (C8), (C9)
- **Figure:** Fig. 3, Fig. 4
- **Appendix:** App. C
- **Page:** 23–25; Figs. 3–4 on p. 6 (arXiv v2)
- **Source version:** arXiv:2507.22054v2 (2026-06-04); same in v1
- **Source URL:** https://arxiv.org/abs/2507.22054v2
- **Short quote:** —
- **Source statement (paraphrase):** All numerics use U(θ) = ∏_i e^(−iθ_i X_i) with global observables: ⊗Z_i (QNG, NN initialisation, rescaled PS) or a sum of four global-Z strings (CVaR). Parameters are initialised uniformly at random for QNG, CVaR and the rescaled parameter-shift rule; the classical-NN method initialises them from a neural network's output.
- **Mathematical expression:** U(θ) = ∏_{i=1}^n e^(−iθ_i X_i); H = ⊗_i Z_i
- **Assumptions:** —
- **Scope:** Parity-type costs on the product circuit, not the |0ⁿ⟩ fidelity.
- **Relation to Stage 7:** Same circuit family as Stage 3–7, different cost. This is Stage 6's parity benchmark setup (STAGE6.md §19 notes it). Stage 7's fidelity cost and the Loschmidt/SWAP readouts are not used in the source's numerics.
- **Does NOT establish:** Fidelity-cost results on this circuit.

<a id="B-08c"></a>
### B-08c · App. C setup read against the Stage 3–7 fidelity landscape
- **Paper:** B
- **Candidate:** CTX
- **Matrix rows:** 27
- **Claim category:** Concentration of the product single-qubit-rotation fidelity landscape
- **Classification:** PARTIAL / RELATED
- **Section:** App. C; IV. Practical step-by-step guidelines
- **Subsection:** Further details of numerical simulation; Subtlety regarding the choice of POVM
- **Equation:** (C1), (C8)
- **Figure:** Fig. 3, Fig. 4
- **Appendix:** App. C
- **Page:** 9–10, 23–25; Figs. 3–4 on p. 6 (arXiv v2)
- **Source version:** arXiv:2507.22054v2 (2026-06-04)
- **Source URL:** https://arxiv.org/abs/2507.22054v2
- **Short quote:** —
- **Source statement (paraphrase):** The source's numerics use the same product circuit (one layer of single-qubit X rotations; random initialisation, uniform for QNG, CVaR and rescaled parameter shift), but with global Pauli-Z costs, citing the globality-induced barren plateau. Its fidelity discussion (Loschmidt vs SWAP) assumes a 2-design, not the product family. The concentration of the product-circuit |0ⁿ⟩ fidelity is not stated.
- **Mathematical expression:** —
- **Assumptions:** —
- **Scope:** Same circuit, different cost; or same cost, different (2-design) ensemble.
- **Relation to Stage 7:** Paper A's Proposition 3 (A-03) is the located explicit statement for the product-rotation fidelity family.
- **Does NOT establish:** Concentration statistics of the product-circuit fidelity or its gradients.

<a id="B-09a"></a>
### B-09a · §IV "Subtlety regarding the choice of POVM", Fig. 5: Loschmidt vs SWAP fixed distributions
- **Paper:** B
- **Candidate:** F-B; F-C; C3
- **Matrix rows:** 1, 2, 3, 14a
- **Claim category:** Same fidelity objective, two POVMs: zero estimate (LE) vs 50/50 outcomes (SWAP)
- **Classification:** EXPLICIT
- **Section:** IV. Practical step-by-step guidelines
- **Subsection:** Subtlety regarding the choice of POVM
- **Equation:** —
- **Figure:** Fig. 5 (a) Loschmidt echo test, (b) SWAP test
- **Appendix:** App. B (Theorem 2 invoked)
- **Page:** 9 (arXiv v2)
- **Source version:** arXiv:2507.22054v2 (2026-06-04); **absent from v1**; published-version status unverified
- **Source URL:** https://arxiv.org/abs/2507.22054v2
- **Short quote:** —
- **Source statement (paraphrase):** For learning |φ⟩ = U₀|0⟩ with infidelity loss and a variational state forming a 2-design: the Loschmidt POVM {|0⟩⟨0|^⊗n, 𝟙 − |0⟩⟨0|^⊗n} has fixed distribution (0, 1), so with probability exponentially close to 1 the estimated fidelity is zero. The SWAP-test POVM {|0⟩⟨0|, |1⟩⟨1|} on the ancilla has fixed distribution (1/2, 1/2). Both approaches are said to suffer from concentration (the fidelity concentrates under the 2-design assumption).
- **Mathematical expression:** P_fixed^(LE) = (0, 1); P_fixed^(SWAP) = (1/2, 1/2)
- **Assumptions:** |ψ(θ)⟩ = U(θ)|0⟩ forms a unitary 2-design; Theorem 2.
- **Scope:** Fidelity/loss estimates; analytic, no numerics; no gradients.
- **Relation to Stage 7:** Stage 7 lists the zero-vs-50/50 fidelity-level contrast as prior work (§3, §19). The source states it for a 2-design ensemble; Stage 7's landscape is the product (non-2-design) family.
- **Does NOT establish:** Gradient-level comparison, exponents, zero/sign probabilities.

<a id="B-09b"></a>
### B-09b · §IV, Loschmidt fixed distribution (0, 1), read against gradient/vector zero probabilities
- **Paper:** B
- **Candidate:** C1; C4
- **Matrix rows:** 9, 19
- **Claim category:** Zero probability of the gradient component / full gradient vector
- **Classification:** PARTIAL / RELATED
- **Section:** IV. Practical step-by-step guidelines
- **Subsection:** Subtlety regarding the choice of POVM
- **Equation:** —
- **Figure:** Fig. 5 (a)
- **Appendix:** —
- **Page:** 9 (arXiv v2)
- **Source version:** arXiv:2507.22054v2 (2026-06-04); absent from v1
- **Source URL:** https://arxiv.org/abs/2507.22054v2
- **Short quote:** —
- **Source statement (paraphrase):** Only the fidelity estimate is said to be zero with high probability. Nothing is said about the parameter-shift difference of two such estimates or about all components simultaneously.
- **Mathematical expression:** —
- **Assumptions:** As B-09a.
- **Scope:** Fidelity level.
- **Relation to Stage 7:** Both shifted Loschmidt estimates being zero makes the gradient component zero. Stage 5–7 compute the exact P(ĝ = 0) (including equal nonzero counts) and the full-vector zero probability.
- **Does NOT establish:** Gradient-level or vector-level zero probabilities or their values.

<a id="B-10a"></a>
### B-10a · §IV: single-shot estimator variances Var^(SWAP) = 1 − F², Var^(LE) = F(1 − F)
- **Paper:** B
- **Candidate:** F-D; C3
- **Matrix rows:** 4
- **Claim category:** Different estimator variances for Loschmidt vs SWAP (same fidelity)
- **Classification:** EXPLICIT
- **Section:** IV. Practical step-by-step guidelines
- **Subsection:** Subtlety regarding the choice of POVM
- **Equation:** Inline, unnumbered (text after Eq. (12))
- **Figure:** —
- **Appendix:** —
- **Page:** 10 (arXiv v2)
- **Source version:** arXiv:2507.22054v2 (2026-06-04); **absent from v1**; published-version status unverified
- **Source URL:** https://arxiv.org/abs/2507.22054v2
- **Short quote:** "These distinct forms of variance result in different statistical behaviors"
- **Source statement (paraphrase):** The landscape variance is exponentially small regardless of the measurement scheme, but the estimator variances differ: 1 − F_θ² for the SWAP test and F_θ(1 − F_θ) for the Loschmidt echo test.
- **Mathematical expression:** Var^(SWAP)_ρ(θ)[ℓ̂(θ)] = 1 − F_θ²; Var^(LE)_ρ(θ)[ℓ̂(θ)] = F_θ(1 − F_θ)
- **Assumptions:** Per-shot variance (the source writes Var_ρ(θ)[ℓ̂(θ)] without a 1/N factor); pure variational state |ψ(θ)⟩ and target |φ⟩ as in the source's example; the landscape-variance clause relies on the stated assumption that |ψ(θ)⟩ = U(θ)|0⟩ forms a unitary 2-design (p. 9).
- **Scope:** Fidelity (loss) estimate, not its parameter-shift gradient.
- **Relation to Stage 7:** These are the per-shot variances Stage 7 uses (Var F̂ = F(1 − F)/M for LE; (1 − F²)/M for SWAP). Stage 7's gradient variances follow by the parameter-shift difference with independent batches (B5_MATH_COMPARISON §6.2).
- **Does NOT establish:** The gradient-estimator variances, SNRs, zero/sign probabilities, or the exponents they imply.

<a id="B-10b"></a>
### B-10b · §IV, Eq. (12): resolution ratio ε_N and the exponential-shot requirement
- **Paper:** B
- **Candidate:** F-E
- **Matrix rows:** 5
- **Claim category:** Exponential measurement burden via the estimator-variance / landscape-variance ratio
- **Classification:** EXPLICIT
- **Section:** IV. Practical step-by-step guidelines
- **Subsection:** Subtlety regarding the choice of POVM
- **Equation:** (12)
- **Figure:** —
- **Appendix:** —
- **Page:** 9–10 (arXiv v2)
- **Source version:** arXiv:2507.22054v2 (2026-06-04); absent from v1
- **Source URL:** https://arxiv.org/abs/2507.22054v2
- **Short quote:** —
- **Source statement (paraphrase):** Define ε_N(α) = Var_ρ(α)[ℓ̂(α)]/(N Var_α[ℓ(α)]). Sufficient resolution needs ε_N ≲ 1, which under exponential concentration requires N ∈ Ω(exp(n)). The source offers this as a rule-of-thumb diagnostic.
- **Mathematical expression:** ε_N(α) = Var_ρ(α)[ℓ̂(α)] / (N Var_α[ℓ(α)]); require ε_N ≲ 1 ⇒ N ∈ Ω(exp(n))
- **Assumptions:** Quantities that are expectation values; exponential concentration of Var_α[ℓ].
- **Scope:** General, applied to fidelity.
- **Relation to Stage 7:** Stage 7 uses a different, per-θ criterion (SNR ≥ ρ, P_correct ≥ q at fixed θ; median over θ). The two criteria give different exponents on the heavy-tailed product landscape (B5_MATH_COMPARISON §6.3).
- **Does NOT establish:** Readout-specific exponents (B-10c) or gradient-level criteria.

<a id="B-10c"></a>
### B-10c · §IV, Eq. (12) with the two variances: readout-dependent exponents implied but not computed
- **Paper:** B
- **Candidate:** C3
- **Matrix rows:** 15, 16, 18
- **Claim category:** Measurement-dependent shot-complexity exponent
- **Classification:** PARTIAL / IMPLIED
- **Section:** IV. Practical step-by-step guidelines
- **Subsection:** Subtlety regarding the choice of POVM
- **Equation:** (12) + inline variances
- **Figure:** —
- **Appendix:** —
- **Page:** 9–10 (arXiv v2)
- **Source version:** arXiv:2507.22054v2 (2026-06-04); absent from v1
- **Source URL:** https://arxiv.org/abs/2507.22054v2
- **Short quote:** "affect the degree to which the problem suffers from shot noise"
- **Source statement (paraphrase):** Changing the POVM does not change the landscape variance and is not expected to remove concentration. Different estimators instead change how badly shot noise affects the problem. The source supplies ε_N and the two variances but does not combine them into shot counts or exponents for either readout.
- **Mathematical expression:** (our algebra, not in source) 2-design: N_LE ~ F/Var_α ~ 2ⁿ, N_SWAP ~ 1/Var_α ~ 4ⁿ (fidelity level); per-point SNR: N_SWAP/N_LE ≈ 1/F (fidelity) or ≈ 2/A (gradient)
- **Assumptions:** For the implication: 2-design statistics (E F = 2^(−n), Var_α F ≈ 2^(−2n)) or a per-point SNR criterion; independent shots; for the gradient version, parameter shift with independent batches and the Stage 5 log-typical A.
- **Scope:** The source is at the fidelity level under a 2-design. Stage 7 is at the gradient level on the product landscape.
- **Relation to Stage 7:** Algebra in B5_MATH_COMPARISON §6.2–§6.3. The implication holds at the fidelity level under the source's 2-design assumption (2ⁿ vs 4ⁿ). Reaching Stage 7's 4ⁿ vs 16ⁿ additionally needs the parameter-shift step, a per-θ criterion, and the product-landscape log-typical A. With the source's own ε_N criterion on the product landscape, the result differs: (2/3)ⁿ at the log-typical F, or (4/3)ⁿ at the mean F, vs (8/3)ⁿ.
- **Does NOT establish:** Any exponent value stated by the authors; gradient-level statements; same-θ paired comparison.

<a id="B-11"></a>
### B-11 · §IV: purity example and the open question about entangled POVMs
- **Paper:** B
- **Candidate:** CTX; F-D
- **Matrix rows:** —
- **Claim category:** POVM choice and multi-copy measurements under concentration
- **Classification:** EXPLICIT
- **Section:** IV. Practical step-by-step guidelines
- **Subsection:** Subtlety regarding the choice of POVM
- **Equation:** —
- **Figure:** —
- **Appendix:** —
- **Page:** 10 (arXiv v2)
- **Source version:** arXiv:2507.22054v2 (2026-06-04); absent from v1
- **Source URL:** https://arxiv.org/abs/2507.22054v2
- **Short quote:** —
- **Source statement (paraphrase):** Two-copy SWAP measurements estimate purity efficiently without concentration. Under a 2-design the landscape variance stays exponentially small, so the advantage vanishes. Whether some entangled POVM could change the indistinguishability conclusion is left open; analysis must use the POVMs employed in practice.
- **Mathematical expression:** —
- **Assumptions:** 2-design.
- **Scope:** Conceptual.
- **Relation to Stage 7:** Positioning context for the readout-dependence question (A4 / B2 belong to the principal and later tracks).
- **Does NOT establish:** Readout exponents.

<a id="B-12"></a>
### B-12 · §IV "Subtlety regarding measure-first-estimate-later approaches", Eq. (11)
- **Paper:** B
- **Candidate:** F-H
- **Matrix rows:** 26
- **Claim category:** Classical-shadow / measure-first schemes are POVM post-processing
- **Classification:** EXPLICIT
- **Section:** IV. Practical step-by-step guidelines
- **Subsection:** Subtlety regarding measure-first-estimate-later approaches
- **Equation:** (11)
- **Figure:** —
- **Appendix:** —
- **Page:** 8–9 (arXiv v2)
- **Source version:** arXiv:2507.22054v2 (2026-06-04); **absent from v1**
- **Source URL:** https://arxiv.org/abs/2507.22054v2
- **Short quote:** —
- **Source statement (paraphrase):** Post-processing a classical representation built from measurements (e.g. ρ̂_Z) to evaluate an observable is equivalent to implementing that observable's POVM. Concentration results therefore apply however the POVM is realised.
- **Mathematical expression:** ρ̂_Z = (1/N) Σ_i |b_i⟩⟨b_i|
- **Assumptions:** —
- **Scope:** Classical shadows / measure-first schemes.
- **Relation to Stage 7:** Background (post-processing limits).
- **Does NOT establish:** —

<a id="B-13"></a>
### B-13 · §II examples, §IV guidelines, §V: QNG, CVaR, NN-assisted initialisation, rescaled parameter shift
- **Paper:** B
- **Candidate:** F-H; F-G
- **Matrix rows:** 26
- **Claim category:** Optimisation/post-processing strategies do not overcome outcome-probability concentration
- **Classification:** EXPLICIT
- **Section:** II. Framework; IV. Practical step-by-step guidelines; V. Discussion
- **Subsection:** Guideline steps 1–3 and following text
- **Equation:** (6), (7)
- **Figure:** Fig. 4
- **Appendix:** App. C
- **Page:** 4–5, 8, 10; Fig. 4 on p. 6 (arXiv v2)
- **Source version:** arXiv:2507.22054v2 (2026-06-04); also in v1
- **Source URL:** https://arxiv.org/abs/2507.22054v2
- **Short quote:** —
- **Source statement (paraphrase):** By the guidelines (polynomial POVM + concentrated outcome probabilities ⇒ no information), these methods still suffer from exponential concentration, though they may help training in other ways.
- **Mathematical expression:** —
- **Assumptions:** Polynomial POVMs; concentration.
- **Scope:** Named methods.
- **Relation to Stage 7:** Background.
- **Does NOT establish:** —

<a id="B-14"></a>
### B-14 · §I Introduction: the source's own description of the random-walk corollary
- **Paper:** B
- **Candidate:** F-G; C4
- **Matrix rows:** 8
- **Claim category:** Source's own positioning of its random-walk result (self-description)
- **Classification:** EXPLICIT
- **Section:** I. Introduction
- **Subsection:** Final paragraph
- **Equation:** —
- **Figure:** —
- **Appendix:** —
- **Page:** 2 (arXiv v2)
- **Source version:** arXiv:2507.22054v2 (2026-06-04); also in v1
- **Source URL:** https://arxiv.org/abs/2507.22054v2
- **Short quote:** "it had not been proven as far we are aware" (source's own wording about its Corollary 2)
- **Source statement (paraphrase):** The authors say random-walk behaviour of vanilla GD on barren plateaus under practical shot budgets was mentioned in passing earlier (their Ref. [38], Arrasmith et al. 2021) and prove it via Corollary 2 by taking the post-processing to be a gradient calculation.
- **Mathematical expression:** —
- **Assumptions:** —
- **Scope:** Self-description, quoted as the source's own claim; B5 makes no determination.
- **Relation to Stage 7:** Locates the random-walk fact in prior work (this source and its Ref. [38]).
- **Does NOT establish:** Anything beyond the source's own positioning.

<a id="B-15"></a>
### B-15 · §V Discussion and App. D (Eqs. D1–D13): scope limits (exponentially many outcomes)
- **Paper:** B
- **Candidate:** CTX
- **Matrix rows:** —
- **Claim category:** Exponentially-many-element POVMs are outside the framework (counterexample)
- **Classification:** EXPLICIT
- **Section:** V. Discussion; App. D
- **Subsection:** App. D 1 Counter example; App. D 2 Indistinguishable example
- **Equation:** (D1)–(D13)
- **Figure:** —
- **Appendix:** App. D
- **Page:** 10–11; 25–27 (arXiv v2)
- **Source version:** arXiv:2507.22054v2 (2026-06-04); also in v1 (appendix title differs)
- **Source URL:** https://arxiv.org/abs/2507.22054v2
- **Short quote:** —
- **Source statement (paraphrase):** With exponentially many outcomes, each probability can be exponentially close to its concentration point while the distributions remain distinguishable. A counterexample with an odd/even test is given.
- **Mathematical expression:** —
- **Assumptions:** —
- **Scope:** Scope boundary of the framework.
- **Relation to Stage 7:** Stage 7's two readouts are two-outcome POVMs, inside the framework.
- **Does NOT establish:** —

<a id="B-16"></a>
### B-16 · App. A: Lemmas 1–2, Proposition 1, Definitions 2–4 (Eqs. A1–A23)
- **Paper:** B
- **Candidate:** F-F
- **Matrix rows:** 6
- **Claim category:** Hypothesis-testing tools for statistical indistinguishability
- **Classification:** EXPLICIT
- **Section:** App. A
- **Subsection:** A 1 One sample; A 2 Many samples; A 3 Statistical indistinguishability
- **Equation:** (A1), (A12), (A19), (A23)
- **Figure:** —
- **Appendix:** App. A
- **Page:** 15–18 (arXiv v2)
- **Source version:** arXiv:2507.22054v2 (2026-06-04); v1 numbering of App. A differs (no Product-distribution definition; A-equations end at (A20))
- **Source URL:** https://arxiv.org/abs/2507.22054v2
- **Short quote:** —
- **Source statement (paraphrase):** One-sample optimal success = 1/2 + ‖P − P′‖₁/4. For product distributions ‖P − P′‖₁ ≤ Σ‖P_i − P′_i‖₁. With N samples, success ≤ 1/2 + N‖P₀ − P′₀‖₁/4. ε-statistical indistinguishability is defined for distributions and outputs.
- **Mathematical expression:** Pr[right decision] ≤ 1/2 + N‖P₀ − P′₀‖₁/4
- **Assumptions:** Equal priors; independent samples.
- **Scope:** General.
- **Relation to Stage 7:** As A-11. These 1-norm bounds do not separate the two readouts' exponents (B5_MATH_COMPARISON §6.4).
- **Does NOT establish:** Hellinger/variance-based exponents.

### Paper B — not-located search log

All entries below: arXiv v2 read in full (27 pp., main text + Appendices A–D, including figure captions), the
LaTeX source of v2 and v1 read, and text searches of both versions. Common search terms: sign, direction,
conditional, tie, binomial, Skellam, cosine, angle, inner product, norm, exponent, 4^n, 16^n, SNR,
signal-to-noise, sample complexity, gradient variance.

<a id="B-NL13"></a>
### B-NL13 · Conditional Loschmidt sign law
- **Paper:** B
- **Candidate:** C2
- **Matrix rows:** 13
- **Claim category:** P(correct sign | nonzero estimate) for the Loschmidt readout
- **Classification:** NOT LOCATED
- **Section:** Whole paper (main text + Appendices A–D)
- **Subsection:** —
- **Equation:** —
- **Figure:** —
- **Appendix:** —
- **Page:** 1–27 (arXiv v2); 1–24 (arXiv v1)
- **Source version:** arXiv:2507.22054v2 (2026-06-04) and v1 (2025-07-29)
- **Source URL:** https://arxiv.org/abs/2507.22054v2
- **Short quote:** —
- **Source statement (paraphrase):** Not located in the reviewed versions after searches for: conditional, given nonzero, sign, direction, single count, (1+|sin|)/2, success probability. The Loschmidt statements are the fixed distribution (0, 1) (B-09a), the per-shot variance F(1 − F) (B-10a) and the overlap-test POVM for kernels (§II p. 5). None concerns signs.
- **Mathematical expression:** —
- **Assumptions:** —
- **Scope:** As above.
- **Relation to Stage 7:** —
- **Does NOT establish:** Absence from the wider literature.

<a id="B-NL14b"></a>
### B-NL14b · Same-objective Loschmidt vs SWAP comparison at the parameter-shift-gradient level
- **Paper:** B
- **Candidate:** C3
- **Matrix rows:** 14b
- **Claim category:** Same-θ, same-gradient comparison of two readouts
- **Classification:** NOT LOCATED
- **Section:** Whole paper (main text + Appendices A–D)
- **Subsection:** —
- **Equation:** —
- **Figure:** —
- **Appendix:** —
- **Page:** 1–27 (arXiv v2)
- **Source version:** arXiv:2507.22054v2 (2026-06-04) and v1 (2025-07-29)
- **Source URL:** https://arxiv.org/abs/2507.22054v2
- **Short quote:** —
- **Source statement (paraphrase):** Not located in the reviewed versions after searches for: gradient + Loschmidt, gradient + SWAP, parameter shift + fidelity, same landscape. The Loschmidt/SWAP comparison is at the fidelity level under a 2-design (B-09a, B-10a). The parameter-shift analysis (B-06a) is for Pauli observables.
- **Mathematical expression:** —
- **Assumptions:** —
- **Scope:** As above.
- **Relation to Stage 7:** —
- **Does NOT establish:** Absence from the wider literature.

<a id="B-NL17"></a>
### B-NL17 · Explicit 4ⁿ vs 16ⁿ comparison
- **Paper:** B
- **Candidate:** C3
- **Matrix rows:** 17
- **Claim category:** Explicit pair of shot exponents 4ⁿ (LE) and 16ⁿ (SWAP)
- **Classification:** NOT LOCATED
- **Section:** Whole paper (main text + Appendices A–D)
- **Subsection:** —
- **Equation:** —
- **Figure:** —
- **Appendix:** —
- **Page:** 1–27 (arXiv v2)
- **Source version:** arXiv:2507.22054v2 (2026-06-04) and v1 (2025-07-29)
- **Source URL:** https://arxiv.org/abs/2507.22054v2
- **Short quote:** —
- **Source statement (paraphrase):** Not located in the reviewed versions after searches for: 4^n, 16^n, 4n, 16n, 2^{2n}, exponent, decades per qubit. The only base-specific shot budgets are the numerical regimes 10 × n, 2ⁿ (Fig. 4) and 150 vs 2¹⁵ (Fig. 3). The general statement is N ∈ Ω(exp(n)).
- **Mathematical expression:** —
- **Assumptions:** —
- **Scope:** As above.
- **Relation to Stage 7:** —
- **Does NOT establish:** Absence from the wider literature.

---

## Version differences — Aghaei Saem et al., arXiv v1 (2025-07-29) → v2 (2026-06-04)

Method: paragraph-level diff of the two LaTeX sources (140 vs 159 paragraphs, 21 changed blocks) plus a comparison
of equation/figure numbering in both PDFs. Only differences relevant to B5 are listed.

| Id | Change | Where (v2) | Relevance |
|---|---|---|---|
| V-01 | **Inserted:** "Subtlety regarding the choice of POVM": Loschmidt vs SWAP fixed distributions, Eq. (12) ε_N, both estimator variances, purity example, Fig. 5. | §IV, pp. 9–10 | All B-09/B-10/B-11 evidence exists **only in v2**. |
| V-02 | **Inserted:** "Subtlety regarding measure-first-estimate-later approaches", Eq. (11). | §IV, pp. 8–9 | B-12 is v2-only. |
| V-03 | **Inserted:** global Pauli-Z parity POVM example {Π₊, Π₋}. | §II, p. 5 | B-02 is v2-only. |
| V-04 | Fig. 3 caption now defines the plotted quantity (1/N_p)‖θ^(t) − θ^(0)‖₁ over initialisations. | Fig. 3, p. 6 | Clarifies B-07a/B-07b; the random-walk comparison itself is in both versions. |
| V-05 | Corollary 2 (informal) and Corollary 4 (formal) wording: v1 "results in a random walk … the updated parameters … follow"; v2 "is statistically indistinguishable from a random walk … with high probability at least 1 − c … are statistically indistinguishable from the update rule". The sentence after (B26) also changed (v1: the update "does not incorporate information about the current parameter values"; v2: indistinguishable from a parameter-independent random variable "with probability at least exponentially close to 1"). | §III p. 7; App. B pp. 21, 23 | Same mathematical content (B-06a); v2 makes the indistinguishability framing and probability qualifier explicit. |
| V-06 | Definition 1 adds "(which are drawn from some distribution D over a certain domain)" and a following remark on D-dependence. | §III, pp. 5–6 | B-03; no change to the bound. |
| V-07 | App. A: added Definition 2 (Product distribution); "Statistical indistinguishability" definitions renamed "ε-statistical…" and renumbered (v1 Defs. 2–3 → v2 Defs. 3–4; v1 App. A equations end at (A20), v2 at (A23)). | App. A | Cite App. A numbers with the version. App. B–D equation numbers are identical in v1 and v2. |
| V-08 | QNG paragraph expanded (block-diagonal QGT); Appendix D title changed ("exponential support" → "exponentially many elements"). | §II p. 4; App. D | Not candidate-relevant. |

Main-text equations (1)–(10) and Figs. 1–4 have the same numbers in both versions. Eqs. (11)–(12) and Fig. 5 exist
only in v2.

## Open items (not resolved in B5)

1. **Published-version check for paper B.** The IOP version of record (online 2026-01-30) could not be retrieved
   (bot protection on every IOP URL; no repository copy listed by OpenAlex). Because arXiv v2 post-dates
   publication, it most likely reflects the accepted text, but this is **unverified**. In particular, whether the
   version of record contains the v2-only Loschmidt/SWAP passage (V-01) and with which section/equation numbers
   must be checked before any of B-02, B-09, B-10, B-11 or B-12 is cited in a paper.
2. Paper A's arXiv versions were not compared (not required by the brief). All A locators are to the published
   article and its published SI.
