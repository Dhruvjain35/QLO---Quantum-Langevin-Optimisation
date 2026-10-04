# B5 prior-work matrix

Evidence packet for the principal novelty audit (A1). The matrix records what two prior papers **explicitly
contain**. It does not decide novelty.

- **Thanasilp 2024** = Thanasilp, Wang, Cerezo, Holmes, Nat. Commun. 15, 5200 (2024). Published article + SI.
- **Aghaei Saem 2026** = Aghaei Saem, Tafreshi, Holmes, Thanasilp, QST 11, 015049 (2026). Reviewed as
  **arXiv:2507.22054v2**: the published IOP text was not retrievable. A dagger (†) marks evidence present only in
  arXiv v2 (absent from v1); its presence in the version of record is unverified.

Labels:
- `EXPLICIT` — stated in the source.
- `PARTIAL / RELATED` — a related but non-identical statement.
- `PARTIAL / IMPLIED` — follows from the source's formulas only after algebra the source does not perform (candidate 3 rows).
- `NOT LOCATED` — not found in the reviewed version after the logged searches. This is **not** a statement about the
  wider literature.

Every cell links to its entries in [`B5_EVIDENCE_LEDGER.md`](B5_EVIDENCE_LEDGER.md). Rows 1–25 are the rows required
by the B5 brief (row 14 is split into 14a/14b). Rows 26–28 are added known-territory/context rows.

| # | Candidate / claim | Thanasilp 2024 | Aghaei Saem 2026 | Notes |
|---|---|---|---|---|
| 1 | Fidelity exponential concentration | EXPLICIT — [A-01](B5_EVIDENCE_LEDGER.md#A-01), [A-02](B5_EVIDENCE_LEDGER.md#A-02), [A-03](B5_EVIDENCE_LEDGER.md#A-03), [A-04](B5_EVIDENCE_LEDGER.md#A-04) | EXPLICIT — [B-03](B5_EVIDENCE_LEDGER.md#B-03), [B-09a](B5_EVIDENCE_LEDGER.md#B-09a)† | A: kernel (fidelity) concentration incl. the product single-qubit-rotation family (Prop. 3). B: general outcome-probability concentration; fidelity-specific only under a 2-design (†). |
| 2 | Loschmidt estimator collapse toward zero | EXPLICIT — [A-05](B5_EVIDENCE_LEDGER.md#A-05), [A-06a](B5_EVIDENCE_LEDGER.md#A-06a) | EXPLICIT — [B-09a](B5_EVIDENCE_LEDGER.md#B-09a)† | Fidelity-estimate level in both. A gives the exact single-estimate factor (1 − s)^N. |
| 3 | SWAP estimator random / data-independent limit | EXPLICIT — [A-08](B5_EVIDENCE_LEDGER.md#A-08) | EXPLICIT — [B-09a](B5_EVIDENCE_LEDGER.md#B-09a)† | Fidelity-estimate level; ±1 equiprobable null in both. |
| 4 | Different Loschmidt / SWAP estimator variances | PARTIAL / RELATED — [A-08b](B5_EVIDENCE_LEDGER.md#A-08b) | EXPLICIT — [B-10a](B5_EVIDENCE_LEDGER.md#B-10a)† | B states the per-shot fidelity-estimator variances F(1−F) and 1−F². Gradient-estimator variances are not located in either. A gives only the outcome models. |
| 5 | Exponential finite-shot measurement burden | EXPLICIT — [A-10a](B5_EVIDENCE_LEDGER.md#A-10a), [A-12a](B5_EVIDENCE_LEDGER.md#A-12a), [A-13](B5_EVIDENCE_LEDGER.md#A-13) | EXPLICIT — [B-04](B5_EVIDENCE_LEDGER.md#B-04), [B-08a](B5_EVIDENCE_LEDGER.md#B-08a), [B-10b](B5_EVIDENCE_LEDGER.md#B-10b)† | Both state exponential (Ω(b^(2n)), Ω(2ⁿ) numerically, Ω(exp n)) shot needs. Neither ties the exponent to the readout. |
| 6 | Outcome-distribution distinguishability framework | EXPLICIT — [A-11](B5_EVIDENCE_LEDGER.md#A-11) | EXPLICIT — [B-03](B5_EVIDENCE_LEDGER.md#B-03), [B-04](B5_EVIDENCE_LEDGER.md#B-04), [B-16](B5_EVIDENCE_LEDGER.md#B-16) | A: binary outcome distributions (kernel tests). B: general POVMs. Both use 1-norm hypothesis-testing bounds. |
| 7 | Parameter-shift finite-shot analysis | NOT LOCATED — [A-NL07](B5_EVIDENCE_LEDGER.md#A-NL07) | EXPLICIT — [B-01](B5_EVIDENCE_LEDGER.md#B-01) | B: generic finite-shot parameter-shift update (Eq. 6), analysed via indistinguishability (Cor. 2/4). No readout-specific estimator statistics. |
| 8 | Parameter-shift random-walk behaviour | NOT LOCATED — [A-NL08](B5_EVIDENCE_LEDGER.md#A-NL08) | EXPLICIT — [B-06a](B5_EVIDENCE_LEDGER.md#B-06a), [B-07a](B5_EVIDENCE_LEDGER.md#B-07a), [B-14](B5_EVIDENCE_LEDGER.md#B-14) | B: Corollary 2/4 (Pauli observables, ±1 outcomes) + 15-qubit global-Z numerics. B-14 is the source's own positioning (quoted as such). |
| 9 | Exact gradient P_zero | PARTIAL / RELATED — [A-06b](B5_EVIDENCE_LEDGER.md#A-06b) | PARTIAL / RELATED — [B-09b](B5_EVIDENCE_LEDGER.md#B-09b)† | Both only at the fidelity-estimate level. The difference-of-counts P(ĝ = 0) (incl. equal nonzero counts) is not located. |
| 10 | Exact gradient P_correct | NOT LOCATED — [A-NL10](B5_EVIDENCE_LEDGER.md#A-NL10) | PARTIAL / RELATED — [B-06c](B5_EVIDENCE_LEDGER.md#B-06c) | B: indistinguishability from a sign-symmetric null bears on the sign, but no probability is stated. |
| 11 | Exact gradient P_wrong | NOT LOCATED — [A-NL11](B5_EVIDENCE_LEDGER.md#A-NL11) | PARTIAL / RELATED — [B-06c](B5_EVIDENCE_LEDGER.md#B-06c) | As row 10. |
| 12 | Exact difference-of-binomials gradient law | NOT LOCATED — [A-NL12](B5_EVIDENCE_LEDGER.md#A-NL12) | PARTIAL / RELATED — [B-06b](B5_EVIDENCE_LEDGER.md#B-06b) | B writes the update as a difference of two means of ±1 outcomes (B18, B26). Its distribution and tie probability are not derived. Loschmidt (0/1-outcome) form not written. |
| 13 | Conditional Loschmidt sign law P(correct \| ĝ ≠ 0) → (1+\|sin θ_k\|)/2 | NOT LOCATED — [A-NL13](B5_EVIDENCE_LEDGER.md#A-NL13) | NOT LOCATED — [B-NL13](B5_EVIDENCE_LEDGER.md#B-NL13) | Not located after concept searches (conditional, sign, direction, single count). |
| 14a | Same-objective Loschmidt vs SWAP comparison — fidelity / kernel-estimate level | EXPLICIT — [A-09a](B5_EVIDENCE_LEDGER.md#A-09a), [A-10a](B5_EVIDENCE_LEDGER.md#A-10a) | EXPLICIT — [B-09a](B5_EVIDENCE_LEDGER.md#B-09a)† | A: same kernels/data, analytic + numerical (incl. product-R_y kernels, n = 5–40). B: same fidelity objective, 2-design, analytic only. |
| 14b | Same-objective Loschmidt vs SWAP comparison — parameter-shift-gradient level, same θ | NOT LOCATED — [A-NL14b](B5_EVIDENCE_LEDGER.md#A-NL14b) | NOT LOCATED — [B-NL14b](B5_EVIDENCE_LEDGER.md#B-NL14b) | — |
| 15 | Gradient shot-complexity exponent — Loschmidt | PARTIAL / RELATED — [A-10b](B5_EVIDENCE_LEDGER.md#A-10b) | PARTIAL / IMPLIED — [B-10c](B5_EVIDENCE_LEDGER.md#B-10c)† | A: kernel-level numerics, "N ∈ Ω(2ⁿ)" for a fixed non-zero fraction. B: variance + Eq. (12) imply readout dependence only after algebra (see B5_MATH_COMPARISON §6). Neither is gradient level. |
| 16 | Gradient shot-complexity exponent — SWAP | PARTIAL / RELATED — [A-10b](B5_EVIDENCE_LEDGER.md#A-10b) | PARTIAL / IMPLIED — [B-10c](B5_EVIDENCE_LEDGER.md#B-10c)† | A: "at least exponentially" (kernel level, numerical). B: as row 15. |
| 17 | Explicit 4ⁿ vs 16ⁿ comparison | NOT LOCATED — [A-NL17](B5_EVIDENCE_LEDGER.md#A-NL17) | NOT LOCATED — [B-NL17](B5_EVIDENCE_LEDGER.md#B-NL17) | No 4ⁿ, 16ⁿ, or fitted-exponent statement in either reviewed source. |
| 18 | Measurement scheme changing the gradient-resolution exponent | NOT LOCATED — [A-NL18](B5_EVIDENCE_LEDGER.md#A-NL18) | PARTIAL / IMPLIED — [B-10c](B5_EVIDENCE_LEDGER.md#B-10c)† | A's resolution bound (SI Supp. Prop. 5) is range-based and the same for both readouts. Its Gram-matrix corollary assumes per-outcome fluctuations "stay constant", which holds for ±1 but not for 0/1 Loschmidt outcomes at small F; the source does not discuss this (A-12a). B says estimators differ in shot-noise impact but computes no exponent. |
| 19 | Full-gradient vector exactly-zero probability | PARTIAL / RELATED — [A-07](B5_EVIDENCE_LEDGER.md#A-07) | PARTIAL / RELATED — [B-09b](B5_EVIDENCE_LEDGER.md#B-09b)† | A: joint all-zero probability of an estimated Gram matrix (product form). B: fidelity-level zero only. |
| 20 | Estimated/exact gradient cosine similarity | NOT LOCATED — [A-NL20](B5_EVIDENCE_LEDGER.md#A-NL20) | PARTIAL / RELATED — [B-06d](B5_EVIDENCE_LEDGER.md#B-06d) | B: random-walk indistinguishability is related to, but not the same statement as, cosine ≈ 0. No cosine is computed. |
| 21 | Gradient norm inflation | NOT LOCATED — [A-NL21](B5_EVIDENCE_LEDGER.md#A-NL21) | PARTIAL / RELATED — [B-06e](B5_EVIDENCE_LEDGER.md#B-06e) | B: null update has parameter-independent magnitude (B15). No norm comparison is made. |
| 22 | Probability ĝ·g > 0 | NOT LOCATED — [A-NL22](B5_EVIDENCE_LEDGER.md#A-NL22) | PARTIAL / RELATED — [B-06d](B5_EVIDENCE_LEDGER.md#B-06d) | As row 20. |
| 23 | Component sign accuracy | NOT LOCATED — [A-NL23](B5_EVIDENCE_LEDGER.md#A-NL23) | PARTIAL / RELATED — [B-06d](B5_EVIDENCE_LEDGER.md#B-06d) | As row 20. |
| 24 | Matched-start finite-shot optimisation trajectories | NOT LOCATED — [A-NL24](B5_EVIDENCE_LEDGER.md#A-NL24) | PARTIAL / RELATED — [B-07c](B5_EVIDENCE_LEDGER.md#B-07c) | B compares shot budgets (Figs. 3–4); identical starts are not stated. |
| 25 | Signal-free random-walk control | PARTIAL / RELATED — [A-09b](B5_EVIDENCE_LEDGER.md#A-09b) | EXPLICIT — [B-07b](B5_EVIDENCE_LEDGER.md#B-07b) | B: "Random Walk" reference for displacement mean/variance (construction unstated). A: random-matrix control for a kernel model (not trajectories). |
| 26 | Limits of classical post-processing / mitigation (fact H) | EXPLICIT — [A-14](B5_EVIDENCE_LEDGER.md#A-14), [A-15](B5_EVIDENCE_LEDGER.md#A-15), [A-16](B5_EVIDENCE_LEDGER.md#A-16) | EXPLICIT — [B-05](B5_EVIDENCE_LEDGER.md#B-05), [B-12](B5_EVIDENCE_LEDGER.md#B-12)†, [B-13](B5_EVIDENCE_LEDGER.md#B-13) | Post-processing maps, multi-copy processing, error mitigation, classical shadows. |
| 27 | Concentration of the product single-qubit-rotation fidelity landscape (Stage 3–7 family) | EXPLICIT — [A-03](B5_EVIDENCE_LEDGER.md#A-03) | PARTIAL / RELATED — [B-08c](B5_EVIDENCE_LEDGER.md#B-08c) | A: Var ≤ E[κ²] = (3/8)ⁿ, mean 2^(−n), for y-rotations; the same per-factor law holds for Stage 3–7's x-rotations (our observation). Log-typical value not stated. |
| 28 | Global-Z (parity) parameter-shift training numerics on a single rotation layer (Stage 6 relation) | NOT LOCATED — [A-NL28](B5_EVIDENCE_LEDGER.md#A-NL28) | EXPLICIT — [B-07a](B5_EVIDENCE_LEDGER.md#B-07a), [B-08b](B5_EVIDENCE_LEDGER.md#B-08b) | Same setup as Stage 6's parity benchmark (already noted in STAGE6.md §19). |

## Reading guide by candidate

- **Candidate 1** (exact finite-shot gradient zero/sign probabilities): rows 9–12.
- **Candidate 2** (conditional Loschmidt sign law): row 13.
- **Candidate 3** (same-landscape readout-dependent gradient shot scaling): rows 4, 14a–18. Algebra is in
  [`B5_MATH_COMPARISON.md`](B5_MATH_COMPARISON.md) §6.
- **Candidate 4** (full-vector / trajectory consequences): rows 8, 19–25.
- **Known prior territory** (facts A–H of the brief): rows 1–6, 8, 26.

`NOT LOCATED` means only: not found in the reviewed version after the searches logged in the ledger.
This matrix is evidence for the principal novelty audit and does not make the novelty decision.
