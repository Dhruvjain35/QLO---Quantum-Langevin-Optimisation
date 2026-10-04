# B5 potential issues

Issues noticed while reading the two prior papers. **Nothing here was fixed.** No Stage 1–7 file was modified.

## Outcome of the numerical cross-check

**No numerical or code issue found.** The Stage 7 outcome models and per-shot variances match the sources exactly:

| Stage 7 (code / STAGE7.md) | Source |
|---|---|
| Loschmidt: one shot ~ Bernoulli(F); `estimators.fidelity_estimate("loschmidt", K, M) = K/M`; `fidelity_estimate_variance` = F(1 − F)/M | A main Eq. (12) and p. 4; A SI Eq. (17); B §IV Var^(LE) = F(1 − F) per shot (B-10a†) |
| SWAP: ancilla +1 w.p. q = (1 + F)/2; `fidelity_estimate("swap", K, M) = 2K/M − 1`; variance (1 − F²)/M | A main p. 4 (p₊ = 1/2 + κ/2); A SI Eq. (46); B §IV Var^(SWAP) = 1 − F² per shot (B-10a†) |
| Fidelity-level Loschmidt zero probability (1 − F)^M (`analysis.prior_work_replication`: `p0 = exp(M log1p(−F))`) | A SI Eqs. (26)–(27), conditional factor (1 − s)^N |

The structural statements Stage 7 lists as prior work (§3, §18, §19) are present in the sources:
- Loschmidt estimates collapse to zero: A Prop. 1, SI Supp. Prop. 2; B §IV†.
- SWAP estimates become data-independent: A Prop. 2, SI Supp. Lemma 3; B §IV†.
- Global-Pauli parameter-shift GD is a random walk: B Corollaries 2/4.

## Items for the principal audit (documentation / wording, not code)

### PI-1 · Version provenance of the Loschmidt/SWAP statements attributed to Aghaei Saem et al.
- **Suspected issue:**
  - STAGE7.md §3, §19 and §22 attribute to *Aghaei Saem et al., QST 11, 015049 (2026)*:
    - that Loschmidt and SWAP give different fixed distributions for the same fidelity;
    - that different POVMs for the same quantity have different estimator variances.
  - In the arXiv record these statements appear in **v2 (2026-06-04)**, §IV "Subtlety regarding the choice of
    POVM", and are **absent from v1 (2025-07-29)**.
  - The published (IOP) version of record could not be retrieved in B5, so whether it contains this passage, and
    under which section/equation numbers, is unverified.
- **Source evidence:** B5_EVIDENCE_LEDGER.md entries B-09a, B-10a, B-10b, B-10c; version table V-01.
- **Why it might matter:** A paper citing the QST article for these points needs the version-of-record locator.
  Citing arXiv v1 for them would be wrong.
- **Stage 7 location potentially affected:** `STAGE7.md` §3, §19, §22 (prose only). No code.

### PI-2 · Arithmetic-mean vs log-typical scales when positioning against Thanasilp et al.
- **Suspected issue:**
  - Paper A states its shot requirements in terms of the mean kernel value μ = 2^(−n) (e.g. "N ∈ Ω(2ⁿ)" for a fixed
    non-zero fraction of Loschmidt estimates, SI Fig. 4a).
  - Stage 5/7 express medians through the log-typical A = 4^(−(n−1)) (≈ 4ⁿ for Loschmidt).
  - Both are consistent (Ω(2ⁿ) is a lower bound), but a side-by-side quotation could look like a conflict unless
    the two scales are named explicitly.
- **Source evidence:** A-10a, A-10b; B5_MATH_COMPARISON §1, §6.3.
- **Why it might matter:** Positioning text (A4) comparing exponents.
- **Stage 7 location potentially affected:** `STAGE5.md` §7, `STAGE7.md` §11–§12 wording when cited against A. No
  code.

### PI-3 · Probabilistic qualifier on the random-walk statement
- **Suspected issue:**
  - Paper B's Corollaries 2/4 hold "with high probability" over the random initialisation (probability ≥ 1 − c,
    c ∈ O(exp(−n))).
  - STAGE7.md §17 and the README Stage 7 paragraph state that SWAP-driven GD at n ≥ 10 "is indistinguishable from a
    signal-free random walk" without that qualifier.
  - Yet Stage 7's own table shows non-zero success for SWAP at n = 10: P(F_final ≥ 0.5) between 0.04 and 0.14
    across variants, against 0 for the random walk.
- **Source evidence:** B-06a (assumptions and scope); B-14 (source's own positioning).
- **Why it might matter:** The prior-work statement is probabilistic over initialisations. Matching its scope avoids
  over-generalising Stage 7's diagnostic. (This point was examined numerically on the separate B3 branch; B5 does not
  use or depend on that work.)
- **Stage 7 location potentially affected:** `STAGE7.md` §17, §18, §20 F, §21; `README.md` Stage 7 paragraph (prose
  only). Function `qlo.stage7.analysis.optimization_diagnostic` is **not** suspected of any error.

This file records observations for the principal audit and does not make the novelty decision.
