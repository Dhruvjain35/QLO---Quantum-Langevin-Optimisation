# QLO principal track: working document

Baseline: Stage 7 as recorded in the September 2026 update (commit ba51866, 180 tests). Prepared 4 October 2026.

This file covers the framing decision and items A1 to A4 from the work allocation. A5 (the update to QuIP Network) is kept out of this public repo. The A2 numerics are in `a2_numerics/`, run from inside that folder with `python general_scaling.py 4000` and `python condsign_check.py`.

**Two things to know before reading.**

1. **Stage 7 is not on GitHub.** The public repo ends at Stage 6 (merge commit 87c4e06). Commit ba51866 exists only on your machine. All Stage 7 numbers below come from your PDF, not from code I could run. Push it before anyone else (the intern, QuIP Network, a referee) is pointed at the repo.
2. **Prior work quotes need a line check.** I read Thanasilp et al. and Aghaei Saem et al. through the arXiv HTML pages using a summarising fetch tool. The section and equation pointers below are a map, not a citation. Open the PDFs and confirm each one before it goes into the paper. This is also exactly what B5 (the intern's extraction table) is for.

---

## 0. One page summary

| Item | Determination | What it means for you |
|---|---|---|
| Framing (Section 2) | **Recommend: new paper, new title. Retire the Langevin apparatus to its own track.** Keep one paragraph that uses Langevin language to explain the Stage 4 failure. | Your call. Section 1 sets out both options. |
| A1, C1 exact zero and sign probabilities | **Refinement.** Mathematically elementary. Not a headline. | Keep as a tool, not a claim. |
| A1, C2 conditional sign law | **New as far as I can find.** Generalises beyond the product circuit (A2). | One of the two results the paper can stand on. |
| A1, C3 same landscape 4^n vs 16^n | **Refinement.** The per shot variances that imply it are already printed in Aghaei Saem et al. They do not state the shot ratio or apply it to gradients. | Claim the exact identity and the "exponent doubles" rule, not the discovery that SWAP costs more. |
| A1, C4 full vector alignment and trajectories | **Mostly already published in qualitative form** (Aghaei Saem et al., Corollary 2 and Fig. 3, random walk trajectories on the same RX circuit). | Supporting evidence only. |
| A1, new C5 from A2 | **New candidate.** An exact identity for the SWAP to Loschmidt shot ratio with no product assumption, and the slope rule that follows. | The other result the paper can stand on. |
| A2 generalisation | **Positive so far.** The doubling of the exponent survives on entangling circuits; the bases 4 and 16 do not. | See Section 3 for regime of validity. |
| A3 normalisation | **Use shots (circuit executions) per gradient component as the primary unit.** Report qubit shots as a secondary row. The exponent claim does not depend on the choice. | One paragraph for the methods section is drafted. |
| A4 positioning | Three paragraphs drafted, ready to paste. | Section 5. |

**Overall A1 verdict.** The paper exists, but in a narrower form than "SWAP needs 16^n shots". What survives is: an exact, model free identity that ties the measurement scheme to the shot exponent (the exponent doubles), plus the conditional sign law showing that the two schemes fail in ways that carry different information. Both are checkable in a few lines, which is good for a referee.

---

## 1. The framing decision

### The two options

**Option A. Reconnect the Langevin apparatus.** Keep the QLO title and present the measurement result as the reason the original mechanism fails, then rebuild the effective temperature and escape time analysis on top of it.

**Option B. Separate paper.** Retitle around measurement and gradient resolvability. Move the effective temperature derivation, the Eyring Kramers and narrow escape analysis, the Riemannian metric extension and the temperature matching conjecture to a separate track.

### Recommendation: Option B, with one bridge paragraph

Reasons:

1. **The Langevin story now has close prior work of its own.** Kaminishi, Mori, Sugawara and Yamamoto, *Impact of Measurement Noise on Escaping Saddles in Variational Quantum Algorithms*, arXiv:2406.09780 (2024) already model shot noise with an SDE whose temperature is set by η/N_s, and show that more noise shortens saddle escape time. So the original QLO hypothesis was partly anticipated, while your Stage 4 result (no help on a real barren plateau) is the useful contrast. A Langevin framed paper would have to argue against them on their own ground. A measurement framed paper only needs to cite them.
2. **Option A asks a referee to accept an apparatus that the paper's own data says is unused.** There is no optimizer in Stage 7.
3. **The measurement result is cleaner on its own.** Same landscape, same exact gradient, one variable changed, closed forms matching to three decimals.

**The bridge paragraph (keeps the Stage 4 result and the original motivation without the apparatus).** In Langevin language the per step noise variance plays the role of temperature. For the parameter shift estimator:

* Loschmidt: $\mathrm{Var}(\hat g_k) = [F_+(1-F_+) + F_-(1-F_-)]/(4M) \approx (F_+ + F_-)/(4M)$. The temperature falls with the fidelity scale, at the same exponential rate as the drift. The system freezes.
* SWAP: $\mathrm{Var}(\hat g_k) = [2 - F_+^2 - F_-^2]/(4M) \approx 1/(2M)$. The temperature stays fixed while the drift $g_k \propto F$ vanishes. The system diffuses with no bias.

So the reason noise cannot help in a barren plateau is that the noise is either tied to the vanishing signal (it switches off with it) or untied from it (it carries no direction). Neither gives the "noise near a weak but present drift" regime that Kaminishi et al. exploit at saddles. This is a few lines, it is honest about the falsification, and it does not need the Eyring Kramers machinery.

### Title ideas (plain, no claims of novelty)

* *Same landscape, different shots: how the fidelity measurement sets the exponent of finite shot gradient cost*
* *Silence or noise: two finite shot failure modes of parameter shift gradients under fidelity concentration*
* *The measurement doubles the exponent: Loschmidt echo versus SWAP test gradients on a barren plateau*

---

## 2. A1 novelty audit

### What the two adjacent papers say (map for the line check)

**Thanasilp, Wang, Cerezo, Holmes, *Exponential concentration in quantum kernel methods*, Nat. Commun. 15, 5200 (2024), arXiv:2208.11060.**

* Section II.2 sets up the Loschmidt echo test (assign +1 to the all zero bitstring) and the SWAP test (±1 outcomes with $p_+ = 1/2 + \kappa/2$).
* Proposition 1 (around Eq. 13): with polynomially many shots and an exponentially concentrated kernel, the Loschmidt estimated Gram matrix equals the identity with probability $1 - O(c^{-n})$. The mechanism stated in the text is $(1-\mu)^N \approx 1 - N\mu$.
* Proposition 2 (Eq. 14): the SWAP estimated Gram matrix is statistically indistinguishable from one built from fair ±1 coin flips.
* Scope: kernel **values**. No parameter shift gradients, no exact zero probability theorem, no sign probabilities, no comparison of shot counts between the two tests.

**Aghaei Saem, Tafreshi, Holmes, Thanasilp, *Pitfalls when tackling the exponential concentration of parameterized quantum models*, Quantum Sci. Technol. 11, 015049 (2026), arXiv:2507.22054.**

* Section IV ("Subtlety regarding the choice of POVM", Fig. 5): Loschmidt POVM $\{|0\rangle\langle0|^{\otimes n}, 1 - |0\rangle\langle0|^{\otimes n}\}$ concentrates to (0, 1) and the estimated fidelity is zero; SWAP POVM concentrates to (1/2, 1/2).
* Same section: **per shot variances are stated, $\mathrm{Var}^{\mathrm{SWAP}} = 1 - F^2$ and $\mathrm{Var}^{\mathrm{LE}} = F(1-F)$**, with the remark that the different forms give different statistical behaviour. As far as the fetch shows, they do **not** turn this into a shot ratio, and do not apply it to gradients.
* Eq. 12: the resolution ratio $\epsilon_N = \mathrm{Var}_\rho[\hat\ell] / (N\,\mathrm{Var}_\alpha[\ell])$, needing $\epsilon_N \lesssim 1$, hence exponential N under concentration.
* Eq. 6 and Corollary 2: parameter shift updates under concentrated outcome probabilities are indistinguishable from parameter independent random variables, so the trajectory is a random walk.
* Fig. 3: **15 qubits, one layer of single qubit X rotations, global Z cost**. This is the Stage 6 parity benchmark circuit. Trajectories with 150 shots look like a random walk (PCA plot, displacement mean and variance).

### Determination per candidate contribution

**C1. Exact finite shot parameter shift zero and sign probabilities. Refinement.**
The estimator is a difference of two independent binomials, so $P(\hat g = 0) = \sum_r \mathrm{Bin}(r;M,F_+)\,\mathrm{Bin}(r;M,F_-)$ and its Poisson limit $e^{-M(F_++F_-)} I_0(2M\sqrt{F_+F_-})$ is the zero mass of a Skellam distribution. Neither paper writes these down for gradients, but a referee will see them as textbook probability applied to a known setting. Use them as the exact machinery that makes everything else checkable. Do not list them as a contribution on their own.

**C2. Conditional sign law. New as far as I can find.**
Neither paper has it, and searches turned up nothing that does. It also says something the prior work does not: the Loschmidt estimator, when it fires, is *informative*. Aghaei Saem et al. Corollary 2 says the update is indistinguishable from a parameter independent variable with high probability. That stays true, because the informative events are the rare ones. The conditional law shows where the information sits: in the rare nonzero events, with sign accuracy around 0.85, not spread thinly over every shot. It also generalises beyond the product circuit (Section 3, Result 2), which makes it a statement about the measurement and not about the toy circuit.
*Risk:* it is a short derivation (Poisson thinning). Present it as a clean exact law with a clear message, not as hard mathematics.

**C3. Same landscape 4^n against 16^n. Refinement, and the most exposed.**
From $\mathrm{Var}^{\mathrm{LE}} = F(1-F)$ and $\mathrm{Var}^{\mathrm{SWAP}} = 1-F^2$ printed in Aghaei Saem et al., resolving a quantity of size F needs about $1/F$ shots with the Loschmidt echo and about $1/F^2$ with the SWAP test. A referee can derive your 4^n against 16^n from their Section IV in two lines. What they did not do: state the ratio, show that it carries over exactly to parameter shift **gradients**, isolate it on a single landscape, or confirm it to three decimals. So the honest claim is *"we make exact and verify for gradients what their variance expressions imply"*, citing their Section IV directly. Do not claim the scaling gap as a discovery.

**C4. Full gradient vector alignment and trajectory consequences. Mostly already published (qualitative).**
Aghaei Saem et al. show random walk trajectories on the same RX circuit with global Z (Fig. 3) and give the random walk corollary. Your median cosine numbers (0.74 when nonzero for Loschmidt, 0.009 for SWAP) and the "74 percent of vectors exactly zero" figure are sharper, and the frozen versus wandering contrast is worth one figure. But it is supporting evidence, not a contribution.

**C5 (new, from A2). Model free shot ratio identity and the exponent doubling rule.** See Section 3. This is the strongest candidate because it is (a) exact, (b) not tied to the product circuit, (c) not stated in either paper, and (d) gives a prediction that B1 can test.

### Other papers found during the search

| Paper | Relevance |
|---|---|
| Kaminishi et al., arXiv:2406.09780 (2024) | Shot noise as Langevin temperature η/N_s helps escape saddles. Central to the framing decision; cite in Discussion. |
| Arrasmith et al., Quantum Sci. Technol. 7, 045015 (2022), arXiv:2104.05868 | Barren plateaus, cost concentration and narrow gorges are equivalent. Needed for A4. |
| *The Cost of Certainty: Shot Budgets in Quantum Program Testing*, arXiv:2510.22418 | Reports that the inverse (Loschmidt type) test is the most sample efficient and the SWAP test costs about a factor of two more, in the near identity regime of program testing. This is the opposite end of the fidelity range from yours (F near 1, not F near 0). Worth one citation to show the gap is regime dependent: with both shifted fidelities equal to F, the identity in Section 3 gives a shot ratio of $(1+F)/F$, which is 2 at F = 1 and grows like $1/F$ as F goes to 0. **Read before citing.** |
| Kang, arXiv:2605.01319 (2026) | Sign organisation of gradient *terms* in barren plateaus. Not about finite shot estimators. No overlap, but the word "sign" might make a referee ask; one line distinguishes it. |
| Li et al., arXiv:2607.11095 (2026) | Measurement cost of gradient based attacks on quantum classifiers. No overlap. |

**Searches that came up empty** (worth repeating with a library database before submission): exact zero probability of a finite shot parameter shift gradient; conditional sign probability of a nonzero finite shot gradient; SWAP versus Loschmidt shot ratio for gradients.

---

## 3. A2 generalisation beyond the RX product circuit

**Short answer.** The measurement dependent scaling survives entanglement, in a precise form. The bases 4 and 16 are properties of the RX product circuit. What is general is that **the SWAP exponent equals the Loschmidt exponent plus the decay rate of the fidelity scale**, and that this **doubles** the exponent whenever the relative gradient does not itself shrink with n. That held on every circuit tested.

### Setting (no product form)

Any n qubit circuit $|\psi(\theta)\rangle = U(\theta)|0^n\rangle$, any pure target $|\phi\rangle = V|0^n\rangle$, fidelity $F(\theta) = |\langle\phi|\psi(\theta)\rangle|^2$, cost $C = 1 - F$. The parameter $\theta_k$ enters through a gate $e^{-i\theta_k P/2}$ with P a Pauli string, so the two term shift rule is exact. Write
$F_\pm = F(\theta \pm \tfrac{\pi}{2} e_k)$, $S = F_+ + F_-$, $g_k = (F_- - F_+)/2$, and the **relative gradient** $r = (F_- - F_+)/S \in [-1, 1]$. One independent batch of M shots at each shift.

* **Loschmidt echo:** run $V^\dagger U(\theta_\pm)$, count all zero outcomes, $K_\pm \sim \mathrm{Bin}(M, F_\pm)$, $\hat g = (K_- - K_+)/(2M)$.
* **SWAP test:** ancilla reads 0 with probability $(1+F_\pm)/2$, $K_\pm \sim \mathrm{Bin}(M, (1+F_\pm)/2)$, $\hat g = (K_- - K_+)/M$.

Both are unbiased for $g_k$.

### Result 1. Exact shot ratio identity (any circuit)

$$\mathrm{Var}_{\mathrm{LE}}(\hat g_k) = \frac{F_+(1-F_+) + F_-(1-F_-)}{4M}, \qquad \mathrm{Var}_{\mathrm{SW}}(\hat g_k) = \frac{2 - F_+^2 - F_-^2}{4M}.$$

So the shots needed for $\mathrm{SNR} \ge \rho$ are

$$M_{\mathrm{LE}} = \rho^2\,\frac{F_+(1-F_+) + F_-(1-F_-)}{(F_- - F_+)^2}, \qquad M_{\mathrm{SW}} = \rho^2\,\frac{2 - F_+^2 - F_-^2}{(F_- - F_+)^2},$$

and their ratio is

$$\boxed{\;R = \frac{M_{\mathrm{SW}}}{M_{\mathrm{LE}}} = \frac{2 - F_+^2 - F_-^2}{S - F_+^2 - F_-^2} = \frac{2}{S}\,\bigl(1 + O(S)\bigr)\;}$$

*Proof.* Two independent binomials; the per shot variances are those of Aghaei Saem et al. Section IV, applied at the two shifts. The gradient cancels in the ratio. ∎

What it says:

* The ratio depends **only on the sum of the shifted fidelities**. Not on the gradient, not on $\rho$, not on the circuit's structure.
* $R \ge 1$ always (since $S \le 2$), with equality only when both shifted states equal the target. So for this estimator the SWAP test never needs fewer shots.
* Because $\log R$ is a decreasing function of S alone, the **median** of $\log R$ over initialisations is exactly $\log(2/\mathrm{median}\,S)$ up to the $O(S)$ term. No approximation about medians of sums is needed for this part.

**Slope rule.** In the concentrated regime $M_{\mathrm{LE}} \approx \rho^2/(r^2 S)$ and $M_{\mathrm{SW}} \approx 2\rho^2/(r^2 S^2)$. If the typical (median) $\log_{10} S$ falls at rate $b_S$ per qubit and the typical $\log_{10}|r|$ falls at rate $b_r$, then

$$b_{\mathrm{LE}} = b_S + 2b_r, \qquad b_{\mathrm{SW}} = 2b_S + 2b_r, \qquad b_{\mathrm{SW}} - b_{\mathrm{LE}} = b_S.$$

The exponent exactly doubles when $b_r = 0$, that is when the relative gradient does not concentrate.

### Result 2. General conditional sign law (any circuit)

As $M S \to 0$, a nonzero Loschmidt estimate almost surely comes from a single count, and that count sits in the minus batch with probability $F_-/S$. The sign is right when the count is in the batch with the larger fidelity. Hence

$$P\bigl(\operatorname{sign}\hat g_k = \operatorname{sign} g_k \,\big|\, \hat g_k \ne 0\bigr) \;\longrightarrow\; \frac{1 + |r|}{2}.$$

On the RX product circuit $r = \sin\theta_k$, which is the Stage 7 law $(1 + |\sin\theta_k|)/2$. **The constant 0.854 is specific to the product circuit** (it is the median of $(1 + |\sin\theta|)/2$ for uniform $\theta$). On entangling circuits the median is about 0.73 to 0.76. The law is general; its median is not.

### Numerical check (`a2_numerics/`)

Standalone NumPy statevector code, checked against dense matrices to $6\times10^{-18}$ and against finite differences to $8\times10^{-13}$. A Monte Carlo run of 200,000 SWAP estimates matched the variance formula to 0.5 percent. 4,000 random initialisations per n, $\theta \sim U[-\pi, \pi]$, shifting the first rotation on qubit 0. Hardware efficient ansatz = RY then RZ on every qubit, then a CNOT ring (the repo's layout). Fidelity to $|0^n\rangle$.

| Circuit | n range (fit) | $b_S$ | $b_r$ | $b_{\mathrm{LE}}$ | $b_{\mathrm{SW}}$ | ratio | $b_{\mathrm{LE}} + b_S$ (predicted $b_{\mathrm{SW}}$) | $R^2$ (LE, SW) |
|---|---|---|---|---|---|---|---|---|
| RX product (regression arm) | 4 to 20 | 0.595 | 0.000 | 0.595 | 1.189 | 1.998 | 1.190 | 0.9999, 0.9999 |
| HEA, depth 2 | 4 to 12 | 0.376 | 0.000 | 0.387 | 0.759 | 1.963 | 0.762 | 0.9998, 0.9999 |
| HEA, depth n | 6 to 12 | 0.302 | 0.001 | 0.304 | 0.606 | 1.994 | 0.606 | 0.9996, 1.0000 |

Reading the table:

* **The regression arm reproduces Stage 7** (0.595 and 1.189 against your 0.600 to 0.607 and 1.197 to 1.203, from 4,000 rather than 100,000 samples). The harness is consistent with yours.
* **Deep circuits:** $b_S = 0.302 \approx \log_{10} 2$, as expected when the state behaves like a random state and $F \sim 2^{-n}$. So the bases become about **2 for Loschmidt and 4 for SWAP**, not 4 and 16.
* **Shallow circuits:** an intermediate rate, about $10^{0.387} \approx 2.4$ against $10^{0.759} \approx 5.7$ per qubit.
* In all three, the predicted SWAP slope $b_{\mathrm{LE}} + b_S$ matches the fitted one within 0.003. The ratio of slopes is 1.96 to 2.00. The small shortfall in the shallow case is the gap between $b_{\mathrm{LE}}$ and $b_S$ (0.011), which is a median effect at these sizes, not a change in the rule.
* The relative gradient has no trend in any family (median $\log_{10}|r|$ near $-0.15$, $-0.26$, $-0.30$, flat in n). That is why the exponent doubles.
* **Conditional sign law** (`condsign_results.json`, exact binomial sums, 400 points per setting, M = 16 and 1024, six settings): the gap to $(1+|r|)/2$ grows in proportion to $MS$. It stays at most 0.0005 for $MS < 0.01$, at most 0.005 for $MS < 0.1$, and at most 0.045 for $MS < 1$, on every circuit. Medians at M = 16: product n = 20, 0.852 against law 0.852; depth 2, n = 12, 0.761 against 0.760; depth n, n = 12, 0.733 against 0.733.

(Note on files: the `condsign_deep_max_abs_err` column in `general_scaling_results.json` came from an earlier version of the sign routine that lost precision when $MS < 10^{-8}$. It is wrong for the product circuit at n ≥ 10. Use `condsign_results.json`, which sums both tails directly.)

### Regime of validity (state it in the paper)

Holds exactly: pure states, ideal noiseless measurement, the ancilla SWAP test (the destructive SWAP test has the same outcome statistics for pure states), gates with Pauli generators, one independent batch per shift. The ratio identity needs nothing else. The "exponent doubles" statement needs in addition that the typical relative gradient does not decay with n ($b_r = 0$).

Not covered: depolarising or other noise (B4), shared or correlated shot allocation, adaptive shot schedules, generators with more than two eigenvalues (multi term shift rules change the constants but not the structure), parameters other than the first layer (the tests shifted only the first rotation; a parameter late in a deep circuit could in principle have a concentrating $r$), and costs that are not fidelities. No proof is given that $b_S = \log_{10} 2$ for deep circuits; it matches the random state picture but is only observed here for n up to 12.

### Hand off to B1

The intern's sweeps are now a direct test of a stated prediction: for every ansatz family and every parameter position, $b_{\mathrm{SW}} = b_{\mathrm{LE}} + b_S$, and the slope ratio is 2 exactly when $b_r = 0$. Ask for $b_S$ and $b_r$ in the B1 table (Section 6). Extending n beyond 12 and adding confidence intervals (B3) is what turns this table into a figure.

---

## 4. A3 normalisation of the comparison

### The question

The SWAP test uses 2n + 1 qubits and the Loschmidt echo uses n. Is "4^n against 16^n" fair per shot, per qubit, or per unit of time?

### The answer to put in the paper

**Primary unit: shots, meaning circuit executions, per gradient component, with one batch of M shots at each of the two shifts.** This is the quantity the parameter shift rule consumes, the quantity hardware providers bill, and the quantity in Aghaei Saem et al.'s $\epsilon_N$. Every number in Stage 5 to 7 is already in this unit.

**The exponent claim does not depend on the choice.** Any fair conversion (qubits, gate count, circuit depth, wall clock time) multiplies the shot count by a factor polynomial in n. That changes the intercept and adds a slowly vanishing term to the slope, never the base of the exponential. From the Stage 7 medians (SNR ≥ 1):

| Unit | Loschmidt n = 10 | Loschmidt n = 20 | Slope | SWAP n = 10 | SWAP n = 20 | Slope | SWAP / Loschmidt slope |
|---|---|---|---|---|---|---|---|
| log10 shots | 5.70 | 11.72 | 0.602 | 11.11 | 23.14 | 1.203 | 2.00 |
| log10 qubit shots (n and 2n + 1) | 6.70 | 13.02 | 0.632 | 12.43 | 24.75 | 1.232 | 1.95 |

At n = 20 the SWAP test needs about $10^{11.4}$ times more shots. The extra qubits change that by a factor of about 2. No polynomial overhead can close an exponential gap.

### Two points a referee may raise, with answers

1. **"The two tests do not need the same resources."** True, and this cuts both ways. The Loschmidt echo needs the *inverse* of the target preparation circuit, so its depth is roughly doubled. The SWAP test needs only a *copy* of the target state, plus n controlled SWAPs (or a constant depth destructive version) and twice the width. When the target is a state you can only prepare and not invert (data from another device, an unknown state), the SWAP test is the only option, and the paper should say plainly that the comparison applies when both are available.
2. **"Per qubit shot is the fair unit."** Report it as a second row (table above). It moves the slope ratio from 2.00 to 1.95 at n = 10 to 20, and the gap tends back to 2 as n grows.

### Methods paragraph (ready to paste)

> We count cost in shots, meaning executions of the measured circuit, with M shots at each of the two parameter shifts of a given component. The SWAP test circuit acts on 2n + 1 qubits and the Loschmidt echo circuit on n qubits with roughly twice the depth, so any per qubit, per gate or wall clock normalisation multiplies our shot counts by a factor polynomial in n. Such factors shift intercepts but cannot change the base of an exponential, and we report per qubit shot slopes alongside per shot slopes to show this directly. The comparison assumes that both the inverse of the target preparation (needed for the echo) and an independent copy of the target state (needed for the SWAP test) are available.

---

## 5. A4 positioning against cost concentration

### Draft related work paragraphs (ready to paste, adjust citations)

> Barren plateaus, exponential concentration of the cost about its mean, and exponentially narrow gorges are equivalent for parameterised quantum circuits [Arrasmith et al. 2022], and it is well understood that concentration forces an exponential number of measurement shots to resolve the cost or its gradient [Cerezo et al. 2021; Wang et al. 2021; Larocca et al. 2025]. Our work takes this exponential cost as given. Its subject is the constant in the exponent and what sets it. We show that this constant is not a property of the landscape alone: on a single landscape with an identical exact gradient, changing only the measurement that estimates the fidelity doubles the exponent of the shot cost.
>
> Closest to our work, Thanasilp et al. [2024] show for quantum kernels that, under concentration and with polynomially many shots, Loschmidt echo estimates collapse to zero while SWAP test estimates become indistinguishable from fair coin flips, and Aghaei Saem et al. [2026] place this in a general framework of concentrated measurement outcome probabilities, giving the per shot variances $F(1-F)$ and $1-F^2$ of the two fidelity estimators and showing that finite shot parameter shift training then resembles a random walk. We build directly on these results. We move from the value of the fidelity to its parameter shift gradient, derive the exact finite shot distribution of the gradient estimator under both measurements, and obtain an exact identity for the ratio of the shots they need, which shows that the SWAP to Loschmidt ratio is set by the shifted fidelities alone and that the shot exponent doubles whenever the relative gradient does not itself concentrate. We also show that the two failure modes differ in information content: a nonzero Loschmidt gradient keeps a sign accuracy of $(1 + |r|)/2$ in the deep plateau, while the SWAP sign accuracy tends to one half.
>
> Shot noise has also been studied as a possible resource, acting as an effective temperature that speeds escape from saddle points [Kaminishi et al. 2024]. Our results explain why this mechanism does not carry over to barren plateaus: under the Loschmidt echo the noise vanishes together with the gradient, and under the SWAP test it persists without carrying direction.

**What these paragraphs do and do not claim.** They claim the constant in the exponent, the gradient level exact identity, and the conditional sign law. They give Thanasilp et al. and Aghaei Saem et al. credit for both failure modes and for the variance expressions. They do not claim that SWAP being worse is new.

---

## 6. What changes downstream

**For you (principal), in order.**

1. Push Stage 7 to GitHub.
2. Make the framing decision (Section 1). Everything below assumes Option B.
3. Line check the A1 map against the two PDFs, especially Aghaei Saem et al. Section IV. If their text *does* state the shot ratio or apply it to gradients, C3 drops to "already published" and C5 carries more weight. The paper still stands on C2 and C5.
4. Turn Section 3 into a proper proposition with proof (the identity is two lines; the slope rule needs the stated regime).
5. Send A5 once the compute numbers are filled in.

**For the intern, updated acceptance checks** (no framing or novelty content):

* **B1.** Besides the existing regression check, report for every ansatz family the fitted slope of the median of $\log_{10}(F_+ + F_-)$ and of the median of $\log_{10}|r|$. The prediction to test is $b_{\mathrm{SWAP}} = b_{\mathrm{LE}} - b_S$ (Section 3). Do not tell them what you expect the deep circuit bases to be until they have run it.
* **B2.** The destructive SWAP test has the same outcome statistics as the ancilla SWAP test for pure states, so it should reproduce the SWAP column within error. A disagreement means a bug, not a finding. The Hadamard test estimates an amplitude, not the fidelity, so the gradient estimator has a different form; ask for its exact variance to be derived and checked against simulation before any slopes are fitted.
* **B5.** Add Section IV and Fig. 5 of Aghaei Saem et al. and Proposition 1 and 2 of Thanasilp et al. as required rows.
