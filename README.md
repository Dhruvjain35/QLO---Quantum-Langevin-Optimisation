# qlo-research — Stage 1: experimental foundation

## Research question

Can finite-shot measurement noise — normally treated as a limitation of
variational-quantum-algorithm (VQA) optimization — act as *useful* stochastic
exploration, e.g. Langevin-style diffusion, in regimes where gradients are
otherwise vanishingly small ("barren plateaus")?

**Stage 1 (this repository state) does not attempt to answer that.** It only
builds and verifies the minimum trustworthy foundation on which that question
can later be studied: a reproducible ansatz, controlled local/global costs,
verified gradients, deterministic seeding, and a pipeline smoke test.

## What is established vs. hypothesized vs. not claimed

**Established (background, not tested here):**
- Finite-shot quantum measurements produce statistical estimator noise in
  expectation values and therefore in gradient estimates.
- Some VQA architecture/cost-function combinations exhibit barren plateaus
  (gradient variance that decays exponentially with qubit count).
- Stochastic / Langevin dynamics produce exploratory behaviour in optimization.

**Hypotheses to test (in later stages, not here):**
- Shot noise may improve exploration in a *verified* barren-plateau regime.
- The covariance structure of shot noise may generate a useful effective
  diffusion for the parameters.

**Not claimed (and must not be assumed by anything in this repo):**
- That shot noise is fluctuation–dissipation-theorem (FDT) noise.
- That any "QLO" method escapes plateaus in polynomial time.
- That QLO is superior to existing optimizers.
- Any quantum advantage.
- That the ansatz here universally exhibits barren plateaus, or that there is a
  universal qubit-count threshold for them.

## Layout

```
src/qlo/
  circuits/hardware_efficient.py   HardwareEfficientAnsatz, param_shape, n_params
  costs/observables.py             global Z⊗…⊗Z and local (1/n)ΣZ_i observables
  gradients/exact.py               analytic cost + autograd / parameter-shift / finite-diff gradients
  gradients/finite_shot.py         finite-shot parameter-shift estimator (Stage 2)
  analysis/bp_variance.py          Var_theta[dC/dtheta_k] across initializations (the real BP statistic)
  analysis/shot_theory.py          analytic shot variance for the ±1 global observable (validation only)
  analysis/statistics.py           replicate summaries, log-log slope fit
  benchmarks/cerezo2021.py         controlled RX-product BP benchmark (closed forms + PennyLane) (Stage 3)
  experiments/bp_smoke.py          tiny pipeline smoke experiment -> results/*.csv
  experiments/shot_noise.py        Stage 2 shot-noise characterization -> results/stage2/
  experiments/controlled_bp.py     Stage 3 controlled barren-plateau scaling -> results/stage3/
  stage4/                          Stage 4: binomial projector shots, update rules A-E, harness, seeds, stats, figures
  experiments/stage4.py            Stage 4 driver (crosscheck | snr | tune | evaluate | analyze | figures | all) -> results/stage4/
  utils/seeding.py                 make_rng / random_params (explicit RNG, no global state)
tests/                             133 pytest tests
configs/bp_smoke_tiny.yaml         the smoke configuration, for reference
results/                           CSV outputs (smoke only)
STATUS.md                          Stage 1 status report
STAGE2.md                          Stage 2 report (finite-shot noise characterization)
STAGE3.md                          Stage 3 report (controlled barren-plateau benchmark)
STAGE4.md                          Stage 4 report (finite-shot stochastic escape test — negative for shot noise)
```

## Circuit

Hardware-efficient ansatz, `depth` layers of

    RY(θ[l,q,0]) RZ(θ[l,q,1]) on every qubit q, then CNOT entangler

- parameter tensor shape: `(depth, n_qubits, 2)`
- parameter count: `p(n, depth) = 2 · n · depth`
- entangler: `"ring"` (CNOT(q,q+1) for q<n−1, then CNOT(n−1,0)) or `"chain"`
  (open nearest-neighbour line)

## Cost functions

| name     | observable                    |
|----------|-------------------------------|
| `global` | `Z_0 Z_1 … Z_{n-1}`           |
| `local`  | `(1/n) Σ_i Z_i`               |

Both are exposed through `qlo.get_observable(name, n_qubits)` and consumed by
`qlo.make_cost(ansatz, name, diff_method)`.

## Gradients

Analytic `default.qubit` (`shots=None`) only. Three paths, all returning the
parameter-tensor shape:

- `gradient_autograd` — backprop through the state-vector simulator (reference)
- `gradient_parameter_shift` — two-term shift rule; the estimator that will be
  reused unchanged on finite-shot devices later
- `gradient_finite_difference` — central differences on the scalar cost;
  independent of PennyLane's gradient machinery

Tests require autograd and parameter-shift to agree to `1e-10` and autograd and
finite-difference to `1e-6`.

## Setup

```
uv venv --python 3.12 .venv
uv pip install --python .venv/bin/python -e ".[dev]"
```

(Or any Python ≥ 3.12 venv + `pip install -e ".[dev]"`.)

## Run

```
.venv/bin/python -m pytest -v
.venv/bin/python -m qlo.experiments.bp_smoke                 # global cost, tiny config
.venv/bin/python -m qlo.experiments.bp_smoke --cost local

.venv/bin/python -m qlo.experiments.shot_noise                # Stage 2, ~25 s
.venv/bin/python -m qlo.experiments.shot_noise --k 2 3 8 --out-dir results/stage2_nonzero_k
.venv/bin/python -m qlo.experiments.controlled_bp             # Stage 3, ~60 s
.venv/bin/python -m qlo.experiments.controlled_bp --n-init 1000000 --no-pennylane --out-dir results/stage3_n1e6
.venv/bin/python -m qlo.experiments.stage4 all                   # Stage 4, ~15 min on 9 cores
```

The smoke experiment is a pipeline test only. Its numbers are not evidence
about barren plateaus. Stage 2 (see `STAGE2.md`) characterizes finite-shot
*estimator noise* only. Stage 3 (see `STAGE3.md`) reproduces a controlled
literature barren plateau on a separate benchmark; the hardware-efficient ansatz
has structural zero gradients with the global cost and is not that benchmark.
Stage 4 (see `STAGE4.md`) tests whether finite-shot stochasticity helps escape that
verified plateau: it does not (null at matched learning rate; indistinguishable from
covariance-matched Gaussian noise; large classical diffusion helps only for n ≤ 8).
