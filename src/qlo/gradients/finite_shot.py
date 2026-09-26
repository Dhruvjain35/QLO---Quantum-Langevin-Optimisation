"""Finite-shot parameter-shift gradient estimator.

For flat parameter index ``k`` and ``M`` shots per circuit evaluation:

    g_hat_k(theta) = 0.5 * [ C_hat(theta + (pi/2) e_k)  -  C_hat(theta - (pi/2) e_k) ]

where each ``C_hat`` is the sample mean of ``M`` independent measurement outcomes
of the cost observable on ``default.qubit``.

Randomness model
----------------
* One ``numpy.random.SeedSequence`` is built from ``master_seed``.
* **Every** circuit evaluation (each ``+pi/2`` call, each ``-pi/2`` call, every
  replicate, every parameter index) spawns a fresh child ``SeedSequence`` and
  gets its own ``numpy.random.Generator`` handed to a fresh ``default.qubit``
  device. Children of one ``SeedSequence`` are statistically independent streams.
* Spawning is sequential and deterministic, so the whole replicate sequence is
  reproducible from ``master_seed`` alone, and the spawn keys are recorded in
  ``self.spawn_keys`` so tests can audit that no stream was reused.
* The device is created with an explicit ``seed=`` — never PennyLane's
  ``seed="global"`` default — so the legacy global NumPy state is not involved.

This module deliberately does **not** import the analytic reference module.
The estimator never sees the analytic gradient; that is only used offline for
validation (see ``qlo.analysis``).
"""

from __future__ import annotations

import numpy as np
import pennylane as qml

from qlo.circuits.hardware_efficient import HardwareEfficientAnsatz
from qlo.costs.observables import get_observable
from qlo.utils.indexing import shifted_params

SHIFT = np.pi / 2.0


class FiniteShotParameterShift:
    """Stateful (only in its RNG stream) finite-shot parameter-shift estimator."""

    def __init__(
        self,
        ansatz: HardwareEfficientAnsatz,
        cost: str,
        shots: int,
        master_seed,
    ) -> None:
        if int(shots) < 1 or int(shots) != shots:
            raise ValueError(f"shots must be a positive integer, got {shots!r}")
        self.ansatz = ansatz
        self.cost = cost
        self.shots = int(shots)
        self.master_seed = master_seed
        self._observable = get_observable(cost, ansatz.n_qubits)
        self._seed_seq = np.random.SeedSequence(master_seed)
        self.spawn_keys: list[tuple[int, ...]] = []  # audit trail, one per circuit evaluation
        self.n_evaluations = 0

    # ----- randomness -------------------------------------------------------------
    def _fresh_generator(self) -> np.random.Generator:
        child = self._seed_seq.spawn(1)[0]
        self.spawn_keys.append(tuple(child.spawn_key))
        return np.random.default_rng(child)

    # ----- single M-shot cost estimate --------------------------------------------
    def cost_estimate(self, params) -> float:
        """One ``M``-shot estimate ``C_hat(params)`` with fresh simulator randomness."""
        params = np.asarray(params, dtype=np.float64)
        if params.shape != self.ansatz.param_shape:
            raise ValueError(f"params shape {params.shape} != {self.ansatz.param_shape}")
        dev = qml.device("default.qubit", wires=self.ansatz.n_qubits, seed=self._fresh_generator())
        ansatz, observable = self.ansatz, self._observable

        @qml.set_shots(shots=self.shots)
        @qml.qnode(dev, diff_method=None)
        def circuit(p):
            ansatz.apply(p)
            return qml.expval(observable)

        self.n_evaluations += 1
        return float(circuit(params))

    # ----- gradient estimates ------------------------------------------------------
    def gradient_component(self, params, k: int) -> float:
        """``g_hat_k`` from two independent ``M``-shot evaluations (+ then -)."""
        c_plus = self.cost_estimate(shifted_params(params, k, +SHIFT))
        c_minus = self.cost_estimate(shifted_params(params, k, -SHIFT))
        return 0.5 * (c_plus - c_minus)

    def gradient(self, params) -> np.ndarray:
        """Full estimate, shaped like ``params``; ``2p`` independent evaluations."""
        params = np.asarray(params, dtype=np.float64)
        flat = np.array([self.gradient_component(params, k) for k in range(params.size)])
        return flat.reshape(params.shape)

    def replicate_component(self, params, k: int, n_replicates: int) -> np.ndarray:
        """``n_replicates`` independent draws of ``g_hat_k``; shape ``(n_replicates,)``."""
        return np.array([self.gradient_component(params, k) for _ in range(int(n_replicates))])

    def replicate_gradient(self, params, n_replicates: int) -> np.ndarray:
        """``n_replicates`` independent full-gradient draws; shape ``(n_replicates, p)``."""
        return np.stack([self.gradient(params).ravel() for _ in range(int(n_replicates))])
