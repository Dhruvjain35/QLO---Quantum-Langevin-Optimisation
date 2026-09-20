"""Exact finite-shot estimator for the projector cost ``C = 1 - |<0^n|psi>|^2``.

Measuring ``|0^n><0^n|`` ``M`` times on the ideal simulator yields ``K ~ Binomial(M, F(theta))``
successes, and ``C_hat = 1 - K/M``. Drawing ``K`` directly from the Binomial is statistically
identical to sampling every shot from the state vector (verified against PennyLane in
``qlo.experiments.stage4`` ``crosscheck``), and is what makes the optimization study cheap.

Parameter-shift gradient (independent ``M``-shot batches for the + and - shift of every k):

    g_hat_k = 1/2 [ C_hat(theta + pi/2 e_k) - C_hat(theta - pi/2 e_k) ]

Conditional (given theta) noise variance, from Binomial variance ``F(1-F)/M`` per evaluation:

    Var(g_hat_k | theta) = [ F_+(1-F_+) + F_-(1-F_-) ] / (4M)

With one independent batch per (k, shift) the conditional covariance is diagonal BY
CONSTRUCTION of this measurement strategy; that is not a universal property.

Randomness: one ``numpy.random.Generator`` per estimator instance, built from an explicit
``SeedSequence`` child (see ``qlo.stage4.seeds``). Every call draws fresh, independent
Binomial variates from that stream (a vector draw with a vector of probabilities gives
independent variates), so the + shift, the - shift, different k, different iterations and
different replicates never share a draw, and everything is reproducible from the master seed.
The estimator never touches the exact gradient.
"""

from __future__ import annotations

import numpy as np

from qlo.stage4.landscape import fidelity, shifted_fidelities


def conditional_shot_variance(theta: np.ndarray, shots: int) -> np.ndarray:
    """``[F_+(1-F_+) + F_-(1-F_-)] / (4M)`` for every k. Analytic; validation/matched-noise use."""
    f_plus, f_minus = shifted_fidelities(theta)
    return (f_plus * (1.0 - f_plus) + f_minus * (1.0 - f_minus)) / (4.0 * int(shots))


class BinomialProjectorShots:
    """Finite-shot projector cost + parameter-shift gradient with resource accounting."""

    def __init__(self, n_qubits: int, shots: int, seed) -> None:
        if int(shots) < 1:
            raise ValueError("shots must be >= 1")
        self.n_qubits = int(n_qubits)
        self.shots = int(shots)
        self.rng = np.random.default_rng(seed)
        self.circuit_evaluations = 0
        self.total_shots = 0

    def cost_estimate(self, theta: np.ndarray) -> float:
        """One ``M``-shot estimate ``1 - K/M``, ``K ~ Binomial(M, F(theta))``."""
        p = float(np.clip(fidelity(theta), 0.0, 1.0))
        k = self.rng.binomial(self.shots, p)
        self.circuit_evaluations += 1
        self.total_shots += self.shots
        return 1.0 - k / self.shots

    def gradient(self, theta: np.ndarray) -> np.ndarray:
        """Full parameter-shift estimate; ``2n`` independent evaluations, ``2nM`` shots."""
        theta = np.asarray(theta, dtype=np.float64)
        if theta.ndim != 1 or theta.size != self.n_qubits:
            raise ValueError(f"theta must have shape ({self.n_qubits},)")
        f_plus, f_minus = shifted_fidelities(theta)
        k_plus = self.rng.binomial(self.shots, f_plus)    # n independent draws
        k_minus = self.rng.binomial(self.shots, f_minus)  # n more independent draws
        self.circuit_evaluations += 2 * self.n_qubits
        self.total_shots += 2 * self.n_qubits * self.shots
        c_plus = 1.0 - k_plus / self.shots
        c_minus = 1.0 - k_minus / self.shots
        return 0.5 * (c_plus - c_minus)
