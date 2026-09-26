"""Explicit PennyLane circuits for both measurement schemes, compared with the direct binomial model.

LOSCHMIDT: n wires, RX(theta_j) on each, computational-basis samples; a shot "succeeds" iff all bits are 0
           (the projector |0^n><0^n|), so K ~ Bin(M, F).
SWAP:      wire 0 = ancilla, wires 1..n = |psi(theta)>, wires n+1..2n = target |0^n> (left in |0>).
           H(anc); CSWAP(anc, a_j, b_j) for every j; H(anc); sample Z(anc). P(Z = +1) = (1 + F)/2.
"""

from __future__ import annotations

import numpy as np
import pennylane as qml


def make_loschmidt(n: int, M: int, seed):
    dev = qml.device("default.qubit", wires=n, seed=seed)

    @qml.set_shots(shots=M)
    @qml.qnode(dev, diff_method=None)
    def circuit(theta):
        for j in range(n):
            qml.RX(theta[j], wires=j)
        return qml.sample(wires=range(n))

    def count(theta) -> int:
        bits = np.asarray(circuit(theta)).reshape(M, n)
        return int(np.sum(np.all(bits == 0, axis=1)))

    return count


def make_swap(n: int, M: int, seed):
    dev = qml.device("default.qubit", wires=2 * n + 1, seed=seed)

    @qml.set_shots(shots=M)
    @qml.qnode(dev, diff_method=None)
    def circuit(theta):
        qml.Hadamard(wires=0)
        for j in range(n):
            qml.RX(theta[j], wires=1 + j)
        for j in range(n):
            qml.CSWAP(wires=[0, 1 + j, 1 + n + j])
        qml.Hadamard(wires=0)
        return qml.sample(qml.PauliZ(0))

    def count(theta) -> int:
        z = np.asarray(circuit(theta)).ravel()
        return int(np.sum(z > 0))  # number of +1 ancilla outcomes

    return count


def swap_expval(n: int, theta) -> float:
    """Analytic (shots=None) ancilla <Z> of the SWAP-test circuit; must equal F = prod cos^2(theta_j/2)."""
    dev = qml.device("default.qubit", wires=2 * n + 1)

    @qml.qnode(dev)
    def circuit(t):
        qml.Hadamard(wires=0)
        for j in range(n):
            qml.RX(t[j], wires=1 + j)
        for j in range(n):
            qml.CSWAP(wires=[0, 1 + j, 1 + n + j])
        qml.Hadamard(wires=0)
        return qml.expval(qml.PauliZ(0))

    return float(circuit(np.asarray(theta, dtype=float)))
