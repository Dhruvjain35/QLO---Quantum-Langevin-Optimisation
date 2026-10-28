"""Hadamard-test circuits for the overlap a = <phi|psi> with |phi> = V|0^n>, |psi> = U|0^n>, W = V^dagger U.

Real quadrature:      anc H, controlled-W, anc H, measure Z_anc       E[Z] = Re <0^n|W|0^n> = x
Imaginary quadrature: anc H, controlled-W, anc S^dagger, anc H         E[Z] = Im <0^n|W|0^n> = y
Simulated two ways: an explicit numpy statevector of the (n+1)-qubit circuit and PennyLane qml.ctrl circuits.
"""

from __future__ import annotations

import numpy as np
import pennylane as qml

H = np.array([[1, 1], [1, -1]], dtype=np.complex128) / np.sqrt(2.0)
SDG = np.diag([1.0, -1j]).astype(np.complex128)


def hadamard_expectation(W: np.ndarray, imag: bool) -> float:
    """E[Z_anc] of the explicit (n+1)-qubit Hadamard test for a dense unitary W (ancilla = qubit 0)."""
    d = W.shape[0]
    zero = np.zeros(d, dtype=np.complex128)
    zero[0] = 1.0
    state = np.concatenate([zero, zero]) / np.sqrt(2.0)                  # after H on anc: (|0>|0^n> + |1>|0^n>)/sqrt2
    state[d:] = W @ state[d:]                                            # controlled-W
    a0, a1 = state[:d], state[d:]
    if imag:
        a1 = SDG[1, 1] * a1
    b0, b1 = (a0 + a1) / np.sqrt(2.0), (a0 - a1) / np.sqrt(2.0)          # final H on anc
    return float(np.sum(np.abs(b0) ** 2) - np.sum(np.abs(b1) ** 2))


def hea_ops(params, n: int, depth: int) -> None:
    for layer in range(depth):
        for q in range(n):
            qml.RY(params[layer, q, 0], wires=q)
            qml.RZ(params[layer, q, 1], wires=q)
        if n > 1:
            for q in range(n):
                if n == 2 and q == 1:
                    qml.CNOT(wires=[1, 0])
                elif q < n - 1:
                    qml.CNOT(wires=[q, q + 1])
                else:
                    qml.CNOT(wires=[n - 1, 0])


def hea_unitary(params, n: int, depth: int) -> np.ndarray:
    return np.asarray(qml.matrix(lambda p: hea_ops(p, n, depth), wire_order=range(n))(params))


def pennylane_hadamard(params, n: int, depth: int, imag: bool, shots=None):
    """PennyLane Hadamard test with qml.ctrl on the HEA (ancilla = wire n). Returns <Z_anc> (or a sample mean)."""
    dev = qml.device("default.qubit", wires=n + 1)
    anc = n

    def circuit():
        qml.Hadamard(wires=anc)
        qml.ctrl(hea_ops, control=anc)(params, n, depth)
        if imag:
            qml.adjoint(qml.S)(wires=anc)
        qml.Hadamard(wires=anc)
        return qml.expval(qml.PauliZ(anc))

    qn = qml.QNode(circuit, dev)
    if shots is None:
        return float(qn())
    return float(qml.set_shots(qn, shots=shots)())
