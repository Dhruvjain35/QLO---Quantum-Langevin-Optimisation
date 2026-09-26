"""Global and local cost observables.

Both are returned as PennyLane operators so a single QNode can compute
``qml.expval(observable)`` and the cost is ``C(theta) = <psi(theta)| O |psi(theta)>``.

* global : ``O_global = Z_0 Z_1 ... Z_{n-1}``           (range ``[-1, 1]``)
* local  : ``O_local  = (1/n) * sum_i Z_i``            (range ``[-1, 1]``)

These are the standard "global vs local" pair used to compare cost-function
locality. Including both does not presuppose that either exhibits a plateau here.
"""

from __future__ import annotations

from functools import reduce

import pennylane as qml

COST_NAMES = ("global", "local")


def global_z_observable(n_qubits: int):
    """Tensor product ``Z_0 ⊗ Z_1 ⊗ ... ⊗ Z_{n-1}``."""
    if n_qubits < 1:
        raise ValueError("n_qubits must be >= 1")
    return reduce(lambda a, b: a @ b, (qml.PauliZ(i) for i in range(n_qubits)))


def local_z_observable(n_qubits: int):
    """Mean single-qubit Z: ``(1/n) * sum_i Z_i`` as one Hamiltonian operator."""
    if n_qubits < 1:
        raise ValueError("n_qubits must be >= 1")
    coeffs = [1.0 / n_qubits] * n_qubits
    ops = [qml.PauliZ(i) for i in range(n_qubits)]
    return qml.Hamiltonian(coeffs, ops)


def get_observable(name: str, n_qubits: int):
    """Uniform accessor: ``name`` in ``COST_NAMES``."""
    if name == "global":
        return global_z_observable(n_qubits)
    if name == "local":
        return local_z_observable(n_qubits)
    raise ValueError(f"unknown cost {name!r}; expected one of {COST_NAMES}")
