"""Analytic (shots=None) cost and gradient evaluation on ``default.qubit``.

Three gradient paths are exposed, all returning an array with the *same shape as
the parameter tensor*:

* :func:`gradient_autograd`        — PennyLane's autograd interface with
  ``diff_method="backprop"`` (reverse-mode differentiation through the
  state-vector simulation). This is the "analytic reference".
* :func:`gradient_parameter_shift` — ``diff_method="parameter-shift"``. Each
  partial derivative is the difference of two shifted circuit evaluations. This
  is the estimator that carries over unchanged to finite-shot devices later.
* :func:`gradient_finite_difference` — plain central finite differences built on
  the scalar cost. Independent of PennyLane's gradient machinery; used as a third
  sanity check in tests.

All QNodes are built with ``shots=None``; no finite-shot logic lives here.
"""

from __future__ import annotations

import numpy as np
import pennylane as qml
from pennylane import numpy as pnp

from qlo.circuits.hardware_efficient import HardwareEfficientAnsatz
from qlo.costs.observables import get_observable

_DIFF_METHODS = {"autograd": "backprop", "parameter-shift": "parameter-shift"}


def make_cost(ansatz: HardwareEfficientAnsatz, cost: str, diff_method: str = "autograd"):
    """Build a QNode ``C(theta) -> float`` on analytic ``default.qubit``.

    ``diff_method`` selects how PennyLane will differentiate the QNode if asked:
    ``"autograd"`` (backprop) or ``"parameter-shift"``.
    """
    if diff_method not in _DIFF_METHODS:
        raise ValueError(f"diff_method must be one of {tuple(_DIFF_METHODS)}")
    observable = get_observable(cost, ansatz.n_qubits)
    dev = qml.device("default.qubit", wires=ansatz.n_qubits, shots=None)

    @qml.qnode(dev, interface="autograd", diff_method=_DIFF_METHODS[diff_method])
    def circuit(params):
        ansatz.apply(params)
        return qml.expval(observable)

    return circuit


def _as_trainable(params) -> pnp.ndarray:
    return pnp.array(np.asarray(params, dtype=np.float64), requires_grad=True)


def cost_value(ansatz: HardwareEfficientAnsatz, cost: str, params) -> float:
    """Scalar analytic cost ``C(theta)``."""
    qnode = make_cost(ansatz, cost, diff_method="autograd")
    return float(qnode(_as_trainable(params)))


def gradient_autograd(ansatz: HardwareEfficientAnsatz, cost: str, params) -> np.ndarray:
    """dC/dtheta via backprop through the state-vector simulator."""
    qnode = make_cost(ansatz, cost, diff_method="autograd")
    grad = qml.grad(qnode)(_as_trainable(params))
    return np.asarray(grad, dtype=np.float64).reshape(ansatz.param_shape)


def gradient_parameter_shift(ansatz: HardwareEfficientAnsatz, cost: str, params) -> np.ndarray:
    """dC/dtheta via the two-term parameter-shift rule (analytic device)."""
    qnode = make_cost(ansatz, cost, diff_method="parameter-shift")
    grad = qml.grad(qnode)(_as_trainable(params))
    return np.asarray(grad, dtype=np.float64).reshape(ansatz.param_shape)


def gradient_finite_difference(
    ansatz: HardwareEfficientAnsatz, cost: str, params, h: float = 1e-5
) -> np.ndarray:
    """Central finite-difference gradient built directly on :func:`cost_value`.

    Deliberately does not use PennyLane's gradient transforms so it is an
    independent check on the two PennyLane paths.
    """
    qnode = make_cost(ansatz, cost, diff_method="autograd")
    base = np.asarray(params, dtype=np.float64).copy()
    grad = np.zeros_like(base)
    for idx in np.ndindex(base.shape):
        plus = base.copy()
        minus = base.copy()
        plus[idx] += h
        minus[idx] -= h
        grad[idx] = (float(qnode(plus)) - float(qnode(minus))) / (2.0 * h)
    return grad
