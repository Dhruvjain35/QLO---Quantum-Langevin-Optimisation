"""Controlled state-preparation benchmark (Cerezo et al., Nat. Commun. 12, 1791 (2021)).

Circuit:      V(theta) = ⊗_j RX(theta_j) acting on |0>^{⊗n};  p = n parameters.
Initializer:  theta_j ~ iid Uniform[-pi, pi].
Target:       |0>^{⊗n}.

Global cost   C_G = 1 - |<0^n|psi>|^2          = 1 - prod_j cos^2(theta_j/2)
              dC_G/dtheta_k = (1/2) sin(theta_k) prod_{j!=k} cos^2(theta_j/2)
              E_theta[dC_G/dtheta_k] = 0,   Var_theta[dC_G/dtheta_k] = (1/8)(3/8)^(n-1)

Local cost    C_L = 1 - (1/n) sum_j P(qubit j = 0) = 1 - (1/n) sum_j cos^2(theta_j/2)
              dC_L/dtheta_k = sin(theta_k) / (2n)
              E_theta[dC_L/dtheta_k] = 0,   Var_theta[dC_L/dtheta_k] = 1/(8 n^2)

Two independent implementations live here on purpose:
  * ``closed_form_*`` / ``theoretical_*`` — pure NumPy, no PennyLane. External
    validation reference and the engine of the large Monte Carlo sweep.
  * ``Cerezo2021Circuit`` — the actual PennyLane circuit (analytic default.qubit),
    used to verify that the closed forms describe what the simulator computes.

Nothing here is quantum-hardware execution; the Monte Carlo sweep is exact
arithmetic on a formula that has been checked against the simulator.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pennylane as qml
from pennylane import numpy as pnp

COST_NAMES = ("global", "local")
SHIFT = np.pi / 2.0

# --------------------------------------------------------------------------------------
# Parameter sampling
# --------------------------------------------------------------------------------------
def sample_params(n_qubits: int, rng: np.random.Generator, n_samples: int | None = None) -> np.ndarray:
    """theta ~ Uniform[-pi, pi]; shape ``(n,)`` or ``(n_samples, n)``."""
    if n_qubits < 1:
        raise ValueError("n_qubits must be >= 1")
    shape = (n_qubits,) if n_samples is None else (int(n_samples), n_qubits)
    return rng.uniform(-np.pi, np.pi, size=shape).astype(np.float64)


# --------------------------------------------------------------------------------------
# Closed forms (NumPy only). All accept theta of shape (..., n) and vectorize over
# leading axes.
# --------------------------------------------------------------------------------------
def _cos2_half(theta: np.ndarray) -> np.ndarray:
    return np.cos(np.asarray(theta, dtype=np.float64) / 2.0) ** 2


def closed_form_global_cost(theta) -> np.ndarray:
    """``1 - prod_j cos^2(theta_j/2)``."""
    return 1.0 - np.prod(_cos2_half(theta), axis=-1)


def closed_form_local_cost(theta) -> np.ndarray:
    """``1 - mean_j cos^2(theta_j/2)``."""
    return 1.0 - np.mean(_cos2_half(theta), axis=-1)


def closed_form_global_gradient(theta) -> np.ndarray:
    """Full gradient ``(..., n)``: ``(1/2) sin(theta_k) prod_{j!=k} cos^2(theta_j/2)``.

    The product over ``j != k`` is formed explicitly (no division by cos^2(theta_k/2)),
    so a parameter at exactly ±pi does not produce 0/0.
    """
    theta = np.asarray(theta, dtype=np.float64)
    n = theta.shape[-1]
    c2 = _cos2_half(theta)
    out = np.empty_like(theta)
    for k in range(n):
        others = np.prod(np.delete(c2, k, axis=-1), axis=-1) if n > 1 else np.ones(theta.shape[:-1])
        out[..., k] = 0.5 * np.sin(theta[..., k]) * others
    return out


def closed_form_local_gradient(theta) -> np.ndarray:
    """Full gradient ``(..., n)``: ``sin(theta_k) / (2n)``."""
    theta = np.asarray(theta, dtype=np.float64)
    return np.sin(theta) / (2.0 * theta.shape[-1])


def closed_form_gradient_component(cost: str, theta, k: int) -> np.ndarray:
    """Single component ``dC/dtheta_k``, vectorized over leading axes (used by the sweep)."""
    theta = np.asarray(theta, dtype=np.float64)
    n = theta.shape[-1]
    if not 0 <= k < n:
        raise IndexError(f"k={k} out of range for n={n}")
    if cost == "global":
        c2 = _cos2_half(theta)
        others = np.prod(np.delete(c2, k, axis=-1), axis=-1) if n > 1 else np.ones(theta.shape[:-1])
        return 0.5 * np.sin(theta[..., k]) * others
    if cost == "local":
        return np.sin(theta[..., k]) / (2.0 * n)
    raise ValueError(f"cost must be one of {COST_NAMES}")


def closed_form_global_gradient_component_logdomain(theta, k: int) -> tuple[np.ndarray, np.ndarray]:
    """Numerical-stability diagnostic: ``(sign, log|dC_G/dtheta_k|)`` computed in the log domain.

    ``log|g| = log(0.5) + log|sin theta_k| + sum_{j!=k} log cos^2(theta_j/2)``.
    Returns ``-inf`` where the derivative is exactly zero. Does not change the estimand.
    """
    theta = np.asarray(theta, dtype=np.float64)
    with np.errstate(divide="ignore"):
        logc2 = np.log(_cos2_half(theta))
        s = np.sin(theta[..., k])
        log_abs = np.log(0.5) + np.log(np.abs(s)) + np.sum(np.delete(logc2, k, axis=-1), axis=-1)
    return np.sign(s), log_abs


def theoretical_global_variance(n_qubits: int) -> float:
    """``(1/8) (3/8)^(n-1)`` — Var over theta ~ U[-pi,pi]^n of dC_G/dtheta_k (any k)."""
    if n_qubits < 1:
        raise ValueError("n_qubits must be >= 1")
    return (1.0 / 8.0) * (3.0 / 8.0) ** (n_qubits - 1)


def theoretical_local_variance(n_qubits: int) -> float:
    """``1 / (8 n^2)`` — Var over theta ~ U[-pi,pi]^n of dC_L/dtheta_k (any k)."""
    if n_qubits < 1:
        raise ValueError("n_qubits must be >= 1")
    return 1.0 / (8.0 * n_qubits**2)


def theoretical_variance(cost: str, n_qubits: int) -> float:
    if cost == "global":
        return theoretical_global_variance(n_qubits)
    if cost == "local":
        return theoretical_local_variance(n_qubits)
    raise ValueError(f"cost must be one of {COST_NAMES}")


THEORETICAL_GLOBAL_SEMILOG_SLOPE = float(np.log(3.0 / 8.0))  # d log V_G / d n
THEORETICAL_LOCAL_LOGLOG_SLOPE = -2.0  # d log V_L / d log n


# --------------------------------------------------------------------------------------
# PennyLane implementation (analytic default.qubit, shots=None)
# --------------------------------------------------------------------------------------
@dataclass(frozen=True)
class Cerezo2021Circuit:
    """``⊗_j RX(theta_j)`` on ``|0>^n`` with the two Cerezo-2021 costs, on PennyLane."""

    n_qubits: int

    def __post_init__(self) -> None:
        if self.n_qubits < 1:
            raise ValueError("n_qubits must be >= 1")

    @property
    def n_params(self) -> int:
        return self.n_qubits

    def apply(self, theta) -> None:
        """Queue the gates. ``theta[..., j]`` supports PennyLane parameter broadcasting."""
        for j in range(self.n_qubits):
            qml.RX(theta[..., j], wires=j)

    def _observable(self, cost: str):
        n = self.n_qubits
        if cost == "global":
            return qml.Projector(np.zeros(n, dtype=int), wires=range(n))
        if cost == "local":
            return qml.Hamiltonian([1.0 / n] * n, [qml.Projector([0], wires=j) for j in range(n)])
        raise ValueError(f"cost must be one of {COST_NAMES}")

    def make_cost(self, cost: str, diff_method: str = "backprop"):
        """Return ``C(theta)`` as a differentiable Python callable (``1 - <observable>``)."""
        if diff_method not in ("backprop", "parameter-shift"):
            raise ValueError("diff_method must be 'backprop' or 'parameter-shift'")
        dev = qml.device("default.qubit", wires=self.n_qubits, shots=None)
        observable = self._observable(cost)

        @qml.qnode(dev, interface="autograd", diff_method=diff_method)
        def fidelity_like(theta):
            self.apply(theta)
            return qml.expval(observable)

        def cost_fn(theta):
            return 1.0 - fidelity_like(theta)

        return cost_fn

    def cost(self, cost: str, theta) -> float:
        return float(self.make_cost(cost)(np.asarray(theta, dtype=np.float64)))

    def gradient(self, cost: str, theta, diff_method: str = "backprop") -> np.ndarray:
        """Full gradient ``(n,)`` via PennyLane (``backprop`` or ``parameter-shift``)."""
        fn = self.make_cost(cost, diff_method)
        t = pnp.array(np.asarray(theta, dtype=np.float64), requires_grad=True)
        return np.asarray(qml.grad(fn)(t), dtype=np.float64).reshape(self.n_qubits)

    def gradient_component_parameter_shift_batched(self, cost: str, Theta: np.ndarray, k: int) -> np.ndarray:
        """``dC/dtheta_k`` for a batch ``Theta`` of shape ``(N, n)`` by the two-term shift rule,
        using two broadcast analytic circuit executions. Pure simulator arithmetic, no closed form."""
        Theta = np.asarray(Theta, dtype=np.float64)
        if Theta.ndim != 2 or Theta.shape[1] != self.n_qubits:
            raise ValueError("Theta must be (N, n)")
        fn = self.make_cost(cost)
        plus, minus = Theta.copy(), Theta.copy()
        plus[:, k] += SHIFT
        minus[:, k] -= SHIFT
        return 0.5 * (np.asarray(fn(plus), dtype=np.float64) - np.asarray(fn(minus), dtype=np.float64))


# --------------------------------------------------------------------------------------
# Exact higher moments: how hard is it to *estimate* the variance by Monte Carlo?
# --------------------------------------------------------------------------------------
# theta ~ U[-pi, pi]:  E[sin^2] = 1/2, E[sin^4] = 3/8, E[cos^4(theta/2)] = 3/8, E[cos^8(theta/2)] = 35/128.
_E_SIN4 = 3.0 / 8.0
_E_COS8_HALF = 35.0 / 128.0


def theoretical_fourth_moment(cost: str, n_qubits: int) -> float:
    """``E_theta[(dC/dtheta_k)^4]`` (any k).

    global: (3/128) (35/128)^(n-1)      local: 3 / (128 n^4)
    """
    if n_qubits < 1:
        raise ValueError("n_qubits must be >= 1")
    if cost == "global":
        return (_E_SIN4 / 16.0) * _E_COS8_HALF ** (n_qubits - 1)
    if cost == "local":
        return _E_SIN4 / (16.0 * n_qubits**4)
    raise ValueError(f"cost must be one of {COST_NAMES}")


def theoretical_variance_estimator_rel_se(cost: str, n_qubits: int, n_samples: int) -> float:
    """Predicted relative standard error of the sample variance (ddof=1) of dC/dtheta_k
    from ``n_samples`` iid theta draws: ``sqrt((E[g^4]/Var^2 - (N-3)/(N-1)) / N)``.

    For the global cost ``E[g^4]/Var^2 = (3/2)(35/18)^(n-1)`` grows exponentially, so
    a fixed sample size resolves the variance only up to some n. For the local cost the
    ratio is the constant 3/2.
    """
    N = int(n_samples)
    if N < 4:
        raise ValueError("n_samples must be >= 4")
    var = theoretical_variance(cost, n_qubits)
    kurt = theoretical_fourth_moment(cost, n_qubits) / var**2
    return float(np.sqrt((kurt - (N - 3) / (N - 1)) / N))


def factorized_global_second_moment_estimate(theta: np.ndarray, k: int) -> float:
    """Diagnostic estimator of ``E[(dC_G/dtheta_k)^2]`` that exploits THIS benchmark's product
    structure: ``E[g_k^2] = (1/4) E[sin^2 theta_k] prod_{j!=k} E[cos^4(theta_j/2)]``, estimated as
    the product of independent per-factor sample means over ``theta`` of shape ``(N, n)``.

    Same estimand as ``Var_theta[dC_G/dtheta_k]`` (the mean is exactly 0 by symmetry), but a
    different estimator with polynomially bounded relative error (~sqrt(0.94 n / N)). It is
    NOT available for generic circuits, which have no such factorization; it is used only to
    check the theory where the generic sample variance is unresolvable.
    """
    theta = np.asarray(theta, dtype=np.float64)
    if theta.ndim != 2:
        raise ValueError("theta must be (N, n)")
    c4 = np.mean(np.cos(theta / 2.0) ** 4, axis=0)
    s2 = float(np.mean(np.sin(theta[:, k]) ** 2))
    return 0.25 * s2 * float(np.prod(np.delete(c4, k)))
