"""Global Pauli-parity benchmark (new in Stage 6).

Circuit (identical to Stages 3-5): |psi(theta)> = prod_j RX(theta_j) |0^n>, theta_j ~ iid U[-pi, pi].

    mu_Z(theta) = <Z_1 ... Z_n> = prod_j cos(theta_j)
    C_Z(theta)  = [1 - mu_Z(theta)] / 2                        in [0, 1]

For component k, with B_k = prod_{j!=k} cos(theta_j) and s_k = sin(theta_k):

    g_{Z,k} = dC_Z/dtheta_k = (1/2) s_k B_k
    E_theta[g_{Z,k}] = 0
    Var_theta[g_{Z,k}] = (1/4) E[sin^2] prod_{j!=k} E[cos^2] = (1/4)(1/2)(1/2)^{n-1} = 2^{-(n+2)}

so this is an exponential barren plateau, like the projector cost but with a different
constant (projector: (1/8)(3/8)^{n-1}).

Typical (log) scale: E[log|cos theta|] = E[log|sin theta|] = -log 2 for theta ~ U[-pi, pi],
hence E[log|s_k B_k|] = -n log 2 and the log-typical gradient magnitude is 2^{-(n+1)},
while the RMS magnitude is sqrt(Var) = 2^{-(n+2)/2} = 2^{-n/2-1}. These differ
exponentially: |g| is heavy-tailed on the log scale (see ``log_moments``).
"""

from __future__ import annotations

import numpy as np

LOG2 = float(np.log(2.0))
E_LOG_ABS_COS = -LOG2  # E_theta[log|cos theta|], theta ~ U[-pi, pi]
E_LOG_ABS_SIN = -LOG2  # E_theta[log|sin theta|]
VAR_LOG_ABS_COS = float(np.pi**2 / 12.0)  # Var[log|cos theta|] = pi^2/12 (quadrature: 0.822467)


def mu_z(theta: np.ndarray) -> np.ndarray:
    """``<Z^{otimes n}> = prod_j cos(theta_j)``; theta of shape (..., n)."""
    return np.prod(np.cos(np.asarray(theta, dtype=np.float64)), axis=-1)


def cost(theta: np.ndarray) -> np.ndarray:
    """``C_Z = [1 - mu_Z] / 2`` in [0, 1]."""
    return (1.0 - mu_z(theta)) / 2.0


def b_k(theta: np.ndarray, k: int) -> np.ndarray:
    """``B_k = prod_{j!=k} cos(theta_j)``."""
    c = np.cos(np.asarray(theta, dtype=np.float64))
    return np.prod(np.delete(c, k, axis=-1), axis=-1)


def log_abs_b_k(theta: np.ndarray, k: int) -> np.ndarray:
    """``log|B_k|`` as a sum of logs (no underflow at large n)."""
    with np.errstate(divide="ignore"):
        la = np.log(np.abs(np.cos(np.asarray(theta, dtype=np.float64))))
    return np.sum(np.delete(la, k, axis=-1), axis=-1)


def s_k(theta: np.ndarray, k: int) -> np.ndarray:
    return np.sin(np.asarray(theta, dtype=np.float64)[..., k])


def exact_gradient(theta: np.ndarray) -> np.ndarray:
    """Full gradient ``dC_Z/dtheta = (1/2) sin(theta_k) prod_{j!=k} cos(theta_j)`` (1-D theta).

    Uses prefix/suffix products (no division), so a parameter at exactly +-pi/2 is fine.
    """
    theta = np.asarray(theta, dtype=np.float64)
    c = np.cos(theta)
    n = c.size
    prefix, suffix = np.ones(n + 1), np.ones(n + 1)
    prefix[1:] = np.cumprod(c)
    suffix[:-1] = np.cumprod(c[::-1])[::-1]
    others = prefix[:n] * suffix[1:]
    return 0.5 * np.sin(theta) * others


def exact_gradient_k(s, B) -> np.ndarray:
    """``g_{Z,k} = s B / 2``."""
    return np.asarray(s, dtype=np.float64) * np.asarray(B, dtype=np.float64) / 2.0


def shifted_mu(s, B) -> tuple[np.ndarray, np.ndarray]:
    """``mu(theta +- pi/2 e_k) = -+ s_k B_k``. (cos(theta_k +- pi/2) = -+ sin theta_k.)"""
    sB = np.asarray(s, dtype=np.float64) * np.asarray(B, dtype=np.float64)
    return -sB, sB


def shifted_probabilities(s, B) -> tuple[np.ndarray, np.ndarray]:
    """Bernoulli cost probabilities at the two shifts: ``p_+- = (1 +- s B)/2``, both near 1/2."""
    sB = np.clip(np.asarray(s, dtype=np.float64) * np.asarray(B, dtype=np.float64), -1.0, 1.0)
    return (1.0 + sB) / 2.0, (1.0 - sB) / 2.0


def theoretical_gradient_variance(n_qubits: int) -> float:
    """``Var_theta[g_{Z,k}] = 2^{-(n+2)}``."""
    if n_qubits < 1:
        raise ValueError("n_qubits must be >= 1")
    return 2.0 ** (-(n_qubits + 2))


def log_moments(n_qubits: int) -> dict:
    """Typical vs RMS scales for ``|g| = |s B|/2`` (they differ exponentially)."""
    return {
        "E_log_abs_sB": -n_qubits * LOG2,
        "geometric_abs_sB": 2.0 ** (-n_qubits),
        "geometric_abs_g": 2.0 ** (-(n_qubits + 1)),
        "var_log_abs_sB": n_qubits * VAR_LOG_ABS_COS,
        "rms_g": float(np.sqrt(theoretical_gradient_variance(n_qubits))),
        "mean_abs_sB": (2.0 / np.pi) ** n_qubits,  # E|sin| = E|cos| = 2/pi
    }
