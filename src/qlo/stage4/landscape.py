"""Exact landscape of the Stage 3 global benchmark (NumPy only, no PennyLane).

    F(theta) = prod_j cos^2(theta_j/2)         C = 1 - F
    dC/dtheta_k = (1/2) sin(theta_k) prod_{j!=k} cos^2(theta_j/2)

``others(theta)[k] = prod_{j!=k} cos^2(theta_j/2)`` is formed with prefix/suffix products
(no division), so a parameter at exactly ±pi never gives 0/0.
"""

from __future__ import annotations

import numpy as np

SHIFT = np.pi / 2.0


def cos2_half(theta: np.ndarray) -> np.ndarray:
    return np.cos(np.asarray(theta, dtype=np.float64) / 2.0) ** 2


def fidelity(theta: np.ndarray) -> float:
    """``F = prod_j cos^2(theta_j/2)`` (scalar for a 1-D theta; vectorized over leading axes)."""
    return np.prod(cos2_half(theta), axis=-1)


def cost(theta: np.ndarray) -> float:
    return 1.0 - fidelity(theta)


def others_product(theta: np.ndarray) -> np.ndarray:
    """``prod_{j!=k} cos^2(theta_j/2)`` for every k, via prefix/suffix products. 1-D theta only."""
    c2 = cos2_half(theta)
    n = c2.size
    prefix = np.ones(n + 1)
    suffix = np.ones(n + 1)
    prefix[1:] = np.cumprod(c2)
    suffix[:-1] = np.cumprod(c2[::-1])[::-1]
    return prefix[:n] * suffix[1:]


def exact_gradient(theta: np.ndarray) -> np.ndarray:
    """``dC/dtheta`` for a 1-D theta."""
    theta = np.asarray(theta, dtype=np.float64)
    return 0.5 * np.sin(theta) * others_product(theta)


def shifted_fidelities(theta: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """``F(theta + pi/2 e_k)`` and ``F(theta - pi/2 e_k)`` for every k (each shape ``(n,)``)."""
    theta = np.asarray(theta, dtype=np.float64)
    oth = others_product(theta)
    f_plus = np.cos((theta + SHIFT) / 2.0) ** 2 * oth
    f_minus = np.cos((theta - SHIFT) / 2.0) ** 2 * oth
    return np.clip(f_plus, 0.0, 1.0), np.clip(f_minus, 0.0, 1.0)


def wrap_angles(theta: np.ndarray) -> np.ndarray:
    """Canonicalize to ``[-pi, pi)``: ``((theta + pi) mod 2pi) - pi``. Exact symmetry of F."""
    return (np.asarray(theta, dtype=np.float64) + np.pi) % (2.0 * np.pi) - np.pi
