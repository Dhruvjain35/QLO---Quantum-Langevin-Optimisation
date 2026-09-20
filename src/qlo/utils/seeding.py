"""Centralized seeding.

All randomness in this package flows through an explicit ``numpy.random.Generator``
created by :func:`make_rng`. Nothing here touches the legacy global ``np.random``
state, and no module-level RNG is kept, so there is no hidden mutable state.

For the analytic (``shots=None``) reference path there is no simulator randomness.
When finite-shot devices are introduced later, the same seed should be passed to
the PennyLane device (``qml.device(..., seed=seed)``) so that the parameter draw and
the sampling noise are both reproducible from one integer.
"""

from __future__ import annotations

import numpy as np


def make_rng(seed: int | None) -> np.random.Generator:
    """Return a fresh, independent ``numpy.random.Generator``.

    ``seed=None`` yields an unseeded (non-reproducible) generator; this is allowed
    but every experiment entry point in this package should pass an explicit seed.
    """
    return np.random.default_rng(seed)


def random_params(
    shape: tuple[int, ...],
    rng: np.random.Generator,
    low: float = 0.0,
    high: float = 2.0 * np.pi,
) -> np.ndarray:
    """Draw a parameter tensor uniformly from ``[low, high)`` with the given shape.

    Uniform-on-``[0, 2π)`` is the conventional "random initialization" used in the
    barren-plateau literature; it is a choice, not a claim about which
    initializations induce plateaus.
    """
    return rng.uniform(low=low, high=high, size=shape).astype(np.float64)
