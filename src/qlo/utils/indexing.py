"""Flat <-> tensor parameter indexing helpers.

The ansatz stores parameters as ``(depth, n_qubits, 2)``. Estimators and
barren-plateau statistics address a single component by a *flat* index
``k in range(p)`` in C (row-major) order, i.e. ``params.ravel()[k]``.
These helpers are thin, tested wrappers around NumPy so the same convention is
used everywhere.
"""

from __future__ import annotations

import numpy as np


def flat_to_multi(k: int, shape: tuple[int, ...]) -> tuple[int, ...]:
    """Multi-index of flat index ``k`` in row-major order."""
    size = int(np.prod(shape))
    if not 0 <= k < size:
        raise IndexError(f"flat index {k} out of range for shape {shape} (size {size})")
    return tuple(int(i) for i in np.unravel_index(k, shape))


def multi_to_flat(idx: tuple[int, ...], shape: tuple[int, ...]) -> int:
    """Flat row-major index of multi-index ``idx``."""
    return int(np.ravel_multi_index(idx, shape))


def shifted_params(params: np.ndarray, k: int, delta: float) -> np.ndarray:
    """Copy of ``params`` with flat component ``k`` shifted by ``delta``."""
    out = np.array(params, dtype=np.float64, copy=True)
    out[flat_to_multi(k, out.shape)] += delta
    return out
