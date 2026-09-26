"""The standard barren-plateau diagnostic: ``Var_theta[ dC/dtheta_k ]``.

For a fixed flat parameter index ``k``, draw ``n_init`` independent random
initializations ``theta^(1..n_init)``, evaluate the exact gradient component
``g_k(theta^(i))`` at each, and take the variance **across initializations**
(``ddof=1``). Doing this for every ``k`` gives a length-``p`` vector; summary
statistics *across k* are then reported separately.

This is distinct from — and must not be confused with — the within-gradient
statistic ``within_gradient_entry_var`` in ``bp_smoke`` (variance over the
entries of a single gradient vector at one theta), which is not a BP diagnostic.
Nothing here pools every entry of every gradient into one variance.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np

from qlo.circuits.hardware_efficient import HardwareEfficientAnsatz
from qlo.gradients.exact import gradient_autograd
from qlo.utils.seeding import make_rng


def gradient_samples_across_inits(
    ansatz: HardwareEfficientAnsatz,
    cost: str,
    n_init: int,
    seed,
    gradient_fn: Callable = gradient_autograd,
) -> np.ndarray:
    """Matrix ``G`` of shape ``(n_init, p)``: row ``i`` is the flattened exact gradient at ``theta^(i)``."""
    if n_init < 2:
        raise ValueError("n_init must be >= 2 to form a variance across initializations")
    rng = make_rng(seed)
    rows = []
    for _ in range(int(n_init)):
        theta = ansatz.init_params(rng)
        rows.append(np.asarray(gradient_fn(ansatz, cost, theta), dtype=np.float64).ravel())
    return np.stack(rows)


def bp_variance_per_parameter(samples: np.ndarray) -> np.ndarray:
    """``Var_theta[g_k]`` for every ``k``; input ``(n_init, p)``, output ``(p,)``. ``ddof=1``."""
    samples = np.asarray(samples, dtype=np.float64)
    if samples.ndim != 2 or samples.shape[0] < 2:
        raise ValueError("samples must be (n_init >= 2, p)")
    return np.var(samples, axis=0, ddof=1)


def bp_variance_summary(per_k: np.ndarray) -> dict[str, float]:
    """Summary across parameter indices of the per-k variances (second step, kept separate)."""
    v = np.asarray(per_k, dtype=np.float64)
    return {
        "n_params": int(v.size),
        "mean_var_k": float(np.mean(v)),
        "median_var_k": float(np.median(v)),
        "min_var_k": float(np.min(v)),
        "max_var_k": float(np.max(v)),
        "std_var_k": float(np.std(v, ddof=1)) if v.size > 1 else 0.0,
    }


def barren_plateau_variance(
    ansatz: HardwareEfficientAnsatz,
    cost: str,
    n_init: int,
    seed,
    k: int | None = None,
    gradient_fn: Callable = gradient_autograd,
) -> dict:
    """Convenience wrapper.

    ``k=None``: returns ``{"per_k": (p,) array, "summary": {...}, "samples": (n_init, p)}``.
    ``k=int`` : returns ``{"k": k, "var_k": float, "samples_k": (n_init,)}``.
    """
    G = gradient_samples_across_inits(ansatz, cost, n_init, seed, gradient_fn)
    if k is not None:
        col = G[:, int(k)]
        return {"k": int(k), "n_init": int(n_init), "var_k": float(np.var(col, ddof=1)), "samples_k": col}
    per_k = bp_variance_per_parameter(G)
    return {"n_init": int(n_init), "per_k": per_k, "summary": bp_variance_summary(per_k), "samples": G}
