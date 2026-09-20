"""Analytic finite-shot variance for the parameter-shift estimator — GLOBAL PAULI COST ONLY.

For the observable ``O = Z_0 Z_1 ... Z_{n-1}`` every single-shot outcome is ±1, so
an ``M``-shot sample mean ``C_hat`` of a state with exact expectation ``C`` has

    Var(C_hat) = (1 - C^2) / M .

With independent ``+pi/2`` and ``-pi/2`` evaluations,

    Var(g_hat_k) = 1/4 [ Var(C_hat_+) + Var(C_hat_-) ]
                 = [ 2 - C_+^2 - C_-^2 ] / (4 M) ,

where ``C_±`` are the EXACT expectation values at the shifted parameters
(computed here with the analytic Stage 1 reference — this module is validation
code and is allowed to use it).

This derivation relies on ±1 single-shot outcomes. It is NOT valid for the
local cost ``(1/n) sum_i Z_i`` (whose single-shot value is a mean of n ±1's and
whose variance involves the Z_i Z_j correlations) or for arbitrary observables,
so any other ``cost`` is rejected rather than silently mis-applied.
"""

from __future__ import annotations

import numpy as np

from qlo.circuits.hardware_efficient import HardwareEfficientAnsatz
from qlo.gradients.exact import cost_value
from qlo.gradients.finite_shot import SHIFT
from qlo.utils.indexing import shifted_params

_SUPPORTED = ("global",)


def single_shot_variance_pm1(c_exact: float, shots: int) -> float:
    """``(1 - C^2)/M`` for a ±1-valued observable with exact expectation ``C``."""
    if not -1.0 - 1e-12 <= c_exact <= 1.0 + 1e-12:
        raise ValueError(f"C must lie in [-1, 1] for a ±1 observable, got {c_exact}")
    return max(0.0, 1.0 - float(c_exact) ** 2) / int(shots)


def theoretical_parameter_shift_variance(
    ansatz: HardwareEfficientAnsatz, params, k: int, shots: int, cost: str = "global"
) -> float:
    """``Var(g_hat_k) = [2 - C_+^2 - C_-^2] / (4M)`` using exact shifted expectations."""
    if cost not in _SUPPORTED:
        raise ValueError(
            f"analytic shot variance is derived only for the ±1-valued global Pauli product; got cost={cost!r}"
        )
    if int(shots) < 1:
        raise ValueError("shots must be >= 1")
    c_plus = cost_value(ansatz, cost, shifted_params(params, k, +SHIFT))
    c_minus = cost_value(ansatz, cost, shifted_params(params, k, -SHIFT))
    return 0.25 * (single_shot_variance_pm1(c_plus, shots) + single_shot_variance_pm1(c_minus, shots))


def theoretical_parameter_shift_variance_vector(
    ansatz: HardwareEfficientAnsatz, params, shots: int, cost: str = "global"
) -> np.ndarray:
    """Length-``p`` vector of ``theoretical_parameter_shift_variance`` over all flat ``k``."""
    p = int(np.asarray(params).size)
    return np.array([theoretical_parameter_shift_variance(ansatz, params, k, shots, cost) for k in range(p)])
