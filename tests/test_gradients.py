"""Gradient validation: three mutually independent paths must agree.

  autograd (backprop)  vs  parameter-shift : tight, analytic-mode tolerance
  autograd (backprop)  vs  finite-difference : looser (O(h^2) truncation)
"""

import numpy as np
import pytest

from qlo import (
    COST_NAMES,
    HardwareEfficientAnsatz,
    cost_value,
    gradient_autograd,
    gradient_finite_difference,
    gradient_parameter_shift,
    make_rng,
)

ATOL_PS = 1e-10  # analytic mode: both are exact up to float rounding
ATOL_FD = 1e-6   # central differences with h=1e-5


@pytest.mark.parametrize("cost", COST_NAMES)
def test_gradient_shape_matches_params(small_ansatz, small_params, cost):
    for fn in (gradient_autograd, gradient_parameter_shift, gradient_finite_difference):
        g = fn(small_ansatz, cost, small_params)
        assert g.shape == small_params.shape == small_ansatz.param_shape
        assert np.all(np.isfinite(g))


@pytest.mark.parametrize("cost", COST_NAMES)
@pytest.mark.parametrize("seed", [0, 1, 2])
def test_autograd_matches_parameter_shift_n2_d1(cost, seed):
    a = HardwareEfficientAnsatz(2, 1)
    p = a.init_params(make_rng(seed))
    g_ag = gradient_autograd(a, cost, p)
    g_ps = gradient_parameter_shift(a, cost, p)
    diff = np.max(np.abs(g_ag - g_ps))
    assert diff < ATOL_PS, f"autograd vs parameter-shift max|diff|={diff:.3e}"


@pytest.mark.parametrize("cost", COST_NAMES)
def test_autograd_matches_finite_difference_n2_d1(cost):
    a = HardwareEfficientAnsatz(2, 1)
    p = a.init_params(make_rng(42))
    g_ag = gradient_autograd(a, cost, p)
    g_fd = gradient_finite_difference(a, cost, p)
    diff = np.max(np.abs(g_ag - g_fd))
    assert diff < ATOL_FD, f"autograd vs finite-difference max|diff|={diff:.3e}"


@pytest.mark.parametrize("cost", COST_NAMES)
@pytest.mark.parametrize("n,depth,entangler", [(3, 2, "ring"), (3, 2, "chain"), (4, 3, "ring")])
def test_three_paths_agree_larger_circuits(cost, n, depth, entangler):
    a = HardwareEfficientAnsatz(n, depth, entangler)
    p = a.init_params(make_rng(n * 100 + depth))
    g_ag = gradient_autograd(a, cost, p)
    g_ps = gradient_parameter_shift(a, cost, p)
    g_fd = gradient_finite_difference(a, cost, p)
    assert np.max(np.abs(g_ag - g_ps)) < ATOL_PS
    assert np.max(np.abs(g_ag - g_fd)) < ATOL_FD


def test_gradient_is_nontrivial_for_small_circuit(small_ansatz, small_params):
    # Guard against a silently-zero gradient path passing the agreement tests.
    g = gradient_autograd(small_ansatz, "global", small_params)
    assert np.linalg.norm(g) > 1e-3


def test_manual_parameter_shift_matches_autograd(small_ansatz, small_params):
    """Hand-rolled two-term shift rule on the scalar cost, no PennyLane gradient code."""
    p = small_params
    manual = np.zeros_like(p)
    for idx in np.ndindex(p.shape):
        plus, minus = p.copy(), p.copy()
        plus[idx] += np.pi / 2
        minus[idx] -= np.pi / 2
        manual[idx] = 0.5 * (cost_value(small_ansatz, "global", plus) - cost_value(small_ansatz, "global", minus))
    g_ag = gradient_autograd(small_ansatz, "global", p)
    assert np.max(np.abs(manual - g_ag)) < ATOL_PS
