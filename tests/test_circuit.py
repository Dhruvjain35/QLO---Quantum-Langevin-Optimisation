import numpy as np
import pennylane as qml
import pytest

from qlo import HardwareEfficientAnsatz, cost_value, make_rng, n_params, param_shape


@pytest.mark.parametrize("n,depth", [(1, 1), (2, 1), (3, 2), (4, 3)])
def test_param_shape_and_count(n, depth):
    a = HardwareEfficientAnsatz(n, depth)
    assert a.param_shape == (depth, n, 2) == param_shape(n, depth)
    assert a.n_params == 2 * n * depth == n_params(n, depth)
    assert a.init_params(make_rng(0)).shape == a.param_shape


def test_entangling_topology():
    assert HardwareEfficientAnsatz(4, 1, "ring").entangling_pairs() == [(0, 1), (1, 2), (2, 3), (3, 0)]
    assert HardwareEfficientAnsatz(4, 1, "chain").entangling_pairs() == [(0, 1), (1, 2), (2, 3)]
    assert HardwareEfficientAnsatz(1, 1, "ring").entangling_pairs() == []


def test_gate_sequence_matches_spec():
    a = HardwareEfficientAnsatz(2, 1, "ring")
    params = np.zeros(a.param_shape)
    with qml.tape.QuantumTape() as tape:
        a.apply(params)
    names = [op.name for op in tape.operations]
    assert names == ["RY", "RZ", "RY", "RZ", "CNOT", "CNOT"]


def test_invalid_config_rejected():
    with pytest.raises(ValueError):
        HardwareEfficientAnsatz(0, 1)
    with pytest.raises(ValueError):
        HardwareEfficientAnsatz(2, 0)
    with pytest.raises(ValueError):
        HardwareEfficientAnsatz(2, 1, "star")


def test_wrong_param_shape_rejected(small_ansatz):
    with pytest.raises(ValueError):
        cost_value(small_ansatz, "global", np.zeros((2, 2, 2)))


@pytest.mark.parametrize("cost", ["global", "local"])
def test_cost_is_finite_scalar_in_range(small_ansatz, small_params, cost):
    v = cost_value(small_ansatz, cost, small_params)
    assert isinstance(v, float)
    assert np.isfinite(v)
    assert -1.0 - 1e-12 <= v <= 1.0 + 1e-12


def test_zero_params_gives_all_zero_state():
    # theta = 0 leaves |0...0>, CNOTs act trivially, so <Z_i> = 1 for all i.
    a = HardwareEfficientAnsatz(3, 2)
    z = np.zeros(a.param_shape)
    assert cost_value(a, "global", z) == pytest.approx(1.0, abs=1e-12)
    assert cost_value(a, "local", z) == pytest.approx(1.0, abs=1e-12)
