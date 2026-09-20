import numpy as np
import pennylane as qml
import pytest

from qlo import COST_NAMES, HardwareEfficientAnsatz, cost_value, get_observable, global_z_observable, local_z_observable, make_rng


def test_cost_names_exposed():
    assert COST_NAMES == ("global", "local")
    for name in COST_NAMES:
        get_observable(name, 3)
    with pytest.raises(ValueError):
        get_observable("nonlocal", 3)


def test_global_observable_matrix_is_z_tensor_product():
    n = 3
    z = np.diag([1.0, -1.0])
    expected = np.kron(np.kron(z, z), z)
    got = qml.matrix(global_z_observable(n), wire_order=range(n))
    assert np.allclose(got, expected)


def test_local_observable_matrix_is_mean_z():
    n = 3
    z = np.diag([1.0, -1.0])
    eye = np.eye(2)
    terms = [np.kron(np.kron(z, eye), eye), np.kron(np.kron(eye, z), eye), np.kron(np.kron(eye, eye), z)]
    expected = sum(terms) / n
    got = qml.matrix(local_z_observable(n), wire_order=range(n))
    assert np.allclose(got, expected)


def test_local_cost_equals_mean_of_single_qubit_z():
    """Independent check: build <Z_i> one qubit at a time and average."""
    a = HardwareEfficientAnsatz(3, 2)
    params = a.init_params(make_rng(7))
    dev = qml.device("default.qubit", wires=3, shots=None)

    @qml.qnode(dev)
    def single_z(p, i):
        a.apply(p)
        return qml.expval(qml.PauliZ(i))

    manual = np.mean([float(single_z(params, i)) for i in range(3)])
    assert cost_value(a, "local", params) == pytest.approx(manual, abs=1e-12)


def test_global_cost_equals_parity_expectation():
    """Independent check: <Z...Z> from probabilities = sum_b (-1)^{popcount(b)} P(b)."""
    a = HardwareEfficientAnsatz(3, 2)
    params = a.init_params(make_rng(8))
    dev = qml.device("default.qubit", wires=3, shots=None)

    @qml.qnode(dev)
    def probs(p):
        a.apply(p)
        return qml.probs(wires=range(3))

    pr = np.asarray(probs(params))
    parity = np.array([(-1) ** bin(b).count("1") for b in range(2**3)])
    assert cost_value(a, "global", params) == pytest.approx(float(parity @ pr), abs=1e-12)


def test_costs_bounded_over_random_draws():
    a = HardwareEfficientAnsatz(4, 2)
    rng = make_rng(99)
    for _ in range(5):
        p = a.init_params(rng)
        for c in COST_NAMES:
            v = cost_value(a, c, p)
            assert np.isfinite(v) and -1.0 - 1e-12 <= v <= 1.0 + 1e-12
