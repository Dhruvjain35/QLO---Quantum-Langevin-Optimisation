import numpy as np
import pytest

from qlo import HardwareEfficientAnsatz, make_rng


@pytest.fixture
def small_ansatz() -> HardwareEfficientAnsatz:
    return HardwareEfficientAnsatz(n_qubits=2, depth=1)


@pytest.fixture
def small_params(small_ansatz) -> np.ndarray:
    return small_ansatz.init_params(make_rng(1234))
