"""Update rules. Each object exposes ``step(theta) -> theta_new`` and a ``resources`` dict.

A  ExactGD           theta - eta g(theta)                                    (analytic oracle)
B  ShotSGD           theta - eta g_hat(theta; M)      finite-shot binomial parameter shift
C  MatchedGaussian   theta - eta [g(theta) + eps],  eps_k ~ N(0, v_k(theta, M))   (analytic oracle)
D  ConstantLangevin  theta - eta g(theta) + sigma z,  z ~ N(0, I);  D_coef = sigma^2 / (2 eta)  (analytic oracle)
E  NoiseOnly         theta - eta eps,  eps_k ~ N(0, v_k(theta, M))                (diagnostic, not an optimizer)

Every rule wraps theta to [-pi, pi) after the update (exact symmetry of the benchmark).
No rule receives anything from the offline evaluator; B sees only finite-shot estimates.
Only B's resources are hardware-meaningful; A, C, D, E count oracle calls for bookkeeping.
"""

from __future__ import annotations

import numpy as np

from qlo.stage4.landscape import exact_gradient, wrap_angles
from qlo.stage4.shots import BinomialProjectorShots, conditional_shot_variance

METHODS = ("exact_gd", "shot_sgd", "matched_gaussian", "constant_langevin", "noise_only")
HARDWARE_REALIZABLE = {"shot_sgd"}


class ExactGD:
    name = "exact_gd"

    def __init__(self, eta: float) -> None:
        self.eta = float(eta)
        self.oracle_gradient_calls = 0

    def step(self, theta: np.ndarray) -> np.ndarray:
        self.oracle_gradient_calls += 1
        return wrap_angles(theta - self.eta * exact_gradient(theta))

    @property
    def resources(self) -> dict:
        return {"circuit_evaluations": 0, "shots": 0, "oracle_gradient_calls": self.oracle_gradient_calls}


class ShotSGD:
    name = "shot_sgd"

    def __init__(self, n_qubits: int, eta: float, shots: int, seed) -> None:
        self.eta = float(eta)
        self.shots = int(shots)
        self.estimator = BinomialProjectorShots(n_qubits, shots, seed)

    def step(self, theta: np.ndarray) -> np.ndarray:
        return wrap_angles(theta - self.eta * self.estimator.gradient(theta))

    @property
    def resources(self) -> dict:
        return {"circuit_evaluations": self.estimator.circuit_evaluations, "shots": self.estimator.total_shots,
                "oracle_gradient_calls": 0}


class MatchedGaussian:
    name = "matched_gaussian"

    def __init__(self, eta: float, shots: int, seed) -> None:
        self.eta = float(eta)
        self.shots = int(shots)
        self.rng = np.random.default_rng(seed)
        self.oracle_gradient_calls = 0

    def step(self, theta: np.ndarray) -> np.ndarray:
        v = conditional_shot_variance(theta, self.shots)
        eps = self.rng.normal(0.0, np.sqrt(v))
        self.oracle_gradient_calls += 1
        return wrap_angles(theta - self.eta * (exact_gradient(theta) + eps))

    @property
    def resources(self) -> dict:
        return {"circuit_evaluations": 0, "shots": 0, "oracle_gradient_calls": self.oracle_gradient_calls}


class ConstantLangevin:
    name = "constant_langevin"

    def __init__(self, eta: float, sigma_step: float, seed) -> None:
        self.eta = float(eta)
        self.sigma_step = float(sigma_step)
        self.rng = np.random.default_rng(seed)
        self.oracle_gradient_calls = 0

    @property
    def diffusion_coefficient(self) -> float:
        """Conventional ``D = sigma_step^2 / (2 eta)``. No temperature/FDT relation is assumed.
        ``eta = 0`` is pure diffusion (diagnostic control); D is then undefined (nan)."""
        return self.sigma_step**2 / (2.0 * self.eta) if self.eta > 0 else float("nan")

    def step(self, theta: np.ndarray) -> np.ndarray:
        self.oracle_gradient_calls += 1
        z = self.rng.standard_normal(theta.size)
        return wrap_angles(theta - self.eta * exact_gradient(theta) + self.sigma_step * z)

    @property
    def resources(self) -> dict:
        return {"circuit_evaluations": 0, "shots": 0, "oracle_gradient_calls": self.oracle_gradient_calls}


class NoiseOnly:
    name = "noise_only"

    def __init__(self, eta: float, shots: int, seed) -> None:
        self.eta = float(eta)
        self.shots = int(shots)
        self.rng = np.random.default_rng(seed)
        self.oracle_gradient_calls = 0  # it calls the variance oracle, not the gradient

    def step(self, theta: np.ndarray) -> np.ndarray:
        v = conditional_shot_variance(theta, self.shots)
        return wrap_angles(theta - self.eta * self.rng.normal(0.0, np.sqrt(v)))

    @property
    def resources(self) -> dict:
        return {"circuit_evaluations": 0, "shots": 0, "oracle_gradient_calls": 0}


def make_method(name: str, n_qubits: int, hp: dict, seed=None):
    """Factory. ``hp`` holds eta / shots / sigma_step as needed."""
    if name == "exact_gd":
        return ExactGD(hp["eta"])
    if name == "shot_sgd":
        return ShotSGD(n_qubits, hp["eta"], hp["shots"], seed)
    if name == "matched_gaussian":
        return MatchedGaussian(hp["eta"], hp["shots"], seed)
    if name == "constant_langevin":
        return ConstantLangevin(hp["eta"], hp["sigma_step"], seed)
    if name == "noise_only":
        return NoiseOnly(hp["eta"], hp["shots"], seed)
    raise ValueError(f"unknown method {name!r}; expected one of {METHODS}")
