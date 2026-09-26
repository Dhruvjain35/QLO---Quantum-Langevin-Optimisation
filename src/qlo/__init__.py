"""qlo: Stage 1 experimental foundation for QLO / barren-plateau research.

Public API (re-exported for convenience):
    - HardwareEfficientAnsatz, n_params, param_shape
    - global_z_observable, local_z_observable, get_observable, COST_NAMES
    - make_cost, cost_value, gradient_autograd, gradient_parameter_shift, gradient_finite_difference
    - make_rng, random_params
"""

from qlo.circuits.hardware_efficient import HardwareEfficientAnsatz, n_params, param_shape
from qlo.costs.observables import COST_NAMES, get_observable, global_z_observable, local_z_observable
from qlo.gradients.exact import (
    cost_value,
    gradient_autograd,
    gradient_finite_difference,
    gradient_parameter_shift,
    make_cost,
)
from qlo.utils.seeding import make_rng, random_params

__all__ = [
    "HardwareEfficientAnsatz",
    "n_params",
    "param_shape",
    "COST_NAMES",
    "get_observable",
    "global_z_observable",
    "local_z_observable",
    "make_cost",
    "cost_value",
    "gradient_autograd",
    "gradient_parameter_shift",
    "gradient_finite_difference",
    "make_rng",
    "random_params",
]
