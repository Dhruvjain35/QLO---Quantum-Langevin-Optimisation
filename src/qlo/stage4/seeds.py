"""Seed policy for Stage 4.

* TUNING start seeds 0-39 and EVALUATION start seeds 1000-1099 are disjoint; evaluation
  seeds are never used to select hyperparameters.
* The initial theta is drawn iid U[-pi, pi] from ``start_rng(start_seed, n)``.
* Every stochastic trajectory gets its own noise stream from
  ``noise_seed(method, n, start_seed, replicate, config)``; methods use distinct ids so
  their draws are never shared, while the same (start_seed, replicate) pairs replicates
  across methods by id.
"""

from __future__ import annotations

import numpy as np

STAGE = 4
TUNING_START_SEEDS = tuple(range(0, 40))
EVALUATION_START_SEEDS = tuple(range(1000, 1100))
REPRESENTATIVE_START_SEEDS = (1000, 1001, 1002)  # predetermined, for trajectory figures/data

METHOD_IDS = {"exact_gd": 1, "shot_sgd": 2, "matched_gaussian": 3, "constant_langevin": 4, "noise_only": 5}


def start_rng(start_seed: int, n_qubits: int) -> np.random.Generator:
    return np.random.default_rng(np.random.SeedSequence([STAGE, 100, int(start_seed), int(n_qubits)]))


def initial_theta(start_seed: int, n_qubits: int) -> np.ndarray:
    return start_rng(start_seed, n_qubits).uniform(-np.pi, np.pi, size=int(n_qubits))


def noise_seed(method: str, n_qubits: int, start_seed: int, replicate: int, config_id: int = 0) -> np.random.SeedSequence:
    return np.random.SeedSequence([STAGE, 200, METHOD_IDS[method], int(n_qubits), int(start_seed), int(replicate), int(config_id)])
