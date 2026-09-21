"""Vectorized sampling of (log A_k, s_k) over theta ~ U[-pi, pi]^n, in log space."""

from __future__ import annotations

import numpy as np

from qlo.stage5.theory import LOG2, log_a_k, s_k


def sample_theta(n_qubits: int, n_samples: int, seed) -> np.ndarray:
    rng = np.random.default_rng(np.random.SeedSequence([5, 100, int(n_qubits), int(seed)]))
    return rng.uniform(-np.pi, np.pi, size=(int(n_samples), int(n_qubits)))


def sample_log_a_s(n_qubits: int, n_samples: int, seed, k: int = 0, chunk: int = 20_000) -> tuple[np.ndarray, np.ndarray]:
    """Returns (log A_k, s_k) for n_samples theta draws without holding all theta at once."""
    logA, s = np.empty(n_samples), np.empty(n_samples)
    rng = np.random.default_rng(np.random.SeedSequence([5, 100, int(n_qubits), int(seed)]))
    for i in range(0, n_samples, chunk):
        m = min(chunk, n_samples - i)
        th = rng.uniform(-np.pi, np.pi, size=(m, n_qubits))
        logA[i:i + m] = log_a_k(th, k)
        s[i:i + m] = s_k(th, k)
    return logA, s


def theory_log_a_stats(n_qubits: int) -> dict:
    """Exact moments of log A_k: mean -2(n-1) log 2; Var = (n-1) Var[log cos^2(theta/2)] = (n-1) pi^2/3."""
    return {"mean_logA": -2.0 * (n_qubits - 1) * LOG2, "var_logA": (n_qubits - 1) * np.pi**2 / 3.0,
            "mean_A": 0.5 ** (n_qubits - 1), "geometric_A": 4.0 ** (-(n_qubits - 1))}
