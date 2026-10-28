"""Per-theta required shots for all four estimators on ONE Stage 7 RX-product theta array (exactly paired).

Loschmidt / ancilla SWAP: the frozen Stage 7 functions. Destructive SWAP: its own +-1 model, E[Z] built from explicit
2-qubit Bell-circuit probabilities on every qubit pair (product states). Hadamard: complex shifted amplitudes, exact
variance, exact integer inversion. Every estimator sees the same theta, the same F_pm and the same g = A s / 2.
"""

from __future__ import annotations

import numpy as np

from qlo.b2 import destructive_swap as DS
from qlo.b2 import hadamard_theory as HT
from qlo.stage5.theory import log_a_k, s_k
from qlo.stage7 import analysis as AN
from qlo.stage7 import required_shots as R

LN10 = float(np.log(10.0))


def ceil_log10(l) -> np.ndarray:
    """Stage 7 integer budget convention: log10 max(ceil M, 1) for log10 M < 15, else max(log10 M, 0)."""
    l = np.asarray(l, dtype=np.float64)
    small = l < 15.0
    out = np.maximum(l, 0.0)
    out[small] = np.log10(np.maximum(np.ceil(np.power(10.0, l[small])), 1.0))
    return out


def _rx_states(theta: np.ndarray) -> np.ndarray:
    """Single-qubit states RX(theta)|0> = [cos(t/2), -i sin(t/2)], shape theta.shape + (2,)."""
    return np.stack([np.cos(theta / 2.0) + 0j, -1j * np.sin(theta / 2.0)], axis=-1)


def shifted_quantities(theta: np.ndarray, k: int = 0) -> dict:
    """F_pm (closed form), destructive-SWAP E[Z_pm] (Bell circuits) and complex amplitudes a_pm (x_pm, y_pm)."""
    logA, s = log_a_k(theta, k), s_k(theta, k)
    A = np.exp(logA)
    out = {"logA": logA, "s": s, "delta": A * s, "F_plus": A * (1 - s) / 2, "F_minus": A * (1 + s) / 2}
    zero = np.array([1.0 + 0j, 0j])
    for tag, sh in (("plus", np.pi / 2), ("minus", -np.pi / 2)):
        th = theta.copy()
        th[:, k] += sh
        psi = _rx_states(th)                                            # (N, n, 2)
        e_pairs = DS.bell_pair_expectation(psi, np.broadcast_to(zero, psi.shape))
        with np.errstate(divide="ignore"):
            out[f"log_e_{tag}"] = np.sum(np.log(e_pairs), axis=1)       # E[Z] = prod over pairs (all >= 0 here)
        amp = psi[..., 0]                                               # <0|psi_j> per qubit, complex
        with np.errstate(divide="ignore"):
            log_mag = np.sum(np.log(np.abs(amp)), axis=1)
        phase = np.sum(np.angle(amp), axis=1)
        mag = np.exp(log_mag)
        out[f"x_{tag}"], out[f"y_{tag}"] = mag * np.cos(phase), mag * np.sin(phase)
    return out


def destructive_swap_log10_m(q: dict, rho: float) -> np.ndarray:
    e_p, e_m = np.exp(q["log_e_plus"]), np.exp(q["log_e_minus"])
    with np.errstate(divide="ignore"):
        l = 2 * np.log10(rho) + np.log10(2.0 - e_p**2 - e_m**2) - 2 * np.log10(np.abs(e_m - e_p))
    return ceil_log10(l)


def required_all(n: int, n_samples: int, seed: int, stream: int, rho: float) -> dict:
    theta = AN.sample_theta(n, n_samples, seed, stream)
    q = shifted_quantities(theta)
    m = {"loschmidt": R.snr_shots("loschmidt", q["logA"], q["s"], rho),
         "swap_ancilla": R.snr_shots("swap", q["logA"], q["s"], rho),
         "swap_destructive": destructive_swap_log10_m(q, rho),
         "hadamard_u": np.log10(HT.required_m_ht(q["x_plus"], q["y_plus"], q["x_minus"], q["y_minus"], q["delta"], rho))}
    return {"theta_q": q, "log10_m_internal": m}
