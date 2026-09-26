"""Loschmidt (projector) and SWAP-test finite-shot parameter-shift estimators on ONE fidelity landscape.

Common landscape (Stage 3-5, frozen): for component k, A = prod_{j!=k} cos^2(theta_j/2), s = sin(theta_k),

    F_+ = F(theta + pi/2 e_k) = A(1 - s)/2,     F_- = A(1 + s)/2,     g = dC/dtheta_k = A s / 2.

LOSCHMIDT. Each shot is Bernoulli(F). K_+- ~ Bin(M, F_+-), F_hat = K/M, C_hat = 1 - F_hat,
    g_hat_LE = (C_hat_+ - C_hat_-)/2 = (K_- - K_+)/(2M)
    Var = [F_+(1-F_+) + F_-(1-F_-)]/(4M) = [A - A^2(1+s^2)/2]/(4M)          (Stage 5)
    SNR^2 = M A s^2 / [1 - A(1+s^2)/2]

SWAP. Ancilla Z = +1 with probability q = (1 + F)/2 (E[Z] = F for two pure states). K counts +1 outcomes,
F_hat = 2K/M - 1, C_hat = 1 - F_hat = 2 - 2K/M, so
    g_hat_SWAP = (C_hat_+ - C_hat_-)/2 = (K_- - K_+)/M                          (sign checked in tests)
    E[g_hat] = q_- - q_+ = (F_- - F_+)/2 = A s / 2 = g
    Var(F_hat) = 4 q(1-q)/M = (1 - F^2)/M
    Var(g_hat) = [q_+(1-q_+) + q_-(1-q_-)]/M = [2 - F_+^2 - F_-^2]/(4M) = [2 - A^2(1+s^2)/2]/(4M)
    SNR^2 = M A^2 s^2 / [2 - A^2(1+s^2)/2]

Both are  g_hat = c (K_pos - K_neg)/M  with independent binomial counts:
    loschmidt  p_pos = F_-,  p_neg = F_+,  c = 1/2
    swap       p_pos = q_-,  p_neg = q_+,  c = 1
and in both cases  c (p_pos - p_neg) = g  exactly.

Precision: for SWAP both q are 1/2 + O(A). Once A|s| < ~1e-16 the float difference q_- - q_+ is lost next to
1/2, so every SWAP quantity that depends on the DIFFERENCE is computed from the exact difference
``delta = A s / 2`` (``pair_delta``) or in log space, never by subtracting the two probabilities.
"""

from __future__ import annotations

import numpy as np

from qlo.stage5.theory import f_plus_minus

SCHEMES = ("loschmidt", "swap")
SCALE = {"loschmidt": 0.5, "swap": 1.0}
LN10 = float(np.log(10.0))


# ---- outcome probabilities -----------------------------------------------------------------------
def shifted_fidelities(A, s):
    return f_plus_minus(A, s)


def outcome_probabilities(scheme: str, A, s):
    """``(p_plus_shift, p_minus_shift)``: success probability of ONE shot at theta +- pi/2 e_k.
    Loschmidt: F_+-; SWAP: q_+- = (1 + F_+-)/2."""
    fp, fm = shifted_fidelities(A, s)
    if scheme == "loschmidt":
        return fp, fm
    if scheme == "swap":
        return (1.0 + fp) / 2.0, (1.0 + fm) / 2.0
    raise ValueError(scheme)


def pair(scheme: str, A, s):
    """``(p_pos, p_neg, c)`` with ``g_hat = c (K_pos - K_neg)/M`` and ``c (p_pos - p_neg) = g``."""
    pp, pm = outcome_probabilities(scheme, A, s)
    return pm, pp, SCALE[scheme]  # K_pos is the MINUS-shift count for both schemes


def pair_delta(scheme: str, A, s):
    """Exact ``p_pos - p_neg`` without cancellation: Loschmidt ``A s``; SWAP ``A s / 2``."""
    A, s = np.asarray(A, dtype=np.float64), np.asarray(s, dtype=np.float64)
    return A * s if scheme == "loschmidt" else A * s / 2.0


# ---- estimators from counts (built from the C_hat definitions, not from the simplified forms) --------
def fidelity_estimate(scheme: str, K, M):
    K, M = np.asarray(K, dtype=np.float64), np.asarray(M, dtype=np.float64)
    return K / M if scheme == "loschmidt" else 2.0 * K / M - 1.0


def gradient_from_counts(scheme: str, k_plus, k_minus, M):
    """``(C_hat(theta + pi/2 e_k) - C_hat(theta - pi/2 e_k))/2`` with ``C_hat = 1 - F_hat``."""
    c_plus = 1.0 - fidelity_estimate(scheme, k_plus, M)
    c_minus = 1.0 - fidelity_estimate(scheme, k_minus, M)
    return 0.5 * (c_plus - c_minus)


def gradient_from_counts_simplified(scheme: str, k_plus, k_minus, M):
    """Simplified forms: LE ``(K_- - K_+)/(2M)``, SWAP ``(K_- - K_+)/M``."""
    k_plus, k_minus, M = (np.asarray(x, dtype=np.float64) for x in (k_plus, k_minus, M))
    return SCALE[scheme] * (k_minus - k_plus) / M


# ---- moments --------------------------------------------------------------------------------------
def fidelity_estimate_variance(scheme: str, F, M):
    F, M = np.asarray(F, dtype=np.float64), np.asarray(M, dtype=np.float64)
    return F * (1.0 - F) / M if scheme == "loschmidt" else (1.0 - F**2) / M


def gradient_variance(scheme: str, A, s, M):
    """Closed forms: LE ``[A - A^2(1+s^2)/2]/(4M)``; SWAP ``[2 - A^2(1+s^2)/2]/(4M)``."""
    A, s, M = (np.asarray(x, dtype=np.float64) for x in (A, s, M))
    sq = A**2 * (1.0 + s**2) / 2.0  # = F_+^2 + F_-^2
    return (A - sq) / (4.0 * M) if scheme == "loschmidt" else (2.0 - sq) / (4.0 * M)


def gradient_variance_generic(scheme: str, A, s, M):
    """From the outcome probabilities: ``c^2 [p_pos(1-p_pos) + p_neg(1-p_neg)]/M``."""
    p_pos, p_neg, c = pair(scheme, A, s)
    return c**2 * (p_pos * (1 - p_pos) + p_neg * (1 - p_neg)) / np.asarray(M, dtype=np.float64)


def snr_squared(scheme: str, A, s, M):
    """LE ``M A s^2/[1 - A(1+s^2)/2]``; SWAP ``M A^2 s^2/[2 - A^2(1+s^2)/2]``. 0 at s = 0 or A = 0."""
    A, s, M = (np.asarray(x, dtype=np.float64) for x in (A, s, M))
    if scheme == "loschmidt":
        num, den = M * A * s**2, 1.0 - A * (1.0 + s**2) / 2.0
    else:
        num, den = M * A**2 * s**2, 2.0 - A**2 * (1.0 + s**2) / 2.0
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(den > 0, num / den, np.where(num > 0, np.inf, 0.0))


def log10_shots_for_snr(scheme: str, logA, s, rho):
    """Log-space continuous M for SNR = rho (integer budget = ceil). Stable for A down to exp(-700)."""
    logA, s, rho = (np.asarray(x, dtype=np.float64) for x in (logA, s, rho))
    A = np.exp(logA)
    with np.errstate(divide="ignore", invalid="ignore"):
        if scheme == "loschmidt":
            return 2 * np.log10(rho) + np.log10(1.0 - A * (1.0 + s**2) / 2.0) - logA / LN10 - np.log10(s**2)
        return 2 * np.log10(rho) + np.log10(2.0 - A**2 * (1.0 + s**2) / 2.0) - 2 * logA / LN10 - np.log10(s**2)


# ---- resource accounting (Task 16) -----------------------------------------------------------------
def resource_accounting(scheme: str, n_qubits: int, shots_per_circuit: int, n_components: int = 1) -> dict:
    """Abstract per-gradient accounting. Shot count is NOT physical cost: the SWAP test needs an ancilla, a
    second n-qubit register holding the target state and n controlled-SWAPs per shot. No gate-level claims."""
    shots = 2 * shots_per_circuit * n_components  # two shifted circuits per component
    if scheme == "loschmidt":
        return {"scheme": scheme, "measurement_shots": shots, "variational_state_preparations": shots,
                "target_state_preparations": 0, "state_copies_per_shot": 1, "total_state_copies": shots,
                "qubits_per_shot": n_qubits, "controlled_swaps_per_shot": 0,
                "note": "target |0^n> is the computational-basis all-zero outcome: overlap read directly"}
    if scheme == "swap":
        return {"scheme": scheme, "measurement_shots": shots, "variational_state_preparations": shots,
                "target_state_preparations": shots, "state_copies_per_shot": 2, "total_state_copies": 2 * shots,
                "qubits_per_shot": 2 * n_qubits + 1, "controlled_swaps_per_shot": n_qubits,
                "note": "standard ancilla SWAP test: ancilla + variational register + target register"}
    raise ValueError(scheme)
