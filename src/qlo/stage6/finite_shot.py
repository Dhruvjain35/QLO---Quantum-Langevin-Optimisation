"""Finite-shot parameter-shift estimator for the global parity cost.

A single measurement of ``Z^{otimes n}`` yields X in {-1,+1} with E[X] = mu_Z, equivalently
a Bernoulli cost outcome Y = (1-X)/2 with P(Y=1) = C_Z. At the two shifts of theta_k,
mu(theta +- pi/2 e_k) = -+ s_k B_k, so

    p_+ = C_Z(theta + pi/2 e_k) = (1 + s_k B_k)/2        p_- = (1 - s_k B_k)/2

With K_+ ~ Bin(M, p_+), K_- ~ Bin(M, p_-) independent (one batch per shift),

    g_hat_{Z,k} = (K_+ - K_-)/(2M)                        (exact identity)
    E[g_hat] = (p_+ - p_-)/2 = s_k B_k / 2 = g_{Z,k}      (unbiased)
    Var(g_hat|theta) = [p_+(1-p_+) + p_-(1-p_-)]/(4M) = [1 - (s_k B_k)^2] / (8M)
                     = [1 - 4 g^2] / (8M)
    SNR_Z^2 = g^2/Var = 2 M (s_k B_k)^2 / [1 - (s_k B_k)^2]
    M_SNR(rho) = rho^2 [1 - (sB)^2] / (2 (sB)^2)

Contrast with the projector benchmark, where the two shifted probabilities are
F_+- = A(1 -+ s)/2, i.e. BOTH near 0 (of order A ~ 4^{-(n-1)}), whereas here both are
near 1/2. Same estimator algebra, different location of the Bernoulli parameters.
"""

from __future__ import annotations

import numpy as np

from qlo.stage6.parity_benchmark import shifted_probabilities


def estimator_from_counts(k_plus, k_minus, M) -> np.ndarray:
    """``g_hat_{Z,k} = (K_+ - K_-)/(2M)``; support ``{j/(2M) : j = -M..M}``."""
    return (np.asarray(k_plus, dtype=np.float64) - np.asarray(k_minus, dtype=np.float64)) / (2.0 * np.asarray(M, dtype=np.float64))


def shot_variance(sB, M) -> np.ndarray:
    """``[1 - (sB)^2] / (8M)``."""
    sB = np.clip(np.asarray(sB, dtype=np.float64), -1.0, 1.0)
    return (1.0 - sB**2) / (8.0 * np.asarray(M, dtype=np.float64))


def shot_variance_from_probabilities(p_plus, p_minus, M) -> np.ndarray:
    """``[p_+(1-p_+) + p_-(1-p_-)]/(4M)`` (the generic parameter-shift form)."""
    pp, pm = np.asarray(p_plus, dtype=np.float64), np.asarray(p_minus, dtype=np.float64)
    return (pp * (1 - pp) + pm * (1 - pm)) / (4.0 * np.asarray(M, dtype=np.float64))


def snr_squared(sB, M) -> np.ndarray:
    """``2 M (sB)^2 / [1 - (sB)^2]``. 0 when sB = 0; inf when |sB| = 1 (deterministic outcomes)."""
    sB = np.clip(np.asarray(sB, dtype=np.float64), -1.0, 1.0)
    M = np.asarray(M, dtype=np.float64)
    num = 2.0 * M * sB**2
    den = 1.0 - sB**2
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(den > 0, num / den, np.where(num > 0, np.inf, 0.0))


def shots_for_snr(sB, rho) -> np.ndarray:
    """``M = rho^2 [1 - (sB)^2] / (2 (sB)^2)``; inf at sB = 0; 0 at |sB| = 1."""
    sB = np.clip(np.asarray(sB, dtype=np.float64), -1.0, 1.0)
    rho = np.asarray(rho, dtype=np.float64)
    num = rho**2 * (1.0 - sB**2)
    den = 2.0 * sB**2
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(den > 0, num / den, np.where(num > 0, np.inf, 0.0))


def log10_shots_for_snr(log_abs_sB, rho) -> np.ndarray:
    """Log-space form for exponentially small |sB|: ``2log10 rho + log10(1-(sB)^2) - log10 2 - 2 log10|sB|``."""
    l = np.asarray(log_abs_sB, dtype=np.float64)
    sB2 = np.exp(2 * l)
    with np.errstate(divide="ignore", invalid="ignore"):
        return 2 * np.log10(rho) + np.log10(np.maximum(1.0 - sB2, 0.0)) - np.log10(2.0) - 2 * l / np.log(10.0)


def sample_counts(p_plus, p_minus, M, rng) -> tuple[np.ndarray, np.ndarray]:
    """Draw ``(K_+, K_-)`` with independent batches."""
    M = int(M)
    return rng.binomial(M, np.asarray(p_plus)), rng.binomial(M, np.asarray(p_minus))


def sample_gradient(s, B, M, rng) -> np.ndarray:
    """One finite-shot parity gradient component estimate."""
    pp, pm = shifted_probabilities(s, B)
    kp, km = sample_counts(pp, pm, M, rng)
    return estimator_from_counts(kp, km, M)
