"""Ancilla and destructive SWAP tests as EXPLICIT circuits on two copies (numpy statevector, qubit 0 = MSB).

Ancilla SWAP (2n + 1 qubits): H(anc), CSWAP(anc, A_i, B_i) for i = 0..n-1, H(anc), measure anc.
    P(anc = 0) = (1 + |<phi|psi>|^2)/2, single-shot Z = +1 on anc = 0.
Destructive SWAP (2n qubits, no ancilla): for every pair (A_i, B_i): CNOT(A_i -> B_i), H(A_i), measure both.
    Pair outcome (a_i, b_i) = (1, 1) is the singlet (SWAP eigenvalue -1); Z = prod_i (-1)^(a_i b_i).
    E[Z] = Tr[SWAP (rho_A x rho_B)] = |<phi|psi>|^2 for pure states.
Nothing here calls a fidelity formula: every probability comes from simulating the circuit.
"""

from __future__ import annotations

import numpy as np

H = np.array([[1, 1], [1, -1]], dtype=np.complex128) / np.sqrt(2.0)


def _apply_1q(state: np.ndarray, q: int, u: np.ndarray) -> np.ndarray:
    return np.moveaxis(np.tensordot(u, state, axes=([1], [q])), 0, q)


def _cnot(state: np.ndarray, c: int, t: int) -> np.ndarray:
    s = state.copy()
    idx = [slice(None)] * s.ndim
    idx[c] = 1
    sub = s[tuple(idx)]
    t_axis = t - (1 if t > c else 0)
    s[tuple(idx)] = np.flip(sub, axis=t_axis)
    return s


def _cswap(state: np.ndarray, c: int, a: int, b: int) -> np.ndarray:
    s = state.copy()
    idx = [slice(None)] * s.ndim
    idx[c] = 1
    sub = s[tuple(idx)]
    aa, bb = a - (1 if a > c else 0), b - (1 if b > c else 0)
    s[tuple(idx)] = np.swapaxes(sub, aa, bb)
    return s


def ancilla_swap_p_plus(psi: np.ndarray, phi: np.ndarray) -> float:
    """P(Z = +1) of the ancilla SWAP test, by explicit circuit simulation."""
    n = int(np.log2(psi.size))
    state = np.kron(np.array([1.0, 0.0], dtype=np.complex128), np.kron(psi, phi)).reshape((2,) * (2 * n + 1))
    state = _apply_1q(state, 0, H)
    for i in range(n):
        state = _cswap(state, 0, 1 + i, 1 + n + i)
    state = _apply_1q(state, 0, H)
    return float(np.sum(np.abs(state[0]) ** 2))


def destructive_swap_outcomes(psi: np.ndarray, phi: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(probabilities of all 2^(2n) bit strings, Z value of each string) for the Bell-basis destructive SWAP."""
    n = int(np.log2(psi.size))
    state = np.kron(psi, phi).reshape((2,) * (2 * n))
    for i in range(n):
        state = _cnot(state, i, n + i)
        state = _apply_1q(state, i, H)
    probs = np.abs(state.reshape(-1)) ** 2
    bits = (np.arange(probs.size)[:, None] >> np.arange(2 * n - 1, -1, -1)[None, :]) & 1
    singlets = np.sum(bits[:, :n] & bits[:, n:], axis=1)
    z = np.where(singlets % 2 == 0, 1, -1)
    return probs, z


def destructive_swap_p_plus(psi: np.ndarray, phi: np.ndarray) -> float:
    p, z = destructive_swap_outcomes(psi, phi)
    return float(np.sum(p[z == 1]))


def bell_pair_expectation(psi1: np.ndarray, phi1: np.ndarray) -> np.ndarray:
    """E[(-1)^(a b)] of the destructive SWAP on ONE qubit pair, vectorized over rows: psi1, phi1 shape (..., 2).

    The 2-qubit Bell circuit is simulated explicitly (CNOT then H on the first qubit); used for product states,
    where the pair outcomes are independent and E[Z] is the product over pairs."""
    st = np.einsum("...i,...j->...ij", psi1, phi1)                     # (..., 2, 2)  [A, B]
    st = st.copy()
    st[..., 1, :] = st[..., 1, ::-1]                                     # CNOT(A -> B)
    st = np.einsum("ai,...ij->...aj", H, st)                             # H on A
    p = np.abs(st) ** 2
    return p[..., 0, 0] + p[..., 0, 1] + p[..., 1, 0] - p[..., 1, 1]


def sample_z(rng, p_plus: float, M: int, reps: int) -> np.ndarray:
    """Sample means of M single-shot Z values, reps times."""
    return 2.0 * rng.binomial(M, p_plus, size=reps) / M - 1.0


def sample_destructive_shots(rng, psi, phi, M: int, reps: int) -> np.ndarray:
    """Finite-shot destructive SWAP: draw full 2n-bit measurement records from the circuit distribution, map to Z."""
    p, z = destructive_swap_outcomes(psi, phi)
    draws = rng.choice(p.size, size=(reps, M), p=p / p.sum())
    return z[draws].mean(axis=1)
