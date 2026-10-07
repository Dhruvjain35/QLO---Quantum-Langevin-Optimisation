"""Batched numpy statevector simulator for the two B1 entangling families.

Conventions (identical to PennyLane; checked against it in the tests and in the structural-zero audit):
    wire 0 is the most significant bit of the computational-basis index;
    RX(x) = [[c, -i s], [-i s, c]],  RY(x) = [[c, -s], [s, c]],  RZ(x) = diag(e^{-ix/2}, e^{ix/2}),  c, s = cos, sin(x/2);
    every rotation is R_P(x) = exp(-i x P / 2), so the ordinary two-term shift rule at +-pi/2 is exact.

Families (parameter tensor theta[l, q, j], shape (depth, n, 2); target |0^n>):
    hea_ring      the repo ``HardwareEfficientAnsatz(n, depth, "ring")``: per layer RY(theta[l,q,0]) RZ(theta[l,q,1])
                  on every qubit, then CNOT(0,1), CNOT(1,2), ..., CNOT(n-2,n-1), CNOT(n-1,0) in that order.
    rxry_czbrick  the B1 fallback family (no second ansatz is defined in the repo or the principal work): per layer
                  RX(theta[l,q,0]) RY(theta[l,q,1]) on every qubit, then open-boundary CZ brickwork,
                  even l: CZ(0,1), CZ(2,3), ...;  odd l: CZ(1,2), CZ(3,4), ...

Shifted fidelities from one forward and one backward pass. For the parameter gate R_P(x) at position (l, 0, j):
    U = U_after R_P(x) U_before,   |a> = U_before |0^n>,   |b> = U_after^dagger |0^n>,
    alpha = <b|a>,   beta = <b|P|a>,   <0^n|U|0^n> = cos(x/2) alpha - i sin(x/2) beta,
so F(theta +- pi/2 e_k) follows from (alpha, beta) without re-simulation. The shift rule is built into this
representation, so the audit checks it against INDEPENDENT full re-simulation at theta +- pi/2 e_k, an analytic
derivative (generator insertion) and PennyLane autograd on the repo class.
All positions used by B1 act on qubit 0; gates on other qubits of the same layer commute with them.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np

FAMILIES = ("hea_ring", "rxry_czbrick")
AXES = {"hea_ring": ("Y", "Z"), "rxry_czbrick": ("X", "Y")}      # rotation axis of theta[l,q,0], theta[l,q,1]
FAMILY_LABEL = {"hea_ring": "repo HEA (RY,RZ + CNOT ring)", "rxry_czbrick": "RX,RY + CZ brickwork (B1 fallback)"}
MAX_STATE_ELEMENTS = 1 << 22                                       # rows x 2^n per state array (64 MiB complex128)


def batch_rows(n: int, max_elements: int = MAX_STATE_ELEMENTS) -> int:
    return max(1, max_elements >> n)


def rotation(axis: str, x) -> np.ndarray:
    """Per-row 2x2 rotation matrices exp(-i x P / 2), shape x.shape + (2, 2)."""
    x = np.asarray(x, dtype=np.float64)
    c, s = np.cos(x / 2.0), np.sin(x / 2.0)
    m = np.zeros(x.shape + (2, 2), dtype=np.complex128)
    if axis == "X":
        m[..., 0, 0], m[..., 0, 1], m[..., 1, 0], m[..., 1, 1] = c, -1j * s, -1j * s, c
    elif axis == "Y":
        m[..., 0, 0], m[..., 0, 1], m[..., 1, 0], m[..., 1, 1] = c, -s, s, c
    elif axis == "Z":
        m[..., 0, 0], m[..., 1, 1] = np.exp(-0.5j * x), np.exp(0.5j * x)
    else:
        raise ValueError(axis)
    return m


def dagger(u: np.ndarray) -> np.ndarray:
    return np.conj(np.swapaxes(u, -1, -2))


def apply_1q(psi: np.ndarray, n: int, q: int, u: np.ndarray) -> np.ndarray:
    """In place: psi (B, 2^n) <- u_q psi with one 2x2 matrix per row, u (B, 2, 2)."""
    rows = psi.shape[0]
    v = psi.reshape(rows, 1 << q, 2, 1 << (n - q - 1))
    a0, a1 = v[:, :, 0, :], v[:, :, 1, :]
    u = u.reshape(rows, 2, 2, 1, 1)
    new0 = u[:, 0, 0] * a0 + u[:, 0, 1] * a1
    v[:, :, 1, :] = u[:, 1, 0] * a0 + u[:, 1, 1] * a1
    v[:, :, 0, :] = new0
    return psi


def apply_pauli0(psi: np.ndarray, n: int, axis: str) -> np.ndarray:
    """In place: the generator P on qubit 0 (the most significant bit)."""
    h = 1 << (n - 1)
    a0, a1 = psi[:, :h].copy(), psi[:, h:].copy()
    if axis == "X":
        psi[:, :h], psi[:, h:] = a1, a0
    elif axis == "Y":
        psi[:, :h], psi[:, h:] = -1j * a1, 1j * a0
    elif axis == "Z":
        psi[:, h:] = -a1
    else:
        raise ValueError(axis)
    return psi


def _bit(idx: np.ndarray, n: int, q: int) -> np.ndarray:
    return (idx >> (n - 1 - q)) & 1


def ring_pairs(n: int) -> list[tuple[int, int]]:
    """Same list as ``HardwareEfficientAnsatz.entangling_pairs`` for entangler="ring"."""
    if n == 1:
        return []
    return [(q, q + 1) for q in range(n - 1)] + [(n - 1, 0)]


def brick_pairs(n: int, layer: int) -> list[tuple[int, int]]:
    return [(a, a + 1) for a in range(layer % 2, n - 1, 2)]


@lru_cache(maxsize=None)
def _cnot_ring_gather(n: int, inverse: bool) -> np.ndarray:
    """Index g with (E psi)[j] = psi[g[j]] for the ordered CNOT ring E (or E^dagger if ``inverse``)."""
    idx = np.arange(1 << n)
    pairs = ring_pairs(n)
    for c, t in (pairs if inverse else pairs[::-1]):
        idx = np.where(_bit(idx, n, c) == 1, idx ^ (1 << (n - 1 - t)), idx)
    idx.setflags(write=False)
    return idx


@lru_cache(maxsize=None)
def _cz_phase(n: int, parity: int) -> np.ndarray:
    idx = np.arange(1 << n)
    ph = np.ones(1 << n)
    for a, b in brick_pairs(n, parity):
        ph[(_bit(idx, n, a) & _bit(idx, n, b)) == 1] *= -1.0
    ph.setflags(write=False)
    return ph


def entangle(psi: np.ndarray, family: str, n: int, layer: int, inverse: bool = False) -> np.ndarray:
    if family == "hea_ring":
        return psi[:, _cnot_ring_gather(n, inverse)] if n > 1 else psi
    if family == "rxry_czbrick":
        psi *= _cz_phase(n, layer % 2)          # diagonal and self-inverse
        return psi
    raise ValueError(family)


def _layer_unitaries(theta: np.ndarray, family: str, layer: int, q: int) -> tuple[np.ndarray, np.ndarray]:
    ax0, ax1 = AXES[family]
    return rotation(ax0, theta[:, layer, q, 0]), rotation(ax1, theta[:, layer, q, 1])


def forward(theta: np.ndarray, family: str, n: int, depth: int, save=(), insert=None):
    """Final state for theta (B, depth, n, 2), plus {position: state just before that gate} for ``save``.

    Positions are (layer, 0, j). ``insert`` = a position after whose gate the generator P is applied (used only for
    the analytic derivative). Within a layer the qubit-0 gates are applied after the other qubits' rotations, which
    commute with them."""
    theta = np.asarray(theta, dtype=np.float64)
    if theta.shape[1:] != (depth, n, 2):
        raise ValueError(f"theta shape {theta.shape[1:]} != {(depth, n, 2)}")
    save = set(save)
    special = {p[0] for p in save} | ({insert[0]} if insert is not None else set())
    psi = np.zeros((theta.shape[0], 1 << n), dtype=np.complex128)
    psi[:, 0] = 1.0
    saved = {}
    axes = AXES[family]
    for layer in range(depth):
        for q in range(1, n):
            r0, r1 = _layer_unitaries(theta, family, layer, q)
            apply_1q(psi, n, q, r1 @ r0)
        if layer in special:
            for j in (0, 1):
                if (layer, 0, j) in save:
                    saved[(layer, 0, j)] = psi.copy()
                apply_1q(psi, n, 0, rotation(axes[j], theta[:, layer, 0, j]))
                if insert == (layer, 0, j):
                    apply_pauli0(psi, n, axes[j])
        else:
            r0, r1 = _layer_unitaries(theta, family, layer, 0)
            apply_1q(psi, n, 0, r1 @ r0)
        psi = entangle(psi, family, n, layer)
    return psi, saved


def _bracket(b: np.ndarray, a: np.ndarray, n: int, axis: str) -> tuple[np.ndarray, np.ndarray]:
    """alpha = <b|a>, beta = <b|P_0|a> row by row."""
    h = 1 << (n - 1)
    b0, b1 = np.conj(b[:, :h]), np.conj(b[:, h:])
    a0, a1 = a[:, :h], a[:, h:]
    s00, s11 = np.einsum("ij,ij->i", b0, a0), np.einsum("ij,ij->i", b1, a1)
    s01, s10 = np.einsum("ij,ij->i", b0, a1), np.einsum("ij,ij->i", b1, a0)
    alpha = s00 + s11
    beta = {"X": s01 + s10, "Y": -1j * s01 + 1j * s10, "Z": s00 - s11}[axis]
    return alpha, beta


def brackets(theta: np.ndarray, family: str, n: int, depth: int, positions) -> tuple[np.ndarray, dict]:
    """(<0^n|U(theta)|0^n>, {position: (alpha, beta)}) from one forward and one backward pass."""
    positions = set(positions)
    psi, saved = forward(theta, family, n, depth, save=positions)
    amp = psi[:, 0].copy()
    del psi
    axes = AXES[family]
    phi = np.zeros((theta.shape[0], 1 << n), dtype=np.complex128)
    phi[:, 0] = 1.0
    out = {}
    lmin = min(p[0] for p in positions)
    for layer in range(depth - 1, lmin - 1, -1):
        phi = entangle(phi, family, n, layer, inverse=True)
        if any(p[0] == layer for p in positions):
            for j in (1, 0):
                if (layer, 0, j) in positions:
                    out[(layer, 0, j)] = _bracket(phi, saved.pop((layer, 0, j)), n, axes[j])
                apply_1q(phi, n, 0, dagger(rotation(axes[j], theta[:, layer, 0, j])))
        else:
            r0, r1 = _layer_unitaries(theta, family, layer, 0)
            apply_1q(phi, n, 0, dagger(r1 @ r0))
        if layer > lmin:
            for q in range(1, n):
                r0, r1 = _layer_unitaries(theta, family, layer, q)
                apply_1q(phi, n, q, dagger(r1 @ r0))
    return amp, out


def fidelity_from_brackets(alpha: np.ndarray, beta: np.ndarray, x) -> np.ndarray:
    """F at parameter value x: |cos(x/2) alpha - i sin(x/2) beta|^2 (non-negative by construction)."""
    x = np.asarray(x, dtype=np.float64)
    amp = np.cos(x / 2.0) * alpha - 1j * np.sin(x / 2.0) * beta
    return amp.real**2 + amp.imag**2


def shifted_fidelities(theta: np.ndarray, family: str, n: int, depth: int, positions) -> tuple[np.ndarray, dict]:
    """F(theta) from the forward pass and, per position, (x, F_plus, F_minus, F(theta) rebuilt from the brackets)."""
    amp, br = brackets(theta, family, n, depth, positions)
    out = {}
    for p, (alpha, beta) in br.items():
        x = theta[:, p[0], p[1], p[2]]
        out[p] = {"x": x, "F_plus": fidelity_from_brackets(alpha, beta, x + np.pi / 2),
                  "F_minus": fidelity_from_brackets(alpha, beta, x - np.pi / 2),
                  "F_bracket": fidelity_from_brackets(alpha, beta, x)}
    return amp.real**2 + amp.imag**2, out


def fidelity(theta: np.ndarray, family: str, n: int, depth: int) -> np.ndarray:
    """F(theta) = |<0^n|U(theta)|0^n>|^2 by plain forward simulation (no brackets)."""
    amp = forward(theta, family, n, depth)[0][:, 0]
    return amp.real**2 + amp.imag**2


def analytic_cost_gradient(theta: np.ndarray, family: str, n: int, depth: int, position) -> np.ndarray:
    """dC/dx = -dF/dx at ``position`` WITHOUT the shift rule: dR/dx = (-i/2) P R(x), so
    d<0|U|0>/dx = (-i/2) <0|U_after P R(x) U_before|0> (generator inserted after the gate)."""
    amp = forward(theta, family, n, depth)[0][:, 0]
    amp_p = forward(theta, family, n, depth, insert=tuple(position))[0][:, 0]
    dF = 2.0 * np.real(np.conj(amp) * (-0.5j) * amp_p)
    return -dF
