"""Access-model and resource accounting (qualitative where the cost depends on the hardware)."""

from __future__ import annotations

import pandas as pd


def resource_table(n: int | None = None) -> pd.DataFrame:
    nn = "n" if n is None else str(n)
    rows = [
        {"estimator": "loschmidt", "measures": "projector |0^n><0^n| after V^dagger U", "state_copies_per_shot": 1,
         "qubits_per_shot": nn, "needs": "target inverse V^dagger", "entangling_ops_added": "none beyond V^dagger U",
         "executions_per_gradient": "2M (M per shift)", "unbiased_for_F": True, "var_F_hat": "F(1-F)/M",
         "deep_F_var": "~F/M (vanishes with the signal)", "relative_depth": "U then V^dagger"},
        {"estimator": "swap_ancilla", "measures": "ancilla Z after H, CSWAP x n, H", "state_copies_per_shot": 2,
         "qubits_per_shot": f"2{nn}+1" if n is None else str(2 * n + 1), "needs": "two copies + ancilla", "entangling_ops_added": "n controlled-SWAPs (3-qubit gates)",
         "executions_per_gradient": "2M (M per shift)", "unbiased_for_F": True, "var_F_hat": "(1-F^2)/M",
         "deep_F_var": "~1/M (does not vanish)", "relative_depth": "deepest: n Fredkin gates"},
        {"estimator": "swap_destructive", "measures": "Bell-basis measurement of every copy pair; Z = prod (-1)^(a_i b_i)",
         "state_copies_per_shot": 2, "qubits_per_shot": f"2{nn}" if n is None else str(2 * n), "needs": "two copies, no ancilla",
         "entangling_ops_added": "n CNOTs (one per pair) + n H, depth 2", "executions_per_gradient": "2M (M per shift)",
         "unbiased_for_F": True, "var_F_hat": "(1-F^2)/M", "deep_F_var": "~1/M (identical to ancilla SWAP)",
         "relative_depth": "shallow: constant depth 2 after state preparation"},
        {"estimator": "hadamard_u", "measures": "ancilla Z after H, controlled-(V^dagger U), [S^dagger], H; real and imaginary quadratures",
         "state_copies_per_shot": 1, "qubits_per_shot": f"{nn}+1" if n is None else str(n + 1), "needs": "coherent controlled-(V^dagger U)",
         "entangling_ops_added": "every gate of V^dagger U controlled by the ancilla", "executions_per_gradient": "4M (M per quadrature per shift)",
         "unbiased_for_F": True, "var_F_hat": "[2(1-x^4)+4(M-2)x^2(1-x^2)]/[M(M-1)] + same in y",
         "deep_F_var": "~4/M^2 + 4F/M", "relative_depth": "deepest per circuit: controlled version of every gate"},
    ]
    return pd.DataFrame(rows)
