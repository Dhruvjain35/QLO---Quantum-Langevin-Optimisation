"""Minimal barren-plateau *pipeline* smoke experiment.

This is NOT a publication experiment and NOT evidence for or against a barren
plateau. It exists to exercise: seeded initialization -> analytic gradient ->
flatten -> summary statistics -> tidy DataFrame -> CSV.

Per (n_qubits, init) row it records the cost value and the following statistics
of the flattened gradient vector g ∈ R^p:

    mean_grad         = mean(g)
    mean_abs_grad     = mean(|g|)
    within_gradient_entry_var = var(g)   (population variance over the p entries
                                          of ONE gradient vector at ONE theta —
                                          NOT a barren-plateau diagnostic; see
                                          ``qlo.analysis.bp_variance`` for
                                          Var_theta[dC/dtheta_k] across inits)
    grad_norm         = ||g||_2
    max_abs_grad      = max(|g|)
    first_partial     = value of g[0]     (per-init sample; its variance across
                                           inits, ``bp_var_first_partial_across_inits``,
                                           is computed in ``aggregate``)

Usage (defaults are the tiny smoke configuration):

    python -m qlo.experiments.bp_smoke
    python -m qlo.experiments.bp_smoke --qubits 2 4 --depth 1 --n-init 2 --cost local --seed 3
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from qlo.circuits.hardware_efficient import HardwareEfficientAnsatz
from qlo.costs.observables import COST_NAMES
from qlo.gradients.exact import cost_value, gradient_autograd, gradient_parameter_shift
from qlo.utils.seeding import make_rng

_GRADIENT_FNS = {"autograd": gradient_autograd, "parameter-shift": gradient_parameter_shift}

RESULTS_DIR = Path(__file__).resolve().parents[3] / "results"


@dataclass(frozen=True)
class SmokeConfig:
    qubits: tuple[int, ...] = (2, 4, 6)
    depth: int = 2
    n_init: int = 3
    cost: str = "global"
    seed: int = 0
    entangler: str = "ring"
    gradient: str = "autograd"

    def __post_init__(self) -> None:
        if self.cost not in COST_NAMES:
            raise ValueError(f"cost must be one of {COST_NAMES}")
        if self.gradient not in _GRADIENT_FNS:
            raise ValueError(f"gradient must be one of {tuple(_GRADIENT_FNS)}")
        if self.n_init < 1:
            raise ValueError("n_init must be >= 1")


def summarize_gradient(grad: np.ndarray) -> dict[str, float]:
    """Summary statistics of one flattened gradient vector."""
    g = np.asarray(grad, dtype=np.float64).ravel()
    return {
        "n_params": int(g.size),
        "mean_grad": float(np.mean(g)),
        "mean_abs_grad": float(np.mean(np.abs(g))),
        "within_gradient_entry_var": float(np.var(g)),
        "grad_norm": float(np.linalg.norm(g)),
        "max_abs_grad": float(np.max(np.abs(g))),
        "first_partial": float(g[0]),
    }


def run_bp_smoke(config: SmokeConfig) -> pd.DataFrame:
    """One row per (n_qubits, init). Fully determined by ``config`` (incl. seed)."""
    grad_fn = _GRADIENT_FNS[config.gradient]
    rows: list[dict] = []
    for n in config.qubits:
        ansatz = HardwareEfficientAnsatz(n_qubits=n, depth=config.depth, entangler=config.entangler)
        # One RNG stream per qubit count, derived from (seed, n) so adding/removing a
        # qubit count does not change the draws for the others.
        rng = make_rng([config.seed, n])
        for init in range(config.n_init):
            params = ansatz.init_params(rng)
            value = cost_value(ansatz, config.cost, params)
            grad = grad_fn(ansatz, config.cost, params)
            rows.append(
                {
                    "n_qubits": n,
                    "depth": config.depth,
                    "cost": config.cost,
                    "entangler": config.entangler,
                    "gradient": config.gradient,
                    "seed": config.seed,
                    "init": init,
                    "cost_value": value,
                    **summarize_gradient(grad),
                }
            )
    return pd.DataFrame(rows)


def aggregate(df: pd.DataFrame) -> pd.DataFrame:
    """Per-qubit-count aggregate across inits (mean of per-init stats + Var across inits of g[0], ddof=1)."""
    return (
        df.groupby(["n_qubits", "depth", "cost"], as_index=False)
        .agg(
            n_init=("init", "count"),
            n_params=("n_params", "first"),
            mean_grad=("mean_grad", "mean"),
            mean_abs_grad=("mean_abs_grad", "mean"),
            within_gradient_entry_var_mean=("within_gradient_entry_var", "mean"),
            grad_norm_mean=("grad_norm", "mean"),
            grad_norm_std=("grad_norm", "std"),
            bp_var_first_partial_across_inits=("first_partial", lambda s: float(np.var(s, ddof=1)) if len(s) > 1 else float("nan")),
        )
    )


def _parse_args(argv: list[str] | None = None) -> SmokeConfig:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--qubits", type=int, nargs="+", default=[2, 4, 6])
    p.add_argument("--depth", type=int, default=2)
    p.add_argument("--n-init", type=int, default=3)
    p.add_argument("--cost", choices=COST_NAMES, default="global")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--entangler", choices=("ring", "chain"), default="ring")
    p.add_argument("--gradient", choices=tuple(_GRADIENT_FNS), default="autograd")
    p.add_argument("--out", type=Path, default=None, help="CSV path (default: results/bp_smoke_<cost>_seed<seed>.csv)")
    a = p.parse_args(argv)
    cfg = SmokeConfig(
        qubits=tuple(a.qubits), depth=a.depth, n_init=a.n_init, cost=a.cost,
        seed=a.seed, entangler=a.entangler, gradient=a.gradient,
    )
    return cfg, a.out


def main(argv: list[str] | None = None) -> Path:
    cfg, out = _parse_args(argv)
    df = run_bp_smoke(cfg)
    out = out or RESULTS_DIR / f"bp_smoke_{cfg.cost}_seed{cfg.seed}.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    agg_out = out.with_name(out.stem + "_agg.csv")
    aggregate(df).to_csv(agg_out, index=False)
    with pd.option_context("display.width", 200, "display.max_columns", 30, "display.float_format", "{:.6g}".format):
        print(f"config: {cfg}")
        print(df.to_string(index=False))
        print("\naggregate per qubit count:")
        print(aggregate(df).to_string(index=False))
    print(f"\nwrote {out}\nwrote {agg_out}")
    return out


if __name__ == "__main__":
    main()
