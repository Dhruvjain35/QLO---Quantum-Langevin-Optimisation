import numpy as np
import pandas as pd
import pytest

from qlo import HardwareEfficientAnsatz, cost_value, gradient_autograd, make_rng, random_params
from qlo.experiments.bp_smoke import SmokeConfig, run_bp_smoke


def test_same_seed_same_params():
    a = HardwareEfficientAnsatz(3, 2)
    p1 = a.init_params(make_rng(2024))
    p2 = a.init_params(make_rng(2024))
    assert np.array_equal(p1, p2)


def test_different_seed_different_params():
    a = HardwareEfficientAnsatz(3, 2)
    assert not np.array_equal(a.init_params(make_rng(1)), a.init_params(make_rng(2)))


def test_random_params_range_and_dtype():
    p = random_params((2, 3, 2), make_rng(0))
    assert p.dtype == np.float64
    assert np.all(p >= 0.0) and np.all(p < 2 * np.pi)


def test_rng_does_not_touch_global_numpy_state():
    np.random.seed(0)
    before = np.random.get_state()[1].copy()
    HardwareEfficientAnsatz(2, 1).init_params(make_rng(5))
    after = np.random.get_state()[1]
    assert np.array_equal(before, after)


def test_same_seed_same_cost_and_gradient():
    a = HardwareEfficientAnsatz(3, 2)
    outs = []
    for _ in range(2):
        p = a.init_params(make_rng(77))
        outs.append((cost_value(a, "local", p), gradient_autograd(a, "local", p)))
    assert outs[0][0] == outs[1][0]
    assert np.array_equal(outs[0][1], outs[1][1])


# ---- smoke experiment ---------------------------------------------------------

TINY = SmokeConfig(qubits=(2, 3), depth=1, n_init=2, cost="global", seed=0)


def test_smoke_runs_and_has_expected_rows():
    df = run_bp_smoke(TINY)
    assert isinstance(df, pd.DataFrame)
    assert len(df) == len(TINY.qubits) * TINY.n_init
    assert set(df["n_qubits"]) == set(TINY.qubits)
    assert (df["n_params"] == 2 * df["n_qubits"] * df["depth"]).all()


def test_smoke_outputs_are_finite():
    df = run_bp_smoke(TINY)
    numeric = df.select_dtypes(include=[np.number])
    assert np.isfinite(numeric.to_numpy()).all(), "NaN/Inf in smoke outputs"


def test_smoke_is_deterministic_under_seed():
    d1 = run_bp_smoke(TINY)
    d2 = run_bp_smoke(TINY)
    pd.testing.assert_frame_equal(d1, d2)


def test_smoke_changes_with_seed():
    d1 = run_bp_smoke(TINY)
    d2 = run_bp_smoke(SmokeConfig(**{**TINY.__dict__, "seed": 1}))
    assert not np.allclose(d1["cost_value"], d2["cost_value"])


def test_smoke_local_and_parameter_shift_paths():
    cfg = SmokeConfig(qubits=(2,), depth=1, n_init=1, cost="local", seed=0, gradient="parameter-shift")
    df = run_bp_smoke(cfg)
    assert len(df) == 1 and df.loc[0, "cost"] == "local"


def test_smoke_rejects_bad_config():
    with pytest.raises(ValueError):
        SmokeConfig(cost="bogus")
    with pytest.raises(ValueError):
        SmokeConfig(n_init=0)


def test_smoke_cli_writes_csv(tmp_path):
    from qlo.experiments.bp_smoke import main

    out = tmp_path / "smoke.csv"
    main(["--qubits", "2", "--depth", "1", "--n-init", "1", "--out", str(out)])
    assert out.exists()
    assert (tmp_path / "smoke_agg.csv").exists()
    assert len(pd.read_csv(out)) == 1
