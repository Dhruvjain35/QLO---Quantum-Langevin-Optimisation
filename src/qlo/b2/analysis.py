"""B2 validation, fixed-budget failure modes and the secondary entangling spot check."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import binom

from qlo.b2 import destructive_swap as DS
from qlo.b2 import hadamard as HD
from qlo.b2 import hadamard_theory as HT
from qlo.b2 import required_shots as RS
from qlo.stage6.directional import direction_probabilities
from qlo.stage7 import analysis as AN
from qlo.stage7.estimators import gradient_variance, pair


def _rand_state(rng, n):
    v = rng.normal(size=2**n) + 1j * rng.normal(size=2**n)
    return v / np.linalg.norm(v)


# ---- validation -------------------------------------------------------------------------------------------------
def destructive_swap_validation(cfg, seed: int = 0) -> pd.DataFrame:
    rng = np.random.default_rng(np.random.SeedSequence([0xB2, 1, seed]))
    rows = []
    for n in (1, 2, 3, 4):
        for trial in range(3):
            psi, phi = _rand_state(rng, n), _rand_state(rng, n)
            F = abs(np.vdot(phi, psi)) ** 2
            pa, pd_ = DS.ancilla_swap_p_plus(psi, phi), DS.destructive_swap_p_plus(psi, phi)
            p, z = DS.destructive_swap_outcomes(psi, phi)
            ez = float(np.sum(p * z))
            base = {"n": n, "trial": trial, "F_exact": F, "p_plus_ancilla_circuit": pa, "p_plus_destructive_circuit": pd_,
                    "E_Z_destructive": ez, "Var_Z_destructive": float(np.sum(p * z**2) - ez**2), "Var_model_1_minus_F2": 1 - F**2,
                    "abs_err_ancilla": abs(pa - (1 + F) / 2), "abs_err_destructive": abs(pd_ - (1 + F) / 2),
                    "abs_err_ancilla_vs_destructive": abs(pa - pd_)}
            for M in cfg.validation_shots:
                reps = 4000 if n <= 3 else 2000
                means = DS.sample_destructive_shots(rng, psi, phi, M, reps)          # full 2n-bit records
                anc = DS.sample_z(rng, pa, M, reps)                                  # ancilla circuit outcome
                se = np.sqrt((1 - F**2) / M / reps)
                rows.append({**base, "M": M, "reps": reps, "dswap_sample_mean": float(means.mean()),
                             "dswap_mean_z": float((means.mean() - F) / se),
                             "dswap_sample_var_x_M": float(means.var(ddof=1) * M), "model_var_x_M": 1 - F**2,
                             "ancilla_sample_mean": float(anc.mean()), "ancilla_mean_z": float((anc.mean() - F) / se),
                             "two_sample_z": float((means.mean() - anc.mean()) / (se * np.sqrt(2)))})
    return pd.DataFrame(rows)


def hadamard_circuit_validation(seed: int = 0) -> pd.DataFrame:
    rng = np.random.default_rng(np.random.SeedSequence([0xB2, 2, seed]))
    rows = []
    for n in (1, 2, 3, 4):
        for trial in range(4):
            params = rng.uniform(-np.pi, np.pi, (2, n, 2))
            U = HD.hea_unitary(params, n, 2)
            a = U[0, 0]
            x_np, y_np = HD.hadamard_expectation(U, False), HD.hadamard_expectation(U, True)
            x_pl, y_pl = HD.pennylane_hadamard(params, n, 2, False), HD.pennylane_hadamard(params, n, 2, True)
            F = abs(a) ** 2
            rows.append({"n": n, "trial": trial, "re_a": a.real, "im_a": a.imag, "F_exact": F, "x_numpy": x_np, "y_numpy": y_np,
                         "x_pennylane": x_pl, "y_pennylane": y_pl,
                         "max_abs_err": max(abs(x_np - a.real), abs(y_np - a.imag), abs(x_pl - a.real), abs(y_pl - a.imag)),
                         "abs_err_x2_plus_y2_vs_F": abs(x_pl**2 + y_pl**2 - F)})
    return pd.DataFrame(rows)


def hadamard_variance_validation(cfg, seed: int = 0) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(np.random.SeedSequence([0xB2, 3, seed]))
    rows, pmf_rows = [], []
    cases = [(0.0, 0.0), (0.31, -0.42), (0.9, 0.1), (-0.05, 0.02), (0.6, 0.6)]
    for x, y in cases:
        F = x * x + y * y
        for M in cfg.hadamard_validation_m:
            r = {"x": x, "y": y, "F": F, "M": M, "var_formula": float(HT.var_f_ht(x, y, M))}
            if M <= cfg.hadamard_exact_validation_max_m:
                f, p = HT.f_ht_pmf(x, y, M)
                mean = float(np.sum(p * f))
                r.update({"exact_pmf_sum": float(p.sum()), "exact_mean": mean, "exact_bias": mean - F,
                          "exact_var": float(np.sum(p * (f - mean) ** 2)), "exact_var_err": float(np.sum(p * (f - mean) ** 2) - r["var_formula"]),
                          "exact_support_size": int(f.size), "exact_min": float(f.min()), "exact_max": float(f.max()),
                          "exact_P_F_lt_0": float(p[f < 0].sum()), "exact_P_F_gt_1": float(p[f > 1].sum())})
                if M <= 8:
                    pmf_rows += [{"x": x, "y": y, "M": M, "F_hat": float(v), "prob": float(w)} for v, w in zip(f, p)]
            draws = HT.sample_f_ht(rng, x, y, M, cfg.validation_mc_reps)
            se = np.sqrt(r["var_formula"] / cfg.validation_mc_reps)
            r.update({"mc_reps": cfg.validation_mc_reps, "mc_mean": float(draws.mean()),
                      "mc_mean_z": float((draws.mean() - F) / se) if se > 0 else 0.0, "mc_var": float(draws.var(ddof=1)),
                      "mc_var_rel_err": float(draws.var(ddof=1) / r["var_formula"] - 1), "mc_P_F_lt_0": float(np.mean(draws < 0)),
                      "mc_P_F_gt_1": float(np.mean(draws > 1)), "naive_plugin_bias_formula": (2 - F) / M})
            rows.append(r)
    # gradient: exact small-M difference distribution (mean and variance of g_hat)
    for xp, xm in ((0.2, 0.35), (0.05, -0.01)):
        for M in (2, 3, 4, 8, 16):
            fp, pp = HT.f_ht_pmf(xp, 0.1, M)
            fm, pm = HT.f_ht_pmf(xm, -0.2, M)
            g = (fm[None, :] - fp[:, None]) / 2.0
            w = pp[:, None] * pm[None, :]
            gm = float(np.sum(w * g))
            gtrue = ((xm**2 + 0.04) - (xp**2 + 0.01)) / 2
            rows.append({"x": xp, "y": 0.1, "F": xp**2 + 0.01, "M": M, "gradient_case": f"x-={xm}, y-=-0.2",
                         "exact_grad_mean": gm, "exact_grad_bias": gm - gtrue, "exact_grad_var": float(np.sum(w * (g - gm) ** 2)),
                         "grad_var_formula": float(HT.var_g_ht(xp, 0.1, xm, -0.2, M))})
    return pd.DataFrame(rows), pd.DataFrame(pmf_rows)


# ---- fixed budgets ----------------------------------------------------------------------------------------------
def fixed_budget(cfg) -> pd.DataFrame:
    rows = []
    rng = np.random.default_rng(np.random.SeedSequence([0xB2, 4, cfg.mc_seed]))
    for n in cfg.fixed_n:
        theta = AN.sample_theta(n, cfg.fixed_samples, cfg.stage7_seed, cfg.fixed_stream)
        q = RS.shifted_quantities(theta)
        A, s, g = np.exp(q["logA"]), q["s"], q["delta"] / 2
        e_p, e_m = np.exp(q["log_e_plus"]), np.exp(q["log_e_minus"])
        for T in cfg.budgets_total_executions:
            M2, M4 = T // 2, T // 4
            for est in ("loschmidt", "swap_ancilla", "swap_destructive"):
                if est == "swap_destructive":
                    p_pos, p_neg = (1 + e_m) / 2, (1 + e_p) / 2
                    var = (2 - e_p**2 - e_m**2) / (4 * M2)
                else:
                    sc = "loschmidt" if est == "loschmidt" else "swap"
                    p_pos, p_neg, _ = pair(sc, A, s)
                    var = gradient_variance(sc, A, s, M2)
                pr = direction_probabilities(p_pos, p_neg, float(M2))
                nz = 1 - pr["p_zero"]
                with np.errstate(divide="ignore", invalid="ignore"):
                    cond = np.where(nz > 0, pr["p_correct"] / nz, np.nan)
                rows.append(_frow(n, T, est, M2, pr["p_zero"], pr["p_correct"], pr["p_wrong"], cond, np.abs(g) / np.sqrt(var), var, "exact"))
            # Hadamard on the first fixed_samples_hadamard theta
            h = cfg.fixed_samples_hadamard
            xp, yp, xm, ym = (q[k][:h] for k in ("x_plus", "y_plus", "x_minus", "y_minus"))
            sg = np.sign(g[:h])
            var = HT.var_g_ht(xp, yp, xm, ym, M4)
            if M4 <= cfg.hadamard_exact_max_m:
                res = [HT.gradient_sign_exact(xp[i], yp[i], xm[i], ym[i], M4) for i in range(h)]
                pz = np.array([r["p_zero"] for r in res])
                pc = np.array([r["p_pos"] if sg[i] > 0 else r["p_neg"] for i, r in enumerate(res)])
                pw = 1 - pz - pc
                method = "exact"
            else:
                pz, pc, pw = np.empty(h), np.empty(h), np.empty(h)
                for i in range(h):
                    d = (HT.sample_f_ht(rng, xm[i], ym[i], M4, cfg.hadamard_mc_reps) - HT.sample_f_ht(rng, xp[i], yp[i], M4, cfg.hadamard_mc_reps))
                    pz[i], pc[i], pw[i] = np.mean(d == 0), np.mean(np.sign(d) == sg[i]), np.mean(np.sign(d) == -sg[i])
                method = f"monte_carlo_{cfg.hadamard_mc_reps}"
            with np.errstate(divide="ignore", invalid="ignore"):
                cond = np.where(1 - pz > 0, pc / (1 - pz), np.nan)
            rows.append(_frow(n, T, "hadamard_u", M4, pz, pc, pw, cond, np.abs(g[:h]) / np.sqrt(var), var, method))
    return pd.DataFrame(rows)


def _frow(n, T, est, M, pz, pc, pw, cond, snr, var, method) -> dict:
    return {"n": n, "total_executions": T, "estimator": est, "M_internal": M, "method": method, "theta_count": int(np.size(pz)),
            "median_P_zero": float(np.median(pz)), "median_P_correct": float(np.median(pc)), "median_P_wrong": float(np.median(pw)),
            "median_P_correct_given_nonzero": float(np.nanmedian(cond)), "mean_P_correct": float(np.mean(pc)),
            "median_SNR": float(np.median(snr)), "median_log10_MSE": float(np.median(np.log10(var)))}


def hadamard_distribution_examples(cfg) -> pd.DataFrame:
    """Exact pmf of the Hadamard gradient estimate at the median-S theta of n = 10, for small budgets."""
    theta = AN.sample_theta(10, cfg.fixed_samples, cfg.stage7_seed, cfg.fixed_stream)
    q = RS.shifted_quantities(theta)
    i = int(np.argsort(q["F_plus"] + q["F_minus"])[cfg.fixed_samples // 2])
    rows = []
    for T in (32, 128, 512):                             # exact joint pmf; T = 2048 (M = 512) needs a ~66k x 66k grid
        M = T // 4
        fp, pp = HT.f_ht_pmf(q["x_plus"][i], q["y_plus"][i], M)
        fm, pm = HT.f_ht_pmf(q["x_minus"][i], q["y_minus"][i], M)
        g = ((fm[None, :] - fp[:, None]) / 2).ravel()
        w = (pp[:, None] * pm[None, :]).ravel()
        vals, inv = np.unique(np.round(g, 15), return_inverse=True)
        pw = np.bincount(inv, weights=w)
        rows += [{"total_executions": T, "M_internal": M, "g_true": float(q["delta"][i] / 2), "g_hat": float(v), "prob": float(p)}
                 for v, p in zip(vals, pw) if p > 1e-12]
    return pd.DataFrame(rows)


# ---- secondary entangling spot check -----------------------------------------------------------------------------
def _hea_amplitude(theta: np.ndarray, n: int) -> np.ndarray:
    """<0^n|U(theta)|0^n> for the repo-layout HEA, depth n, batch of theta (B, n, n, 2); numpy statevector."""
    B = theta.shape[0]
    psi = np.zeros((B,) + (2,) * n, dtype=np.complex128)
    psi[(slice(None),) + (0,) * n] = 1
    pairs = [(q, q + 1) for q in range(n - 1)] + [(n - 1, 0)]
    for l in range(n):
        for q in range(n):
            ry, rz = theta[:, l, q, 0], theta[:, l, q, 1]
            c, s = np.cos(ry / 2), np.sin(ry / 2)
            m = np.stack([np.stack([c, -s], -1), np.stack([s, c], -1)], -2).astype(np.complex128)
            m = m * np.stack([np.exp(-0.5j * rz), np.exp(0.5j * rz)], -1)[..., :, None]
            psi = np.moveaxis(np.einsum("bij,bj...->bi...", m, np.moveaxis(psi, 1 + q, 1)), 1, 1 + q)
        for c_, t in pairs:
            idx = [slice(None)] * (n + 1)
            idx[1 + c_] = 1
            sub = psi[tuple(idx)]
            psi[tuple(idx)] = np.flip(sub, axis=1 + t - (1 if t > c_ else 0))
    return psi[(slice(None),) + (0,) * n]


def entangling_spotcheck(cfg) -> pd.DataFrame:
    rows = []
    for n in cfg.spot_n:
        rng = np.random.default_rng(np.random.SeedSequence([0xB2, 5, cfg.spot_seed, n]))
        theta = rng.uniform(-np.pi, np.pi, (cfg.spot_samples, n, n, 2))
        amps = {}
        for tag, sh in (("plus", np.pi / 2), ("minus", -np.pi / 2)):
            t = theta.copy()
            t[:, 0, 0, 0] += sh
            amps[tag] = np.concatenate([_hea_amplitude(t[i:i + 250], n) for i in range(0, len(t), 250)])
        Fp, Fm = np.abs(amps["plus"]) ** 2, np.abs(amps["minus"]) ** 2
        delta = Fm - Fp
        with np.errstate(divide="ignore"):
            l_le = np.log10(Fp * (1 - Fp) + Fm * (1 - Fm)) - 2 * np.log10(np.abs(delta))
            l_sw = np.log10(2 - Fp**2 - Fm**2) - 2 * np.log10(np.abs(delta))
        m_ht = HT.required_m_ht(amps["plus"].real, amps["plus"].imag, amps["minus"].real, amps["minus"].imag, delta, 1.0)
        tot = {"loschmidt": np.log10(2) + RS.ceil_log10(l_le), "swap_ancilla": np.log10(2) + RS.ceil_log10(l_sw),
               "swap_destructive": np.log10(2) + RS.ceil_log10(l_sw), "hadamard_u": np.log10(4 * m_ht)}
        for est, v in tot.items():
            rows.append({"n": n, "estimator": est, "median_log10_total_executions": float(np.median(v)),
                         "median_abs_imag_over_abs_amp": float(np.median(np.abs(amps["plus"].imag) / np.abs(amps["plus"]))),
                         "note": "destructive SWAP uses the ancilla +-1 model here (validated identical); SECONDARY ENTANGLING SPOT CHECK"})
    return pd.DataFrame(rows)
