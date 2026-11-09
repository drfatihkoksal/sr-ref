"""Power simulation for H1 (analysis plan 6.11), before the freeze. Label-free.

Scenario grid: baseline AUC of M1 (sex + 6 features) in {0.86, 0.90, 0.945}; effect on the d′
scale of the equal-variance binormal model, d′ = √2 Φ^-1(AUC), full (Nakasone, Δd′ = 0.368), half
(0.184) or none; correlation ρ in {0.6, 0.8, 0.9} between the M1 and M3 features.

Data per simulated set: n = 310 (the valid Ningbo records, docs/decisions.md item 23) with the
database class balance (77/334 LVOT: 239 RVOT, 71 LVOT); sex drawn per class with the published
female rates (RVOT 192/257, LVOT 38/77); M1 features x ~ N(δ u y, I_6), u = 1/√6; M3 features
z = ρ x + √(1 - ρ²) e, e ~ N(γ u y, I_6). δ is set so that the population AUC of the Bayes score of
(sex, x) equals the baseline AUC, γ so that that of (sex, z) equals the target AUC (Monte Carlo
calibration, 2 x 10^6 draws, bisection).

Fitting: the plan's nested cross-validation (6.5): 10 x 10 repeated stratified folds (SEED), inner
stratified 5-fold choice of C among 13 log-spaced values in [1e-4, 1e2] by log-loss, refit, no class
weights; L2 logistic regression solved by Newton's method with the scikit-learn objective
(½||w||² + C Σ log-loss, intercept not penalised); features standardised on the training folds.

Decision rule (docs/decisions.md item 27): Nadeau-Bengio corrected resampled t-test on the 100
per-fold mean loss differences, t = mean / sqrt((1/100 + 1/9) s²), df 99, one-sided p < 0.025.
Power = share of simulated sets in which H1 is supported; the 'none' rows give the type I error.
Also reported: the earlier rule (upper bound of the 95% CI of ΔLL < 0, 200 patient bootstrap
resamples of the per-patient loss difference averaged over repeats), whose type I error exceeds
0.025 because the bootstrap ignores refitting.

    python -m src.power [--n-sets 500]   -> docs/wp2_power.md; ARTIFACTS/wp2/power.parquet
"""
import argparse
import itertools
import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from scipy.stats import norm, t as tdist
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold

from src.config import ARTIFACTS, DOCS, SEED

N_RVOT, N_LVOT = 239, 71
P_FEMALE = {0: 192 / 257, 1: 38 / 77}
BASELINES = [0.86, 0.90, 0.945]
EFFECTS = {"full": 0.368, "half": 0.184, "none": 0.0}
RHOS = [0.6, 0.8, 0.9]
C_GRID = np.logspace(-4, 2, 13)
N_BOOT_SIM = 200
U = np.ones(6) / np.sqrt(6)


def dprime(auc):
    return np.sqrt(2) * norm.ppf(auc)


def auc_from_dprime(d):
    return norm.cdf(d / np.sqrt(2))


def population_auc(shift, rng, n=2_000_000):
    """AUC of the Bayes score shift·t + log LR(sex), t ~ N(shift·y, 1)."""
    y = np.r_[np.zeros(n), np.ones(n)]
    t = rng.normal(shift * y, 1.0)
    fem = rng.random(2 * n) < np.where(y == 1, P_FEMALE[1], P_FEMALE[0])
    llr_sex = np.where(fem, np.log(P_FEMALE[1] / P_FEMALE[0]), np.log((1 - P_FEMALE[1]) / (1 - P_FEMALE[0])))
    s = shift * t + llr_sex
    pos, neg = s[y == 1], s[y == 0]
    rank = pd.Series(np.r_[neg, pos]).rank().to_numpy()
    return (rank[n:].sum() - n * (n + 1) / 2) / (n * n)


_CAL = {}


def calibrate(target_auc, rng):
    key = round(float(target_auc), 10)
    if key in _CAL:
        return _CAL[key]
    _CAL[key] = _calibrate(target_auc, rng)
    return _CAL[key]


def _calibrate(target_auc, rng):
    lo, hi = 0.0, 6.0
    for _ in range(40):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if population_auc(mid, rng, 400_000) < target_auc else (lo, mid)
    return (lo + hi) / 2


def fit_lr(X, y, C, w0=None, iters=50):
    """L2 logistic regression, scikit-learn objective; X includes no intercept column."""
    n, p = X.shape
    A = np.c_[np.ones(n), X]
    w = np.zeros(p + 1) if w0 is None else w0.copy()
    reg = np.r_[0.0, np.full(p, 1.0 / C)]
    for _ in range(iters):
        z = A @ w
        pr = 1 / (1 + np.exp(-z))
        g = A.T @ (pr - y) + reg * w
        H = (A * (pr * (1 - pr))[:, None]).T @ A + np.diag(reg)
        step = np.linalg.solve(H, g)
        w -= step
        if np.abs(step).max() < 1e-8:
            break
    return w


def predict(w, X):
    return 1 / (1 + np.exp(-(w[0] + X @ w[1:])))


def logloss(y, p):
    p = np.clip(p, 1e-12, 1 - 1e-12)
    return -(y * np.log(p) + (1 - y) * np.log(1 - p))


def oof_loss(X, y, outer):
    """Per-patient OOF log-loss averaged over repeats, and the per-fold losses (nested CV of plan
    6.5)."""
    n = len(y)
    loss = np.zeros(n)
    folds = []
    for k, (tr, te) in enumerate(outer):
        mu, sd = X[tr].mean(0), X[tr].std(0) + 1e-12
        Xtr, Xte = (X[tr] - mu) / sd, (X[te] - mu) / sd
        inner = StratifiedKFold(5, shuffle=True, random_state=SEED + k).split(Xtr, y[tr])
        cv = np.zeros(len(C_GRID))
        for itr, ite in inner:
            w = None
            for c_i, C in enumerate(C_GRID):
                w = fit_lr(Xtr[itr], y[tr][itr], C, w)
                cv[c_i] += logloss(y[tr][ite], predict(w, Xtr[ite])).sum()
        w = fit_lr(Xtr, y[tr], C_GRID[int(np.argmin(cv))])
        lte = logloss(y[te], predict(w, Xte))
        loss[te] += lte
        folds.append(lte)
    return loss / (len(outer) / 10), folds


def nadeau_bengio(f1, f3, ratio=1 / 9):
    """One-sided p for mean(ΔLL) < 0 from per-fold losses of two models on the same folds."""
    d = np.array([(b - a).mean() for a, b in zip(f1, f3)])
    t = d.mean() / np.sqrt((1 / len(d) + ratio) * d.var(ddof=1))
    return float(tdist.cdf(t, len(d) - 1))


def one_set(args):
    delta, gamma, rho, seed = args
    rng = np.random.default_rng(seed)
    y = np.r_[np.zeros(N_RVOT), np.ones(N_LVOT)]
    female = rng.random(len(y)) < np.where(y == 1, P_FEMALE[1], P_FEMALE[0])
    x = rng.normal(0, 1, (len(y), 6)) + delta * np.outer(y, U)
    e = rng.normal(0, 1, (len(y), 6)) + gamma * np.outer(y, U)
    z = rho * x + np.sqrt(1 - rho ** 2) * e
    outer = list(RepeatedStratifiedKFold(n_splits=10, n_repeats=10, random_state=SEED).split(x, y))
    l1, f1 = oof_loss(np.c_[female, x], y, outer)
    l3, f3 = oof_loss(np.c_[female, z], y, outer)
    d = l3 - l1
    boots = [d[rng.integers(0, len(d), len(d))].mean() for _ in range(N_BOOT_SIM)]
    return {"dll": d.mean(), "upper": float(np.percentile(boots, 97.5)), "p_nb": nadeau_bengio(f1, f3)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-sets", type=int, default=500)
    n_sets = ap.parse_args().n_sets
    rng = np.random.default_rng(SEED)
    cells, jobs = [], []
    for b_i, (base, (eff, dd), rho) in enumerate(itertools.product(BASELINES, EFFECTS.items(), RHOS)):
        target = auc_from_dprime(dprime(base) + dd)
        delta, dz = calibrate(base, rng), calibrate(target, rng)
        gamma = (dz - rho * delta) / np.sqrt(1 - rho ** 2)
        cells.append({"baseline": base, "effect": eff, "rho": rho, "target_auc": target, "delta": delta, "gamma": gamma})
        jobs += [(delta, gamma, rho, SEED * 1000 + b_i * n_sets + s) for s in range(n_sets)]
    with ProcessPoolExecutor(max_workers=max(os.cpu_count() - 4, 1)) as ex:
        res = list(ex.map(one_set, jobs, chunksize=4))
    rows = []
    for c_i, c in enumerate(cells):
        r = pd.DataFrame(res[c_i * n_sets:(c_i + 1) * n_sets])
        rows.append({**c, "power": float((r.p_nb < 0.025).mean()), "power_boot": float((r.upper < 0).mean()),
                     "median_dll": float(r.dll.median())})
    tab = pd.DataFrame(rows)
    tab.to_parquet(ARTIFACTS / "wp2" / "power.parquet")
    lines = ["# WP2: power simulation for H1 (plan 6.11)", "",
             f"Generated by `python -m src.power --n-sets {n_sets}`; definitions in the module docstring. n = 310",
             "(239 RVOT, 71 LVOT), the plan's nested CV. Decision rule: Nadeau-Bengio, one-sided p < 0.025. Power",
             "= share of sets with H1 supported; the 'none' rows give the type I error (nominal 0.025). In brackets:",
             "the earlier bootstrap rule (upper 95% bound of ΔLL < 0).", "",
             "| Baseline AUC (M1) | Effect | Target AUC (M3) | ρ = 0.6 | ρ = 0.8 | ρ = 0.9 | median ΔLL (ρ = 0.8) |",
             "|---|---|---|---|---|---|---|"]
    for (base, eff), g in tab.groupby(["baseline", "effect"], sort=False):
        pw = {r.rho: f"{r.power:.2f} [{r.power_boot:.2f}]" for r in g.itertuples()}
        mid = g[g.rho == 0.8].median_dll.iloc[0]
        lines.append(f"| {base} | {eff} | {g.target_auc.iloc[0]:.3f} | {pw[0.6]} | {pw[0.8]} | {pw[0.9]} | {mid:+.4f} |")
    plan = tab[(tab.baseline == 0.945) & (tab.effect == "half")]
    lines += ["", "**Planning scenario (baseline 0.945, half effect):** power "
              + ", ".join(f"{r.power:.2f} at ρ = {r.rho}" for r in plan.itertuples()) + ". "
              + ("Below 0.80: recorded as \"H1 may be inconclusive\" (plan 6.11); the study is not stopped."
                 if (plan.power < 0.8).any() else "At least 0.80 for every ρ.")]
    (DOCS / "wp2_power.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
