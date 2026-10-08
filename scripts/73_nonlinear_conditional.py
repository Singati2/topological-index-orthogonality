"""Conditional (partial) NON-LINEAR association test: does candidate z carry
target information beyond the 30-index baseline X under non-linear dependence?

Responds to the referee point that scripts/61 (dCor vs Pearson among baseline
indices) never tests the actual question, and that dCor and |r| at a common
0.95 threshold are not comparable scales.

Method (double-ML / cross-fitted residualization):
  * 5-fold KFold (shuffle, seed). In each fold, fit
    RandomForestRegressor(300 trees, min_samples_leaf=5) of y on X and of z on X
    on the training folds, predict the held-out fold.
    ry = y - yhat_oof, rz = z - zhat_oof.  For BBBP (y in {0,1}) the RF
    regressor estimates P(y=1|X), so ry is a probability residual.
  * rho_nl = Pearson(rz, ry) with Fisher-z 95% CI (se = 1/sqrt(n-3)).
  * dCor(rz, ry) (Szekely-Rizzo V-statistic, double-centred distance matrices,
    O(n^2) memory) with a 199-permutation p-value. For the dCor part ONLY,
    datasets with n > 2000 are subsampled to 2000 molecules (fixed seed);
    rho_nl uses all molecules.
  * linear_pcor: OLS residualization of y and z on [1, X] (in-sample), Pearson.
  * Reference cases: REF_SCI_inspan (sum over edges 1/sqrt(du+dv); certified
    exactly in span([1,X]) by scripts/70), REF_noise (N(0,1) column),
    REF_leak (z = y + N(0, sd(y))).
  * BH adjustment of perm_p across all real candidates x datasets
    (reference rows are excluded from the BH family and reported with bh_q=NaN).

Outputs:
  results/nonlinear_conditional_tests.csv
  results/nonlinear_conditional_report.md
"""
from __future__ import annotations
import math
import os
import sys
import time

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold

THIS_DIR = os.path.abspath(os.path.dirname(__file__))
PROJECT = os.path.abspath(os.path.join(THIS_DIR, ".."))
sys.path.insert(0, PROJECT)

from src.load_data import load
from src.mol_to_graph import smiles_to_graph
from src.standard_indices import compute_all
from src.novel_candidates import CANDIDATE_INDICES

SEED = 20261005
DATASETS = ["esol", "freesolv", "lipophilicity", "bbbp"]
N_FOLDS = 5
N_TREES = 300
MIN_LEAF = 5
N_PERM = 199
DCOR_MAX_N = 2000
RES = os.path.join(PROJECT, "results")


# ---------------------------------------------------------------- helpers
def inverse_nirmala(G):
    deg = dict(G.degree())
    return float(sum(1.0 / math.sqrt(deg[u] + deg[v]) for u, v in G.edges()))


def impute(v):
    v = np.asarray(v, dtype=float)
    if np.any(~np.isfinite(v)):
        m = np.nanmean(np.where(np.isfinite(v), v, np.nan))
        v = np.where(np.isfinite(v), v, m)
    return v


def rf_oof(X, t, seed):
    """Cross-fitted out-of-fold RF predictions of t from X."""
    pred = np.empty_like(t, dtype=float)
    kf = KFold(N_FOLDS, shuffle=True, random_state=seed)
    for tr, te in kf.split(X):
        rf = RandomForestRegressor(n_estimators=N_TREES, min_samples_leaf=MIN_LEAF,
                                   random_state=seed, n_jobs=-1)
        rf.fit(X[tr], t[tr])
        pred[te] = rf.predict(X[te])
    return pred


def ols_resid(Xc, t):
    beta, *_ = np.linalg.lstsq(Xc, t, rcond=None)
    return t - Xc @ beta


def pearson(a, b):
    a = a - a.mean(); b = b - b.mean()
    d = math.sqrt(float(a @ a) * float(b @ b))
    return float(a @ b) / d if d > 0 else float("nan")


def fisher_ci(r, n):
    if not np.isfinite(r):
        return float("nan"), float("nan")
    z = np.arctanh(np.clip(r, -0.999999, 0.999999))
    se = 1.0 / math.sqrt(n - 3)
    return float(np.tanh(z - 1.959964 * se)), float(np.tanh(z + 1.959964 * se))


def _dcentre(x):
    D = np.abs(x[:, None] - x[None, :])
    return D - D.mean(0)[None, :] - D.mean(1)[:, None] + D.mean()


def dcor_perm(x, y, rng, n_perm=N_PERM):
    """dCor (V-statistic) of 1-D x, y + permutation p-value."""
    A = _dcentre(x); B = _dcentre(y)
    dvx = (A * A).mean(); dvy = (B * B).mean()
    denom = math.sqrt(dvx * dvy)
    if denom <= 0:
        return float("nan"), float("nan")
    stat = (A * B).mean()
    dc = math.sqrt(max(stat, 0.0) / denom)
    n = len(x); ge = 0
    for _ in range(n_perm):
        p = rng.permutation(n)
        if (A * B[np.ix_(p, p)]).mean() >= stat - 1e-15:
            ge += 1
    return dc, (1 + ge) / (1 + n_perm)


def bh(p):
    p = np.asarray(p, float); out = np.full_like(p, np.nan)
    ok = np.isfinite(p); pv = p[ok]; m = len(pv)
    if m == 0:
        return out
    o = np.argsort(pv); q = pv[o] * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    res = np.empty(m); res[o] = np.minimum(q, 1.0)
    out[ok] = res
    return out


# ---------------------------------------------------------------- main
def run_dataset(name):
    t0 = time.time()
    df = load(name)
    graphs, ys = [], []
    for s, t in zip(df["smiles"], df["target"]):
        G = smiles_to_graph(s)
        if G is None:
            continue
        graphs.append(G); ys.append(t)
    y = np.asarray(ys, float); n = len(y)
    X = np.column_stack([impute(c) for c in pd.DataFrame([compute_all(G) for G in graphs]).values.T])
    Xs = (X - X.mean(0)) / np.where(X.std(0) > 0, X.std(0), 1.0)
    Xc = np.column_stack([np.ones(n), Xs])
    print(f"[{name}] n={n}, baseline built ({time.time()-t0:.0f}s)", flush=True)

    cands = {}
    for cname, fn in CANDIDATE_INDICES.items():
        vals = []
        for G in graphs:
            try:
                vals.append(fn(G))
            except Exception:
                vals.append(float("nan"))
        cands[cname] = impute(vals)
    rng_ref = np.random.default_rng(SEED)
    cands["REF_SCI_inspan"] = impute([inverse_nirmala(G) for G in graphs])
    cands["REF_noise"] = rng_ref.standard_normal(n)
    cands["REF_leak"] = y + rng_ref.normal(0.0, y.std(), n)

    ry = y - rf_oof(Xs, y, SEED)
    ry_lin = ols_resid(Xc, y)
    sub_rng = np.random.default_rng(SEED)
    sub = np.sort(sub_rng.choice(n, DCOR_MAX_N, replace=False)) if n > DCOR_MAX_N else np.arange(n)

    rows = []
    for cname, z in cands.items():
        if np.std(z) < 1e-12:
            rows.append(dict(dataset=name, candidate=cname, n=n, n_dcor=len(sub),
                             linear_pcor=np.nan, rho_nl=np.nan, ci_lo=np.nan, ci_hi=np.nan,
                             dcor=np.nan, perm_p=np.nan, note="constant"))
            continue
        zs = (z - z.mean()) / z.std()
        rz = zs - rf_oof(Xs, zs, SEED)
        rz_lin = ols_resid(Xc, zs)
        rel_lin = float(np.linalg.norm(rz_lin) / np.linalg.norm(zs))
        lin = pearson(rz_lin, ry_lin) if rel_lin > 1e-6 else 0.0
        r = pearson(rz, ry)
        lo, hi = fisher_ci(r, n)
        dc, pp = dcor_perm(rz[sub], ry[sub], np.random.default_rng(SEED))
        rows.append(dict(dataset=name, candidate=cname, n=n, n_dcor=len(sub),
                         linear_pcor=lin, rho_nl=r, ci_lo=lo, ci_hi=hi, dcor=dc, perm_p=pp,
                         lin_rel_resid=rel_lin,
                         rf_rz_var_frac=float(rz.var()),
                         note="in_span(linear)" if rel_lin <= 1e-6 else ""))
        print(f"  {cname:15s} lin={lin:+.3f} rho_nl={r:+.3f} dcor={dc:.3f} p={pp:.3f}", flush=True)
    print(f"[{name}] done ({time.time()-t0:.0f}s)", flush=True)
    return rows


def main():
    rows = []
    for d in DATASETS:
        rows += run_dataset(d)
    out = pd.DataFrame(rows)
    is_ref = out["candidate"].str.startswith("REF_")
    out["bh_q"] = np.nan
    out.loc[~is_ref, "bh_q"] = bh(out.loc[~is_ref, "perm_p"].values)
    cols = ["dataset", "candidate", "linear_pcor", "rho_nl", "ci_lo", "ci_hi", "dcor",
            "perm_p", "bh_q", "n", "n_dcor", "lin_rel_resid", "rf_rz_var_frac", "note"]
    out = out[cols]
    out.to_csv(os.path.join(RES, "nonlinear_conditional_tests.csv"), index=False)

    # ------------------------------------------------------------ report
    min_p = 1.0 / (N_PERM + 1)
    L = ["# Conditional non-linear association test (scripts/73_nonlinear_conditional.py)\n",
         "Question: does candidate z carry target information beyond the 30-index baseline X "
         "under NON-LINEAR dependence?\n",
         "## Method\n",
         f"- Cross-fitted (double-ML) residualization: {N_FOLDS}-fold KFold; within each fold "
         f"RandomForestRegressor({N_TREES} trees, min_samples_leaf={MIN_LEAF}, seed={SEED}) of y on X "
         "and of z on X fitted on training folds, predicted on the held-out fold; ry = y - yhat, rz = z - zhat "
         "(z standardized first). BBBP: y in {0,1}, RF regressor gives a probability residual.",
         "- rho_nl = Pearson(rz, ry), Fisher-z 95% CI with se = 1/sqrt(n-3), all molecules.",
         f"- dcor = distance correlation(rz, ry) (V-statistic, from scratch, O(n^2) memory), "
         f"permutation p-value with {N_PERM} permutations (minimum attainable p = {min_p:.3f}). "
         f"**For the dCor part only, Lipophilicity and BBBP (n > {DCOR_MAX_N}) are subsampled to "
         f"{DCOR_MAX_N} molecules with a fixed seed**; rho_nl uses all molecules.",
         "- linear_pcor = Pearson of OLS residuals of y and z on [1, X] (in-sample).",
         "- BH adjustment across all real candidates x datasets (reference rows excluded from the family).",
         "- Reference rows: REF_SCI_inspan (sum 1/sqrt(du+dv); exactly in linear span of X), "
         "REF_noise (N(0,1) column; calibration), REF_leak (z = y + N(0, sd(y)); positive control).\n",
         "Caveat: a candidate's non-linear residual rz is exactly zero only if the RF reproduces z perfectly; "
         "RF cannot extrapolate/represent exact linear identities, so an in-span index has rz != 0 "
         "(nuisance error). Cross-fitting makes this error independent of the held-out y noise, so "
         "rho_nl stays centred near 0 unless the nuisance error itself tracks ry.\n"]
    pd.set_option("display.width", 200)
    fmt = out.copy()
    for c in ["linear_pcor", "rho_nl", "ci_lo", "ci_hi", "dcor", "lin_rel_resid", "rf_rz_var_frac"]:
        fmt[c] = fmt[c].map(lambda v: f"{v:+.3f}" if pd.notna(v) else "nan")
    fmt["perm_p"] = fmt["perm_p"].map(lambda v: f"{v:.3f}" if pd.notna(v) else "nan")
    fmt["bh_q"] = fmt["bh_q"].map(lambda v: f"{v:.3f}" if pd.notna(v) else "-")
    show = ["candidate", "linear_pcor", "rho_nl", "ci_lo", "ci_hi", "dcor", "perm_p", "bh_q", "note"]
    for d in DATASETS:
        sub = fmt[fmt["dataset"] == d]
        L.append(f"## {d} (n = {out.loc[out.dataset == d, 'n'].iloc[0]}, "
                 f"n_dcor = {out.loc[out.dataset == d, 'n_dcor'].iloc[0]})\n")
        L.append("```\n" + sub[show].to_string(index=False) + "\n```\n")
    real = out[~is_ref]
    flag = real[(real.bh_q < 0.05) | (real.rho_nl.abs() >= 0.10)]
    L.append("## Flagged real candidates (bh_q < 0.05 or |rho_nl| >= 0.10)\n")
    L.append("```\n" + (flag[["dataset", "candidate", "linear_pcor", "rho_nl", "ci_lo", "ci_hi",
                              "dcor", "perm_p", "bh_q"]].round(3).to_string(index=False)
                        if len(flag) else "(none)") + "\n```\n")
    lin_hit = real.linear_pcor.abs() >= 0.10
    nl_hit = real.rho_nl.abs() >= 0.10
    L.append("## Linear vs non-linear screen agreement (|.| >= 0.10)\n")
    L.append(f"- both: {int((lin_hit & nl_hit).sum())}; linear only: {int((lin_hit & ~nl_hit).sum())}; "
             f"non-linear only: {int((~lin_hit & nl_hit).sum())}; neither: {int((~lin_hit & ~nl_hit).sum())}; "
             f"BH q < 0.05 (dCor perm): {int((real.bh_q < 0.05).sum())} of {int(real.perm_p.notna().sum())}.\n")
    with open(os.path.join(RES, "nonlinear_conditional_report.md"), "w") as fh:
        fh.write("\n".join(L))
    print("\n".join(L))


if __name__ == "__main__":
    main()
