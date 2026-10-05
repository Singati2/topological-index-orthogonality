"""Benchmark v2: leakage-controlled, curated, correctly-tested downstream benchmark.

Replaces the downstream benchmark of scripts/11_ml_benchmark.py and
scripts/18_pipeline_then_lasso.py, addressing referee findings:

  (1) exact duplicates: SMILES canonicalized with RDKit (largest fragment
      kept), exact canonical duplicates merged (regression: averaged if the
      label range is within DUP_TOL, else dropped; BBBP: conflicting labels
      dropped).
  (2) graph-level leakage: many molecules share an identical hydrogen-
      suppressed graph (identical 30-index feature vector). Groups are the
      unlabeled Weisfeiler-Lehman hash (iterations=4) of that graph, and CV
      is GroupKFold (StratifiedGroupKFold for BBBP), 5 repeats x 5 folds.
  (3) fold dependence: Nadeau-Bengio corrected repeated k-fold t-test
      (variance factor 1/(k r) + n_test/n_train, df = k r - 1).
  (4) degenerate folds: if combined pruning keeps zero features the fold is
      scored with an intercept-only model and counted, not dropped.
  (5) ceiling: RDKit 2D descriptor RandomForest.

All preprocessing (standardization, pruning, pcor filter, LASSO/L1 penalty
selection via grouped inner CV, descriptor column filtering / median
imputation) is fit inside the outer training fold only.

Pruning functions are copied verbatim from scripts/18_pipeline_then_lasso.py
(which are themselves byte-identical copies of scripts/11_ml_benchmark.py).

Outputs (results/):
  benchmark_v2_perfold.csv, benchmark_v2_summary.csv, benchmark_v2_tests.csv,
  benchmark_v2_curation.csv, benchmark_v2_report.md
"""
from __future__ import annotations

import os
import sys
import time
import warnings

import numpy as np
import pandas as pd
import networkx as nx

THIS_DIR = os.path.abspath(os.path.dirname(__file__))
PROJECT = os.path.abspath(os.path.join(THIS_DIR, ".."))
sys.path.insert(0, PROJECT)

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors
from scipy import stats
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.linear_model import LassoCV, LinearRegression, LogisticRegressionCV
from sklearn.metrics import mean_squared_error, roc_auc_score
from sklearn.model_selection import GroupKFold, StratifiedGroupKFold
from sklearn.preprocessing import StandardScaler

from src.load_data import load
from src.mol_to_graph import smiles_to_graph
from src.standard_indices import compute_all

RDLogger.DisableLog("rdApp.*")
warnings.filterwarnings("ignore")

SEED = 20261005
N_REPEATS = 5
N_FOLDS = 5
PAIRWISE_THRESHOLD = 0.95
PCOR_THRESHOLD = 0.10
DUP_TOL = 0.5   # max label range (target units) for merging regression duplicates
WL_ITER = 4

DATASETS = [("esol", "regression"), ("freesolv", "regression"),
            ("lipophilicity", "regression"), ("bbbp", "classification")]

MODELS = ["a_full30_lasso", "b_pairwise_lasso", "c_combined_lasso",
          "d_full30_rf", "e_pairwise_rf", "f_rdkit2d_rf",
          "g_mij_rf", "g_mij_lasso"]
COMPARISONS = [("b_pairwise_lasso", "a_full30_lasso"),
               ("c_combined_lasso", "a_full30_lasso"),
               ("e_pairwise_rf", "d_full30_rf"),
               ("g_mij_rf", "d_full30_rf"),
               ("f_rdkit2d_rf", "d_full30_rf")]

MIJ_PAIRS = [(i, j) for i in range(1, 5) for j in range(i, 5)]
MIJ_COLS = [f"m{i}{j}" for i, j in MIJ_PAIRS]
# Descriptors that segfault in this environment (rdkit 2022.09 compiled against
# numpy 1.x, running under numpy 2.x: they call numpy-returning C routines such as
# GetDistanceMatrix). Probed one by one; excluded from the ceiling and reported.
RDKIT_EXCLUDED = {"MaxEStateIndex", "MinEStateIndex", "MaxAbsEStateIndex",
                  "MinAbsEStateIndex", "BalabanJ", "BertzCT", "Ipc"}
RDKIT_EXCLUDED |= {f"EState_VSA{i}" for i in range(1, 12)}
RDKIT_EXCLUDED |= {f"VSA_EState{i}" for i in range(1, 11)}
RDKIT_DESC = [(n, f) for n, f in Descriptors.descList if n not in RDKIT_EXCLUDED]
RDKIT_NAMES = [n for n, _ in RDKIT_DESC]


# ---- pruning: verbatim from scripts/18_pipeline_then_lasso.py ----------------
def select_pairwise_pruned(X: pd.DataFrame, threshold: float):
    kept = []
    for col in X.columns:
        col_vals = X[col].values
        if np.std(col_vals) < 1e-12:
            continue
        if not kept:
            kept.append(col)
            continue
        kept_mat = X[kept].values
        corrs = np.abs(np.corrcoef(col_vals, kept_mat, rowvar=False)[0, 1:])
        if np.nanmax(corrs) < threshold:
            kept.append(col)
    return kept


def select_combined_pruned(X: pd.DataFrame, y: np.ndarray,
                           pair_threshold: float, pcor_threshold: float):
    pair_kept = select_pairwise_pruned(X, pair_threshold)
    if len(pair_kept) <= 1:
        return pair_kept
    final_kept = []
    for col in pair_kept:
        others = [c for c in pair_kept if c != col]
        if not others:
            final_kept.append(col)
            continue
        Xo = X[others].values
        z = X[col].values
        if np.std(z) < 1e-12:
            continue
        z_res = z - LinearRegression().fit(Xo, z).predict(Xo)
        y_res = y - LinearRegression().fit(Xo, y).predict(Xo)
        if np.std(z_res) < 1e-12 or np.std(y_res) < 1e-12:
            continue
        pcor = float(np.corrcoef(z_res, y_res)[0, 1])
        if abs(pcor) >= pcor_threshold:
            final_kept.append(col)
    return final_kept


# ---- curation -----------------------------------------------------------------
def canonicalize(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None, None, False
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    multi = len(frags) > 1
    if multi:
        mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())
    can = Chem.MolToSmiles(mol)
    mol2 = Chem.MolFromSmiles(can)
    return can, mol2, multi


def mij_counts(G):
    deg = dict(G.degree())
    c = {k: 0 for k in MIJ_COLS}
    for u, v in G.edges():
        a, b = sorted((min(deg[u], 4), min(deg[v], 4)))
        c[f"m{a}{b}"] += 1
    return c


def rdkit_desc(mol):
    out = []
    for _, fn in RDKIT_DESC:
        try:
            out.append(float(fn(mol)))
        except Exception:
            out.append(np.nan)
    return out


def curate(name, task, features=True):
    raw = load(name)
    n_raw = len(raw)
    recs, n_invalid, n_multi = [], 0, 0
    for smi, y in zip(raw["smiles"], raw["target"]):
        can, mol, multi = canonicalize(str(smi))
        if can is None or mol is None or pd.isna(y):
            n_invalid += 1
            continue
        n_multi += int(multi)
        recs.append((can, float(y)))
    df = pd.DataFrame(recs, columns=["can", "target"])
    n_valid = len(df)

    # exact canonical duplicates
    grp = df.groupby("can")["target"]
    agg = grp.agg(["count", "mean", "min", "max"]).reset_index()
    dup_groups = agg[agg["count"] > 1]
    n_dup_groups = len(dup_groups)
    n_dup_rows = int(dup_groups["count"].sum())
    if task == "regression":
        conflict = (agg["max"] - agg["min"]) > DUP_TOL
    else:
        conflict = agg["min"] != agg["max"]
    n_conflict_groups = int((conflict & (agg["count"] > 1)).sum())
    agg = agg[~conflict]
    df = agg[["can", "mean"]].rename(columns={"mean": "target"}).reset_index(drop=True)
    n_dedup = len(df)

    # graphs, features
    rows = []
    n_nograph = 0
    n_deg_gt4 = 0
    for can, y in zip(df["can"], df["target"]):
        G = smiles_to_graph(can)
        if G is None:
            n_nograph += 1
            continue
        mol = Chem.MolFromSmiles(can)
        h = nx.weisfeiler_lehman_graph_hash(G, iterations=WL_ITER)
        if max(dict(G.degree()).values()) > 4:
            n_deg_gt4 += 1
        r = {"can": can, "target": y, "group": h}
        if not features:
            rows.append(r)
            continue
        r.update({f"TI_{k}": v for k, v in compute_all(G).items()})
        r.update(mij_counts(G))
        r.update({f"RD_{k}": v for k, v in zip(RDKIT_NAMES, rdkit_desc(mol))})
        rows.append(r)
    D = pd.DataFrame(rows)
    n_final = len(D)
    if task == "classification":
        D["target"] = D["target"].astype(int)

    gsize = D.groupby("group")["group"].transform("size")
    n_groups = D["group"].nunique()
    pct_shared = 100.0 * float((gsize > 1).mean())
    multi = D[gsize > 1]
    if len(multi):
        dev = multi["target"] - multi.groupby("group")["target"].transform("mean")
        dof = len(multi) - multi["group"].nunique()
        within_sd = float(np.sqrt((dev ** 2).sum() / dof)) if dof > 0 else np.nan
    else:
        within_sd = np.nan
    # Graph-only RMSE floor: any model that sees only the H-suppressed graph must
    # predict one value per graph; the best in-sample choice is the group mean.
    ss_within = float(((D["target"] - D.groupby("group")["target"].transform("mean")) ** 2).sum())
    floor_all = float(np.sqrt(ss_within / n_final))             # sqrt(mean within-group var), all mols
    floor_unb = float(np.sqrt(ss_within / (n_final - n_groups))) if n_final > n_groups else np.nan
    overall_sd = float(D["target"].std(ddof=1))
    mixed = np.nan
    if task == "classification":
        g = multi.groupby("group")["target"].nunique()
        mixed = float((g > 1).mean()) if len(g) else np.nan

    cur = {
        "dataset": name, "n_raw": n_raw, "n_invalid_smiles": n_invalid,
        "n_multifragment_stripped": n_multi, "n_valid": n_valid,
        "n_dup_canonical_groups": n_dup_groups, "n_rows_in_dup_groups": n_dup_rows,
        "n_conflicting_dup_groups_dropped": n_conflict_groups,
        "n_after_dedup": n_dedup, "n_no_graph_dropped(<2 heavy atoms)": n_nograph,
        "n_final": n_final, "n_distinct_graphs": n_groups,
        "pct_mols_graph_shared": round(pct_shared, 2),
        "n_mols_max_degree_gt4(mij_capped)": n_deg_gt4,
        "within_graph_target_sd(pooled, shared groups)": round(within_sd, 4),
        "overall_target_sd": round(overall_sd, 4),
        "graph_only_rmse_floor(sqrt SSwithin/N)": round(floor_all, 4),
        "graph_only_rmse_floor_dofcorr(sqrt SSwithin/(N-G))": round(floor_unb, 4),
        "frac_shared_groups_mixed_labels": (round(mixed, 4) if task == "classification" else ""),
        "max_group_size": int(gsize.max()),
    }
    return D, cur


# ---- models -------------------------------------------------------------------
def inner_cv(task, Xtr, ytr, gtr, seed):
    if task == "regression":
        cv = GroupKFold(n_splits=5, shuffle=True, random_state=seed)
        return list(cv.split(Xtr, ytr, gtr))
    cv = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=seed)
    return list(cv.split(Xtr, ytr, gtr))


def fit_linear(task, Xtr, ytr, Xte, gtr, seed):
    """Returns (pred_or_prob, n_active, intercept_only_flag)."""
    if Xtr.shape[1] == 0:
        if task == "regression":
            return np.full(len(Xte), ytr.mean()), 0, True
        return np.full(len(Xte), ytr.mean()), 0, True
    sc = StandardScaler().fit(Xtr)
    Xtr_s, Xte_s = sc.transform(Xtr), sc.transform(Xte)
    splits = inner_cv(task, Xtr_s, ytr, gtr, seed)
    if task == "regression":
        m = LassoCV(cv=splits, max_iter=10000, n_alphas=50).fit(Xtr_s, ytr)
        coef = m.coef_
        pred = m.predict(Xte_s)
    else:
        m = LogisticRegressionCV(cv=splits, l1_ratios=[1.0], solver="liblinear",
                                 Cs=20, max_iter=2000, random_state=SEED,
                                 scoring="roc_auc", class_weight="balanced")
        m.fit(Xtr_s, ytr)
        coef = np.asarray(m.coef_).ravel()
        pred = m.predict_proba(Xte_s)[:, 1]
    return pred, int(np.sum(np.abs(coef) > 1e-10)), False


def fit_rf(task, Xtr, ytr, Xte):
    if task == "regression":
        m = RandomForestRegressor(n_estimators=300, random_state=SEED, n_jobs=-1)
        m.fit(Xtr, ytr)
        return m.predict(Xte)
    m = RandomForestClassifier(n_estimators=300, random_state=SEED, n_jobs=-1)
    m.fit(Xtr, ytr)
    return m.predict_proba(Xte)[:, 1]


def rdkit_prep(Rtr, Rte):
    Rtr = Rtr.astype(float).copy(); Rte = Rte.astype(float).copy()
    Rtr[~np.isfinite(Rtr) | (np.abs(Rtr) > 1e30)] = np.nan
    Rte[~np.isfinite(Rte) | (np.abs(Rte) > 1e30)] = np.nan
    # drop columns with any non-finite value in TRAINING fold, and constant columns
    keep = ~np.isnan(Rtr).any(axis=0)
    keep &= np.nanstd(np.where(keep, Rtr, 0.0), axis=0) > 0
    Rtr, Rte = Rtr[:, keep], Rte[:, keep]
    med = np.median(Rtr, axis=0)
    idx = np.where(np.isnan(Rte))
    Rte[idx] = med[idx[1]]
    return Rtr, Rte, int(keep.sum())


def score(task, y, p):
    if task == "regression":
        return float(np.sqrt(mean_squared_error(y, p)))
    return float(roc_auc_score(y, p))


def run_dataset(name, task, D):
    ti_cols = [c for c in D.columns if c.startswith("TI_")]
    rd_cols = [c for c in D.columns if c.startswith("RD_")]
    X_ti = D[ti_cols]
    X_mij = D[MIJ_COLS].values.astype(float)
    X_rd = D[rd_cols].values
    y = D["target"].values
    groups = D["group"].values
    rows = []
    for r in range(N_REPEATS):
        seed_r = SEED + r
        if task == "regression":
            cv = GroupKFold(n_splits=N_FOLDS, shuffle=True, random_state=seed_r)
            splits = list(cv.split(X_ti, y, groups))
        else:
            cv = StratifiedGroupKFold(n_splits=N_FOLDS, shuffle=True, random_state=seed_r)
            splits = list(cv.split(X_ti, y, groups))
        for f, (tr, te) in enumerate(splits):
            assert not set(groups[tr]) & set(groups[te])
            ytr, yte, gtr = y[tr], y[te], groups[tr]
            Xtr_df = X_ti.iloc[tr].reset_index(drop=True)
            Xte_df = X_ti.iloc[te].reset_index(drop=True)
            pw = select_pairwise_pruned(Xtr_df, PAIRWISE_THRESHOLD)
            cb = select_combined_pruned(Xtr_df, ytr.astype(float),
                                        PAIRWISE_THRESHOLD, PCOR_THRESHOLD)
            base = {"dataset": name, "task": task, "repeat": r, "fold": f,
                    "n_train": len(tr), "n_test": len(te),
                    "metric": "rmse" if task == "regression" else "roc_auc"}

            def add(model, val, n_in, n_active, io=False):
                rows.append({**base, "model": model, "value": val,
                             "n_features_input": n_in, "n_features_active": n_active,
                             "intercept_only": int(io)})

            def sub(cols):
                return Xtr_df[cols].values, Xte_df[cols].values

            for model, cols in (("a_full30_lasso", list(ti_cols)),
                                ("b_pairwise_lasso", pw), ("c_combined_lasso", cb)):
                Xa, Xb = sub(cols) if cols else (np.empty((len(tr), 0)), np.empty((len(te), 0)))
                p, na, io = fit_linear(task, Xa, ytr, Xb, gtr, seed_r)
                add(model, score(task, yte, p), len(cols), na, io)

            p = fit_rf(task, Xtr_df.values, ytr, Xte_df.values)
            add("d_full30_rf", score(task, yte, p), len(ti_cols), len(ti_cols))
            Xa, Xb = sub(pw)
            p = fit_rf(task, Xa, ytr, Xb)
            add("e_pairwise_rf", score(task, yte, p), len(pw), len(pw))
            Rtr, Rte, nk = rdkit_prep(X_rd[tr], X_rd[te])
            p = fit_rf(task, Rtr, ytr, Rte)
            add("f_rdkit2d_rf", score(task, yte, p), nk, nk)
            p = fit_rf(task, X_mij[tr], ytr, X_mij[te])
            add("g_mij_rf", score(task, yte, p), 10, 10)
            p, na, io = fit_linear(task, X_mij[tr], ytr, X_mij[te], gtr, seed_r)
            add("g_mij_lasso", score(task, yte, p), 10, na, io)
        print(f"  [{name}] repeat {r} done ({time.strftime('%H:%M:%S')})", flush=True)
    return pd.DataFrame(rows)


# ---- statistics ----------------------------------------------------------------
def nb_test(diff, ratio, k, r):
    J = k * r
    m = float(np.mean(diff))
    s2 = float(np.var(diff, ddof=1))
    se = np.sqrt((1.0 / J + ratio) * s2)
    dfree = J - 1
    if se == 0:
        return m, se, np.nan, np.nan, m, m, dfree
    t = m / se
    p = 2 * stats.t.sf(abs(t), dfree)
    tc = stats.t.ppf(0.95, dfree)
    return m, se, t, p, m - tc * se, m + tc * se, dfree


def holm(pvals):
    p = np.asarray(pvals, float)
    order = np.argsort(p)
    adj = np.empty_like(p)
    running = 0.0
    n = len(p)
    for i, idx in enumerate(order):
        running = max(running, min(1.0, (n - i) * p[idx]))
        adj[idx] = running
    return adj


def tests_for(per):
    out = []
    for ds, sub in per.groupby("dataset", sort=False):
        task = sub["task"].iloc[0]
        piv = sub.pivot_table(index=["repeat", "fold"], columns="model", values="value")
        nt = sub.drop_duplicates(["repeat", "fold"])
        ratio = float((nt["n_test"] / nt["n_train"]).mean())
        recs = []
        for mdl, ref in COMPARISONS:
            diff = (piv[mdl] - piv[ref]).values
            refmean = float(piv[ref].mean())
            margin = 0.05 * refmean if task == "regression" else 0.05 * (refmean - 0.5)
            m, se, t, p, lo, hi, dfree = nb_test(diff, ratio, N_FOLDS, N_REPEATS)
            # TOST (two one-sided corrected t-tests, alpha=0.05 each)
            if se > 0:
                p_lower = stats.t.sf((m + margin) / se, dfree)
                p_upper = stats.t.cdf((m - margin) / se, dfree)
                p_tost = float(max(p_lower, p_upper))
            else:
                p_tost = 0.0 if abs(m) < margin else 1.0
            equiv = (lo > -margin) and (hi < margin)
            recs.append({
                "dataset": ds, "metric": "rmse" if task == "regression" else "roc_auc",
                "comparison": f"{mdl} - {ref}", "ref_mean": refmean,
                "model_mean": float(piv[mdl].mean()), "mean_diff": m,
                "corrected_se": se, "t": t, "df": dfree, "p_raw": p,
                "ci90_lo": lo, "ci90_hi": hi, "tost_margin": margin,
                "p_tost": p_tost,
                "tost_verdict": "equivalent" if equiv else "not shown equivalent",
                "n_test_over_n_train": ratio,
            })
        adj = holm([r_["p_raw"] for r_ in recs])
        for r_, a in zip(recs, adj):
            r_["p_holm"] = float(a)
        out.extend(recs)
    return pd.DataFrame(out)


def main():
    t0 = time.time()
    res = os.path.join(PROJECT, "results")
    cur_rows, per_all = [], []
    for name, task in DATASETS:
        print(f"--- {name} ---", flush=True)
        D, cur = curate(name, task)
        cur_rows.append(cur)
        print("  curation:", cur, flush=True)
        per_all.append(run_dataset(name, task, D))
    cur = pd.DataFrame(cur_rows)
    per = pd.concat(per_all, ignore_index=True)

    summ = (per.groupby(["dataset", "model"], sort=False)
            .agg(metric=("metric", "first"), mean=("value", "mean"),
                 sd=("value", "std"), mean_n_features=("n_features_active", "mean"),
                 mean_n_features_input=("n_features_input", "mean"),
                 n_intercept_only_folds=("intercept_only", "sum"),
                 n_folds=("value", "size"))
            .reset_index())
    tests = tests_for(per)

    cur.to_csv(os.path.join(res, "benchmark_v2_curation.csv"), index=False)
    per.to_csv(os.path.join(res, "benchmark_v2_perfold.csv"), index=False)
    summ.to_csv(os.path.join(res, "benchmark_v2_summary.csv"), index=False)
    tests.to_csv(os.path.join(res, "benchmark_v2_tests.csv"), index=False)

    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 50)
    with open(os.path.join(res, "benchmark_v2_report.md"), "w") as fh:
        fh.write("# Benchmark v2 (curated, graph-grouped repeated CV)\n\n")
        fh.write("Generated by `scripts/72_benchmark_v2.py`.\n\n")
        fh.write(f"- Seed {SEED}; {N_REPEATS} repeats x {N_FOLDS} folds; groups = unlabeled "
                 f"WL hash (iterations={WL_ITER}) of the hydrogen-suppressed graph; "
                 "GroupKFold (regression) / StratifiedGroupKFold (BBBP), shuffled per repeat.\n")
        fh.write("- Curation: RDKit canonical SMILES, largest fragment; exact canonical "
                 f"duplicates averaged if label range <= {DUP_TOL} (regression) else dropped; "
                 "BBBP conflicting duplicates dropped; molecules with <2 heavy atoms dropped "
                 "(no graph).\n")
        fh.write("- All preprocessing (scaling, pairwise |r|<0.95 pruning, |pcor|>=0.10 "
                 "filter, LASSO/L1 penalty via grouped inner 5-fold CV, RDKit column "
                 "filtering + median imputation) fit on the training fold only. Zero-feature "
                 "combined-pruned folds use an intercept-only model and are counted.\n")
        fh.write(f"- RDKit ceiling uses {len(RDKIT_NAMES)} of {len(Descriptors.descList)} "
                 f"descriptors: {len(RDKIT_EXCLUDED)} EState/BalabanJ/BertzCT/Ipc descriptors "
                 "segfault under the rdkit-2022.09/numpy-2 ABI mismatch and were excluded.\n")
        fh.write("- m_ij: edge counts by sorted end-vertex degree pair, degrees capped at 4.\n")
        fh.write("- Tests: Nadeau-Bengio corrected repeated k-fold t (variance factor "
                 "1/(k r) + n_test/n_train, df = k r - 1); 90% corrected CI; TOST margin = "
                 "5% of reference mean RMSE or 0.05*(AUC_ref - 0.5); Holm across the 5 "
                 "comparisons within each dataset. Differences are model - reference "
                 "(RMSE: positive = worse; AUC: positive = better).\n\n")
        fh.write("## Curation\n\n```\n" + cur.T.to_string() + "\n```\n\n")
        fh.write("## Summary\n\n```\n" + summ.round(4).to_string(index=False) + "\n```\n\n")
        cols = ["dataset", "comparison", "ref_mean", "model_mean", "mean_diff",
                "ci90_lo", "ci90_hi", "tost_margin", "p_tost", "tost_verdict",
                "p_raw", "p_holm"]
        fh.write("## Tests\n\n```\n" + tests[cols].round(4).to_string(index=False) + "\n```\n\n")
        fh.write(f"Runtime: {time.time() - t0:.0f} s\n")
    print(open(os.path.join(res, "benchmark_v2_report.md")).read())


def curation_only():
    res = os.path.join(PROJECT, "results")
    cur = pd.DataFrame([curate(n, t, features=False)[1] for n, t in DATASETS])
    cur.to_csv(os.path.join(res, "benchmark_v2_curation.csv"), index=False)
    print(cur.T.to_string())
    rp = os.path.join(res, "benchmark_v2_report.md")
    if os.path.exists(rp):
        txt = open(rp).read()
        a = txt.index("## Curation"); b = txt.index("## Summary")
        txt = txt[:a] + "## Curation\n\n```\n" + cur.T.to_string() + "\n```\n\n" + txt[b:]
        open(rp, "w").write(txt)


if __name__ == "__main__":
    if "--curation-only" in sys.argv:
        curation_only()
    else:
        main()
