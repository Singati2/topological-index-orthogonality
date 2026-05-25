"""Pipeline-then-LASSO: matched-model-class benchmark for Paper 2 §6.3.2.

Reviewer concern this script resolves: the LASSO/ElasticNet/mRMR/Boruta
comparators in scripts/15_feature_selection_comparison.py are evaluated
under a *linear* model class, while the screen-pipeline rows of
results/ml_benchmark.csv are *RandomForest*. The cross-model-class
comparison is honestly disclaimed in the manuscript, but a stricter
reviewer can still ask: ``the only honest comparison is screen-pipeline-
then-LASSO vs LASSO-only; you didn't run it.''

This script runs exactly that comparison. For each of the four datasets,
under the same KFold seeds as scripts/11_ml_benchmark.py (RNG_SEED=42,
N_FOLDS=5):

  full                LASSO (or L1-Logistic on BBBP) on the full 30 features
  pairwise_pruned     same model on the in-fold pairwise-pruned subset
  combined_pruned     same model on the in-fold combined-pruned subset

All feature selection (pairwise pruning, combined pruning) is done
*inside each training fold* using the exact same utilities from
scripts/11_ml_benchmark.py, so the pruned subsets exactly match the
ones used by the RF benchmark.

Output:
  results/paper2_pipeline_then_lasso.csv

This complements (does not replace) the existing
paper2_feature_selection_comparison.csv. The two together form Table 10
in the manuscript: the existing CSV gives LASSO/ElasticNet/mRMR/Boruta
as *external* feature-selection comparators (cross-model-class); the
new CSV gives the matched-model-class same-model-class pipeline-then-
LASSO numbers.
"""
from __future__ import annotations
import csv
import os
import sys
import warnings

import numpy as np
import pandas as pd

THIS_DIR = os.path.abspath(os.path.dirname(__file__))
PROJECT = os.path.abspath(os.path.join(THIS_DIR, ".."))
sys.path.insert(0, PROJECT)
sys.path.insert(0, THIS_DIR)

from sklearn.linear_model import LassoCV, LogisticRegressionCV, LinearRegression
from sklearn.model_selection import KFold, StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    mean_squared_error, mean_absolute_error, r2_score,
    roc_auc_score, accuracy_score,
)
from scipy import stats as _scipy_stats


# === In-fold feature selectors: BYTE-IDENTICAL to those in scripts/11_ml_benchmark.py ===
# (Inlined rather than imported because script 11 has __main__-only side effects
# that exec/import-as-module triggers; the two functions below are copy-paste
# from scripts/11_ml_benchmark.py:90-140.)

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

warnings.filterwarnings("ignore")

RNG_SEED = 42                  # matches scripts/11_ml_benchmark.py
N_FOLDS = 5
PAIRWISE_THRESHOLD = 0.95
PCOR_THRESHOLD = 0.10

REGRESSION_DATASETS = [
    ("esol",          "logS"),
    ("freesolv",      "dG_hyd"),
    ("lipophilicity", "logD"),
]
CLASSIFICATION_DATASETS = [
    ("bbbp", "BBB_penetration"),
]
CONFIGS = ["full", "pairwise_pruned", "combined_pruned"]


def _load_dataset(name: str):
    path = os.path.join(PROJECT, "results", name, "descriptors_baseline.csv")
    if not os.path.exists(path):
        return None, None
    df = pd.read_csv(path)
    target_col = df.columns[-1]
    baseline_cols = list(df.columns[:30])
    X = df[baseline_cols]
    y = df[target_col].values.astype(float)
    return X, y


def _select_features(X_train_df, y_train, config):
    if config == "full":
        return list(X_train_df.columns)
    if config == "pairwise_pruned":
        return select_pairwise_pruned(X_train_df, PAIRWISE_THRESHOLD)
    if config == "combined_pruned":
        return select_combined_pruned(
            X_train_df, y_train, PAIRWISE_THRESHOLD, PCOR_THRESHOLD,
        )
    raise ValueError(config)


def _fit_eval_regression(Xtr, ytr, Xte, yte):
    if Xtr.shape[1] == 0:
        return {"rmse": float("nan"), "mae": float("nan"), "r2": float("nan"),
                "n_features": 0}
    sc = StandardScaler().fit(Xtr)
    Xtr_s = sc.transform(Xtr); Xte_s = sc.transform(Xte)
    model = LassoCV(cv=5, random_state=RNG_SEED, max_iter=10000, n_alphas=50)
    model.fit(Xtr_s, ytr)
    pred = model.predict(Xte_s)
    nf = int(np.sum(np.abs(model.coef_) > 1e-10))
    return {
        "rmse": float(np.sqrt(mean_squared_error(yte, pred))),
        "mae": float(mean_absolute_error(yte, pred)),
        "r2": float(r2_score(yte, pred)),
        "n_features": nf,
    }


def _fit_eval_classification(Xtr, ytr, Xte, yte):
    if Xtr.shape[1] == 0:
        return {"roc_auc": float("nan"), "accuracy": float("nan"), "n_features": 0}
    sc = StandardScaler().fit(Xtr)
    Xtr_s = sc.transform(Xtr); Xte_s = sc.transform(Xte)
    model = LogisticRegressionCV(
        cv=5, penalty="l1", solver="liblinear",
        Cs=20, max_iter=2000, random_state=RNG_SEED,
        scoring="roc_auc", class_weight="balanced",
    )
    model.fit(Xtr_s, ytr)
    prob = model.predict_proba(Xte_s)[:, 1]
    pred = model.predict(Xte_s)
    coef = np.asarray(model.coef_).ravel()
    nf = int(np.sum(np.abs(coef) > 1e-10))
    return {
        "roc_auc": float(roc_auc_score(yte, prob)),
        "accuracy": float(accuracy_score(yte, pred)),
        "n_features": nf,
    }


def benchmark_dataset(name: str, task: str):
    X, y = _load_dataset(name)
    if X is None:
        print(f"  [{name}] CACHED DESCRIPTORS NOT FOUND; skipping")
        return []
    n = len(y)
    print(f"\n--- {name}  ({task})  n={n} ---")

    if task == "classification":
        y_int = y.astype(int)
        splitter = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=RNG_SEED)
        splits = list(splitter.split(X.values, y_int))
        targets = y_int
    else:
        splitter = KFold(n_splits=N_FOLDS, shuffle=True, random_state=RNG_SEED)
        splits = list(splitter.split(X.values))
        targets = y

    per_cfg = {c: [] for c in CONFIGS}
    n_kept_per_cfg = {c: [] for c in CONFIGS}

    for fold_idx, (tr, te) in enumerate(splits):
        X_train_df = X.iloc[tr].reset_index(drop=True)
        X_test_df = X.iloc[te].reset_index(drop=True)
        y_train = targets[tr]; y_test = targets[te]
        for cfg in CONFIGS:
            kept = _select_features(X_train_df, y_train, cfg)
            n_kept_per_cfg[cfg].append(len(kept))
            Xtr = X_train_df[kept].values if kept else np.empty((len(y_train), 0))
            Xte = X_test_df[kept].values if kept else np.empty((len(y_test), 0))
            if task == "regression":
                m = _fit_eval_regression(Xtr, y_train, Xte, y_test)
            else:
                m = _fit_eval_classification(Xtr, y_train, Xte, y_test)
            m["n_kept_pre_lasso"] = len(kept)
            per_cfg[cfg].append(m)

    rows = []
    primary_metric = "rmse" if task == "regression" else "roc_auc"
    # Per-fold scores used for paired t-tests
    per_fold_primary = {cfg: np.array([v[primary_metric] for v in per_cfg[cfg]],
                                      dtype=float) for cfg in CONFIGS}

    for cfg in CONFIGS:
        vals = per_cfg[cfg]
        rec = {"dataset": name, "task": task, "config": cfg, "n": n,
               "model": "LASSO_CV" if task == "regression" else "L1_Logistic_CV"}
        rec["n_kept_pre_lasso_mean"] = round(float(np.mean(n_kept_per_cfg[cfg])), 2)
        rec["n_kept_pre_lasso_std"]  = round(float(np.std(n_kept_per_cfg[cfg])), 2)
        for metric in (["rmse", "mae", "r2"] if task == "regression" else ["roc_auc", "accuracy"]):
            arr = np.array([v[metric] for v in vals], dtype=float)
            rec[f"{metric}_mean"] = round(float(np.nanmean(arr)), 4)
            rec[f"{metric}_std"]  = round(float(np.nanstd(arr)), 4)
        nfa = np.array([v["n_features"] for v in vals], dtype=float)
        rec["n_features_mean"] = round(float(nfa.mean()), 2)
        rec["n_features_std"]  = round(float(nfa.std()), 2)
        # Paired t-test vs full: H0 = no difference in primary metric
        # across folds. Two-sided, df = N_FOLDS - 1.
        if cfg == "full":
            rec["paired_t_vs_full"] = ""; rec["paired_p_vs_full"] = ""
        else:
            diffs = per_fold_primary[cfg] - per_fold_primary["full"]
            # Use only finite paired diffs (combined_pruned on Lipo can have NaN)
            mask = np.isfinite(diffs)
            if mask.sum() >= 2:
                t_stat, p_val = _scipy_stats.ttest_rel(
                    per_fold_primary[cfg][mask],
                    per_fold_primary["full"][mask],
                )
                rec["paired_t_vs_full"] = round(float(t_stat), 3)
                rec["paired_p_vs_full"] = round(float(p_val), 4)
            else:
                rec["paired_t_vs_full"] = "NaN"
                rec["paired_p_vs_full"] = "NaN"
        rows.append(rec)

    # console summary
    if task == "regression":
        print(f"  {'config':>17s}  {'pre-LASSO':>10s}  {'final':>7s}  {'RMSE (mean±std)':>17s}  {'t vs full':>10s}  {'p vs full':>10s}")
        for r in rows:
            print(f"  {r['config']:>17s}  {r['n_kept_pre_lasso_mean']:>5.1f}      "
                  f"{r['n_features_mean']:>5.1f}    "
                  f"{r['rmse_mean']:>10.4f}±{r['rmse_std']:.4f}  "
                  f"{str(r['paired_t_vs_full']):>10s}  {str(r['paired_p_vs_full']):>10s}")
    else:
        print(f"  {'config':>17s}  {'pre-L1':>7s}  {'final':>7s}  {'AUC (mean±std)':>17s}  {'t vs full':>10s}  {'p vs full':>10s}")
        for r in rows:
            print(f"  {r['config']:>17s}  {r['n_kept_pre_lasso_mean']:>5.1f}    "
                  f"{r['n_features_mean']:>5.1f}    "
                  f"{r['roc_auc_mean']:>10.4f}±{r['roc_auc_std']:.4f}  "
                  f"{str(r['paired_t_vs_full']):>10s}  {str(r['paired_p_vs_full']):>10s}")
    return rows


def main():
    all_rows = []
    for ds, _ in REGRESSION_DATASETS:
        all_rows.extend(benchmark_dataset(ds, "regression"))
    for ds, _ in CLASSIFICATION_DATASETS:
        all_rows.extend(benchmark_dataset(ds, "classification"))

    out_path = os.path.join(PROJECT, "results", "paper2_pipeline_then_lasso.csv")
    if all_rows:
        # union of all keys, preserving order
        keys = []
        for r in all_rows:
            for k in r.keys():
                if k not in keys:
                    keys.append(k)
        with open(out_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            for r in all_rows:
                w.writerow(r)
        print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
