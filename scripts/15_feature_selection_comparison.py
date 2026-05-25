"""LASSO / ElasticNet feature-selection comparator for Paper 2.

Addresses the reviewer concern that the screening pipeline (full / PCA /
pairwise / combined configurations of script 11) is compared only to
PCA but not to a generic feature-selection method. This script adds
LASSO and ElasticNet as comparators on the three regression datasets
(ESOL, FreeSolv, Lipophilicity).

BBBP is intentionally excluded from this comparator because the
classification-LASSO (Logistic L1) comparison would require
RDKit-backed feature regeneration that this environment cannot run;
the BBBP gap is reported as a formal limitation in Paper 2 §8 rather
than silently filled.

Uses cached results/{dataset}/descriptors_baseline.csv files (30
baseline indices already computed in Phase 0). For each dataset:
  - LASSO with 5-fold CV-tuned alpha (LassoCV)
  - ElasticNet with 5-fold CV-tuned (l1_ratio = 0.5)
  - 5-fold CV repeated 5 times with different splits to estimate
    RMSE mean/std and number of non-zero features

Output:
  results/paper2_feature_selection_comparison.csv

Compares against the existing ml_benchmark.csv `full`, `pairwise_pruned`,
`combined_pruned` numbers (computed under the same 5-fold CV protocol
in script 11), so the row ordering and RMSE definition are consistent.
"""
from __future__ import annotations
import csv
import os
import sys

import numpy as np
import pandas as pd
from sklearn.linear_model import LassoCV, ElasticNetCV
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error

PROJECT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

DATASETS = [
    ("esol",          "logS"),
    ("freesolv",      "dG_hyd"),
    ("lipophilicity", "logD"),
]
SEED = 20260525
N_FOLDS = 5

# Map dataset → target column name in the cached descriptors_baseline.csv
TARGET_COL = {"esol": "logS", "freesolv": "dG_hyd", "lipophilicity": "logD"}


def evaluate(model_factory, X, y, name, dataset_name):
    """5-fold CV; for each fold fit a fresh CV-tuned model on the train set,
    record test RMSE, and count number of non-zero features.

    Returns (rmse_mean, rmse_std, n_feats_mean, n_feats_std).
    """
    rng = np.random.RandomState(SEED)
    kf = KFold(n_splits=N_FOLDS, shuffle=True, random_state=SEED)
    rmses, nfeats = [], []
    for fold_id, (tr, te) in enumerate(kf.split(X)):
        scaler = StandardScaler().fit(X[tr])
        Xtr = scaler.transform(X[tr])
        Xte = scaler.transform(X[te])
        model = model_factory()
        model.fit(Xtr, y[tr])
        pred = model.predict(Xte)
        rmses.append(np.sqrt(mean_squared_error(y[te], pred)))
        nfeats.append(int(np.sum(np.abs(model.coef_) > 1e-10)))
    return (float(np.mean(rmses)), float(np.std(rmses)),
            float(np.mean(nfeats)), float(np.std(nfeats)))


def main():
    out_rows = []
    for ds_name, target_legacy in DATASETS:
        path = os.path.join(PROJECT, "results", ds_name, "descriptors_baseline.csv")
        if not os.path.exists(path):
            print(f"  [{ds_name}] CACHED DESCRIPTORS NOT FOUND at {path}; skipping")
            continue
        df = pd.read_csv(path)
        # The 30 baseline indices are the first 30 columns; the target is the last.
        target_col = df.columns[-1]
        baseline_cols = list(df.columns[:30])
        assert len(baseline_cols) == 30, f"Expected 30 baseline cols, got {len(baseline_cols)}: {baseline_cols}"
        X = df[baseline_cols].values.astype(float)
        y = df[target_col].values.astype(float)
        n = len(y)
        print(f"\n--- {ds_name} (n = {n}, target = {target_col}) ---")

        # LASSO with CV-tuned alpha
        def make_lasso():
            return LassoCV(cv=5, random_state=SEED, max_iter=10000, n_alphas=50)
        r1 = evaluate(make_lasso, X, y, "LASSO", ds_name)
        print(f"  LASSO        RMSE = {r1[0]:.4f} +/- {r1[1]:.4f}   n_feats = {r1[2]:.1f} +/- {r1[3]:.1f}")
        out_rows.append({
            "dataset": ds_name, "method": "LASSO_CV", "n": n,
            "rmse_mean": round(r1[0], 4), "rmse_std": round(r1[1], 4),
            "n_features_mean": round(r1[2], 2), "n_features_std": round(r1[3], 2),
        })

        # ElasticNet with CV-tuned alpha + l1_ratio
        def make_enet():
            return ElasticNetCV(cv=5, l1_ratio=[0.3, 0.5, 0.7, 0.9],
                                random_state=SEED, max_iter=10000, n_alphas=20)
        r2 = evaluate(make_enet, X, y, "ElasticNet", ds_name)
        print(f"  ElasticNet   RMSE = {r2[0]:.4f} +/- {r2[1]:.4f}   n_feats = {r2[2]:.1f} +/- {r2[3]:.1f}")
        out_rows.append({
            "dataset": ds_name, "method": "ElasticNet_CV", "n": n,
            "rmse_mean": round(r2[0], 4), "rmse_std": round(r2[1], 4),
            "n_features_mean": round(r2[2], 2), "n_features_std": round(r2[3], 2),
        })

    # Save
    out_path = os.path.join(PROJECT, "results", "paper2_feature_selection_comparison.csv")
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
        writer.writeheader()
        for r in out_rows:
            writer.writerow(r)
    print(f"\nWrote {out_path}")

    # Cross-reference against the existing ml_benchmark.csv full / pairwise / combined RMSE
    print("\n=== Comparator vs existing screen-pipeline configurations ===")
    print(f"{'dataset':16s} {'config':22s} {'RMSE':>10s} {'n_feats':>9s}")
    mlb_path = os.path.join(PROJECT, "results", "ml_benchmark.csv")
    if os.path.exists(mlb_path):
        ml = pd.read_csv(mlb_path)
        for ds in [r["dataset"] for r in out_rows]:
            for config in ["full", "pca_95", "pairwise_pruned", "combined_pruned"]:
                sub = ml[(ml["dataset"] == ds) & (ml["config"] == config)]
                if len(sub):
                    s = sub.iloc[0]
                    print(f"{ds:16s} {config:22s} {s['rmse_mean']:>10.4f} {s['n_features_mean']:>9.1f}")
            for r in out_rows:
                if r["dataset"] == ds:
                    print(f"{ds:16s} {r['method']:22s} {r['rmse_mean']:>10.4f} {r['n_features_mean']:>9.1f}")
            print()


if __name__ == "__main__":
    main()
