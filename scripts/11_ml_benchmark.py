"""ML benchmark: does orthogonality screening change downstream model
performance?

For each dataset (ESOL, FreeSolv, Lipophilicity, BBBP), trains a
RandomForest under four feature configurations across 5-fold CV.
Reports mean +/- std across folds.

Feature configurations
----------------------
  full              all 30 baseline indices
  pca_95            top-k principal components capturing >=95% variance,
                    fit on the training fold only
  pairwise_pruned   greedy pruning at |r| >= 0.95 on the training fold
  combined_pruned   pairwise-pruned PLUS partial-correlation filter,
                    fit on the training fold only

All feature selection is performed inside each CV fold to avoid
test-fold information leakage into the feature-selection step.

For regression targets we report RMSE, MAE, R^2.
For classification targets we report ROC-AUC, accuracy, F1.

The headline claim of the methodology paper is operationalized as:
``combined_pruned'' / ``pairwise_pruned'' should match `full' within
cross-validation noise on the held-out folds, while using fewer features.

Outputs:
  results/ml_benchmark.csv
  results/ml_benchmark.md

This script does not output figures; numerical tables are sufficient
for the methodology paper's claim. Figure generation happens in
scripts/08_make_figures.py::fig9_ml_benchmark from this CSV.
"""
from __future__ import annotations
import math
import os
import sys
import time
from typing import Any

import numpy as np
import pandas as pd

THIS_DIR = os.path.abspath(os.path.dirname(__file__))
PROJECT = os.path.abspath(os.path.join(THIS_DIR, ".."))
sys.path.insert(0, PROJECT)

from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.linear_model import LinearRegression
from sklearn.metrics import (
    accuracy_score, f1_score, mean_absolute_error,
    mean_squared_error, r2_score, roc_auc_score,
)
from sklearn.model_selection import KFold, StratifiedKFold
from sklearn.preprocessing import StandardScaler

from src.load_data import DATASETS, load
from src.mol_to_graph import smiles_to_graph
from src.standard_indices import compute_all

DATASET_NAMES = ["esol", "freesolv", "lipophilicity", "bbbp"]
RNG_SEED = 42
N_FOLDS = 5
PAIRWISE_THRESHOLD = 0.95
PCOR_THRESHOLD = 0.10
PCA_VARIANCE = 0.95
CONFIGS = ["full", "pca_95", "pairwise_pruned", "combined_pruned"]

REGRESSION_METRICS = ("rmse", "mae", "r2")
CLASSIFICATION_METRICS = ("roc_auc", "accuracy", "f1")


def compute_baseline_descriptors(df: pd.DataFrame) -> tuple[pd.DataFrame, np.ndarray]:
    """Compute the 30-index baseline matrix and the target vector for `df`."""
    rows = []
    keep_target = []
    for _, row in df.iterrows():
        G = smiles_to_graph(row["smiles"])
        if G is None:
            continue
        rows.append(compute_all(G))
        keep_target.append(row["target"])
    X = pd.DataFrame(rows)
    y = np.asarray(keep_target, dtype=float)
    return X, y


def select_pairwise_pruned(X: pd.DataFrame, threshold: float) -> list[str]:
    """Greedy redundancy-aware feature selection at |r| < threshold.

    Iterates columns in DataFrame order; keeps each column whose maximum
    absolute correlation with the already-kept columns is below threshold.
    """
    kept: list[str] = []
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


def select_combined_pruned(
    X: pd.DataFrame, y: np.ndarray,
    pair_threshold: float, pcor_threshold: float,
) -> list[str]:
    """Pairwise-pruned set further filtered by |partial corr with y|.

    Both filters are computed on the supplied (X, y) -- which is the
    training fold's slice, not the full dataset.
    """
    pair_kept = select_pairwise_pruned(X, pair_threshold)
    if len(pair_kept) <= 1:
        return pair_kept
    final_kept: list[str] = []
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


def fit_pca_train_only(
    X_train: np.ndarray, X_test: np.ndarray, variance: float,
) -> tuple[np.ndarray, np.ndarray, int]:
    """Fit StandardScaler and PCA on X_train only, transform X_test.

    Returns (X_train_pc, X_test_pc, k) where k is the smallest number of
    PCs whose cumulative explained variance >= `variance` on the training
    fold.
    """
    scaler = StandardScaler().fit(X_train)
    Xs_train = scaler.transform(X_train)
    Xs_test = scaler.transform(X_test)
    pca = PCA().fit(Xs_train)
    cum = np.cumsum(pca.explained_variance_ratio_)
    k = int(np.searchsorted(cum, variance)) + 1
    return pca.transform(Xs_train)[:, :k], pca.transform(Xs_test)[:, :k], k


def nan_metrics_for(task: str) -> dict[str, float]:
    nan = float("nan")
    return {"rmse": nan, "mae": nan, "r2": nan,
            "roc_auc": nan, "accuracy": nan, "f1": nan}


def evaluate(
    X_train: np.ndarray, y_train: np.ndarray,
    X_test: np.ndarray, y_test: np.ndarray,
    task_type: str,
) -> dict[str, float]:
    """Train RandomForest on train, score on test, return task-appropriate metrics."""
    if X_train.shape[1] == 0:
        # Empty feature set; cannot fit a model.
        return nan_metrics_for(task_type)
    if task_type == "regression":
        model = RandomForestRegressor(
            n_estimators=300, random_state=RNG_SEED, n_jobs=-1,
        )
        model.fit(X_train, y_train)
        preds = model.predict(X_test)
        rmse = float(np.sqrt(mean_squared_error(y_test, preds)))
        mae = float(mean_absolute_error(y_test, preds))
        r2 = float(r2_score(y_test, preds))
        return {"rmse": rmse, "mae": mae, "r2": r2,
                "roc_auc": float("nan"), "accuracy": float("nan"), "f1": float("nan")}
    if task_type == "classification":
        model = RandomForestClassifier(
            n_estimators=300, random_state=RNG_SEED, n_jobs=-1,
            class_weight="balanced",
        )
        y_int = y_train.astype(int)
        model.fit(X_train, y_int)
        proba = model.predict_proba(X_test)[:, 1]
        preds = (proba >= 0.5).astype(int)
        y_true = y_test.astype(int)
        return {
            "rmse": float("nan"), "mae": float("nan"), "r2": float("nan"),
            "roc_auc": float(roc_auc_score(y_true, proba)),
            "accuracy": float(accuracy_score(y_true, preds)),
            "f1": float(f1_score(y_true, preds)),
        }
    raise ValueError(f"unknown task_type {task_type!r}")


def aggregate(values: list[float]) -> tuple[float, float]:
    """Return (mean, std) ignoring NaNs; if all NaN, return (NaN, NaN)."""
    arr = np.asarray(values, dtype=float)
    if np.all(np.isnan(arr)):
        return float("nan"), float("nan")
    return float(np.nanmean(arr)), float(np.nanstd(arr, ddof=0))


def benchmark_dataset(dataset: str) -> list[dict[str, Any]]:
    spec = DATASETS[dataset]
    task = spec["task_type"]
    print(f"\n{'='*64}\nML benchmark: {dataset}  ({task})  {N_FOLDS}-fold CV\n{'='*64}")
    df = load(dataset)
    t0 = time.time()
    X, y = compute_baseline_descriptors(df)
    print(f"  Built {len(X)} graphs in {time.time()-t0:.1f}s")
    if task == "classification":
        keep = ~np.isnan(y)
        X = X.loc[keep].reset_index(drop=True)
        y = y[keep]

    if task == "classification":
        splitter = StratifiedKFold(
            n_splits=N_FOLDS, shuffle=True, random_state=RNG_SEED,
        )
        splits = list(splitter.split(X, y.astype(int)))
    else:
        splitter = KFold(
            n_splits=N_FOLDS, shuffle=True, random_state=RNG_SEED,
        )
        splits = list(splitter.split(X))

    per_config: dict[str, dict[str, list]] = {
        cfg: {"n_features": [],
              "rmse": [], "mae": [], "r2": [],
              "roc_auc": [], "accuracy": [], "f1": []}
        for cfg in CONFIGS
    }

    for fold_idx, (train_idx, test_idx) in enumerate(splits):
        X_train_df = X.iloc[train_idx]
        X_test_df = X.iloc[test_idx]
        y_train = y[train_idx]
        y_test = y[test_idx]

        # 1) full
        m = evaluate(X_train_df.values, y_train,
                     X_test_df.values, y_test, task)
        for metric_name, v in m.items():
            per_config["full"][metric_name].append(v)
        per_config["full"]["n_features"].append(X_train_df.shape[1])

        # 2) pca_95 -- fit on train fold only
        Xpc_train, Xpc_test, k = fit_pca_train_only(
            X_train_df.values, X_test_df.values, PCA_VARIANCE,
        )
        m = evaluate(Xpc_train, y_train, Xpc_test, y_test, task)
        for metric_name, v in m.items():
            per_config["pca_95"][metric_name].append(v)
        per_config["pca_95"]["n_features"].append(k)

        # 3) pairwise_pruned -- fit on train fold only
        pair_kept = select_pairwise_pruned(X_train_df, PAIRWISE_THRESHOLD)
        if pair_kept:
            m = evaluate(
                X_train_df[pair_kept].values, y_train,
                X_test_df[pair_kept].values, y_test, task,
            )
        else:
            m = nan_metrics_for(task)
        for metric_name, v in m.items():
            per_config["pairwise_pruned"][metric_name].append(v)
        per_config["pairwise_pruned"]["n_features"].append(len(pair_kept))

        # 4) combined_pruned -- fit on train fold only
        comb_kept = select_combined_pruned(
            X_train_df, y_train, PAIRWISE_THRESHOLD, PCOR_THRESHOLD,
        )
        if comb_kept:
            m = evaluate(
                X_train_df[comb_kept].values, y_train,
                X_test_df[comb_kept].values, y_test, task,
            )
        else:
            # Combined screen kept nothing on this fold; report NaN metrics
            # and 0 features. The aggregator preserves the NaN.
            m = nan_metrics_for(task)
        for metric_name, v in m.items():
            per_config["combined_pruned"][metric_name].append(v)
        per_config["combined_pruned"]["n_features"].append(len(comb_kept))

    rows: list[dict[str, Any]] = []
    for cfg in CONFIGS:
        nf = per_config[cfg]["n_features"]
        nf_mean, nf_std = aggregate(nf)
        nf_min, nf_max = int(np.min(nf)), int(np.max(nf))
        row = {
            "dataset": dataset,
            "task_type": task,
            "config": cfg,
            "n_features_mean": nf_mean,
            "n_features_std":  nf_std,
            "n_features_min":  nf_min,
            "n_features_max":  nf_max,
        }
        for metric_name in REGRESSION_METRICS + CLASSIFICATION_METRICS:
            mean_v, std_v = aggregate(per_config[cfg][metric_name])
            row[f"{metric_name}_mean"] = mean_v
            row[f"{metric_name}_std"] = std_v
        rows.append(row)

        # Pretty print
        if task == "regression":
            print(f"  {cfg:18s}  n_feat={nf_mean:.1f} (range {nf_min}-{nf_max})  "
                  f"RMSE={row['rmse_mean']:.3f}+-{row['rmse_std']:.3f}  "
                  f"R^2={row['r2_mean']:.3f}+-{row['r2_std']:.3f}")
        else:
            print(f"  {cfg:18s}  n_feat={nf_mean:.1f} (range {nf_min}-{nf_max})  "
                  f"ROC-AUC={row['roc_auc_mean']:.3f}+-{row['roc_auc_std']:.3f}  "
                  f"F1={row['f1_mean']:.3f}+-{row['f1_std']:.3f}")
    return rows


def fmt_mean_std(mean: float, std: float, places: int = 3) -> str:
    if math.isnan(mean):
        return "---"
    return f"{mean:.{places}f} +/- {std:.{places}f}"


def fmt_n_features(mean: float, mn: int, mx: int) -> str:
    if mn == mx:
        return f"{mn}"
    return f"{mean:.1f} ({mn}-{mx})"


def main() -> None:
    rows: list[dict[str, Any]] = []
    for name in DATASET_NAMES:
        rows.extend(benchmark_dataset(name))
    df = pd.DataFrame(rows)
    out_csv = os.path.join(PROJECT, "results", "ml_benchmark.csv")
    df.to_csv(out_csv, index=False)

    md_lines: list[str] = []
    md_lines.append("# ML benchmark: does orthogonality screening preserve predictive performance?\n")
    md_lines.append(
        f"Generated by `scripts/11_ml_benchmark.py`. RandomForest "
        f"(300 trees, fixed seed) trained under {N_FOLDS}-fold "
        f"cross-validation (stratified for the classification "
        f"dataset). All feature selection (PCA, pairwise pruning, "
        f"combined pruning) is performed inside each training fold "
        f"to avoid test-fold information leakage. Reported values "
        f"are mean +/- std across the {N_FOLDS} folds.\n")
    md_lines.append(
        f"Feature configurations:\n\n"
        f"- **full** -- all 30 baseline indices\n"
        f"- **pca_95** -- top-k principal components on the training "
        f"fold, cumulative variance >= {PCA_VARIANCE:.0%}\n"
        f"- **pairwise_pruned** -- greedy redundancy filter at "
        f"|r| < {PAIRWISE_THRESHOLD} on the training fold\n"
        f"- **combined_pruned** -- pairwise filter plus "
        f"|partial corr with target| >= {PCOR_THRESHOLD} on the "
        f"training fold\n")
    md_lines.append("## Regression datasets (RMSE / MAE / R^2; lower / lower / higher is better)\n")
    md_lines.append("| Dataset | Config | n_features | RMSE | MAE | R^2 |")
    md_lines.append("|---|---|---|---|---|---|")
    for r in rows:
        if r["task_type"] != "regression":
            continue
        md_lines.append(
            f"| {r['dataset']} | {r['config']} | "
            f"{fmt_n_features(r['n_features_mean'], r['n_features_min'], r['n_features_max'])} | "
            f"{fmt_mean_std(r['rmse_mean'], r['rmse_std'])} | "
            f"{fmt_mean_std(r['mae_mean'], r['mae_std'])} | "
            f"{fmt_mean_std(r['r2_mean'], r['r2_std'])} |"
        )
    md_lines.append("\n## Classification dataset (ROC-AUC / accuracy / F1; higher is better)\n")
    md_lines.append("| Dataset | Config | n_features | ROC-AUC | Accuracy | F1 |")
    md_lines.append("|---|---|---|---|---|---|")
    for r in rows:
        if r["task_type"] != "classification":
            continue
        md_lines.append(
            f"| {r['dataset']} | {r['config']} | "
            f"{fmt_n_features(r['n_features_mean'], r['n_features_min'], r['n_features_max'])} | "
            f"{fmt_mean_std(r['roc_auc_mean'], r['roc_auc_std'])} | "
            f"{fmt_mean_std(r['accuracy_mean'], r['accuracy_std'])} | "
            f"{fmt_mean_std(r['f1_mean'], r['f1_std'])} |"
        )
    md_lines.append("\n## Interpretation\n")
    md_lines.append(
        "Compare the `full` and `pairwise_pruned` rows within each dataset.\n"
        "If the mean +/- std intervals overlap, the screen has reduced the\n"
        "feature count without statistically meaningful loss in predictive\n"
        "performance. The combined-criterion screen is more aggressive; on\n"
        "datasets where the underlying topology-target partial-correlation\n"
        "structure is weak, the screen may retain zero features in some or\n"
        "all folds (reported as '---').\n")
    out_md = os.path.join(PROJECT, "results", "ml_benchmark.md")
    with open(out_md, "w") as f:
        f.write("\n".join(md_lines))

    print()
    print(df.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    print()
    print(f"Written {out_csv}")
    print(f"Written {out_md}")


if __name__ == "__main__":
    main()
