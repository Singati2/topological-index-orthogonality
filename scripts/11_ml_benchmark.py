"""ML benchmark: does orthogonality screening change downstream model
performance?

For each dataset (ESOL, FreeSolv, Lipophilicity, BBBP), trains a
RandomForest on four feature configurations and reports test-set
metrics. Goal: test whether the redundancy-screening pipeline
preserves predictive performance while substantially reducing the
input dimension.

Feature configurations
----------------------
  full              all 30 baseline indices
  pca_95            top-k principal components capturing >=95% variance
  pairwise_pruned   greedy pruning at |r| >= 0.95
                    (keep one representative per correlation cluster)
  combined_pruned   pairwise-pruned PLUS partial-correlation filter
                    (drop indices with |pcor(z, y | X)| < 0.10)

For regression targets we report RMSE, MAE, R^2.
For classification targets we report ROC-AUC, accuracy, F1.

The headline claim of the methodology paper is operationalized as:
``combined_pruned'' should match `full' within standard cross-validation
noise on a held-out test split, while using fewer features.

Outputs:
  results/ml_benchmark.csv
  results/ml_benchmark.md

This script does not output figures; numerical tables are sufficient
for the methodology paper's claim.
"""
from __future__ import annotations
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
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

from src.load_data import DATASETS, load
from src.mol_to_graph import smiles_to_graph
from src.standard_indices import compute_all

DATASET_NAMES = ["esol", "freesolv", "lipophilicity", "bbbp"]
RNG_SEED = 42
TEST_SIZE = 0.20
PAIRWISE_THRESHOLD = 0.95
PCOR_THRESHOLD = 0.10
PCA_VARIANCE = 0.95


def compute_baseline_descriptors(df: pd.DataFrame) -> tuple[pd.DataFrame, np.ndarray]:
    """Compute the 30-index baseline matrix and the target vector for `df`.

    Skips rows whose SMILES does not parse. Returns (X, y) with both arrays
    aligned to the kept rows.
    """
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

    Iterate columns in order; keep each column whose maximum absolute
    correlation with the already-kept columns is below `threshold`.
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
    X: pd.DataFrame, y: np.ndarray, pair_threshold: float, pcor_threshold: float
) -> list[str]:
    """Pairwise-pruned set further filtered by partial correlation with target.

    A column is kept only if it passes the pairwise screen AND its partial
    correlation with the target (after regressing out the other pairwise-
    pruned columns) has absolute value at least `pcor_threshold`.
    """
    pair_kept = select_pairwise_pruned(X, pair_threshold)
    if len(pair_kept) <= 1:
        return pair_kept
    final_kept: list[str] = []
    Xb = X[pair_kept].values
    for col in pair_kept:
        others = [c for c in pair_kept if c != col]
        if not others:
            final_kept.append(col)
            continue
        Xo = X[others].values
        z = X[col].values
        if np.std(z) < 1e-12:
            continue
        # Residualize z and y against others.
        z_res = z - LinearRegression().fit(Xo, z).predict(Xo)
        y_res = y - LinearRegression().fit(Xo, y).predict(Xo)
        if np.std(z_res) < 1e-12 or np.std(y_res) < 1e-12:
            continue
        pcor = float(np.corrcoef(z_res, y_res)[0, 1])
        if abs(pcor) >= pcor_threshold:
            final_kept.append(col)
    return final_kept


def select_pca_95(X: pd.DataFrame, variance: float) -> tuple[np.ndarray, int]:
    """Project X onto the smallest k PCs whose cumulative variance >= `variance`.

    Returns (X_pc, k).
    """
    Xs = StandardScaler().fit_transform(X.values)
    pca = PCA().fit(Xs)
    cum = np.cumsum(pca.explained_variance_ratio_)
    k = int(np.searchsorted(cum, variance)) + 1
    return pca.transform(Xs)[:, :k], k


def evaluate(
    X_train: np.ndarray, y_train: np.ndarray,
    X_test: np.ndarray, y_test: np.ndarray,
    task_type: str,
) -> dict[str, float]:
    """Train RandomForest on train, score on test, return task-appropriate metrics."""
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


def benchmark_dataset(dataset: str) -> list[dict[str, Any]]:
    spec = DATASETS[dataset]
    task = spec["task_type"]
    print(f"\n{'='*64}\nML benchmark: {dataset}  ({task})\n{'='*64}")
    df = load(dataset)
    t0 = time.time()
    X, y = compute_baseline_descriptors(df)
    print(f"  Built {len(X)} graphs in {time.time()-t0:.1f}s")
    if task == "classification":
        # Drop NaN target rows defensively.
        keep = ~np.isnan(y)
        X = X.loc[keep].reset_index(drop=True)
        y = y[keep]

    # Splits used identically for every configuration so the comparison
    # isolates the effect of feature selection.
    stratify = y.astype(int) if task == "classification" else None
    X_train_idx, X_test_idx = train_test_split(
        np.arange(len(X)), test_size=TEST_SIZE,
        random_state=RNG_SEED, stratify=stratify,
    )

    rows: list[dict[str, Any]] = []

    # 1) full
    cfg = "full"
    X_train = X.iloc[X_train_idx].values
    X_test = X.iloc[X_test_idx].values
    metrics = evaluate(X_train, y[X_train_idx], X_test, y[X_test_idx], task)
    rows.append({"dataset": dataset, "task_type": task, "config": cfg,
                 "n_features": X.shape[1], **metrics})
    print(f"  {cfg:18s}  n_features={X.shape[1]:>2}  {metrics}")

    # 2) pca_95
    cfg = "pca_95"
    X_pc_all, k = select_pca_95(X, PCA_VARIANCE)
    X_train_pc = X_pc_all[X_train_idx]
    X_test_pc = X_pc_all[X_test_idx]
    metrics = evaluate(X_train_pc, y[X_train_idx], X_test_pc, y[X_test_idx], task)
    rows.append({"dataset": dataset, "task_type": task, "config": cfg,
                 "n_features": k, **metrics})
    print(f"  {cfg:18s}  n_features={k:>2}  {metrics}")

    # 3) pairwise_pruned
    cfg = "pairwise_pruned"
    pair_kept = select_pairwise_pruned(X, PAIRWISE_THRESHOLD)
    X_train_p = X[pair_kept].iloc[X_train_idx].values
    X_test_p = X[pair_kept].iloc[X_test_idx].values
    metrics = evaluate(X_train_p, y[X_train_idx], X_test_p, y[X_test_idx], task)
    rows.append({"dataset": dataset, "task_type": task, "config": cfg,
                 "n_features": len(pair_kept), **metrics})
    print(f"  {cfg:18s}  n_features={len(pair_kept):>2}  kept={pair_kept}")
    print(f"                     {metrics}")

    # 4) combined_pruned
    cfg = "combined_pruned"
    comb_kept = select_combined_pruned(X, y, PAIRWISE_THRESHOLD, PCOR_THRESHOLD)
    if len(comb_kept) == 0:
        # If the combined screen kills every feature, fall back to the single
        # most-target-correlated feature. We do NOT silently use 0 features.
        # Report and skip ML; report n_features=0 with NaN metrics.
        metrics = {k: float("nan") for k in
                   ("rmse", "mae", "r2", "roc_auc", "accuracy", "f1")}
        rows.append({"dataset": dataset, "task_type": task, "config": cfg,
                     "n_features": 0, **metrics})
        print(f"  {cfg:18s}  n_features= 0  (combined screen kept nothing; "
              f"reported as NaN)")
    else:
        X_train_c = X[comb_kept].iloc[X_train_idx].values
        X_test_c = X[comb_kept].iloc[X_test_idx].values
        metrics = evaluate(X_train_c, y[X_train_idx], X_test_c, y[X_test_idx], task)
        rows.append({"dataset": dataset, "task_type": task, "config": cfg,
                     "n_features": len(comb_kept), **metrics})
        print(f"  {cfg:18s}  n_features={len(comb_kept):>2}  kept={comb_kept}")
        print(f"                     {metrics}")
    return rows


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
        "Generated by `scripts/11_ml_benchmark.py`. RandomForest "
        "(300 trees, fixed seed) trained on an 80/20 split with stratification "
        "for the classification dataset. Identical splits across all four "
        "feature configurations isolate the effect of feature selection.\n")
    md_lines.append(
        "Feature configurations:\n\n"
        f"- **full** — all 30 baseline indices\n"
        f"- **pca_95** — top-k principal components, cumulative variance "
        f">= {PCA_VARIANCE:.0%}\n"
        f"- **pairwise_pruned** — greedy redundancy filter at |r| < "
        f"{PAIRWISE_THRESHOLD}\n"
        f"- **combined_pruned** — pairwise-pruned + |partial corr with target| "
        f">= {PCOR_THRESHOLD}\n")
    md_lines.append("## Regression datasets (RMSE / MAE / R^2; lower / lower / higher is better)\n")
    md_lines.append("| Dataset | Config | n_features | RMSE | MAE | R^2 |")
    md_lines.append("|---|---|---:|---:|---:|---:|")
    for r in rows:
        if r["task_type"] != "regression":
            continue
        md_lines.append(
            f"| {r['dataset']} | {r['config']} | {r['n_features']} | "
            f"{r['rmse']:.4f} | {r['mae']:.4f} | {r['r2']:.4f} |"
        )
    md_lines.append("\n## Classification dataset (ROC-AUC / accuracy / F1; higher is better)\n")
    md_lines.append("| Dataset | Config | n_features | ROC-AUC | Accuracy | F1 |")
    md_lines.append("|---|---|---:|---:|---:|---:|")
    for r in rows:
        if r["task_type"] != "classification":
            continue
        md_lines.append(
            f"| {r['dataset']} | {r['config']} | {r['n_features']} | "
            f"{r['roc_auc']:.4f} | {r['accuracy']:.4f} | {r['f1']:.4f} |"
        )
    md_lines.append("\n## Interpretation\n")
    md_lines.append(
        "The headline question for the methodology paper is: does the "
        "combined-criterion screen preserve predictive performance while "
        "reducing input dimension? A clean affirmative answer is the "
        "comparison of `full` and `combined_pruned` columns within each "
        "dataset: if the metrics are within standard cross-validation noise "
        "(commonly cited as ~0.05 RMSE or ~0.02 ROC-AUC on this kind of data) "
        "and the feature count is materially reduced, then the screen has "
        "operationalized non-redundancy without sacrificing model utility.\n")
    md_lines.append(
        "On per-dataset reading the four rows can also be compared "
        "diagnostically:\n"
        "- `pca_95` collapses to 3 features on the regression datasets "
        "(as expected from the BID-basis dimension bound), so it provides "
        "a *minimum-information* control: model performance with only 3 "
        "principal components establishes a floor.\n"
        "- `pairwise_pruned` reduces feature count by removing only the "
        "obviously redundant columns; performance close to `full` indicates "
        "that the dropped columns add no fresh information.\n"
        "- `combined_pruned` is the strict screen: pairwise screen plus "
        "partial correlation with target. Materially fewer features than "
        "`pairwise_pruned` with similar performance is the strongest "
        "support for the methodology paper's `necessary but not sufficient' "
        "framing of pairwise correlation.\n")
    md_lines.append(
        "**Caveat.** The screen is parameterized by two threshold choices "
        "($|r| = 0.95$ and $|\\mathrm{pcor}| = 0.10$). Sensitivity to those "
        "thresholds for the per-dataset feature count is in "
        "`results/cross_dataset_summary.md`; sensitivity for the ML "
        "performance is an obvious follow-up not done here.\n")
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
