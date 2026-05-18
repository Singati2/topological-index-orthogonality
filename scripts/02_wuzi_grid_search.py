"""Wuzi parameter sweep over (α, β, γ) under the orthogonality screening
pipeline on a chosen MoleculeNet dataset.

Usage:
    python scripts/02_wuzi_grid_search.py --dataset esol
    python scripts/02_wuzi_grid_search.py --dataset freesolv
    python scripts/02_wuzi_grid_search.py --dataset lipophilicity

For each grid point (α, β, γ):
  - max |r| with any baseline index ("kill-test verdict")
  - raw Pearson correlation with the target property
  - partial correlation with target controlling for the 30 baselines
"""
from __future__ import annotations
import argparse
import os
import sys
import time
import numpy as np
import pandas as pd
from itertools import product
from sklearn.linear_model import LinearRegression

THIS_DIR = os.path.abspath(os.path.dirname(__file__))
PROJECT = os.path.abspath(os.path.join(THIS_DIR, ".."))
sys.path.insert(0, PROJECT)

from src.mol_to_graph import smiles_to_graph
from src.standard_indices import compute_all
from src.wuzi_index import wuzi
from src.load_data import load, DATASETS


ALPHA_GRID = [-1.0, -0.5, 0.0, 0.5, 1.0]
BETA_GRID  = [-1.0, -0.5, 0.0, 0.5, 1.0]
GAMMA_GRID = [0.0, 0.5, 1.0, 2.0]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--dataset", required=True, choices=list(DATASETS))
    return p.parse_args()


def main():
    args = parse_args()
    dataset = args.dataset
    out_dir = os.path.join(PROJECT, "results", dataset)
    os.makedirs(out_dir, exist_ok=True)

    t0 = time.time()
    print("=" * 64)
    print(f"Wuzi parameter sweep — dataset = {dataset}")
    print("=" * 64)

    df = load(dataset)
    target_name = df.attrs["target_name"]
    print(f"\n[1/4] Loaded {dataset}: n={len(df)}  target={target_name}")

    print(f"\n[2/4] Building graphs ...")
    graphs, y_vals = [], []
    for _, row in df.iterrows():
        G = smiles_to_graph(row["smiles"])
        if G is None:
            continue
        graphs.append(G)
        y_vals.append(row["target"])
    print(f"      Built {len(graphs)} graphs.")

    print(f"\n[3/4] Computing 30 baseline indices ...")
    baseline = pd.DataFrame([compute_all(G) for G in graphs])
    Xb = baseline.values
    y = np.array(y_vals)
    y_res = y - LinearRegression().fit(Xb, y).predict(Xb)

    total_pts = len(ALPHA_GRID) * len(BETA_GRID) * len(GAMMA_GRID)
    print(f"\n[4/4] Sweeping {total_pts} grid points ...")
    rows = []
    for i, (alpha, beta, gamma) in enumerate(product(ALPHA_GRID, BETA_GRID, GAMMA_GRID)):
        w = np.array([wuzi(G, alpha, beta, gamma) for G in graphs])
        if np.std(w) < 1e-12:
            rows.append({"alpha": alpha, "beta": beta, "gamma": gamma,
                         "max_abs_r_baseline": float("nan"),
                         "most_correlated_baseline": "(constant)",
                         "raw_corr_target": float("nan"),
                         "partial_corr_target": float("nan"),
                         "kill_test": "DEGENERATE"})
            continue
        abs_corrs = {col: abs(float(np.corrcoef(w, baseline[col].values)[0, 1]))
                     for col in baseline.columns
                     if np.std(baseline[col].values) > 1e-12}
        best = max(abs_corrs, key=abs_corrs.get)
        max_r = abs_corrs[best]
        raw_corr = float(np.corrcoef(w, y)[0, 1])
        w_res = w - LinearRegression().fit(Xb, w).predict(Xb)
        partial = 0.0 if np.std(w_res) < 1e-12 else float(np.corrcoef(w_res, y_res)[0, 1])
        rows.append({"alpha": alpha, "beta": beta, "gamma": gamma,
                     "max_abs_r_baseline": max_r,
                     "most_correlated_baseline": best,
                     "raw_corr_target": raw_corr,
                     "partial_corr_target": partial,
                     "kill_test": "PASS" if max_r < 0.95 else "FAIL"})
        if (i + 1) % 25 == 0:
            print(f"      {i+1}/{total_pts}")
    grid = pd.DataFrame(rows)
    grid.to_csv(os.path.join(out_dir, "wuzi_grid.csv"), index=False)

    n_pass = int((grid["kill_test"] == "PASS").sum())
    n_fail = int((grid["kill_test"] == "FAIL").sum())
    n_degen = int((grid["kill_test"] == "DEGENERATE").sum())

    summary = []
    summary.append(f"WUZI PARAMETER SWEEP — dataset = {dataset}  target = {target_name}")
    summary.append(f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    summary.append(f"Grid: alpha={ALPHA_GRID}, beta={BETA_GRID}, gamma={GAMMA_GRID}")
    summary.append(f"Total points: {total_pts}; molecules: {len(graphs)}")
    summary.append("")
    summary.append(f"Kill-test: PASS={n_pass}, FAIL={n_fail}, DEGENERATE={n_degen}")
    summary.append("")
    summary.append("--- 10 best points by |partial corr with target| ---")
    cols = ["alpha","beta","gamma","max_abs_r_baseline","most_correlated_baseline",
            "raw_corr_target","partial_corr_target","kill_test"]
    summary.append(grid.assign(_a=grid["partial_corr_target"].abs())
                       .sort_values("_a", ascending=False)
                       .head(10)[cols].to_string(index=False))
    summary.append("")
    summary.append("--- 10 failing points closest to passing (smallest max|r|) ---")
    fails = grid[grid["kill_test"] == "FAIL"].sort_values("max_abs_r_baseline").head(10)
    summary.append(fails[cols].to_string(index=False))

    txt = "\n".join(summary)
    with open(os.path.join(out_dir, "wuzi_summary.txt"), "w") as f:
        f.write(txt + "\n")
    print()
    print(txt)
    print(f"\nDone in {time.time()-t0:.1f}s. Outputs in {out_dir}/")


if __name__ == "__main__":
    main()
