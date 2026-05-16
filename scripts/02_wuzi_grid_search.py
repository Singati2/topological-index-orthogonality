"""Sweep the Wuzi parameter space (α, β, γ) and characterize each point by:

  - max |r| with the 30 baseline indices on ESOL (the kill-test verdict)
  - raw Pearson correlation with logS (predictive standalone signal)
  - partial correlation with logS controlling for the 30 baselines
    (the strongest claim — orthogonal predictive information)

Output:
  results/wuzi_grid.csv  — one row per (α, β, γ) with the three metrics above
  results/wuzi_passing.csv  — subset with max|r| < 0.95 (orthogonality-passing)
  results/wuzi_summary.txt   — readable narrative of where the family lives
"""
from __future__ import annotations
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
from src.load_data import load_esol


# Grid: 5 x 5 x 4 = 100 points.
ALPHA_GRID = [-1.0, -0.5, 0.0, 0.5, 1.0]
BETA_GRID  = [-1.0, -0.5, 0.0, 0.5, 1.0]
GAMMA_GRID = [0.0, 0.5, 1.0, 2.0]


def main():
    t0 = time.time()
    print("=" * 60)
    print("Wuzi parameter sweep on ESOL")
    print("=" * 60)

    print("\n[1/4] Loading ESOL and building graphs ...")
    df = load_esol()
    graphs, logS_vals = [], []
    for _, row in df.iterrows():
        G = smiles_to_graph(row["smiles"])
        if G is None:
            continue
        graphs.append(G)
        logS_vals.append(row["logS"])
    print(f"      Built {len(graphs)} graphs.")

    print("\n[2/4] Computing 30 baseline indices (for kill-test comparison) ...")
    baseline_rows = [compute_all(G) for G in graphs]
    baseline = pd.DataFrame(baseline_rows)
    y = np.array(logS_vals)
    # Pre-residualize y on baselines once -- used in every partial-corr below.
    y_res = y - LinearRegression().fit(baseline.values, y).predict(baseline.values)

    print("\n[3/4] Sweeping {} parameter combinations ...".format(
        len(ALPHA_GRID) * len(BETA_GRID) * len(GAMMA_GRID)))
    rows = []
    Xb = baseline.values
    for i, (alpha, beta, gamma) in enumerate(product(ALPHA_GRID, BETA_GRID, GAMMA_GRID)):
        w = np.array([wuzi(G, alpha, beta, gamma) for G in graphs])
        if np.std(w) < 1e-12:
            # Degenerate (e.g., constant value across molecules — happens for
            # α=β=γ=0 on graphs of equal edge count, very rare)
            max_r = float("nan")
            best_baseline = "(constant)"
            raw_corr = float("nan")
            partial = float("nan")
        else:
            # Correlations with each baseline
            corrs = {}
            for col in baseline.columns:
                b = baseline[col].values
                if np.std(b) < 1e-12:
                    corrs[col] = 0.0
                else:
                    corrs[col] = float(np.corrcoef(w, b)[0, 1])
            abs_corrs = {k: abs(v) for k, v in corrs.items()}
            best_baseline = max(abs_corrs, key=abs_corrs.get)
            max_r = abs_corrs[best_baseline]
            # Raw corr with target
            raw_corr = float(np.corrcoef(w, y)[0, 1])
            # Partial corr with target controlling for baselines
            w_res = w - LinearRegression().fit(Xb, w).predict(Xb)
            if np.std(w_res) < 1e-12:
                partial = 0.0
            else:
                partial = float(np.corrcoef(w_res, y_res)[0, 1])

        rows.append({
            "alpha": alpha,
            "beta": beta,
            "gamma": gamma,
            "max_abs_r_baseline": max_r,
            "most_correlated_baseline": best_baseline,
            "raw_corr_logS": raw_corr,
            "partial_corr_logS": partial,
            "kill_test": "PASS" if max_r < 0.95 else "FAIL",
        })
        if (i + 1) % 20 == 0:
            print(f"      {i+1}/{len(ALPHA_GRID)*len(BETA_GRID)*len(GAMMA_GRID)}")
    grid = pd.DataFrame(rows)

    print("\n[4/4] Writing results ...")
    grid.to_csv(os.path.join(PROJECT, "results", "wuzi_grid.csv"), index=False)
    passing = grid[grid["kill_test"] == "PASS"].copy()
    passing = passing.sort_values("partial_corr_logS", key=lambda s: s.abs(), ascending=False)
    passing.to_csv(os.path.join(PROJECT, "results", "wuzi_passing.csv"), index=False)

    # Summary narrative
    summary = []
    summary.append("WUZI PARAMETER SWEEP — RESULTS")
    summary.append(f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    summary.append(f"Grid: α∈{ALPHA_GRID}, β∈{BETA_GRID}, γ∈{GAMMA_GRID}")
    summary.append(f"Total grid points: {len(grid)}")
    summary.append(f"Molecules: {len(graphs)}  (ESOL)")
    summary.append("")
    pass_count = (grid["kill_test"] == "PASS").sum()
    fail_count = (grid["kill_test"] == "FAIL").sum()
    summary.append(f"--- Kill-test verdict (max |r| with any baseline < 0.95) ---")
    summary.append(f"PASS: {pass_count} / {len(grid)} ({100*pass_count/len(grid):.0f}%)")
    summary.append(f"FAIL: {fail_count} / {len(grid)} ({100*fail_count/len(grid):.0f}%)")
    summary.append("")
    summary.append("--- Best parameter point by |partial correlation with logS| ---")
    if len(passing) > 0:
        b = passing.iloc[0]
        summary.append(
            f"(α={b['alpha']}, β={b['beta']}, γ={b['gamma']}):"
            f"  max|r|={b['max_abs_r_baseline']:.3f} vs {b['most_correlated_baseline']},"
            f"  raw corr w/ logS={b['raw_corr_logS']:.3f},"
            f"  partial corr={b['partial_corr_logS']:.3f}"
        )
    else:
        summary.append("No parameter point passes the kill-test.")
    summary.append("")
    summary.append("--- Top 10 passing parameter points (by |partial corr|) ---")
    if len(passing) > 0:
        cols = ["alpha","beta","gamma","max_abs_r_baseline","most_correlated_baseline","raw_corr_logS","partial_corr_logS"]
        summary.append(passing[cols].head(10).to_string(index=False))
    summary.append("")
    summary.append("--- 10 failing points with smallest max|r| (closest to passing) ---")
    failing = grid[grid["kill_test"] == "FAIL"].sort_values("max_abs_r_baseline").head(10)
    if len(failing) > 0:
        cols = ["alpha","beta","gamma","max_abs_r_baseline","most_correlated_baseline","raw_corr_logS","partial_corr_logS"]
        summary.append(failing[cols].to_string(index=False))

    txt = "\n".join(summary)
    with open(os.path.join(PROJECT, "results", "wuzi_summary.txt"), "w") as f:
        f.write(txt + "\n")
    print(txt)
    print(f"\nDone in {time.time()-t0:.1f}s. Outputs in {PROJECT}/results/")


if __name__ == "__main__":
    main()
