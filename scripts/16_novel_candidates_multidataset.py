"""Cross-dataset replication of the alternative-family candidate screen.

Generalises scripts/04_novel_candidates_test.py (ESOL-only) to run on
ESOL, FreeSolv, and Lipophilicity. BBBP is deferred because the binary
classification target makes the partial-correlation diagnostic ambiguous.

For each (dataset, candidate index) pair:
  - max |r| with any of the 30 baselines  (pairwise verdict: < 0.95 = PASS)
  - raw Pearson r with the regression target
  - partial correlation with target after residualising the candidate
    on the 30 baselines (combined-diagnostic verdict: pairwise PASS AND
    |partial| >= 0.10)

Output:
  results/novel_candidates_experiment/novel_candidates_multidataset.csv
  results/novel_candidates_experiment/novel_candidates_multidataset_summary.csv

Honest scope: candidate-index implementations live in
src/novel_candidates.py and are dataset-agnostic in principle, but
some candidates may fail or be degenerate on larger / more diverse
chemistry (e.g. Lipophilicity has 4200 molecules covering broader
chemical space than the 1127 ESOL solubility set). Failures are
recorded explicitly rather than discarded.
"""
from __future__ import annotations
import csv
import os
import sys
import time

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression

THIS_DIR = os.path.abspath(os.path.dirname(__file__))
PROJECT = os.path.abspath(os.path.join(THIS_DIR, ".."))
sys.path.insert(0, PROJECT)

from src.mol_to_graph import smiles_to_graph
from src.standard_indices import compute_all
from src.novel_candidates import CANDIDATE_INDICES
from src.load_data import load


FAMILY_OF = {
    "InfoH_deg": "info-theoretic", "InfoH_dist": "info-theoretic",
    "InfoH_spec": "info-theoretic", "InfoH_ecc": "info-theoretic",
    "Bonchev_Id": "info-theoretic",
    "Btw_sum": "centrality", "Cls_sum": "centrality",
    "Eig_sum": "centrality", "Harm_cent": "centrality",
    "Btw_var": "centrality",
    "AlgConn": "spectral-hybrid", "Rho_x_Dia": "spectral-hybrid",
    "EE_x_W": "spectral-hybrid",
    "Triangles": "motif", "MeanClust": "motif",
    "FourCyc": "motif", "CutVerts": "motif",
    "TotEcc": "eccentricity", "Radius": "eccentricity",
    "MeanEcc": "eccentricity",
}


DATASETS = ["esol", "freesolv", "lipophilicity"]

PAIRWISE_THRESHOLD = 0.95
PCOR_THRESHOLD = 0.10


def screen_dataset(dataset_name):
    print(f"\n--- {dataset_name} ---")
    df = load(dataset_name)
    target_col = "target"
    graphs, y_vals = [], []
    for _, row in df.iterrows():
        G = smiles_to_graph(row["smiles"])
        if G is None:
            continue
        graphs.append(G)
        y_vals.append(row[target_col])
    n = len(graphs)
    y = np.array(y_vals, dtype=float)
    print(f"  Built {n} graphs (target: {target_col}).")

    print(f"  Computing 30 baseline indices ...")
    baseline_rows = [compute_all(G) for G in graphs]
    baseline = pd.DataFrame(baseline_rows)
    Xb = baseline.values
    y_res = y - LinearRegression().fit(Xb, y).predict(Xb)

    print(f"  Computing {len(CANDIDATE_INDICES)} candidate indices ...")
    rows = []
    for name, fn in CANDIDATE_INDICES.items():
        t0 = time.time()
        vals = []
        n_fail = 0
        for G in graphs:
            try:
                vals.append(fn(G))
            except Exception:
                vals.append(float("nan"))
                n_fail += 1
        vals = np.array(vals, dtype=float)
        # Impute NaN with column mean
        if np.any(np.isnan(vals)):
            mean = np.nanmean(vals)
            vals = np.where(np.isnan(vals), mean, vals)
        dt = time.time() - t0
        if np.std(vals) < 1e-12:
            rows.append({
                "dataset": dataset_name, "candidate": name, "family": FAMILY_OF[name],
                "max_abs_r_baseline": float("nan"),
                "most_correlated_baseline": "(constant)",
                "raw_corr_target": float("nan"),
                "partial_corr_target": float("nan"),
                "pairwise_verdict": "DEGENERATE",
                "combined_verdict": "DEGENERATE",
                "n_compute_fail": n_fail,
                "compute_sec": round(dt, 2),
            })
            continue
        # max |r| with baselines
        abs_corrs = {col: abs(float(np.corrcoef(vals, baseline[col].values)[0, 1])) for col in baseline.columns if np.std(baseline[col].values) > 1e-12}
        best_baseline = max(abs_corrs, key=abs_corrs.get)
        max_r = abs_corrs[best_baseline]
        raw_corr = float(np.corrcoef(vals, y)[0, 1])
        v_res = vals - LinearRegression().fit(Xb, vals).predict(Xb)
        partial = 0.0 if np.std(v_res) < 1e-12 else float(np.corrcoef(v_res, y_res)[0, 1])
        pairwise_pass = max_r < PAIRWISE_THRESHOLD
        combined_pass = pairwise_pass and abs(partial) >= PCOR_THRESHOLD
        rows.append({
            "dataset": dataset_name, "candidate": name, "family": FAMILY_OF[name],
            "max_abs_r_baseline": round(max_r, 4),
            "most_correlated_baseline": best_baseline,
            "raw_corr_target": round(raw_corr, 4),
            "partial_corr_target": round(partial, 4),
            "pairwise_verdict": "PASS" if pairwise_pass else "FAIL",
            "combined_verdict": "PASS" if combined_pass else "FAIL",
            "n_compute_fail": n_fail,
            "compute_sec": round(dt, 2),
        })
    return rows


def main():
    all_rows = []
    for name in DATASETS:
        try:
            all_rows.extend(screen_dataset(name))
        except Exception as e:
            print(f"  [{name}] FAILED: {e!r}")

    out_dir = os.path.join(PROJECT, "results", "novel_candidates_experiment")
    os.makedirs(out_dir, exist_ok=True)
    df = pd.DataFrame(all_rows)
    detail_csv = os.path.join(out_dir, "novel_candidates_multidataset.csv")
    df.to_csv(detail_csv, index=False)
    print(f"\nWrote {detail_csv}")

    # Summary table: per dataset, pairwise pass count, combined pass count, max |pcor|
    summary_rows = []
    for ds in df["dataset"].unique():
        sub = df[df["dataset"] == ds]
        n = len(sub)
        n_pair = int((sub["pairwise_verdict"] == "PASS").sum())
        n_comb = int((sub["combined_verdict"] == "PASS").sum())
        n_degen = int((sub["pairwise_verdict"] == "DEGENERATE").sum())
        max_pcor = float(sub["partial_corr_target"].abs().max(skipna=True))
        summary_rows.append({
            "dataset": ds, "n_candidates": n,
            "pairwise_pass": n_pair, "combined_pass": n_comb,
            "degenerate": n_degen,
            "max_abs_partial_corr": round(max_pcor, 4),
        })
    summary = pd.DataFrame(summary_rows)
    summary_csv = os.path.join(out_dir, "novel_candidates_multidataset_summary.csv")
    summary.to_csv(summary_csv, index=False)
    print(f"Wrote {summary_csv}")
    print()
    print("Summary:")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
