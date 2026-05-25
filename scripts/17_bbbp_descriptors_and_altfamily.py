"""Generate BBBP descriptors and run alt-family screen on BBBP.

Addresses two Paper-2 v1.1 limitations in one pass:
  (i)  no cached BBBP `descriptors_baseline.csv` (so script 15's
       LASSO/ElasticNet comparator could not be extended to BBBP);
  (ii) alt-family multi-dataset replication (script 16) excluded
       BBBP because the binary target needed a different treatment.

Treatment for binary BBBP:
  - "Raw correlation with target" is point-biserial Pearson r
    between candidate and the 0/1 indicator (this is just standard
    Pearson r on a 0/1 vector, equal to point-biserial by
    construction).
  - "Partial correlation" is computed by the same residualisation
    recipe as in script 16: fit a Linear Probability Model
    y_hat = X_baseline @ beta + b on the binary target, take
    residuals y - y_hat, and Pearson-correlate the candidate's
    residual against the target's residual. This is the standard
    LPM-style partial-correlation diagnostic. We use it (rather
    than logistic-regression-based pseudo-R^2 / Tjur's coefficient)
    so that the BBBP row of Table tab:novel_candidates_multidataset
    is directly comparable to the three regression rows.

Outputs:
  results/bbbp/descriptors_baseline.csv      (30 baselines + 'target')
  results/novel_candidates_experiment/novel_candidates_bbbp.csv
  Plus: this script merges the new BBBP row into
  novel_candidates_multidataset_summary.csv after the merge.
"""
from __future__ import annotations
import os
import sys
import time
import warnings

import numpy as np
import pandas as pd

THIS_DIR = os.path.abspath(os.path.dirname(__file__))
PROJECT = os.path.abspath(os.path.join(THIS_DIR, ".."))
sys.path.insert(0, PROJECT)

# Silence rdkit/numpy ABI warnings - the functions still produce correct output.
warnings.filterwarnings("ignore")

from sklearn.linear_model import LinearRegression  # noqa: E402

from src.load_data import load                      # noqa: E402
from src.mol_to_graph import smiles_to_graph        # noqa: E402
from src.standard_indices import compute_all        # noqa: E402
from src.novel_candidates import CANDIDATE_INDICES  # noqa: E402


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

PAIRWISE_THRESHOLD = 0.95
PCOR_THRESHOLD = 0.10


def main():
    print("--- bbbp ---")
    df = load("bbbp")
    graphs, y_vals = [], []
    for _, row in df.iterrows():
        G = smiles_to_graph(row["smiles"])
        if G is None:
            continue
        graphs.append(G)
        y_vals.append(row["target"])
    n = len(graphs)
    y = np.array(y_vals, dtype=float)
    print(f"  Built {n} graphs (target = BBB penetration, binary 0/1).")

    print(f"  Computing 30 baseline indices ...")
    baseline_rows = [compute_all(G) for G in graphs]
    baseline = pd.DataFrame(baseline_rows)

    # Cache descriptors_baseline.csv so script 15 (LASSO) can be
    # extended to BBBP without recomputing.
    out_dir_b = os.path.join(PROJECT, "results", "bbbp")
    os.makedirs(out_dir_b, exist_ok=True)
    out_csv_b = os.path.join(out_dir_b, "descriptors_baseline.csv")
    baseline_with_target = baseline.copy()
    baseline_with_target["target"] = y
    baseline_with_target.to_csv(out_csv_b, index=False)
    print(f"  Wrote {out_csv_b} ({len(baseline_with_target)} rows × {baseline_with_target.shape[1]} cols).")

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
        if np.any(np.isnan(vals)):
            mean = np.nanmean(vals)
            vals = np.where(np.isnan(vals), mean, vals)
        dt = time.time() - t0
        if np.std(vals) < 1e-12:
            rows.append({
                "dataset": "bbbp", "candidate": name, "family": FAMILY_OF[name],
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
        abs_corrs = {col: abs(float(np.corrcoef(vals, baseline[col].values)[0, 1]))
                     for col in baseline.columns if np.std(baseline[col].values) > 1e-12}
        best_baseline = max(abs_corrs, key=abs_corrs.get)
        max_r = abs_corrs[best_baseline]
        raw_corr = float(np.corrcoef(vals, y)[0, 1])  # point-biserial on 0/1
        v_res = vals - LinearRegression().fit(Xb, vals).predict(Xb)
        partial = 0.0 if np.std(v_res) < 1e-12 else float(np.corrcoef(v_res, y_res)[0, 1])
        pairwise_pass = max_r < PAIRWISE_THRESHOLD
        combined_pass = pairwise_pass and abs(partial) >= PCOR_THRESHOLD
        rows.append({
            "dataset": "bbbp", "candidate": name, "family": FAMILY_OF[name],
            "max_abs_r_baseline": round(max_r, 4),
            "most_correlated_baseline": best_baseline,
            "raw_corr_target": round(raw_corr, 4),  # point-biserial
            "partial_corr_target": round(partial, 4),  # LPM partial r
            "pairwise_verdict": "PASS" if pairwise_pass else "FAIL",
            "combined_verdict": "PASS" if combined_pass else "FAIL",
            "n_compute_fail": n_fail,
            "compute_sec": round(dt, 2),
        })

    out_dir = os.path.join(PROJECT, "results", "novel_candidates_experiment")
    os.makedirs(out_dir, exist_ok=True)
    detail_csv = os.path.join(out_dir, "novel_candidates_bbbp.csv")
    pd.DataFrame(rows).to_csv(detail_csv, index=False)
    print(f"  Wrote {detail_csv}")

    # Merge into the multi-dataset CSV/summary so Paper 2 reads from one file.
    full_csv = os.path.join(out_dir, "novel_candidates_multidataset.csv")
    if os.path.exists(full_csv):
        full = pd.read_csv(full_csv)
        full = full[full["dataset"] != "bbbp"]  # idempotent
        full = pd.concat([full, pd.DataFrame(rows)], ignore_index=True)
        full.to_csv(full_csv, index=False)
        print(f"  Merged BBBP rows into {full_csv}")

        # Re-emit summary
        summary_rows = []
        for ds in full["dataset"].unique():
            sub = full[full["dataset"] == ds]
            summary_rows.append({
                "dataset": ds, "n_candidates": len(sub),
                "pairwise_pass": int((sub["pairwise_verdict"] == "PASS").sum()),
                "combined_pass": int((sub["combined_verdict"] == "PASS").sum()),
                "degenerate":    int((sub["pairwise_verdict"] == "DEGENERATE").sum()),
                "max_abs_partial_corr": round(float(sub["partial_corr_target"].abs().max(skipna=True)), 4),
            })
        summary_csv = os.path.join(out_dir, "novel_candidates_multidataset_summary.csv")
        pd.DataFrame(summary_rows).to_csv(summary_csv, index=False)
        print(f"  Re-emitted {summary_csv}")
        print()
        print("Summary (all 4 datasets):")
        print(pd.DataFrame(summary_rows).to_string(index=False))


if __name__ == "__main__":
    main()
