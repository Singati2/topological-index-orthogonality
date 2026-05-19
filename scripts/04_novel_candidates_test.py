"""Apply the 30-baseline orthogonality screen on ESOL to 20 candidate
indices drawn from 5 alternative design families (information-theoretic,
centrality, spectral-hybrid, motif, eccentricity).

For each candidate compute:
  - max |r| with any baseline (PASS if < 0.95)
  - raw Pearson correlation with logS
  - partial correlation with logS controlling for the 30 baselines

The experiment is purely empirical: it quantifies how much independent
QSPR signal these alternative families add on ESOL relative to the
classical 30-index baseline. It is *not* a claim that the candidate
indices themselves are new — every one of them has prior chemistry
literature (see docs/literature_notes.md).
"""
from __future__ import annotations
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
from src.load_data import load_esol


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


def main():
    t0 = time.time()
    print("=" * 62)
    print("Orthogonality screen on ESOL: do alternative design families")
    print("provide independent signal beyond the 30-index baseline?")
    print("=" * 62)

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

    print("\n[2/4] Computing 30 baseline indices ...")
    baseline_rows = [compute_all(G) for G in graphs]
    baseline = pd.DataFrame(baseline_rows)
    Xb = baseline.values
    y = np.array(logS_vals)
    y_res = y - LinearRegression().fit(Xb, y).predict(Xb)

    print(f"\n[3/4] Computing {len(CANDIDATE_INDICES)} candidate indices ...")
    cand_data = {}
    for name, fn in CANDIDATE_INDICES.items():
        vals = []
        for G in graphs:
            try:
                vals.append(fn(G))
            except Exception as e:
                vals.append(float("nan"))
        cand_data[name] = np.array(vals)
        # Replace NaN with column mean (rare)
        v = cand_data[name]
        if np.any(np.isnan(v)):
            mean = np.nanmean(v)
            v[np.isnan(v)] = mean
            cand_data[name] = v
    print("      Done.")

    print("\n[4/4] Evaluating redundancy screen for each candidate ...")
    rows = []
    for name, vals in cand_data.items():
        if np.std(vals) < 1e-12:
            # Constant on this dataset (rare; e.g., triangle_count = 0 on all)
            rows.append({
                "candidate": name,
                "family": FAMILY_OF[name],
                "max_abs_r_baseline": float("nan"),
                "most_correlated_baseline": "(constant)",
                "raw_corr_logS": float("nan"),
                "partial_corr_logS": float("nan"),
                "screen_verdict": "DEGENERATE",
            })
            continue
        # max |r| with baselines
        corrs = {}
        for col in baseline.columns:
            b = baseline[col].values
            if np.std(b) < 1e-12:
                corrs[col] = 0.0
            else:
                corrs[col] = float(np.corrcoef(vals, b)[0, 1])
        abs_corrs = {k: abs(v) for k, v in corrs.items()}
        best_baseline = max(abs_corrs, key=abs_corrs.get)
        max_r = abs_corrs[best_baseline]
        raw_corr = float(np.corrcoef(vals, y)[0, 1])
        v_res = vals - LinearRegression().fit(Xb, vals).predict(Xb)
        if np.std(v_res) < 1e-12:
            partial = 0.0
        else:
            partial = float(np.corrcoef(v_res, y_res)[0, 1])
        rows.append({
            "candidate": name,
            "family": FAMILY_OF[name],
            "max_abs_r_baseline": max_r,
            "most_correlated_baseline": best_baseline,
            "raw_corr_logS": raw_corr,
            "partial_corr_logS": partial,
            "screen_verdict": "PASS" if max_r < 0.95 else "FAIL",
        })

    df_out = pd.DataFrame(rows).sort_values("partial_corr_logS", key=lambda s: s.abs(), ascending=False).reset_index(drop=True)
    out_dir = os.path.join(PROJECT, "results", "novel_candidates_experiment")
    os.makedirs(out_dir, exist_ok=True)
    df_out.to_csv(os.path.join(out_dir, "novel_candidates.csv"), index=False)

    # Summary
    n_pass = int((df_out["screen_verdict"] == "PASS").sum())
    n_fail = int((df_out["screen_verdict"] == "FAIL").sum())
    n_degen = int((df_out["screen_verdict"] == "DEGENERATE").sum())

    print()
    print("=" * 62)
    print("RESULTS")
    print("=" * 62)
    print(f"Tested: {len(df_out)} candidates across 5 design families")
    print(f"  PASS:       {n_pass} (carry information not in the 30-index baseline)")
    print(f"  FAIL:       {n_fail} (orthogonality ceiling, |r| >= 0.95)")
    print(f"  DEGENERATE: {n_degen} (constant or near-constant on ESOL)")
    print()
    cols = ["candidate", "family", "max_abs_r_baseline", "most_correlated_baseline",
            "raw_corr_logS", "partial_corr_logS", "screen_verdict"]
    print("Full ranking by |partial correlation with logS|:")
    print(df_out[cols].to_string(index=False))
    print()

    # Family-level summary
    print("Family-level summary (pass rate, best partial corr):")
    fam_rows = []
    for fam in df_out["family"].unique():
        sub = df_out[df_out["family"] == fam]
        passing = sub[sub["screen_verdict"] == "PASS"]
        fam_rows.append({
            "family": fam,
            "n_tested": len(sub),
            "n_pass": int((sub["screen_verdict"] == "PASS").sum()),
            "best_partial_corr": float(sub["partial_corr_logS"].abs().max()) if len(sub) else float("nan"),
            "min_max_r": float(sub["max_abs_r_baseline"].min(skipna=True)) if len(sub) else float("nan"),
        })
    fam_df = pd.DataFrame(fam_rows).sort_values("n_pass", ascending=False)
    print(fam_df.to_string(index=False))
    fam_df.to_csv(os.path.join(out_dir, "novel_candidates_by_family.csv"), index=False)

    # Summary
    print()
    print("=" * 62)
    print("SUMMARY")
    print("=" * 62)
    if n_pass >= 1:
        winners = df_out[df_out["screen_verdict"] == "PASS"].sort_values(
            "partial_corr_logS", key=lambda s: s.abs(), ascending=False)
        print(f"{n_pass} candidate(s) pass the |r| < 0.95 screen on ESOL.")
        print("Top candidates by |partial correlation with logS|:")
        print(winners[cols].head(5).to_string(index=False))
        print()
        print("Note: passing the screen is necessary but not sufficient for QSPR")
        print("usefulness. Inspect the partial correlation with logS to gauge")
        print("how much independent target-relevant signal each candidate carries.")
    else:
        print("0 of the 20 candidate indices pass the |r| < 0.95 screen on ESOL.")
        print("Every candidate is highly correlated with at least one classical")
        print("baseline index. This is consistent with the documented overlap")
        print("between alternative-family descriptors and the 30-index baseline")
        print("on drug-like chemistry; it is a single-dataset empirical finding,")
        print("not a general impossibility result.")

    txt_path = os.path.join(out_dir, "novel_candidates_summary.txt")
    with open(txt_path, "w") as f:
        f.write(f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Candidates tested: {len(df_out)} across 5 design families\n")
        f.write(f"PASS: {n_pass}, FAIL: {n_fail}, DEGENERATE: {n_degen}\n\n")
        f.write("Full ranking:\n")
        f.write(df_out[cols].to_string(index=False))
        f.write("\n\nFamily-level summary:\n")
        f.write(fam_df.to_string(index=False))

    print(f"\nDone in {time.time()-t0:.1f}s. Outputs in {PROJECT}/results/")


if __name__ == "__main__":
    main()
