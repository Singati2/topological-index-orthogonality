"""End-to-end run: load ESOL, compute 30 baseline topological indices, write
correlation analysis to results/.

This script establishes the BASELINE redundancy structure of the existing
~30 indices on real drug-like chemistry. When Advik proposes a new index,
add it via the --new-index hook and rerun.

Outputs:
  results/descriptors_baseline.csv      — index values per molecule
  results/correlation_matrix.csv        — 30x30 Pearson correlations
  results/redundant_pairs.csv           — pairs with |r| >= 0.95 (the
                                           noisy zone for any "new" tweak)
  results/pca_variance.csv              — variance explained by PCA components
  results/vif.csv                       — VIF table
  results/partial_corr_logS.csv         — partial correlation with logS target
  results/summary.txt                   — readable summary
"""
from __future__ import annotations
import os
import sys
import time
import numpy as np
import pandas as pd

THIS_DIR = os.path.abspath(os.path.dirname(__file__))
PROJECT = os.path.abspath(os.path.join(THIS_DIR, ".."))
sys.path.insert(0, PROJECT)

from src.mol_to_graph import smiles_to_graph
from src.standard_indices import ALL_INDICES, compute_all
from src.load_data import load_esol
from src import orthogonality as ortho


RESULTS = os.path.join(PROJECT, "results")
os.makedirs(RESULTS, exist_ok=True)


def main():
    print("=" * 60)
    print("Baseline orthogonality analysis — 30 standard indices on ESOL")
    print("=" * 60)

    t0 = time.time()
    df = load_esol()
    print(f"\n[1/4] Loaded ESOL: n={len(df)} molecules")

    # Build graphs; skip unparseable
    print(f"\n[2/4] Building molecular graphs (SMILES -> NetworkX) ...")
    rows = []
    skipped = 0
    for i, row in df.iterrows():
        G = smiles_to_graph(row["smiles"])
        if G is None:
            skipped += 1
            continue
        rows.append({"smiles": row["smiles"], "logS": row["logS"], "G": G})
    print(f"      Built {len(rows)} graphs; skipped {skipped} (single-atom / invalid).")

    # Compute indices
    print(f"\n[3/4] Computing {len(ALL_INDICES)} indices per molecule ...")
    desc_rows = []
    for k, r in enumerate(rows):
        if (k + 1) % 200 == 0:
            print(f"      {k+1}/{len(rows)} ...")
        try:
            d = compute_all(r["G"])
            d["__logS"] = r["logS"]
            desc_rows.append(d)
        except Exception as e:
            print(f"      WARN: failed on row {k}: {e}")
    desc = pd.DataFrame(desc_rows)
    y = desc.pop("__logS")
    print(f"      Descriptor matrix: {desc.shape[0]} rows x {desc.shape[1]} indices")

    # Save raw descriptors
    out = desc.copy()
    out["logS"] = y.values
    out.to_csv(os.path.join(RESULTS, "descriptors_baseline.csv"), index=False)

    # Correlation + redundancy
    print(f"\n[4/4] Running orthogonality analysis ...")
    corr = ortho.correlation_matrix(desc)
    corr.to_csv(os.path.join(RESULTS, "correlation_matrix.csv"))
    pairs = ortho.redundancy_report(corr, threshold=0.95)
    pairs.to_csv(os.path.join(RESULTS, "redundant_pairs.csv"), index=False)
    pca = ortho.pca_variance_explained(desc)
    pca.to_csv(os.path.join(RESULTS, "pca_variance.csv"), index=False)
    vif = ortho.vif(desc)
    vif.to_csv(os.path.join(RESULTS, "vif.csv"), index=False)
    pcy = ortho.partial_corr_with_target(desc, y)
    pcy.to_csv(os.path.join(RESULTS, "partial_corr_logS.csv"), index=False)

    # Summary
    summary = []
    summary.append("BASELINE ORTHOGONALITY REPORT — 30 standard topological indices on ESOL")
    summary.append(f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    summary.append(f"Molecules: {len(desc)}")
    summary.append(f"Indices:   {desc.shape[1]}")
    summary.append("")
    summary.append("--- Redundant pairs (|r| >= 0.95) ---")
    summary.append(f"Count: {len(pairs)}")
    summary.append(pairs.head(20).to_string(index=False))
    summary.append("")
    summary.append("--- PCA: how many components needed for 95% / 99% variance? ---")
    pc95 = int((pca["cumulative"] >= 0.95).idxmax()) + 1
    pc99 = int((pca["cumulative"] >= 0.99).idxmax()) + 1
    summary.append(f"95% variance: {pc95} components (of {desc.shape[1]})")
    summary.append(f"99% variance: {pc99} components (of {desc.shape[1]})")
    summary.append(f"Effective rank of the 30-index block is roughly {pc95} — that's how much")
    summary.append(f"orthogonal information the current zoo actually carries on ESOL.")
    summary.append("")
    summary.append("--- VIF top 10 (highest collinearity) ---")
    summary.append(vif.head(10).to_string(index=False))
    summary.append("")
    summary.append("--- Partial correlation with logS (top 10 |partial r|) ---")
    summary.append(pcy.head(10).to_string(index=False))
    summary.append("")
    summary.append("--- KILL-TEST INSTRUCTIONS for a new candidate index ---")
    summary.append("1. Implement your candidate index as a function f(G) -> float")
    summary.append("   in src/new_index.py with name e.g. 'NEW'.")
    summary.append("2. Add it to the registry and rerun this script.")
    summary.append("3. Inspect correlation_matrix.csv row 'NEW':")
    summary.append("     if max |r| with any baseline >= 0.95 -> redundant, redesign.")
    summary.append("     if max |r| < 0.95 -> proceed to multi-dataset CV phase.")
    summary.append("4. Inspect partial_corr_logS.csv rank of 'NEW':")
    summary.append("     a high rank means it adds information for prediction even")
    summary.append("     controlling for the others — strongest possible case.")

    txt = "\n".join(summary)
    with open(os.path.join(RESULTS, "summary.txt"), "w") as f:
        f.write(txt + "\n")
    print()
    print(txt)
    print(f"\nDone in {time.time()-t0:.1f}s. Outputs in {RESULTS}/")


if __name__ == "__main__":
    main()
