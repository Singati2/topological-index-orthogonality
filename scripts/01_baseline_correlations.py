"""Baseline orthogonality screening for the 30 standard topological indices
on a chosen MoleculeNet dataset (ESOL / FreeSolv / Lipophilicity).

Usage:
    python scripts/01_baseline_correlations.py --dataset esol
    python scripts/01_baseline_correlations.py --dataset freesolv
    python scripts/01_baseline_correlations.py --dataset lipophilicity

Outputs land in results/<dataset>/. The same analysis is run on each dataset
so that cross-dataset redundancy patterns can be compared.
"""
from __future__ import annotations
import argparse
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
from src.load_data import load, DATASETS
from src import orthogonality as ortho


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--dataset", required=True, choices=list(DATASETS),
                   help="Which benchmark to run on.")
    return p.parse_args()


def main():
    args = parse_args()
    dataset = args.dataset
    out_dir = os.path.join(PROJECT, "results", dataset)
    os.makedirs(out_dir, exist_ok=True)

    print("=" * 64)
    print(f"Baseline orthogonality screening — 30 indices — dataset = {dataset}")
    print("=" * 64)

    t0 = time.time()
    df = load(dataset)
    target_name = df.attrs["target_name"]
    print(f"\n[1/4] Loaded {dataset}: n={len(df)}  target={target_name}")

    print(f"\n[2/4] Building molecular graphs ...")
    rows = []
    skipped = 0
    for i, row in df.iterrows():
        G = smiles_to_graph(row["smiles"])
        if G is None:
            skipped += 1
            continue
        rows.append({"smiles": row["smiles"], "target": row["target"], "G": G})
    print(f"      Built {len(rows)} graphs; skipped {skipped}.")

    print(f"\n[3/4] Computing {len(ALL_INDICES)} indices per molecule ...")
    desc_rows = []
    for k, r in enumerate(rows):
        if (k + 1) % 500 == 0:
            print(f"      {k+1}/{len(rows)} ...")
        try:
            d = compute_all(r["G"])
            d["__target"] = r["target"]
            desc_rows.append(d)
        except Exception as e:
            print(f"      WARN: failed on row {k}: {e}")
    desc = pd.DataFrame(desc_rows)
    y = desc.pop("__target")
    print(f"      Descriptor matrix: {desc.shape[0]} rows x {desc.shape[1]} indices")

    out = desc.copy()
    out[target_name] = y.values
    out.to_csv(os.path.join(out_dir, "descriptors_baseline.csv"), index=False)

    print(f"\n[4/4] Orthogonality analysis ...")
    corr = ortho.correlation_matrix(desc)
    corr.to_csv(os.path.join(out_dir, "correlation_matrix.csv"))
    pairs95 = ortho.redundancy_report(corr, threshold=0.95)
    pairs95.to_csv(os.path.join(out_dir, "redundant_pairs_r95.csv"), index=False)
    pairs90 = ortho.redundancy_report(corr, threshold=0.90)
    pairs90.to_csv(os.path.join(out_dir, "redundant_pairs_r90.csv"), index=False)
    pca = ortho.pca_variance_explained(desc)
    pca.to_csv(os.path.join(out_dir, "pca_variance.csv"), index=False)
    vif = ortho.vif(desc)
    vif.to_csv(os.path.join(out_dir, "vif.csv"), index=False)
    pcy = ortho.partial_corr_with_target(desc, y)
    pcy.to_csv(os.path.join(out_dir, f"partial_corr_{target_name}.csv"), index=False)

    summary = []
    summary.append(f"BASELINE ORTHOGONALITY REPORT — dataset = {dataset}")
    summary.append(f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    summary.append(f"Molecules: {len(desc)}  Target: {target_name}")
    summary.append(f"Indices:   {desc.shape[1]}")
    summary.append("")
    summary.append("--- Redundancy ---")
    summary.append(f"Pairs with |r| >= 0.90: {len(pairs90)} / {desc.shape[1]*(desc.shape[1]-1)//2}")
    summary.append(f"Pairs with |r| >= 0.95: {len(pairs95)} / {desc.shape[1]*(desc.shape[1]-1)//2}")
    summary.append("")
    summary.append("--- PCA effective rank ---")
    pc90 = int((pca["cumulative"] >= 0.90).idxmax()) + 1
    pc95 = int((pca["cumulative"] >= 0.95).idxmax()) + 1
    pc99 = int((pca["cumulative"] >= 0.99).idxmax()) + 1
    summary.append(f"Components for 90% variance: {pc90} of {desc.shape[1]}")
    summary.append(f"Components for 95% variance: {pc95} of {desc.shape[1]}")
    summary.append(f"Components for 99% variance: {pc99} of {desc.shape[1]}")
    summary.append("")
    summary.append("--- VIF (top 10 highest) ---")
    summary.append(vif.head(10).to_string(index=False))
    summary.append("")
    summary.append(f"--- Partial correlation with {target_name} (top 10) ---")
    summary.append(pcy.head(10).to_string(index=False))

    txt = "\n".join(summary)
    with open(os.path.join(out_dir, "summary.txt"), "w") as f:
        f.write(txt + "\n")
    print()
    print(txt)
    print(f"\nDone in {time.time()-t0:.1f}s. Outputs in {out_dir}/")


if __name__ == "__main__":
    main()
