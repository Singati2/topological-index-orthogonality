"""§6.4 Degeneracy of Wuzi on trees of order 10 — direct analog of MATCH
95:141-162 Section 5.3 (Figure 5).

Degeneracy of a topological index TI on a class of graphs G_1, ..., G_N
measures how many distinct values TI attains. We report:

    %degeneracy(TI) := (1 - #distinct values / N) × 100

so that 0% means every graph has a unique value, 100% means all graphs
have the same value. The MATCH paper reports ~17–22% on trees of size 10
for the Sombor-family indices. We compute the analog for several Wuzi
parameter triples and for the 30-index baseline, on the same setting.

Output: results/wuzi_degeneracy.csv and results/wuzi_degeneracy.md
"""
from __future__ import annotations
import os
import sys
import numpy as np
import pandas as pd
import networkx as nx

THIS_DIR = os.path.abspath(os.path.dirname(__file__))
PROJECT = os.path.abspath(os.path.join(THIS_DIR, ".."))
sys.path.insert(0, PROJECT)

from src.standard_indices import ALL_INDICES, compute_all
from src.wuzi_index import wuzi

TREE_ORDER = 10
ROUND_DECIMALS = 6


def degeneracy_pct(values: list[float]) -> float:
    rounded = [round(v, ROUND_DECIMALS) for v in values]
    return 100.0 * (1.0 - len(set(rounded)) / len(rounded))


def main():
    print(f"Enumerating non-isomorphic trees of order {TREE_ORDER} ...")
    trees = list(nx.nonisomorphic_trees(TREE_ORDER))
    print(f"Found {len(trees)} non-isomorphic trees.")

    # Baseline: all 30 standard indices
    print("Computing 30 baseline indices on each tree ...")
    baseline_rows = [compute_all(T) for T in trees]
    base_df = pd.DataFrame(baseline_rows)

    rows = []
    for name in base_df.columns:
        rows.append({
            "index": name,
            "family": "baseline",
            "n_distinct": len(set(round(v, ROUND_DECIMALS) for v in base_df[name])),
            "degeneracy_pct": degeneracy_pct(list(base_df[name].values)),
        })

    # Wuzi at canonical & empirically-relevant parameter points
    wuzi_params = [
        (-0.5, 0.0, 2.0),    # ESOL "closest to passing"
        (-1.0, -1.0, 1.0),   # FreeSolv "passing"
        (1.0, 0.0, 0.0),     # = M_2 (sanity check; should match baseline M2)
        (-0.5, 0.0, 0.0),    # = Randic R (should match baseline R)
        (0.5, 0.3, 1.0),     # generic mid-range
        (0.0, 0.0, 1.0),     # pure irregularity exponential
        (0.0, 0.0, 2.0),     # stronger irregularity exponential
    ]
    print(f"Computing Wuzi at {len(wuzi_params)} parameter triples ...")
    for (a, b, g) in wuzi_params:
        vals = [wuzi(T, a, b, g) for T in trees]
        rows.append({
            "index": f"Wuzi(α={a}, β={b}, γ={g})",
            "family": "wuzi",
            "n_distinct": len(set(round(v, ROUND_DECIMALS) for v in vals)),
            "degeneracy_pct": degeneracy_pct(vals),
        })

    df = pd.DataFrame(rows).sort_values("degeneracy_pct").reset_index(drop=True)
    out_csv = os.path.join(PROJECT, "results", "wuzi_degeneracy.csv")
    df.to_csv(out_csv, index=False)

    # Markdown table
    md = [f"# §6.4 Degeneracy on trees of order {TREE_ORDER}\n"]
    md.append(f"Among {len(trees)} non-isomorphic trees of order {TREE_ORDER}, ")
    md.append("we report the percentage of trees that share a value with at least ")
    md.append("one other tree, for each index. Lower = better discrimination.\n")
    md.append("**Wuzi entries are highlighted; baseline entries are for comparison.**\n")
    md.append("| Index | Family | # distinct values | % degeneracy |")
    md.append("|---|---|---:|---:|")
    for _, r in df.iterrows():
        emph = "**" if r["family"] == "wuzi" else ""
        md.append(f"| {emph}{r['index']}{emph} | {r['family']} | "
                  f"{r['n_distinct']} | {r['degeneracy_pct']:.2f}% |")
    md.append(f"\n*Reference: MATCH 95:141-162 Figure 5 reports the analogous "
              "experiment for the Sombor-family indices, observing degeneracies "
              "in the 17–22% range. Values closer to 0% indicate better "
              "discrimination among non-isomorphic structures.*")
    md_path = os.path.join(PROJECT, "results", "wuzi_degeneracy.md")
    with open(md_path, "w") as f:
        f.write("\n".join(md))

    print(f"\n=== Top 10 (least degenerate) ===")
    print(df.head(10).to_string(index=False))
    print(f"\n=== Bottom 10 (most degenerate) ===")
    print(df.tail(10).to_string(index=False))
    print(f"\nWritten {out_csv} and {md_path}")


if __name__ == "__main__":
    main()
