"""§6.5 Structure sensitivity on decanes — analog of MATCH 95:141-162 Section 5.4.

Decanes = C_{10} alkane isomers = trees of order 10 with maximum vertex
degree ≤ 4. There are exactly 75 such isomers.

We report three coefficients following the Furtula–Gutman convention used
in MATCH 95:141-162 Table 2:

    SS(TI)  := std(TI values) / mean(TI values)        (coefficient of variation; "structure sensitivity")
    Abr(TI) := (max − min) / mean(TI values)            ("abruptness")
    SA(TI)  := SS(TI) / Abr(TI)                         (ratio; higher = more gradual)

The MATCH paper reports SS ≈ 0.19–0.21 and Abr ≈ 0.38–0.42 across the
Sombor family on the same decane set. We compute the analog for Wuzi at
representative parameter triples and for the 30-index baseline.

Output: results/structure_sensitivity.csv and .md
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
MAX_DEGREE = 4   # C_10 alkanes: tree with max degree ≤ 4


def metrics(values: list[float]) -> tuple[float, float, float]:
    v = np.array(values, dtype=float)
    m = v.mean()
    if abs(m) < 1e-12:
        return (float("nan"), float("nan"), float("nan"))
    ss = v.std() / abs(m)
    abr = (v.max() - v.min()) / abs(m)
    sa = ss / abr if abs(abr) > 1e-12 else float("nan")
    return float(ss), float(abr), float(sa)


def main():
    print("Enumerating non-isomorphic decanes (trees of order 10 with max degree ≤ 4) ...")
    decanes = [T for T in nx.nonisomorphic_trees(TREE_ORDER)
               if max(d for _, d in T.degree()) <= MAX_DEGREE]
    print(f"Found {len(decanes)} decane structures.")
    # Sanity: the count of constitutional alkane isomers of C_10 H_22 is 75.
    if len(decanes) != 75:
        print(f"WARNING: expected 75 decanes, got {len(decanes)}.")

    print("Computing baseline indices ...")
    baseline_rows = [compute_all(T) for T in decanes]
    base_df = pd.DataFrame(baseline_rows)

    rows = []
    for col in base_df.columns:
        ss, abr, sa = metrics(list(base_df[col].values))
        rows.append({"index": col, "family": "baseline", "SS": ss, "Abr": abr, "SA": sa})

    wuzi_params = [
        (-0.5, 0.0, 2.0),
        (-1.0, -1.0, 1.0),
        (1.0, 0.0, 0.0),    # = M_2
        (-0.5, 0.0, 0.0),   # = Randic R
        (0.5, 0.3, 1.0),
        (0.0, 0.0, 1.0),
        (0.0, 0.0, 2.0),
    ]
    print(f"Computing Wuzi at {len(wuzi_params)} parameter triples ...")
    for (a, b, g) in wuzi_params:
        vals = [wuzi(T, a, b, g) for T in decanes]
        ss, abr, sa = metrics(vals)
        rows.append({"index": f"Wuzi(α={a}, β={b}, γ={g})",
                     "family": "wuzi", "SS": ss, "Abr": abr, "SA": sa})

    df = pd.DataFrame(rows).sort_values("SS", ascending=False).reset_index(drop=True)
    out_csv = os.path.join(PROJECT, "results", "structure_sensitivity.csv")
    df.to_csv(out_csv, index=False)

    md = ["# §6.5 Structure sensitivity on the 75 decane isomers\n"]
    md.append("Following the Furtula–Gutman convention used in MATCH 95:141-162 "
              "Table 2:\n\n"
              "- `SS = std / mean` (coefficient of variation)\n"
              "- `Abr = (max − min) / mean`\n"
              "- `SA = SS / Abr` (higher → more gradual response)\n")
    md.append("\nMATCH 95:141-162 reports SS ≈ 0.19–0.21 across Sombor variants "
              "on the same decane set.\n")
    md.append("| Index | Family | SS | Abr | SA |")
    md.append("|---|---|---:|---:|---:|")
    for _, r in df.iterrows():
        emph = "**" if r["family"] == "wuzi" else ""
        md.append(f"| {emph}{r['index']}{emph} | {r['family']} | "
                  f"{r['SS']:.4f} | {r['Abr']:.4f} | {r['SA']:.4f} |")
    md_path = os.path.join(PROJECT, "results", "structure_sensitivity.md")
    with open(md_path, "w") as f:
        f.write("\n".join(md))

    print(f"\n=== Top 10 by SS (most sensitive) ===")
    print(df.head(10).to_string(index=False))
    print(f"\nWritten {out_csv} and {md_path}")


if __name__ == "__main__":
    main()
