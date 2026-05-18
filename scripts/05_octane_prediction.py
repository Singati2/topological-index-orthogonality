"""§6.1 Octane prediction ability — direct analog of MATCH 95:141-162 Section 5.1.

Compares the Wuzi family (at canonical and at empirically-best parameter
points) against the 30-index baseline on the standard 18 octane isomers
benchmark, for 5 physicochemical properties:

  T_B    — normal boiling point (°C)
  dHf    — enthalpy of formation (kcal/mol)
  dHvap  — enthalpy of vaporization (kcal/mol)
  S      — standard entropy (cal/(mol·K))
  omega  — acentric factor (dimensionless)

Property values: public physical constants from NIST WebBook and standard
references; see `docs/literature_notes.md` for the [REFERENCE NEEDED]
citations to fill before submission. The same table is used in essentially
every Sombor-variant paper, including MATCH 95:141-162 Table 1.

Output: results/octane_prediction.csv and results/octane_prediction.md.
"""
from __future__ import annotations
import os
import sys
import numpy as np
import pandas as pd
from itertools import product

THIS_DIR = os.path.abspath(os.path.dirname(__file__))
PROJECT = os.path.abspath(os.path.join(THIS_DIR, ".."))
sys.path.insert(0, PROJECT)

from src.mol_to_graph import smiles_to_graph
from src.standard_indices import compute_all
from src.wuzi_index import wuzi


# 18 octane isomers (C8H18), SMILES + property values from NIST WebBook
# and standard physical-chemistry references. Values are public physical
# constants.
OCTANES = [
    # (name,                                  SMILES,                       T_B (°C), dHf (kcal/mol), dHvap (kcal/mol), S (cal/mol/K), omega)
    ("n-octane",                              "CCCCCCCC",                   125.6,    -49.82,         9.92,             111.55,        0.398),
    ("2-methylheptane",                       "CC(C)CCCCC",                 117.6,    -51.50,         9.48,             109.84,        0.378),
    ("3-methylheptane",                       "CCC(C)CCCC",                 118.9,    -50.82,         9.52,             111.26,        0.371),
    ("4-methylheptane",                       "CCCC(C)CCC",                 117.7,    -50.69,         9.48,             109.32,        0.372),
    ("2,2-dimethylhexane",                    "CC(C)(C)CCCC",               106.8,    -53.71,         8.92,             103.13,        0.339),
    ("2,3-dimethylhexane",                    "CC(C)C(C)CCC",               115.6,    -51.13,         9.27,             108.02,        0.348),
    ("2,4-dimethylhexane",                    "CC(C)CC(C)CC",               109.4,    -52.44,         9.03,             106.98,        0.344),
    ("2,5-dimethylhexane",                    "CC(C)CCC(C)C",               109.1,    -53.21,         9.05,             105.72,        0.357),
    ("3,3-dimethylhexane",                    "CCCC(C)(C)CC",               111.9,    -52.61,         9.04,             104.74,        0.322),
    ("3,4-dimethylhexane",                    "CCC(C)C(C)CC",               117.7,    -50.91,         9.32,             106.59,        0.340),
    ("2-methyl-3-ethylpentane",               "CCC(CC)C(C)C",               115.6,    -50.48,         9.21,             106.06,        0.330),
    ("3-methyl-3-ethylpentane",               "CCC(C)(CC)CC",               118.3,    -51.38,         9.21,             101.48,        0.302),
    ("3-ethylhexane",                         "CCCCC(CC)C",                 118.5,    -50.40,         9.48,             109.43,        0.362),
    ("2,2,3-trimethylpentane",                "CCC(C)C(C)(C)C",             109.8,    -52.61,         8.88,             101.31,        0.300),
    ("2,2,4-trimethylpentane",                "CC(C)CC(C)(C)C",             99.2,     -53.57,         8.40,             101.81,        0.305),
    ("2,3,3-trimethylpentane",                "CCC(C)(C)C(C)C",             114.8,    -51.73,         9.02,             101.31,        0.291),
    ("2,3,4-trimethylpentane",                "CC(C)C(C)C(C)C",             113.5,    -51.97,         9.01,             102.39,        0.317),
    ("2,2,3,3-tetramethylbutane",             "CC(C)(C)C(C)(C)C",           106.5,    -53.99,         8.41,             93.06,         0.247),
]

PROPS = ["T_B", "dHf", "dHvap", "S", "omega"]


def build_table() -> pd.DataFrame:
    rows = []
    for name, smiles, tb, hf, hv, s, om in OCTANES:
        G = smiles_to_graph(smiles)
        if G is None:
            raise RuntimeError(f"failed to parse SMILES {smiles!r} for {name}")
        rec = {"name": name, "smiles": smiles, "T_B": tb, "dHf": hf,
               "dHvap": hv, "S": s, "omega": om}
        rec.update(compute_all(G))
        # Wuzi at several representative (α, β, γ) for illustration; the
        # canonical-best from our cross-dataset analysis was around
        # (α, β, γ) = (−0.5, 0, 2) on ESOL, (−1, −1, 1) on FreeSolv.
        for (a, b, g) in [(-0.5, 0.0, 2.0), (-1.0, -1.0, 1.0),
                          (1.0, 0.0, 0.0), (-0.5, 0.0, 0.0),
                          (0.5, 0.3, 1.0)]:
            rec[f"Wuzi_a{a}_b{b}_g{g}"] = wuzi(G, a, b, g)
        rows.append(rec)
    return pd.DataFrame(rows)


def main():
    df = build_table()
    out_dir = os.path.join(PROJECT, "results")
    os.makedirs(out_dir, exist_ok=True)
    df.to_csv(os.path.join(out_dir, "octane_descriptors.csv"), index=False)

    # Correlation table — every index column vs every property
    index_cols = [c for c in df.columns if c not in ("name", "smiles", *PROPS)]
    corr_rows = []
    for c in index_cols:
        rec = {"index": c}
        for p in PROPS:
            rec[f"r_{p}"] = float(np.corrcoef(df[c], df[p])[0, 1])
        corr_rows.append(rec)
    corr = pd.DataFrame(corr_rows)
    corr.to_csv(os.path.join(out_dir, "octane_prediction.csv"), index=False)

    # Markdown table (manuscript-ready)
    md = ["# §6.1 Octane prediction ability — correlation of each index with each property\n"]
    md.append("Following the MATCH 95:141-162 Section 5.1 convention. n = 18 octane isomers; threshold for \"predicts\" is |r| > 0.8.\n")
    md.append("| Index | r(T_B) | r(dHf) | r(dHvap) | r(S) | r(omega) |")
    md.append("|---|---:|---:|---:|---:|---:|")
    # Sort by max |r| over properties, descending
    abs_max = corr[[f"r_{p}" for p in PROPS]].abs().max(axis=1)
    order = abs_max.sort_values(ascending=False).index
    for idx in order:
        row = corr.loc[idx]
        md.append(f"| {row['index']} | {row['r_T_B']:.3f} | {row['r_dHf']:.3f} | "
                  f"{row['r_dHvap']:.3f} | {row['r_S']:.3f} | {row['r_omega']:.3f} |")
    md.append("\n*Property values from NIST WebBook and standard physical-chemistry "
              "references; same dataset used in MATCH 95:141-162 Table 1. SMILES "
              "structures embedded in `scripts/05_octane_prediction.py`.*")
    md_path = os.path.join(out_dir, "octane_prediction.md")
    with open(md_path, "w") as f:
        f.write("\n".join(md))

    # Print top hits
    print(f"Octane prediction analysis on n = {len(df)} isomers.")
    print("Top 10 indices by |max r| across properties:")
    top10 = corr.iloc[order].head(10)
    print(top10.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    print(f"\nWritten {md_path}")


if __name__ == "__main__":
    main()
