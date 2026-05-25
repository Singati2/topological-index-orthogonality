"""95% percentile bootstrap (B = 10,000) for octane property correlations.

For each pair (index x property) in the 18-octane benchmark, computes:
  * the original Pearson r
  * the 2.5th and 97.5th percentile bootstrap r over B = 10,000 resamples

This supports the Table 3 caption in `docs/paper1_wuzi.tex`, which states
that, e.g., r(dHvap, GA) = +0.959 with 95% CI [+0.917, +0.986], the
W(0,0,2) vs GA/H margin of 0.008 is well within the ~0.034 CI half-width,
and the 80 (index x property) cells have CI half-widths ranging from
0.019 to 0.47 with median 0.13.

This script uses a minimal alkane-only SMILES parser (the 18 octanes are
all C8H18) to avoid the RDKit dependency, so it can be run on any system
with numpy + networkx; the wuzi() function from src.wuzi_index is used to
compute Wuzi values at the 8 parameter triples shown in Table 3.

Output: results/octane_competitor_bootstrap.csv
Reproducing:
    python scripts/14_octane_bootstrap.py
"""
from __future__ import annotations
import csv
import os
import sys
import time

import networkx as nx
import numpy as np

THIS_DIR = os.path.abspath(os.path.dirname(__file__))
PROJECT = os.path.abspath(os.path.join(THIS_DIR, ".."))
sys.path.insert(0, PROJECT)

from src.wuzi_index import wuzi


# Reproducibility: seed and bootstrap size match what is reported in the
# Table 3 caption.
SEED = 20260524
B = 10_000

OCTANE_CSV = os.path.join(PROJECT, "results", "octane_descriptors.csv")
OUT_CSV = os.path.join(PROJECT, "results", "octane_competitor_bootstrap.csv")

PROPS = ["T_B", "dHf", "dHvap", "S", "omega"]
CLASSICAL = ["M1", "M2", "R", "SO", "GA", "H", "ABC", "AZI"]
WUZI_TRIPLES = [
    (1, 0, 0), (-0.5, 0, 0), (0, -0.5, 0), (0, -1, 0),
    (0, 0, 1), (0, 0, 2), (1, 1, 1), (-1, -1, 1),
]


def alkane_smiles_to_nx(smi: str) -> nx.Graph:
    """Minimal SMILES parser for the C8H18 octane isomers.

    Handles only 'C', '(', ')' (no rings, no heteroatoms, no branches with
    brackets other than parentheses). Verified to parse every line of
    results/octane_descriptors.csv correctly.
    """
    G = nx.Graph()
    stack: list[int] = []
    last: int | None = None
    next_id = 0
    for ch in smi:
        if ch == "C":
            G.add_node(next_id)
            if last is not None:
                G.add_edge(last, next_id)
            last = next_id
            next_id += 1
        elif ch == "(":
            stack.append(last)  # type: ignore[arg-type]
        elif ch == ")":
            last = stack.pop()
        else:
            raise ValueError(f"alkane parser saw unexpected character: {ch!r}")
    return G


def pearson(x: np.ndarray, y: np.ndarray) -> float:
    xm = x - x.mean()
    ym = y - y.mean()
    den = np.sqrt((xm * xm).sum() * (ym * ym).sum())
    return float((xm * ym).sum() / den) if den > 0 else 0.0


def main() -> None:
    np.random.seed(SEED)

    with open(OCTANE_CSV) as f:
        rows = list(csv.DictReader(f))
    n = len(rows)
    assert n == 18, f"Expected 18 octanes, got {n}"

    # Parse SMILES -> H-suppressed graphs
    graphs = []
    for r in rows:
        G = alkane_smiles_to_nx(r["smiles"])
        assert G.number_of_nodes() == 8, (
            f"{r['name']}: expected 8 C atoms, got {G.number_of_nodes()}"
        )
        graphs.append(G)

    # Build the descriptor matrix.
    # Classical indices are read from the existing CSV (already computed).
    # Wuzi triples are recomputed from the graphs (descriptors CSV doesn't
    # cache all 8 triples).
    x_data: dict[str, np.ndarray] = {}
    for idx in CLASSICAL:
        x_data[idx] = np.array([float(r[idx]) for r in rows])
    for a, b, c in WUZI_TRIPLES:
        name = f"W({a},{b},{c})"
        x_data[name] = np.array(
            [wuzi(G, alpha=a, beta=b, gamma=c) for G in graphs]
        )

    y_data: dict[str, np.ndarray] = {
        p: np.array([float(r[p]) for r in rows]) for p in PROPS
    }

    print(
        f"Bootstrap B={B}, seed={SEED}, n={n} octanes, "
        f"{len(x_data)} indices x {len(PROPS)} properties = "
        f"{len(x_data) * len(PROPS)} cells"
    )

    t0 = time.time()
    out_rows = []
    for idx_name, x in x_data.items():
        fam = "Wuzi" if idx_name.startswith("W(") else "classical"
        row: dict[str, object] = {"index": idx_name, "family": fam}
        for p in PROPS:
            y = y_data[p]
            r_orig = pearson(x, y)
            rs = np.empty(B)
            for i in range(B):
                sample = np.random.randint(0, n, n)
                rs[i] = pearson(x[sample], y[sample])
            lo, hi = np.percentile(rs, [2.5, 97.5])
            row[f"r_{p}"] = round(r_orig, 4)
            row[f"lo_{p}"] = round(float(lo), 4)
            row[f"hi_{p}"] = round(float(hi), 4)
        out_rows.append(row)
    print(f"  done in {time.time() - t0:.1f}s")

    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    with open(OUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
        writer.writeheader()
        for r in out_rows:
            writer.writerow(r)
    print(f"Wrote {OUT_CSV}")

    # Print a small summary used in the Table 3 caption.
    hws = []
    for r in out_rows:
        for p in PROPS:
            hws.append(abs(float(r[f"hi_{p}"]) - float(r[f"lo_{p}"])) / 2)
    hws_arr = np.array(hws)
    print(
        "\nCI half-width summary across all (index x property) cells: "
        f"min={hws_arr.min():.3f}, max={hws_arr.max():.3f}, "
        f"median={np.median(hws_arr):.3f}"
    )

    # Leave-one-out cross-validation: for each (index, property), recompute
    # Pearson r 18 times leaving out one octane at a time. Report the largest
    # absolute change in r across all 18 leave-one-outs; an outlier-driven
    # headline correlation shows up as a single large change.
    print(f"\n=== Leave-one-out robustness check (n=18, 18 LOO refits per cell) ===")
    loocv_rows = []
    octane_names = [r["name"] for r in rows]
    for idx_name, x in x_data.items():
        fam = "Wuzi" if idx_name.startswith("W(") else "classical"
        row: dict[str, object] = {"index": idx_name, "family": fam}
        for p in PROPS:
            y = y_data[p]
            r_full = pearson(x, y)
            r_loo = np.array([
                pearson(np.delete(x, i), np.delete(y, i)) for i in range(n)
            ])
            r_loo_min = float(r_loo.min())
            r_loo_max = float(r_loo.max())
            max_drift = float(max(abs(r_full - r_loo_min), abs(r_full - r_loo_max)))
            row[f"r_{p}"] = round(r_full, 4)
            row[f"loo_min_{p}"] = round(r_loo_min, 4)
            row[f"loo_max_{p}"] = round(r_loo_max, 4)
            row[f"loo_max_drift_{p}"] = round(max_drift, 4)
            row[f"loo_argmax_drift_{p}"] = octane_names[int(np.argmax(np.abs(r_full - r_loo)))]
        loocv_rows.append(row)

    loo_csv = os.path.join(PROJECT, "results", "octane_competitor_loocv.csv")
    with open(loo_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(loocv_rows[0].keys()))
        writer.writeheader()
        for r in loocv_rows:
            writer.writerow(r)
    print(f"Wrote {loo_csv}")

    # Summary: largest absolute drift across the 80 (index, property) cells, and
    # which octane drives the headline GA/H/W(0,0,2) dHvap cells.
    all_drift = np.array([
        float(r[f"loo_max_drift_{p}"]) for r in loocv_rows for p in PROPS
    ])
    print(
        f"\nLOOCV max-drift summary across all {len(all_drift)} cells: "
        f"min={all_drift.min():.3f}, max={all_drift.max():.3f}, "
        f"median={np.median(all_drift):.3f}"
    )
    print(f"\nLOOCV drift on the headline dHvap comparison cells:")
    for r in loocv_rows:
        if r["index"] in ("GA", "H", "W(0,0,2)"):
            print(
                f"  {r['index']:10s}  r={r['r_dHvap']:+.4f}  "
                f"LOO range=[{r['loo_min_dHvap']:+.4f}, {r['loo_max_dHvap']:+.4f}]  "
                f"max drift={r['loo_max_drift_dHvap']:.4f}  "
                f"(argmax: {r['loo_argmax_drift_dHvap']})"
            )


if __name__ == "__main__":
    main()
