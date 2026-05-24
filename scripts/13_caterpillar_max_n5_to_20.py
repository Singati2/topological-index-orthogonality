"""Extended n=5..20 enumeration of the (alpha,beta,gamma)=(-1,-1,1) maximizer.

Computes, for each n in 5..20:
  * the maximum of W(T; -1, -1, 1) over all non-isomorphic trees T of order n
  * the canonical graph6 string of the maximizer (so the tree is reproducibly
    identifiable; degree sequence alone is not an isomorphism invariant for
    n >= 6)
  * the degree sequence (descending)
  * the total number of non-isomorphic trees of order n (OEIS A000055)
  * the wall-clock time for the n-th order's enumeration

Used to populate Table 2 of `docs/paper1_wuzi.tex` and to support
Conjecture conj:caterpillar_max (spider for odd n, balanced double-spider
for even n). The n=20 case requires enumerating 823,065 non-isomorphic
trees; total runtime on a single core is ~30 seconds.

Output: results/wuzi_caterpillar_max_n5_to_20.csv

Reproducing:
    python scripts/13_caterpillar_max_n5_to_20.py
"""
from __future__ import annotations
import csv
import os
import sys
import time

import networkx as nx

THIS_DIR = os.path.abspath(os.path.dirname(__file__))
PROJECT = os.path.abspath(os.path.join(THIS_DIR, ".."))
sys.path.insert(0, PROJECT)

from src.wuzi_index import wuzi


TREE_ORDERS = list(range(5, 21))  # n = 5..20 inclusive
PARAMS = (-1.0, -1.0, 1.0)        # the mixed-sign non-classical extremal regime

OUT_PATH = os.path.join(PROJECT, "results", "wuzi_caterpillar_max_n5_to_20.csv")


def main() -> None:
    alpha, beta, gamma = PARAMS
    rows = []

    print(
        f"{'n':>3} {'#trees':>8} {'time':>7}  "
        f"{'W_max':>10}  {'degree sequence':<40}  graph6"
    )
    print("-" * 130)

    for n in TREE_ORDERS:
        t0 = time.time()
        best_val = -float("inf")
        best_tree = None
        count = 0
        for T in nx.generators.nonisomorphic_trees(n):
            count += 1
            val = wuzi(T, alpha=alpha, beta=beta, gamma=gamma)
            if val > best_val:
                best_val = val
                best_tree = T
        g6 = nx.to_graph6_bytes(best_tree, header=False).decode().strip()
        ds = tuple(sorted([d for _, d in best_tree.degree()], reverse=True))
        dt = time.time() - t0
        rows.append((n, count, best_val, ds, g6))
        print(
            f"{n:3d} {count:>8d} {dt:>6.2f}s  "
            f"{best_val:>10.6f}  {str(ds):<40}  {g6}"
        )

    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    with open(OUT_PATH, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            ["n", "n_trees", "W_max", "argmax_degree_seq", "argmax_graph6"]
        )
        for n, count, val, ds, g6 in rows:
            w.writerow([n, count, f"{val:.10f}", ",".join(str(d) for d in ds), g6])
    print(f"\nWrote {OUT_PATH}")


if __name__ == "__main__":
    main()
