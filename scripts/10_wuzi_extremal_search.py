"""Computational extremal-graph search for the Wuzi index family on trees.

Enumerates all non-isomorphic trees of order 5 <= n <= 12 using
networkx.generators.nonisomorphic_trees(n), computes
W(T; alpha, beta, gamma) at each parameter triple, and records the
minimizer and maximizer over the tree class.

Then checks whether the observed minimizer matches the path P_n and
the observed maximizer matches the star S_n at each parameter
triple (Conjecture D.1 of docs/wuzi_bounds_strategy.md).

This is *conjecture-generating* evidence, not proof. The output
should help Advik / Arockiaraj sir target a formal extremal-graph
characterization.

Outputs:
  results/wuzi_extremal_trees.csv
  results/wuzi_extremal_trees.md

Unicyclic and bicyclic enumeration is documented as a TODO below;
networkx does not provide a built-in non-isomorphic-unicyclic
enumerator and a hand-rolled version is left for future work.
"""
from __future__ import annotations
import math
import os
import sys
import time
from typing import Iterable

import networkx as nx
import pandas as pd

THIS_DIR = os.path.abspath(os.path.dirname(__file__))
PROJECT = os.path.abspath(os.path.join(THIS_DIR, ".."))
sys.path.insert(0, PROJECT)

TREE_ORDERS = list(range(5, 13))   # 5..12 inclusive

PARAMETER_TRIPLES: list[tuple[float, float, float]] = [
    (1.0, 0.0, 0.0),       # M_2
    (-0.5, 0.0, 0.0),      # Randic R
    (0.0, -0.5, 0.0),      # Sum-connectivity
    (0.0, -1.0, 0.0),      # Half-harmonic
    (0.0, 0.0, 1.0),       # gamma>0
    (0.0, 0.0, 2.0),
    (1.0, 1.0, 1.0),
    (-1.0, -1.0, 1.0),
]


def wuzi(G: nx.Graph, alpha: float, beta: float, gamma: float) -> float:
    """Compute W(G; alpha, beta, gamma) directly from the definition."""
    total = 0.0
    for u, v in G.edges():
        du, dv = G.degree(u), G.degree(v)
        prod = du * dv
        sm = du + dv
        imb = abs(du - dv) / sm
        total += (prod ** alpha) * (sm ** beta) * math.exp(gamma * imb)
    return total


def canonical_label(G: nx.Graph) -> str:
    """Deterministic isomorphism-invariant label for a tree.

    The 1-Weisfeiler-Lehman test is known to identify all
    non-isomorphic trees (folklore; see e.g. Cai-Furer-Immerman
    1992). For trees of order at most 12 the WL hash is therefore
    a sound isomorphism invariant, not merely a heuristic.
    """
    return nx.weisfeiler_lehman_graph_hash(G)


def degree_sequence(G: nx.Graph) -> tuple[int, ...]:
    """Multiset of vertex degrees, sorted descending."""
    return tuple(sorted((d for _, d in G.degree()), reverse=True))


def is_caterpillar(T: nx.Graph) -> bool:
    """True iff T is a caterpillar tree.

    A caterpillar is a tree whose non-leaf vertices induce a path
    (possibly empty or a single vertex). Equivalently: remove all
    leaves and the remaining graph is a path P_k for some k >= 0.

    For T with at most 2 vertices the property holds vacuously.
    """
    if T.number_of_nodes() <= 2:
        return True
    leaves = [v for v, d in T.degree() if d == 1]
    interior = T.copy()
    interior.remove_nodes_from(leaves)
    n_int = interior.number_of_nodes()
    if n_int == 0:
        # All vertices were leaves -> T was a star K_{1, m} or P_2; caterpillar.
        return True
    if n_int == 1:
        # Single interior vertex; the leaves all attach to it -> star.
        return True
    # interior must be a path: connected, every vertex of degree <= 2.
    if not nx.is_connected(interior):
        return False
    return max(d for _, d in interior.degree()) <= 2


def path_label(n: int) -> str:
    return canonical_label(nx.path_graph(n))


def star_label(n: int) -> str:
    return canonical_label(nx.star_graph(n - 1))   # star_graph(k) has k+1 vertices


def main() -> None:
    t0 = time.time()
    rows = []
    print("=" * 64)
    print("Wuzi extremal-graph search on trees of order 5..12")
    print("=" * 64)

    for n in TREE_ORDERS:
        trees = list(nx.generators.nonisomorphic_trees(n))
        n_trees = len(trees)
        path_lbl = path_label(n)
        star_lbl = star_label(n)
        print(f"\nn = {n}: {n_trees} non-isomorphic trees")
        for (alpha, beta, gamma) in PARAMETER_TRIPLES:
            best_min = math.inf
            best_max = -math.inf
            argmin_idx = -1
            argmax_idx = -1
            for idx, T in enumerate(trees):
                w = wuzi(T, alpha, beta, gamma)
                if w < best_min:
                    best_min = w
                    argmin_idx = idx
                if w > best_max:
                    best_max = w
                    argmax_idx = idx
            argmin_lbl = canonical_label(trees[argmin_idx])
            argmax_lbl = canonical_label(trees[argmax_idx])
            min_is_path = argmin_lbl == path_lbl
            max_is_star = argmax_lbl == star_lbl
            # Also record the dual case (Randic-region sign flip)
            min_is_star = argmin_lbl == star_lbl
            max_is_path = argmax_lbl == path_lbl
            argmin_degseq = degree_sequence(trees[argmin_idx])
            argmax_degseq = degree_sequence(trees[argmax_idx])
            argmin_is_caterpillar = is_caterpillar(trees[argmin_idx])
            argmax_is_caterpillar = is_caterpillar(trees[argmax_idx])
            classical_match = (
                "P_n=min,S_n=max" if (min_is_path and max_is_star) else
                "S_n=min,P_n=max" if (min_is_star and max_is_path) else
                "neither"
            )
            w_path = wuzi(nx.path_graph(n), alpha, beta, gamma)
            w_star = wuzi(nx.star_graph(n - 1), alpha, beta, gamma)
            rows.append({
                "n": n,
                "n_trees": n_trees,
                "alpha": alpha,
                "beta": beta,
                "gamma": gamma,
                "w_min": best_min,
                "w_max": best_max,
                "w_path": w_path,
                "w_star": w_star,
                "min_is_path": min_is_path,
                "max_is_star": max_is_star,
                "min_is_star": min_is_star,
                "max_is_path": max_is_path,
                "classical_match": classical_match,
                "argmin_degree_seq": ",".join(str(d) for d in argmin_degseq),
                "argmax_degree_seq": ",".join(str(d) for d in argmax_degseq),
                "argmin_is_caterpillar": argmin_is_caterpillar,
                "argmax_is_caterpillar": argmax_is_caterpillar,
            })

    df = pd.DataFrame(rows)
    out_csv = os.path.join(PROJECT, "results", "wuzi_extremal_trees.csv")
    df.to_csv(out_csv, index=False)

    md_lines: list[str] = []
    md_lines.append("# Wuzi extremal-graph search on trees of order 5..12\n")
    md_lines.append(f"Generated by `scripts/10_wuzi_extremal_search.py`.\n")
    md_lines.append(
        "For each tree order $n \\in \\{5, \\ldots, 12\\}$ and each "
        "parameter triple $(\\alpha, \\beta, \\gamma)$, the table below "
        "reports the minimum and maximum of $W(T; \\alpha, \\beta, \\gamma)$ "
        "over all non-isomorphic trees $T$ of order $n$, together with "
        "the Wuzi value at the path $P_n$ and the star $S_n$.\n\n"
        "Columns `min_is_path` and `max_is_star` indicate whether the "
        "observed minimizer coincides (via Weisfeiler--Lehman hash) with "
        "$P_n$ and the observed maximizer with $S_n$ respectively. A "
        "consistent `True` across orders is empirical support for "
        "Conjecture D.1 of `docs/wuzi_bounds_strategy.md`; a `False` "
        "indicates the conjecture fails at that parameter triple and "
        "the actual extremal is a different tree.\n"
    )

    md_lines.append("## Detailed table\n")
    md_lines.append(
        "| $n$ | trees | $\\alpha$ | $\\beta$ | $\\gamma$ | $W_{\\min}$ | $W_{\\max}$ "
        "| $W(P_n)$ | $W(S_n)$ | classical match | argmin deg seq | argmax deg seq |"
    )
    md_lines.append("|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|")
    for r in rows:
        md_lines.append(
            f"| {r['n']} | {r['n_trees']} | "
            f"{r['alpha']:g} | {r['beta']:g} | {r['gamma']:g} | "
            f"{r['w_min']:.4f} | {r['w_max']:.4f} | "
            f"{r['w_path']:.4f} | {r['w_star']:.4f} | "
            f"{r['classical_match']} | "
            f"{r['argmin_degree_seq']} | {r['argmax_degree_seq']} |"
        )

    md_lines.append("\n## Summary across $n$\n")
    md_lines.append(
        "For each parameter triple, the fraction of orders for which "
        "the four possible classical extremal patterns hold:\n\n"
        "- $P_n = \\min$, $S_n = \\max$ (Conjecture D.1 stated form)\n"
        "- $S_n = \\min$, $P_n = \\max$ (the sign-flipped dual, expected "
        "in the Randić / sum-connectivity / harmonic sign region)\n"
        "- neither (the extremal trees are neither $P_n$ nor $S_n$)\n"
    )
    summary_rows = []
    for (alpha, beta, gamma) in PARAMETER_TRIPLES:
        sub = df[(df["alpha"] == alpha) & (df["beta"] == beta) & (df["gamma"] == gamma)]
        f_p_min_s_max = float(((sub["min_is_path"]) & (sub["max_is_star"])).mean())
        f_s_min_p_max = float(((sub["min_is_star"]) & (sub["max_is_path"])).mean())
        f_neither = float(((sub["classical_match"]) == "neither").mean())
        summary_rows.append((alpha, beta, gamma, f_p_min_s_max, f_s_min_p_max, f_neither))
    md_lines.append("\n| $\\alpha$ | $\\beta$ | $\\gamma$ | $P_n=\\min, S_n=\\max$ | $S_n=\\min, P_n=\\max$ | neither |")
    md_lines.append("|---:|---:|---:|---:|---:|---:|")
    for (alpha, beta, gamma, fa, fb, fc) in summary_rows:
        md_lines.append(f"| {alpha:g} | {beta:g} | {gamma:g} | {fa:.2f} | {fb:.2f} | {fc:.2f} |")

    md_lines.append(
        "\n## Interpretation\n\n"
        "A row with $P_n = \\min, S_n = \\max$ fraction $= 1.00$ across "
        "all $n$ is empirical support for Conjecture D.1 of "
        "`docs/wuzi_bounds_strategy.md` in its stated $\\alpha,\\beta,\\gamma \\ge 0$ "
        "sign region. The sign-flipped column $S_n = \\min, P_n = \\max$ "
        "captures the dual extremal pattern expected for Randić-style negative "
        "parameter values. The neither column flags the most interesting "
        "case for genuine new mathematics: parameter triples where the "
        "extremal trees are neither $P_n$ nor $S_n$ and a non-classical "
        "extremal structure must be characterized. The columns "
        "`argmin_degree_seq` and `argmax_degree_seq` in the detailed table "
        "record the degree sequence of the observed extremal at each "
        "$(n, \\alpha, \\beta, \\gamma)$, identifying candidate non-classical "
        "extremal trees.\n"
    )
    md_lines.append(
        "## Unicyclic and bicyclic enumeration (TODO)\n\n"
        "`networkx` does not currently expose a built-in enumerator for "
        "non-isomorphic unicyclic or bicyclic graphs. A hand-rolled "
        "enumeration via SageMath, the `nauty` `geng` tool, or a custom "
        "Python implementation is left as a follow-up. For small orders "
        "($n \\le 9$) the counts are tractable: $46$ non-isomorphic "
        "unicyclic graphs at $n = 6$, $215$ at $n = 7$, etc.\n"
    )

    out_md = os.path.join(PROJECT, "results", "wuzi_extremal_trees.md")
    with open(out_md, "w") as f:
        f.write("\n".join(md_lines))

    print()
    print(df.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    print()
    print(f"Written {out_csv}")
    print(f"Written {out_md}")
    print(f"Done in {time.time() - t0:.1f} s.")


if __name__ == "__main__":
    main()
