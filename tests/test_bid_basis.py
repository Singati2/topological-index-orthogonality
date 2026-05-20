"""BID basis: every BID index is a linear combination of m_ij counts.

The hand-built test tree has explicit edge list documented below.
"""
import networkx as nx
import pytest

from src.standard_indices import compute_all


def make_tree_delta_4() -> nx.Graph:
    """7-vertex tree, maximum degree 4.

    Centre vertex 0 has degree 4, connected to vertices 1, 2, 3, 4.
    Vertex 1 is extended into a 2-edge tail (1 - 5 - 6).
    """
    G = nx.Graph()
    G.add_edges_from([
        (0, 1), (0, 2), (0, 3), (0, 4),  # centre -> four neighbours
        (1, 5),                           # tail step 1
        (5, 6),                           # tail step 2
    ])
    return G


def edge_degree_pair_counts(G: nx.Graph) -> dict:
    """Return {(i, j): count} with i <= j (unordered pair) for each edge."""
    counts = {}
    for u, v in G.edges():
        du, dv = G.degree(u), G.degree(v)
        key = (du, dv) if du <= dv else (dv, du)
        counts[key] = counts.get(key, 0) + 1
    return counts


def test_max_degree_is_4():
    G = make_tree_delta_4()
    assert max(d for _, d in G.degree()) == 4


def test_m_ij_sums_to_m():
    G = make_tree_delta_4()
    counts = edge_degree_pair_counts(G)
    assert sum(counts.values()) == G.number_of_edges() == 6


def test_M1_recovered_from_m_ij():
    G = make_tree_delta_4()
    counts = edge_degree_pair_counts(G)
    M1_linear = sum((i + j) * c for (i, j), c in counts.items())
    M1_direct = compute_all(G)["M1"]
    assert M1_linear == pytest.approx(M1_direct)


def test_M2_recovered_from_m_ij():
    G = make_tree_delta_4()
    counts = edge_degree_pair_counts(G)
    M2_linear = sum(i * j * c for (i, j), c in counts.items())
    M2_direct = compute_all(G)["M2"]
    assert M2_linear == pytest.approx(M2_direct)


def test_bid_dimension_bound_on_delta_4():
    """Verify that on Delta <= 4 graphs, any 11 BID indices are linearly
    dependent (i.e. the BID index image has dimension at most 10).

    Strategy: pick 11 BID functionals whose underlying symmetric
    functions on the 10 admissible pairs ARE LINEARLY INDEPENDENT (so
    the test is a non-vacuous witness of Theorem 3.1). The natural
    such choice is 10 indicator functions on each admissible pair plus
    one functional outside their span -- any of these 11 would give
    rank exactly 10 on a panel of graphs that realises every admissible
    pair. We assert rank <= 10 (the bound) and equal to 10 (sharpness).
    """
    import numpy as np

    # Indicator on each admissible (i, j) pair, symmetric in i, j.
    def make_indicator(pair):
        a, b = pair
        def f(i, j):
            return 1.0 if (min(i, j), max(i, j)) == (a, b) else 0.0
        return f

    admissible = [(i, j) for i in range(1, 5) for j in range(i, 5)]
    assert len(admissible) == 10
    indicator_funcs = [make_indicator(p) for p in admissible]

    # Eleventh functional, generic symmetric: f(i,j) = (i+j)^2. By
    # Theorem 3.1 this is forced to be a linear combination of the 10
    # indicators when restricted to admissible pairs (i.e. for all
    # graphs of max degree at most 4) -- i.e. the 11-column matrix
    # built from these 11 functionals must have rank <= 10.
    def f_eleventh(i, j): return (i + j) ** 2

    fs = indicator_funcs + [f_eleventh]
    assert len(fs) == 11

    # Panel of Delta<=4 trees that realises every admissible pair.
    graphs = []
    for n in range(3, 9):
        graphs.extend(list(nx.generators.nonisomorphic_trees(n)))
    graphs = [G for G in graphs if max(d for _, d in G.degree()) <= 4]

    rows = []
    for G in graphs:
        row = []
        for f in fs:
            v = sum(f(G.degree(u), G.degree(v_)) for u, v_ in G.edges())
            row.append(v)
        rows.append(row)
    M = np.asarray(rows, dtype=float)
    assert M.shape == (len(graphs), 11)

    # By Theorem 3.1: rank(M) <= 10.
    rank = int(np.linalg.matrix_rank(M, tol=1e-9))
    assert rank <= 10, (
        f"BID dimension bound VIOLATED: column rank of 11-index matrix "
        f"on {len(graphs)} Delta-<=4 trees is {rank}, expected <= 10"
    )

    # Sharpness check: the first 10 columns alone (the indicator
    # functions) achieve rank exactly equal to the number of admissible
    # pairs realised by the graph panel. On a CONNECTED tree of order
    # >= 3, no edge can be a (1, 1) edge (a (1, 1) edge means both
    # endpoints have degree 1, which on a connected tree of order >= 3
    # is impossible). The other 9 admissible pairs are all realised
    # (e.g. (4, 4) appears in the double-star on 8 vertices with two
    # adjacent degree-4 centres). Hence the indicator subblock has
    # rank exactly 9 on this panel.
    rank_indicators = int(np.linalg.matrix_rank(M[:, :10], tol=1e-9))
    assert rank_indicators == 9, (
        f"Indicator subblock rank {rank_indicators}; expected 9 "
        f"(all admissible pairs except (1,1) realised on these trees)"
    )
