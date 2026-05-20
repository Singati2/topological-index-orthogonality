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

    Strategy: pick 11 distinct symmetric f(i, j) -> R functions; evaluate
    each on a panel of small Delta <= 4 graphs to get an 11-vector of
    index values per graph; stack into a matrix and assert the column rank
    is at most 10.
    """
    import numpy as np

    def f_const(i, j): return 1.0
    def f_sum(i, j): return i + j
    def f_prod(i, j): return i * j
    def f_max(i, j): return max(i, j)
    def f_min(i, j): return min(i, j)
    def f_diff_abs(i, j): return abs(i - j)
    def f_sumsq(i, j): return i * i + j * j
    def f_diffsq(i, j): return (i - j) ** 2
    def f_invsum(i, j): return 1.0 / (i + j)
    def f_invprod(i, j): return 1.0 / (i * j)
    def f_random_symmetric(i, j):
        # Another symmetric BID functional; choose any closed form.
        return (i + j) * abs(i - j) + i * j

    fs = [f_const, f_sum, f_prod, f_max, f_min, f_diff_abs,
          f_sumsq, f_diffsq, f_invsum, f_invprod, f_random_symmetric]
    assert len(fs) == 11, "Need 11 BID functions"

    # Build a panel of Delta <= 4 trees of various orders.
    graphs = []
    for n in range(3, 9):
        graphs.extend(list(nx.generators.nonisomorphic_trees(n)))
    # Keep only graphs with max degree <= 4 (all trees of order <= 8
    # already satisfy this, but be defensive).
    graphs = [G for G in graphs if max(d for _, d in G.degree()) <= 4]

    # For each graph, compute the 11-vector of (sum over edges) BID values.
    rows = []
    for G in graphs:
        row = []
        for f in fs:
            v = sum(f(G.degree(u), G.degree(v_)) for u, v_ in G.edges())
            row.append(v)
        rows.append(row)
    M = np.asarray(rows, dtype=float)
    assert M.shape == (len(graphs), 11)

    # The image of the BID map (graphs -> R^11) sits in a subspace of
    # dimension at most 10 (by Theorem 3.1). Hence the column rank of M
    # is at most 10.
    rank = int(np.linalg.matrix_rank(M, tol=1e-9))
    assert rank <= 10, (
        f"BID dimension bound violated: column rank of 11-index matrix "
        f"on {len(graphs)} Delta-<=4 trees is {rank}, expected <= 10"
    )
