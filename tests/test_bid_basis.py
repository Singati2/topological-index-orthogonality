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
