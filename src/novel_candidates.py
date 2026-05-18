"""Candidate topological indices from five design philosophies that are
*structurally* different from edge-sum-of-degree-functions families.

The Wuzi grid search showed degree-based parametric families cannot escape
orthogonality with the existing 30-index baseline on real chemistry (0/100).
This module tests whether *any* of these alternative design spaces can.

Design philosophies tested:
  1. Information-theoretic    — entropy of degree / distance / spectrum
  2. Centrality-based         — betweenness, closeness, eigenvector, harmonic
  3. Spectral hybrids         — algebraic connectivity, spectral*diameter
  4. Motif / structural       — triangles, clustering, 4-cycles, cut vertices
  5. Eccentricity-based       — total eccentricity, radius, mean eccentricity

Each function takes a NetworkX graph G and returns a float.
"""
from __future__ import annotations
import math
import numpy as np
import networkx as nx
from collections import Counter


def _adj_eigs(G):
    A = nx.adjacency_matrix(G).toarray().astype(float)
    return np.linalg.eigvalsh(A)

def _laplacian_eigs(G):
    L = nx.laplacian_matrix(G).toarray().astype(float)
    return np.linalg.eigvalsh(L)

def _dist_matrix(G):
    return nx.floyd_warshall_numpy(G).astype(int)

def _shannon(counts):
    """Shannon entropy of a count distribution, log base e."""
    total = sum(counts)
    if total == 0:
        return 0.0
    return -sum((c / total) * math.log(c / total) for c in counts if c > 0)


# ---------- 1. Information-theoretic ----------

def shannon_degree_entropy(G):
    degs = [d for _, d in G.degree()]
    return _shannon(list(Counter(degs).values()))

def shannon_distance_entropy(G):
    D = _dist_matrix(G)
    n = G.number_of_nodes()
    iu = np.triu_indices(n, k=1)
    pair_dists = D[iu]
    return _shannon(list(Counter(pair_dists.tolist()).values()))

def spectral_entropy(G):
    eigs = _adj_eigs(G)
    sq = eigs * eigs
    s = sq.sum()
    if s <= 0:
        return 0.0
    p = sq / s
    return -float(np.sum([pi * math.log(pi) for pi in p if pi > 0]))

def eccentricity_entropy(G):
    eccs = list(nx.eccentricity(G).values())
    return _shannon(list(Counter(eccs).values()))

def bonchev_distance_information(G):
    """Bonchev/Trinajstic I_d: information-theoretic distance index.

    I_d = W * log2(W) - Σ k * d_k * log2(d_k), where d_k = number of pairs
    at distance k and W = Wiener.
    """
    D = _dist_matrix(G)
    n = G.number_of_nodes()
    iu = np.triu_indices(n, k=1)
    pair_dists = D[iu].tolist()
    counts = Counter(pair_dists)
    W = sum(k * v for k, v in counts.items())
    if W <= 0:
        return 0.0
    total = W * math.log2(W) - sum(k * v * math.log2(k * v) for k, v in counts.items() if k > 0 and v > 0)
    return total


# ---------- 2. Centrality-based ----------

def betweenness_sum(G):
    return float(sum(nx.betweenness_centrality(G).values()))

def closeness_sum(G):
    return float(sum(nx.closeness_centrality(G).values()))

def eigenvector_centrality_sum(G):
    try:
        ec = nx.eigenvector_centrality_numpy(G)
        return float(sum(ec.values()))
    except Exception:
        return 0.0

def harmonic_centrality_sum(G):
    return float(sum(nx.harmonic_centrality(G).values()))

def betweenness_variance(G):
    vals = list(nx.betweenness_centrality(G).values())
    return float(np.var(vals))


# ---------- 3. Spectral hybrids ----------

def algebraic_connectivity(G):
    """Second-smallest Laplacian eigenvalue λ_2 (Fiedler value).
    Sensitive to bottlenecks / cut-edges — captures global structure
    that no degree-based index sees directly.
    """
    leigs = _laplacian_eigs(G)
    leigs = np.sort(leigs)
    if len(leigs) < 2:
        return 0.0
    return float(leigs[1])

def spectral_times_diameter(G):
    eigs = _adj_eigs(G)
    rho = float(np.max(np.abs(eigs)))
    return rho * nx.diameter(G)

def estrada_times_wiener(G):
    eigs = _adj_eigs(G)
    ee = float(np.sum(np.exp(eigs)))
    D = _dist_matrix(G)
    n = G.number_of_nodes()
    W = float(D.sum() / 2.0)
    # Take log to tame the dynamic range
    return math.log1p(ee) * math.log1p(W)


# ---------- 4. Motif / structural ----------

def triangle_count(G):
    return float(sum(nx.triangles(G).values()) // 3)

def mean_clustering(G):
    return float(nx.average_clustering(G))

def four_cycle_count(G):
    """Number of induced 4-cycles. NetworkX has nx.cycle_basis but for
    counting 4-cycles directly, use a quick combinatorial approach.
    """
    A = nx.adjacency_matrix(G).toarray()
    # (A^4)_{ii} counts closed walks of length 4. Subtract degeneracies:
    # 4-walks through degree edges (vertex back-and-forth).
    A2 = A @ A
    A4 = A2 @ A2
    n = G.number_of_nodes()
    trace_A4 = np.trace(A4)
    # Number of closed 4-walks = trace_A4
    # Each 4-cycle contributes 8 closed walks (2 directions * 4 starting positions)
    # Subtract: 4-walks that backtrack: for each vertex v with degree d_v,
    # closed walks v-a-v-b-v contribute d_v * (2*d_v - 1) walks... cleaner:
    # standard formula: # 4-cycles = (trace(A^4) - 2*Σ d_v(d_v-1) - 2*|E|) / 8
    degs = np.array([d for _, d in G.degree()])
    m = G.number_of_edges()
    correction = 2 * np.sum(degs * (degs - 1)) + 2 * m
    cycles = (trace_A4 - correction) / 8.0
    return float(max(0.0, cycles))

def cut_vertex_count(G):
    return float(len(list(nx.articulation_points(G))))


# ---------- 5. Eccentricity-based ----------

def total_eccentricity(G):
    return float(sum(nx.eccentricity(G).values()))

def radius(G):
    return float(nx.radius(G))

def mean_eccentricity(G):
    eccs = list(nx.eccentricity(G).values())
    return float(np.mean(eccs))


# ---------- registry ----------

CANDIDATE_INDICES = {
    # information-theoretic
    "InfoH_deg":  shannon_degree_entropy,
    "InfoH_dist": shannon_distance_entropy,
    "InfoH_spec": spectral_entropy,
    "InfoH_ecc":  eccentricity_entropy,
    "Bonchev_Id": bonchev_distance_information,
    # centrality
    "Btw_sum":    betweenness_sum,
    "Cls_sum":    closeness_sum,
    "Eig_sum":    eigenvector_centrality_sum,
    "Harm_cent":  harmonic_centrality_sum,
    "Btw_var":    betweenness_variance,
    # spectral hybrids
    "AlgConn":    algebraic_connectivity,
    "Rho_x_Dia":  spectral_times_diameter,
    "EE_x_W":     estrada_times_wiener,
    # motif / structural
    "Triangles":  triangle_count,
    "MeanClust":  mean_clustering,
    "FourCyc":    four_cycle_count,
    "CutVerts":   cut_vertex_count,
    # eccentricity
    "TotEcc":     total_eccentricity,
    "Radius":     radius,
    "MeanEcc":    mean_eccentricity,
}


def _test():
    P5 = nx.path_graph(5)
    K4 = nx.complete_graph(4)
    C6 = nx.cycle_graph(6)
    for name, fn in CANDIDATE_INDICES.items():
        for G in (P5, K4, C6):
            v = fn(G)
            assert isinstance(v, float), f"{name} on {G}: not float ({type(v)})"
            assert not math.isnan(v), f"{name} on {G}: NaN"
    print(f"novel_candidates self-test passed. {len(CANDIDATE_INDICES)} candidates.")


if __name__ == "__main__":
    _test()
