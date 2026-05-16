"""30 standard topological indices on NetworkX simple graphs.

All functions take an undirected NetworkX graph G with integer node labels
(0..n-1) and return a Python float. Distance-based indices use unweighted
shortest paths. Spectral indices use the adjacency-matrix spectrum.

Each function operates on the H-suppressed molecular graph produced by
mol_to_graph.smiles_to_graph.

Indices grouped by family:
  degree-based (18):  M1, M2, mM1, mM2, F, R, SCI, H, GA, AG, ABC, ABS,
                      AZI, SO, SO_red, Alb, Sigma, redM1
  distance-based (8): W, WW, HR, Schultz, Gutman, Mostar, Sz, PI
  spectral (3):       EE, Energy, SpecRad
  other (1):          Balaban_J
"""
from __future__ import annotations
import math
import numpy as np
import networkx as nx


# ---------- helpers ----------

def _degrees(G):
    return dict(G.degree())

def _dist_matrix(G):
    # Floyd-Warshall on unweighted graph -> int distances
    D = nx.floyd_warshall_numpy(G)
    return D.astype(int)

def _adj_eigs(G):
    A = nx.adjacency_matrix(G).toarray().astype(float)
    return np.linalg.eigvalsh(A)


# ---------- degree-based ----------

def first_zagreb(G):
    return float(sum(d*d for _, d in G.degree()))

def second_zagreb(G):
    d = _degrees(G)
    return float(sum(d[u]*d[v] for u, v in G.edges()))

def modified_first_zagreb(G):
    return float(sum(1.0/(d*d) for _, d in G.degree() if d > 0))

def modified_second_zagreb(G):
    d = _degrees(G)
    return float(sum(1.0/(d[u]*d[v]) for u, v in G.edges()))

def forgotten(G):
    return float(sum(d*d*d for _, d in G.degree()))

def randic(G):
    d = _degrees(G)
    return float(sum(1.0/math.sqrt(d[u]*d[v]) for u, v in G.edges()))

def sum_connectivity(G):
    d = _degrees(G)
    return float(sum(1.0/math.sqrt(d[u]+d[v]) for u, v in G.edges()))

def harmonic_index(G):
    d = _degrees(G)
    return float(sum(2.0/(d[u]+d[v]) for u, v in G.edges()))

def geometric_arithmetic(G):
    d = _degrees(G)
    return float(sum(2.0*math.sqrt(d[u]*d[v])/(d[u]+d[v]) for u, v in G.edges()))

def arithmetic_geometric(G):
    d = _degrees(G)
    return float(sum((d[u]+d[v])/(2.0*math.sqrt(d[u]*d[v])) for u, v in G.edges()))

def atom_bond_connectivity(G):
    d = _degrees(G)
    total = 0.0
    for u, v in G.edges():
        s = (d[u]+d[v]-2.0)/(d[u]*d[v])
        if s > 0:
            total += math.sqrt(s)
    return total

def atom_bond_sum_connectivity(G):
    d = _degrees(G)
    total = 0.0
    for u, v in G.edges():
        x = 1.0 - 2.0/(d[u]+d[v])
        if x > 0:
            total += math.sqrt(x)
    return total

def augmented_zagreb(G):
    d = _degrees(G)
    total = 0.0
    for u, v in G.edges():
        denom = d[u]+d[v]-2.0
        if denom > 0:
            total += (d[u]*d[v]/denom)**3
    return total

def sombor(G):
    d = _degrees(G)
    return float(sum(math.sqrt(d[u]**2 + d[v]**2) for u, v in G.edges()))

def reduced_sombor(G):
    d = _degrees(G)
    return float(sum(math.sqrt((d[u]-1)**2 + (d[v]-1)**2) for u, v in G.edges()))

def albertson(G):
    d = _degrees(G)
    return float(sum(abs(d[u]-d[v]) for u, v in G.edges()))

def sigma_index(G):
    d = _degrees(G)
    return float(sum((d[u]-d[v])**2 for u, v in G.edges()))

def reduced_first_zagreb(G):
    return float(sum((d-1)**2 for _, d in G.degree()))


# ---------- distance-based ----------

def wiener(G):
    D = _dist_matrix(G)
    return float(D.sum() / 2.0)

def hyper_wiener(G):
    D = _dist_matrix(G).astype(float)
    n = G.number_of_nodes()
    iu = np.triu_indices(n, k=1)
    d = D[iu]
    return float(np.sum(d*d + d) / 2.0)

def harary(G):
    D = _dist_matrix(G).astype(float)
    n = G.number_of_nodes()
    iu = np.triu_indices(n, k=1)
    d = D[iu]
    return float(np.sum(1.0/d))

def schultz(G):
    d = np.array([G.degree(i) for i in range(G.number_of_nodes())], dtype=float)
    D = _dist_matrix(G).astype(float)
    s = (d[:, None] + d[None, :]) * D
    n = G.number_of_nodes()
    iu = np.triu_indices(n, k=1)
    return float(s[iu].sum())

def gutman_index(G):
    d = np.array([G.degree(i) for i in range(G.number_of_nodes())], dtype=float)
    D = _dist_matrix(G).astype(float)
    p = (d[:, None] * d[None, :]) * D
    n = G.number_of_nodes()
    iu = np.triu_indices(n, k=1)
    return float(p[iu].sum())

def mostar(G):
    D = _dist_matrix(G)
    total = 0.0
    n = G.number_of_nodes()
    for u, v in G.edges():
        # nu = vertices strictly closer to u than to v
        diff = D[u] - D[v]   # diff[w] = d(u,w) - d(v,w)
        nu = int(np.sum(diff < 0))
        nv = int(np.sum(diff > 0))
        total += abs(nu - nv)
    return total

def szeged(G):
    D = _dist_matrix(G)
    total = 0.0
    for u, v in G.edges():
        diff = D[u] - D[v]
        nu = int(np.sum(diff < 0))
        nv = int(np.sum(diff > 0))
        total += nu * nv
    return float(total)

def padmakar_ivan(G):
    """Edge-PI: PI(G) = sum over edges of (m_u(e) + m_v(e)).
    m_u(e) = number of edges parallel to e on u's side (closer to u than to v).
    """
    D = _dist_matrix(G)
    edges = list(G.edges())
    total = 0
    for u, v in edges:
        mu = mv = 0
        for a, b in edges:
            if (a == u and b == v) or (a == v and b == u):
                continue
            a_u = D[a, u] < D[a, v]
            b_u = D[b, u] < D[b, v]
            a_v = D[a, v] < D[a, u]
            b_v = D[b, v] < D[b, u]
            if a_u and b_u:
                mu += 1
            elif a_v and b_v:
                mv += 1
        total += mu + mv
    return float(total)


# ---------- spectral ----------

def estrada(G):
    eigs = _adj_eigs(G)
    return float(np.sum(np.exp(eigs)))

def graph_energy(G):
    eigs = _adj_eigs(G)
    return float(np.sum(np.abs(eigs)))

def spectral_radius(G):
    eigs = _adj_eigs(G)
    return float(np.max(np.abs(eigs)))


# ---------- other ----------

def balaban_j(G):
    n = G.number_of_nodes()
    m = G.number_of_edges()
    mu = m - n + 1   # cyclomatic number (connected graph)
    D = _dist_matrix(G).astype(float)
    s = D.sum(axis=1)
    total = 0.0
    for u, v in G.edges():
        if s[u] > 0 and s[v] > 0:
            total += 1.0 / math.sqrt(s[u] * s[v])
    return float((m / (mu + 1.0)) * total)


# ---------- registry ----------

ALL_INDICES = {
    # degree-based
    "M1": first_zagreb,
    "M2": second_zagreb,
    "mM1": modified_first_zagreb,
    "mM2": modified_second_zagreb,
    "F": forgotten,
    "R": randic,
    "SCI": sum_connectivity,
    "H": harmonic_index,
    "GA": geometric_arithmetic,
    "AG": arithmetic_geometric,
    "ABC": atom_bond_connectivity,
    "ABS": atom_bond_sum_connectivity,
    "AZI": augmented_zagreb,
    "SO": sombor,
    "SO_red": reduced_sombor,
    "Alb": albertson,
    "Sigma": sigma_index,
    "redM1": reduced_first_zagreb,
    # distance-based
    "W": wiener,
    "WW": hyper_wiener,
    "HR": harary,
    "Schultz": schultz,
    "Gutman": gutman_index,
    "Mostar": mostar,
    "Sz": szeged,
    "PI": padmakar_ivan,
    # spectral
    "EE": estrada,
    "Energy": graph_energy,
    "SpecRad": spectral_radius,
    # other
    "Balaban_J": balaban_j,
}


def compute_all(G) -> dict:
    return {name: fn(G) for name, fn in ALL_INDICES.items()}


def _test():
    import networkx as nx
    # Path P4: 4 vertices, 3 edges, degrees [1,2,2,1]
    P4 = nx.path_graph(4)
    vals = compute_all(P4)
    # Known values for P4:
    # M1 = 1+4+4+1 = 10
    # M2 = 1*2 + 2*2 + 2*1 = 8
    # Wiener = 1+2+3 + 1+2 + 1 = 10
    # Randic = 1/sqrt(2) + 1/2 + 1/sqrt(2) = 1.4142...
    assert vals["M1"] == 10.0, vals["M1"]
    assert vals["M2"] == 8.0, vals["M2"]
    assert vals["W"] == 10.0, vals["W"]
    assert abs(vals["R"] - (2/math.sqrt(2) + 0.5)) < 1e-9, vals["R"]
    print(f"  P4 spot-checks pass (M1={vals['M1']}, M2={vals['M2']}, W={vals['W']}, R={vals['R']:.4f})")

    # Cycle C6 (benzene-like graph): n=6, all degrees=2
    C6 = nx.cycle_graph(6)
    vals = compute_all(C6)
    # M1 = 6*4 = 24
    # M2 = 6*4 = 24
    # Wiener of C6: sum of distances = 6*(1+2+3+2+1)/2 = 6*9/2 = 27. Actually for C_n,
    # Wiener = n*floor(n^2/8) for even n. For n=6: 6*4=24. Let me compute by hand:
    # pairs in C6: distances {1,2,3,2,1}, {1,2,3,2}, {1,2,3}, {1,2}, {1}
    # actually for C6 all 15 pairs: each vertex has dists 1,2,3,2,1 -> sum 9 per vertex
    # total over pairs = 6*9/2 = 27
    assert vals["M1"] == 24.0
    assert vals["M2"] == 24.0
    assert vals["W"] == 27.0, vals["W"]
    # Albertson on regular graph = 0; Sigma = 0
    assert vals["Alb"] == 0.0
    assert vals["Sigma"] == 0.0
    print(f"  C6 spot-checks pass (M1={vals['M1']}, W={vals['W']}, Alb={vals['Alb']})")

    print(f"standard_indices self-test passed. {len(ALL_INDICES)} indices computed.")


if __name__ == "__main__":
    _test()
