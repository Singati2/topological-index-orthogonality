"""Closed-form values of the Wuzi parametric index family on standard graph
classes, following the structural template of Proposition 1 in Movahedi,
Gutman, Redžepović & Furtula, MATCH Commun. Math. Comput. Chem. 95:141-162
(2026) — the origin paper for the Diminished Sombor index.

The Wuzi family is defined for a simple graph G = (V, E) and parameters
α, β, γ ∈ ℝ by

    W(G; α, β, γ) = Σ_{uv ∈ E(G)} (d_u · d_v)^α · (d_u + d_v)^β
                                  · exp( γ · |d_u - d_v| / (d_u + d_v) )

Each function in this module returns the closed-form W on a named graph
family. The companion `wuzi_numerical_value(G, α, β, γ)` (from
`src/wuzi_index.py`) lets us verify the closed forms by direct
computation on a NetworkX instance.

All formulas are derived under the convention that the graph is simple,
connected (unless explicitly stated), and contains at least one edge.
"""
from __future__ import annotations
import math


# ----- Regular graphs -----

def W_regular(n: int, k: int, alpha: float, beta: float, gamma: float) -> float:
    """W on a k-regular graph of order n (k ≥ 1, n ≥ k+1, nk even).

    Every edge contributes (k·k)^α · (2k)^β · exp(γ·0) = k^{2α+β} · 2^β.
    Total edges: m = n k / 2.
    """
    m = n * k / 2.0
    return m * (k ** (2 * alpha + beta)) * (2.0 ** beta)


# ----- Complete graph K_n  -----

def W_complete(n: int, alpha: float, beta: float, gamma: float) -> float:
    """W(K_n) — every vertex has degree n-1; (n-1)-regular."""
    if n < 2:
        return 0.0
    return W_regular(n, n - 1, alpha, beta, gamma)


# ----- Cycle C_n  -----

def W_cycle(n: int, alpha: float, beta: float, gamma: float) -> float:
    """W(C_n) — every vertex has degree 2; 2-regular."""
    if n < 3:
        return 0.0
    return W_regular(n, 2, alpha, beta, gamma)


# ----- Path P_n  -----

def W_path(n: int, alpha: float, beta: float, gamma: float) -> float:
    """W(P_n) for n ≥ 2.

    Two pendant edges of type (1, 2) and (n - 3) internal edges of type (2, 2).
        (1, 2)-edge contribution:  (1·2)^α · 3^β · exp(γ · 1/3) = 2^α · 3^β · exp(γ/3)
        (2, 2)-edge contribution:  4^α · 4^β · exp(0)           = 2^{2α + 2β}
    Special case n = 2: a single (1, 1)-edge of contribution 1 · 2^β · 1 = 2^β.
    Special case n = 3: two (1, 2)-edges, no interior.
    """
    if n < 2:
        return 0.0
    if n == 2:
        return (1.0 ** alpha) * (2.0 ** beta) * math.exp(gamma * 0.0)
    pendant = 2.0 * (2.0 ** alpha) * (3.0 ** beta) * math.exp(gamma / 3.0)
    internal = max(n - 3, 0) * (2.0 ** (2 * alpha + 2 * beta))
    return pendant + internal


# ----- Star S_n  -----

def W_star(n: int, alpha: float, beta: float, gamma: float) -> float:
    """W(S_n) where S_n is K_{1, n-1} — center has degree n-1, n-1 leaves of degree 1.

    All n-1 edges are of type (1, n-1):
        contribution per edge = (n-1)^α · n^β · exp(γ · (n-2)/n)
    """
    if n < 2:
        return 0.0
    per_edge = ((n - 1) ** alpha) * (n ** beta) * math.exp(gamma * (n - 2) / n)
    return (n - 1) * per_edge


# ----- Complete bipartite K_{p,q}  -----

def W_complete_bipartite(p: int, q: int, alpha: float, beta: float, gamma: float) -> float:
    """W(K_{p,q}) — p·q edges, every edge has endpoint degrees (q, p) (i.e. {p, q})."""
    if p < 1 or q < 1:
        return 0.0
    per_edge = ((p * q) ** alpha) * ((p + q) ** beta) * math.exp(gamma * abs(p - q) / (p + q))
    return p * q * per_edge


# ----- Hypercube Q_k  -----

def W_hypercube(k: int, alpha: float, beta: float, gamma: float) -> float:
    """W(Q_k) — Q_k is k-regular on 2^k vertices, so |E| = k · 2^{k-1}."""
    if k < 1:
        return 0.0
    n = 2 ** k
    return W_regular(n, k, alpha, beta, gamma)


# ----- Wheel W_n  -----

def W_wheel(n: int, alpha: float, beta: float, gamma: float) -> float:
    """W(W_n) where W_n is C_{n-1} plus one hub connected to all n-1 rim vertices.

    Hub has degree n-1; rim vertices have degree 3.
    Edges:
      - n-1 spokes:  endpoint degrees (3, n-1)
      - n-1 rim edges: endpoint degrees (3, 3)
    """
    if n < 4:
        return 0.0
    spoke = (3 * (n - 1)) ** alpha
    spoke *= (3 + n - 1) ** beta
    spoke *= math.exp(gamma * abs(3 - (n - 1)) / (3 + n - 1))
    rim = (9.0) ** alpha * (6.0) ** beta * math.exp(0.0)
    return (n - 1) * spoke + (n - 1) * rim


# ----- Friendship graph F_p (p triangles sharing one vertex)  -----

def W_friendship(p: int, alpha: float, beta: float, gamma: float) -> float:
    """F_p: p triangles glued at a common vertex.

    Vertices: 2p outer vertices of degree 2; central vertex of degree 2p.
    Edges (3p in total):
      - p outer-outer edges: (2, 2)
      - 2p central-outer edges: (2, 2p)
    """
    if p < 1:
        return 0.0
    outer_outer = p * (4.0 ** alpha) * (4.0 ** beta) * math.exp(0.0)
    central_outer = 2 * p * ((2 * 2 * p) ** alpha) * ((2 + 2 * p) ** beta) \
                          * math.exp(gamma * abs(2 - 2 * p) / (2 + 2 * p))
    return outer_outer + central_outer


# ----- Verification ---------------------------------------------------------

def _verify():
    import networkx as nx
    from src.wuzi_index import wuzi

    test_params = [
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),     # = M_2
        (-0.5, 0.0, 0.0),    # = Randić
        (0.0, -0.5, 0.0),    # = sum-connectivity SCI
        (0.0, -1.0, 0.0),    # = harmonic / 2
        (0.5, 0.3, 1.0),     # generic
    ]

    # K_5
    for a, b, g in test_params:
        G = nx.complete_graph(5)
        lhs = wuzi(G, a, b, g)
        rhs = W_complete(5, a, b, g)
        assert abs(lhs - rhs) < 1e-9, f"K_5 ({a},{b},{g}): {lhs} vs {rhs}"

    # C_6
    for a, b, g in test_params:
        G = nx.cycle_graph(6)
        lhs = wuzi(G, a, b, g)
        rhs = W_cycle(6, a, b, g)
        assert abs(lhs - rhs) < 1e-9, f"C_6: {lhs} vs {rhs}"

    # P_5, P_7
    for n in (5, 7, 10):
        for a, b, g in test_params:
            G = nx.path_graph(n)
            lhs = wuzi(G, a, b, g)
            rhs = W_path(n, a, b, g)
            assert abs(lhs - rhs) < 1e-9, f"P_{n} ({a},{b},{g}): {lhs} vs {rhs}"

    # S_n (star)
    for n in (4, 6, 8):
        for a, b, g in test_params:
            G = nx.star_graph(n - 1)   # nx.star_graph(k) makes K_{1,k}, n=k+1
            lhs = wuzi(G, a, b, g)
            rhs = W_star(n, a, b, g)
            assert abs(lhs - rhs) < 1e-9, f"S_{n}: {lhs} vs {rhs}"

    # K_{3,4}
    for a, b, g in test_params:
        G = nx.complete_bipartite_graph(3, 4)
        lhs = wuzi(G, a, b, g)
        rhs = W_complete_bipartite(3, 4, a, b, g)
        assert abs(lhs - rhs) < 1e-9, f"K_{{3,4}}: {lhs} vs {rhs}"

    # Q_3 (hypercube)
    for a, b, g in test_params:
        G = nx.hypercube_graph(3)
        G = nx.convert_node_labels_to_integers(G)
        lhs = wuzi(G, a, b, g)
        rhs = W_hypercube(3, a, b, g)
        assert abs(lhs - rhs) < 1e-9, f"Q_3: {lhs} vs {rhs}"

    # W_5 (wheel on 5 vertices: hub + C_4)
    for a, b, g in test_params:
        G = nx.wheel_graph(5)
        lhs = wuzi(G, a, b, g)
        rhs = W_wheel(5, a, b, g)
        assert abs(lhs - rhs) < 1e-9, f"W_5 ({a},{b},{g}): {lhs} vs {rhs}"

    # F_3 (friendship)
    for a, b, g in test_params:
        # construct F_3 manually: 3 triangles sharing vertex 0
        G = __import__("networkx").Graph()
        # vertex 0 is central; triangles are (0, 1, 2), (0, 3, 4), (0, 5, 6)
        triangles = [(0, 1, 2), (0, 3, 4), (0, 5, 6)]
        for i, j, k in triangles:
            G.add_edge(i, j); G.add_edge(j, k); G.add_edge(i, k)
        lhs = wuzi(G, a, b, g)
        rhs = W_friendship(3, a, b, g)
        assert abs(lhs - rhs) < 1e-9, f"F_3 ({a},{b},{g}): {lhs} vs {rhs}"

    print("wuzi_analytical: all closed-form identities verified on K_5, C_6, P_{5,7,10}, S_{4,6,8}, K_{3,4}, Q_3, W_5, F_3.")


if __name__ == "__main__":
    import sys, os
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
    _verify()
