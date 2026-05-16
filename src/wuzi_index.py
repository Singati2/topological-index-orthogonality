"""Wuzi parametric index family.

Definition:
    Wuzi(G; α, β, γ) = Σ_{uv ∈ E(G)} (d_u · d_v)^α · (d_u + d_v)^β · exp(γ · |d_u - d_v| / (d_u + d_v))

Three knobs:
    α  controls the degree-product weighting (Randić-family axis)
    β  controls the degree-sum weighting     (sum-connectivity-family axis)
    γ  controls the exponential irregularity weighting (novel)

The exponential-irregularity term exp(γ · |d_u-d_v|/(d_u+d_v)) is the
genuinely novel contribution. Most irregularity-aware indices in the
literature (Albertson, Sigma) use polynomial functions of |d_u-d_v|.
The normalized ratio |d_u-d_v|/(d_u+d_v) ∈ [0, 1) bounds the exponent.

Special cases (with γ=0):
    (α, β) = ( 1, 0)  → Second Zagreb M2
    (α, β) = (-0.5, 0) → Randić R
    (α, β) = ( 0, -0.5) → Sum-connectivity SCI
    (α, β) = ( 0, -1)  → 2 × Harmonic
    (α, β) = ( 0,  0)  → number of edges |E|
"""
from __future__ import annotations
import math
import networkx as nx


def wuzi(G: nx.Graph, alpha: float, beta: float, gamma: float) -> float:
    deg = dict(G.degree())
    total = 0.0
    for u, v in G.edges():
        du, dv = deg[u], deg[v]
        prod = du * dv
        s = du + dv
        diff_ratio = abs(du - dv) / s   # well-defined: s >= 2 always
        term = (prod ** alpha) * (s ** beta) * math.exp(gamma * diff_ratio)
        total += term
    return float(total)


def wuzi_factory(alpha: float, beta: float, gamma: float):
    """Return a closure f(G) -> Wuzi(G; α, β, γ) for use in the index registry."""
    def f(G):
        return wuzi(G, alpha, beta, gamma)
    f.__name__ = f"wuzi_a{alpha}_b{beta}_g{gamma}"
    return f


def special_case_check():
    """Verify the special-case identities to within numerical tolerance."""
    from src.standard_indices import (
        second_zagreb, randic, sum_connectivity, harmonic_index,
    )
    P5 = nx.path_graph(5)         # path, mixed degrees
    K4 = nx.complete_graph(4)     # regular degree 3
    C7 = nx.cycle_graph(7)        # regular degree 2 with odd cycle
    cases = [P5, K4, C7]
    for G in cases:
        # γ=0 zeroes the exponential, so exp factor = 1
        assert abs(wuzi(G, 1, 0, 0) - second_zagreb(G)) < 1e-9, G
        assert abs(wuzi(G, -0.5, 0, 0) - randic(G)) < 1e-9, G
        assert abs(wuzi(G, 0, -0.5, 0) - sum_connectivity(G)) < 1e-9, G
        assert abs(wuzi(G, 0, -1, 0) - harmonic_index(G) / 2.0) < 1e-9, G
        assert abs(wuzi(G, 0, 0, 0) - G.number_of_edges()) < 1e-9, G
    print(f"wuzi_index special-case identities verified on {len(cases)} test graphs.")


if __name__ == "__main__":
    import sys, os
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
    special_case_check()
