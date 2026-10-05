"""One-function span certificate for bond-incident-degree (BID) indices.

    >>> from src.certify import certify
    >>> report = certify(smiles_list, f=lambda a, b: (a*a + b*b + a*b) ** 0.5)   # Euler Sombor
    >>> report["in_span"], report["coefficients"]

A BID index is I_f(G) = sum_{uv in E(G)} f(d_u, d_v).  On graphs of maximum
degree at most Delta it is the linear combination sum_{i<=j} f(i,j) m_ij(G)
of the edge-degree-pair counts m_ij (Gutman 2013; Rada & Cruz 2014; Rada
2026).  ``certify`` turns this into a checkable statement about a baseline:

* data-free: evaluate the baseline indices on the complete bipartite graphs
  K_{i,j}; if the resulting value matrix Phi_B has full rank over the degree
  pairs P_Delta, EVERY BID index is identically a linear combination of the
  baseline on the whole class, and the coefficient vector is returned;
* on data: the realised pair set P_D, the ranks of the count matrix M and of
  the baseline block B, and the relative residual of the candidate on the
  baseline (zero to machine precision when certified).

The default baseline is the 18 degree-based indices of the paper.
"""
from __future__ import annotations

from collections import Counter
from typing import Callable, Iterable, Sequence

import networkx as nx
import numpy as np

from . import orthogonality as ortho
from .standard_indices import ALL_INDICES

BID18 = ["M1", "M2", "mM1", "mM2", "F", "R", "SCI", "H", "GA", "AG",
         "ABC", "ABS", "AZI", "SO", "SO_red", "Alb", "Sigma", "redM1"]


def degree_pairs(delta: int):
    return [(i, j) for i in range(1, delta + 1) for j in range(i, delta + 1)]


def phi_matrix(baseline: Sequence[str], delta: int = 4) -> np.ndarray:
    """Value matrix Phi_B[b, (i,j)] = f_b(i,j), from K_{i,j} (all ij edges of type (i,j))."""
    return np.array([[ALL_INDICES[b](nx.complete_bipartite_graph(i, j)) / (i * j)
                      for (i, j) in degree_pairs(delta)] for b in baseline])


def edge_pair_counts(G: nx.Graph) -> Counter:
    d = dict(G.degree())
    return Counter(tuple(sorted((d[u], d[v]))) for u, v in G.edges())


def certify(graphs_or_smiles: Iterable, f: Callable[[int, int], float] | None = None,
            baseline: Sequence[str] = BID18, delta: int = 4, rtol: float = ortho.SPAN_RTOL) -> dict:
    """Span certificate of ``baseline`` for BID indices, optionally for a candidate f.

    ``graphs_or_smiles``: NetworkX graphs or SMILES strings (RDKit needed for SMILES).
    Returns a dict with keys
      data_free_rank, data_free_certified (rank Phi_B == |P_delta|),
      n_graphs, max_degree, realised_pairs, rank_M, rank_B, rank_BM, certified,
      max_rel_resid_mij  (largest relative residual of an m_ij column on B),
    and, if f is given,
      coefficients (dict baseline index -> coefficient, exact on the class when
      data_free_certified), identity_max_error (on the graphs), rel_resid, in_span.
    """
    graphs = []
    for g in graphs_or_smiles:
        if isinstance(g, str):
            from .mol_to_graph import smiles_to_graph
            g = smiles_to_graph(g)
        if g is not None and g.number_of_edges() > 0:
            graphs.append(g)
    P = degree_pairs(delta)
    Phi = phi_matrix(baseline, delta)
    rank_phi = int(np.linalg.matrix_rank(Phi, tol=1e-9))
    out = dict(data_free_rank=rank_phi, data_free_certified=rank_phi == len(P), n_graphs=len(graphs))
    counts = [edge_pair_counts(g) for g in graphs]
    keys = sorted({k for c in counts for k in c})
    out["max_degree"] = max((k[1] for k in keys), default=0)
    out["realised_pairs"] = keys
    M = np.array([[c.get(k, 0) for k in keys] for c in counts], float)
    B = np.array([[ALL_INDICES[b](g) for b in baseline] for g in graphs])
    rM, rB, rBM = ortho.numerical_rank(M), ortho.numerical_rank(B), ortho.numerical_rank(np.hstack([B, M]))
    res = max((ortho.relative_residual(M[:, j], B) for j in range(M.shape[1])), default=0.0)
    out.update(rank_M=rM, rank_B=rB, rank_BM=rBM, certified=bool(rB == rBM and res < rtol), max_rel_resid_mij=res)
    if f is not None:
        fv = np.array([f(i, j) for (i, j) in P], float)
        c, *_ = np.linalg.lstsq(Phi.T, fv, rcond=None)
        z = np.array([sum(f(a, b) * n for (a, b), n in cnt.items()) for cnt in counts])
        approx = B @ c
        out["coefficients"] = {b: float(c[k]) for k, b in enumerate(baseline)}
        out["identity_max_error"] = float(np.max(np.abs(approx - z)) / max(np.max(np.abs(z)), 1e-300))
        out["rel_resid"] = ortho.relative_residual(z, B)
        out["in_span"] = out["rel_resid"] < rtol
    return out


if __name__ == "__main__":  # worked example: Euler Sombor index on a few molecules
    import math
    smi = ["CCO", "c1ccccc1", "CC(C)(C)O", "C1CCCCC1", "CC(=O)OC1=CC=CC=C1C(=O)O"]
    rep = certify(smi, f=lambda a, b: math.sqrt(a * a + b * b + a * b))
    print({k: v for k, v in rep.items() if k != "coefficients"})
    print({k: round(v, 4) for k, v in rep["coefficients"].items() if abs(v) > 1e-9})
