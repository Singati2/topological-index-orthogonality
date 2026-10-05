"""
Two checks requested by the Feynman review (Paper 2, v2).

1. Exactness of the graph grouping: every group of curated molecules sharing a
   Weisfeiler-Lehman hash is checked with an exact isomorphism test
   (networkx.is_isomorphic) against the group's first member.
2. R^2 bound for BID indices: because every BID index lies in col[1, M]
   (M = edge-degree-pair counts), the in-sample OLS R^2 of y on a single BID
   index cannot exceed the OLS R^2 of y on [1, M].  We report R^2(y ~ m_ij)
   and the largest single-index R^2 over the 18 baseline BID indices, the
   published BID indices and the Loyola grid.

Outputs: results/isomorphism_check.csv, results/bid_r2_bound.csv
"""
from __future__ import annotations

import os
import sys
from itertools import product

import networkx as nx
import numpy as np
import pandas as pd

PROJECT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, PROJECT)
sys.path.insert(0, os.path.join(PROJECT, "scripts"))

from src.load_data import load                              # noqa: E402
from src.mol_to_graph import smiles_to_graph                # noqa: E402
import importlib.util                                       # noqa: E402

spec = importlib.util.spec_from_file_location("cert", os.path.join(PROJECT, "scripts", "70_bid_span_certificate.py"))
cert = importlib.util.module_from_spec(spec); spec.loader.exec_module(cert)
RES = os.path.join(PROJECT, "results")


def r2(y, X):
    A = np.hstack([np.ones((len(y), 1)), X.reshape(len(y), -1)])
    b, *_ = np.linalg.lstsq(A, y, rcond=None)
    return 1 - np.sum((y - A @ b) ** 2) / np.sum((y - y.mean()) ** 2)


def main():
    iso_rows, r2_rows = [], []
    for name in ["esol", "freesolv", "lipophilicity", "bbbp"]:
        cur = pd.read_csv(os.path.join(PROJECT, "data", "curated", f"{name}_curated.csv"))
        graphs = [smiles_to_graph(s) for s in cur.canonical_smiles]
        groups = cur.groupby("graph_wl_hash").indices
        n_groups_multi = n_bad = 0
        for idx in groups.values():
            if len(idx) < 2:
                continue
            n_groups_multi += 1
            g0 = graphs[idx[0]]
            if not all(nx.is_isomorphic(g0, graphs[i]) for i in idx[1:]):
                n_bad += 1
        iso_rows.append(dict(dataset=name, shared_groups=n_groups_multi, non_isomorphic_groups=n_bad))

        df = load(name)
        gs, y = [], []
        for s, t in zip(df["smiles"], df["target"]):
            G = smiles_to_graph(s)
            if G is not None:
                gs.append(G); y.append(float(t))
        y = np.asarray(y)
        pairs = [cert.edge_degree_pairs(G) for G in gs]
        M = np.array([[ps.count(k) for k in cert.P4] for ps in pairs], float)
        fns = {f"baseline:{b}": None for b in cert.BID18}
        from src.standard_indices import ALL_INDICES
        best, best_name = -1, ""
        for b in cert.BID18:
            v = np.array([ALL_INDICES[b](G) for G in gs])
            r = r2(y, v)
            if r > best:
                best, best_name = r, b
        for lab, f in cert.VARIANTS.items():
            v = np.array([cert.bid_value(ps, f) for ps in pairs])
            r = r2(y, v)
            if r > best:
                best, best_name = r, lab
        for a, b_, g in product(cert.ALPHA, cert.BETA, cert.GAMMA):
            v = np.array([cert.bid_value(ps, cert.loyola_f(a, b_, g)) for ps in pairs])
            r = r2(y, v)
            if r > best:
                best, best_name = r, f"LO({a:g},{b_:g},{g:g})"
        r2_rows.append(dict(dataset=name, n=len(y), r2_ols_on_mij=r2(y, M),
                            best_single_bid_r2=best, best_single_bid_index=best_name))
        print(iso_rows[-1], r2_rows[-1])
    pd.DataFrame(iso_rows).to_csv(os.path.join(RES, "isomorphism_check.csv"), index=False)
    pd.DataFrame(r2_rows).to_csv(os.path.join(RES, "bid_r2_bound.csv"), index=False)


if __name__ == "__main__":
    main()
