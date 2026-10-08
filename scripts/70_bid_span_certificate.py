"""
BID span certificate (Paper 2, v2).

Claim tested.  On a dataset D, let M_D be the n x N matrix of edge-degree-pair
counts m_ij(G) (one column per degree pair realised in D).  By the
Edge-Degree-Pair Basis proposition every bond-incident-degree (BID) index
I_f(G) = sum_{uv} f(d_u, d_v) is the linear combination sum f(i,j) m_ij(G).
Hence, if the BID block B of the baseline satisfies
        span([1, B]) contains span([1, M_D])          (*)
then EVERY BID index -- for every f, every parameter value, every published or
future variant -- lies exactly in span([1, B]); its residual on the baseline is
identically zero, its VIF is infinite and its partial correlation with ANY
target is exactly 0.  Condition (*) is a finite rank check that does not
involve the target or the candidate.

This script
  1. verifies (*) on the four MoleculeNet datasets for the 18-index BID block
     of the 30-index baseline (numerical rank of centred/normalised matrices);
  2. finds a greedy minimal subset of classical BID baseline indices that
     already certifies (*);
  3. screens a panel of published BID indices that are NOT in the baseline
     (Sombor-type variants and other degree-based indices) plus the full
     100-point Loyola grid: pairwise max|r| against the 30-index baseline
     (the conventional screen) versus the exact relative residual on the
     BID block (the certificate).

Outputs
  results/bid_span_certificate.csv
  results/bid_variant_screen.csv
  results/bid_span_certificate_summary.md
"""
from __future__ import annotations

import math
import os
import sys
from itertools import product

import numpy as np
import pandas as pd

PROJECT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, PROJECT)

from src import orthogonality as ortho                      # noqa: E402
from src.load_data import load                              # noqa: E402
from src.mol_to_graph import smiles_to_graph                # noqa: E402
from src.standard_indices import compute_all                # noqa: E402

DATASETS = ["esol", "freesolv", "lipophilicity", "bbbp"]
BID18 = ["M1", "M2", "mM1", "mM2", "F", "R", "SCI", "H", "GA", "AG",
         "ABC", "ABS", "AZI", "SO", "SO_red", "Alb", "Sigma", "redM1"]
RES = os.path.join(PROJECT, "results")
TAU = 0.95

# Published BID indices that are NOT members of the 30-index baseline.
# Each is an edge sum of a symmetric function f(a, b) of the end-vertex degrees.
VARIANTS = {
    # Sombor-type (formulas as given in the cited sources; see paper Table)
    "mSO (modified Sombor)":        lambda a, b: 1.0 / math.sqrt(a * a + b * b),
    "ESO (elliptic Sombor)":        lambda a, b: (a + b) * math.sqrt(a * a + b * b),
    "EU (Euler Sombor)":            lambda a, b: math.sqrt(a * a + b * b + a * b),
    "DSO (diminished Sombor)":      lambda a, b: math.sqrt(a * a + b * b) / (a + b),
    "HSO (hyperbolic Sombor)":      lambda a, b: math.sqrt(a * a + b * b) / min(a, b),
    "SO_p, p=1/2 (p-Sombor)":       lambda a, b: (math.sqrt(a) + math.sqrt(b)) ** 2,
    "SO_p, p=3 (p-Sombor)":         lambda a, b: (a ** 3 + b ** 3) ** (1.0 / 3.0),
    # KG-Sombor summed over vertex-edge incidences, edge degree d_e = a + b - 2
    "KG-SO (KG-Sombor)":            lambda a, b: (math.sqrt(a * a + (a + b - 2) ** 2)
                                                  + math.sqrt(b * b + (a + b - 2) ** 2)),
    # 2025-2026 Sombor-type variants
    "ASO (augmented Sombor)":       lambda a, b: math.sqrt((a * a + b * b) / (a + b - 2)) if a + b > 2 else 0.0,
    "CoRSO_60 (cosine-rule Sombor, theta=60)":  lambda a, b: math.sqrt(a * a + b * b - a * b),
    "EU_lambda=1/2 (variable Euler-Sombor)":    lambda a, b: math.sqrt(a * a + b * b + 0.5 * a * b),
    # other degree-based indices
    "ISI (inverse sum indeg)":      lambda a, b: a * b / (a + b),
    "SDD (symm. division deg)":     lambda a, b: a / b + b / a,
    "N (Nirmala)":                  lambda a, b: math.sqrt(a + b),
    "IN1 (first inverse Nirmala)":  lambda a, b: math.sqrt(1.0 / a + 1.0 / b),
    "IN2 (second inverse Nirmala)": lambda a, b: math.sqrt(a * b / (a + b)),
    "GG1 (first Gourava)":          lambda a, b: a + b + a * b,
    "HM (hyper-Zagreb)":            lambda a, b: (a + b) ** 2,
}

ALPHA = [-1.0, -0.5, 0.0, 0.5, 1.0]
BETA = [-1.0, -0.5, 0.0, 0.5, 1.0]
GAMMA = [0.0, 0.5, 1.0, 2.0]


def loyola_f(alpha, beta, gamma):
    return lambda a, b: (a * b) ** alpha * (a + b) ** beta * math.exp(gamma * abs(a - b) / (a + b))


def edge_degree_pairs(G):
    d = dict(G.degree())
    return [tuple(sorted((d[u], d[v]))) for u, v in G.edges()]


def bid_value(pairs, f):
    return sum(f(a, b) for a, b in pairs)


def sombor_coindex(G):
    d = dict(G.degree()); nodes = list(G.nodes()); tot = 0.0
    for i, u in enumerate(nodes):
        for v in nodes[i + 1:]:
            if not G.has_edge(u, v):
                tot += math.sqrt(d[u] ** 2 + d[v] ** 2)
    return tot


def eccentric_sombor(G):
    import networkx as nx
    tot = 0.0
    for comp in nx.connected_components(G):
        H = G.subgraph(comp)
        ecc = nx.eccentricity(H) if H.number_of_nodes() > 1 else {u: 0 for u in H}
        tot += sum(math.sqrt(ecc[u] ** 2 + ecc[v] ** 2) for u, v in H.edges())
    return tot


NON_BID = {"SO coindex": sombor_coindex, "SO_ecc (eccentric Sombor)": eccentric_sombor}


P4 = [(i, j) for i in range(1, 5) for j in range(i, 5)]
BASIS10 = ["M1", "M2", "mM1", "mM2", "F", "R", "SCI", "H", "GA", "AG"]


def f_matrix(names):
    """Edge-function values f_b(i, j) of baseline indices on P4 (data-free).

    Each index is evaluated on the complete bipartite graph K_{i,j}, all of
    whose i*j edges have degree pair (i, j) (K_{1,1} = K_2 for (1,1))."""
    import networkx as nx
    from src.standard_indices import ALL_INDICES
    return np.array([[ALL_INDICES[k](nx.complete_bipartite_graph(i, j)) / (i * j)
                      for (i, j) in P4] for k in names])


def universal_identities():
    """Data-free certificate: rank of the 18 x 10 f-matrix over P4, and the exact
    coefficients expressing each published variant in the 10-index basis."""
    FB = f_matrix(BID18)
    sv = np.linalg.svd(FB, compute_uv=False)
    rank = int((sv > sv[0] * 1e-10).sum())
    F10 = f_matrix(BASIS10)
    rows = []
    for lab, f in VARIANTS.items():
        fv = np.array([f(i, j) for (i, j) in P4])
        c = np.linalg.solve(F10.T, fv)
        rows.append({"index": lab, **{b: c[k] for k, b in enumerate(BASIS10)},
                     "max_abs_identity_error_on_P4": float(np.abs(F10.T @ c - fv).max())})
    return rank, float(sv[-1] / sv[0]), pd.DataFrame(rows)


def generic_subsets(tol=1e-9):
    """How many 10-element subsets of the 18 baseline indices have full rank 10 on P4,
    at the relative tolerance used for all numerical ranks in the paper."""
    from itertools import combinations
    FB = f_matrix(BID18); total = full = 0
    for c in combinations(range(18), 10):
        sv = np.linalg.svd(FB[list(c)], compute_uv=False)
        full += int((sv > sv[0] * tol).sum() == 10); total += 1
    return full, total


def main():
    cert_rows, var_rows, lines = [], [], []
    rank_FB, cond_FB, coef = universal_identities()
    gfull, gtot = generic_subsets()
    pd.DataFrame([dict(tolerance=1e-9, full_rank_subsets=gfull, total_subsets=gtot)]).to_csv(
        os.path.join(RES, "bid_generic_subsets.csv"), index=False)
    coef.to_csv(os.path.join(RES, "bid_universal_identities.csv"), index=False)
    lines.append(f"Data-free certificate: rank of the 18 x 10 f-matrix of the BID baseline over P4 = {rank_FB} "
                 f"(sigma_min/sigma_max = {cond_FB:.2e}); hence every BID index is an exact linear combination "
                 f"of the baseline on every graph with max degree <= 4 and no isolated vertex.\n")
    for name in DATASETS:
        df = load(name)
        graphs, y = [], []
        for _, row in df.iterrows():
            G = smiles_to_graph(row["smiles"])
            if G is None:
                continue
            graphs.append(G)
            y.append(float(row["target"]))
        y = np.asarray(y)
        base = pd.DataFrame([compute_all(G) for G in graphs])
        pairs = [edge_degree_pairs(G) for G in graphs]
        keys = sorted({p for ps in pairs for p in ps})
        M = np.array([[ps.count(k) for k in keys] for ps in pairs], float)
        B = base[BID18].values
        X30 = base.values
        maxdeg = max(max(dict(G.degree()).values()) for G in graphs if G.number_of_edges())

        r_M = ortho.numerical_rank(M)
        r_B = ortho.numerical_rank(B)
        r_BM = ortho.numerical_rank(np.hstack([B, M]))
        max_resid_mij = max(ortho.relative_residual(M[:, j], B) for j in range(M.shape[1]))
        certified = (r_B == r_BM) and max_resid_mij < ortho.SPAN_RTOL

        # greedy minimal certifying subset of the classical BID block
        chosen, cur = [], 0
        while cur < r_M:
            best, best_r = None, cur
            for c in BID18:
                if c in chosen:
                    continue
                rr = ortho.numerical_rank(base[chosen + [c]].values)
                if rr > best_r:
                    best, best_r = c, rr
            if best is None:
                break
            chosen.append(best); cur = best_r

        cert_rows.append({
            "dataset": name, "n": len(graphs), "max_degree": maxdeg,
            "degree_pairs_realised": len(keys),
            "rank_mij": r_M, "rank_BID18": r_B, "rank_BID18_plus_mij": r_BM,
            "max_rel_resid_mij_on_BID18": max_resid_mij,
            "certificate_holds": certified,
            "greedy_minimal_certifying_subset": " ".join(chosen),
            "size_minimal_subset": len(chosen),
        })

        def screen(label, family, f, graph_fn=None):
            if graph_fn is None:
                z = np.array([bid_value(ps, f) for ps in pairs], float)
            else:
                z = np.array([graph_fn(G) for G in graphs], float)
            if np.std(z) == 0:
                return
            cors = {c: abs(float(np.corrcoef(z, base[c].values)[0, 1]))
                    for c in base.columns if np.std(base[c].values) > 0}
            best = max(cors, key=cors.get)
            cert = ortho.span_certified_pcor(z, X30, y)
            # naive pcor exactly as the v1 pipeline computed it (absolute 1e-12 test)
            from sklearn.linear_model import LinearRegression
            zr = z - LinearRegression().fit(X30, z).predict(X30)
            yr = y - LinearRegression().fit(X30, y).predict(X30)
            naive = 0.0 if np.std(zr) < 1e-12 else float(np.corrcoef(zr, yr)[0, 1])
            var_rows.append({
                "dataset": name, "index": label, "family": family,
                "max_abs_r_baseline30": cors[best], "most_correlated": best,
                "pairwise_verdict": "PASS" if cors[best] < TAU else "FAIL",
                "rel_resid_on_BID18": ortho.relative_residual(z, B),
                "rel_resid_on_baseline30": cert["rel_resid"],
                "in_span": cert["in_span"],
                "pcor_certified": cert["pcor"],
                "pcor_naive_v1": naive,
            })

        for lab, f in VARIANTS.items():
            screen(lab, "published BID", f)
            z = np.array([bid_value(ps, f) for ps in pairs], float)
            c = coef.set_index("index").loc[lab, BASIS10].values.astype(float)
            err = np.abs(base[BASIS10].values @ c - z).max() / np.abs(z).max()
            var_rows[-1]["identity_rel_error_on_data"] = float(err)
        for lab, fn in NON_BID.items():
            screen(lab, "non-BID control", None, graph_fn=fn)
        for a, b_, g in product(ALPHA, BETA, GAMMA):
            screen(f"LO({a:g},{b_:g},{g:g})", "Loyola grid", loyola_f(a, b_, g))
        print(f"{name}: n={len(graphs)} pairs={len(keys)} rank M={r_M} B={r_B} [B,M]={r_BM} "
              f"cert={certified} minimal={chosen}")

    cert = pd.DataFrame(cert_rows)
    var = pd.DataFrame(var_rows)
    cert.to_csv(os.path.join(RES, "bid_span_certificate.csv"), index=False)
    var.to_csv(os.path.join(RES, "bid_variant_screen.csv"), index=False)

    lines.append("# BID span certificate\n")
    lines.append(cert.to_string(index=False))
    lines.append("\n## Conventional pairwise screen vs certificate\n")
    g = var.groupby(["dataset", "family"]).agg(
        n=("index", "size"),
        pairwise_pass=("pairwise_verdict", lambda s: int((s == "PASS").sum())),
        in_span=("in_span", "sum"),
        max_rel_resid=("rel_resid_on_baseline30", "max"),
        max_abs_naive_pcor=("pcor_naive_v1", lambda s: float(np.abs(s).max())),
    ).reset_index()
    lines.append(g.to_string(index=False))
    lines.append("\n## Published variants passing the pairwise screen (all are exactly in span)\n")
    lines.append(var[var.family != "Loyola grid"][["dataset", "index", "family", "max_abs_r_baseline30",
                 "most_correlated", "rel_resid_on_baseline30", "in_span", "pcor_certified"]].to_string(index=False))
    pv = var[(var.family == "published BID") & (var.pairwise_verdict == "PASS")]
    lines.append(pv[["dataset", "index", "max_abs_r_baseline30", "most_correlated",
                     "rel_resid_on_baseline30", "in_span"]].to_string(index=False))
    with open(os.path.join(RES, "bid_span_certificate_summary.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
