"""
Null C for the non-linear conditional statistic (Paper 2, v2.2).

Null B (scripts/80) draws columns that are LINEAR functions of the edge-type
counts and therefore of the baseline X.  A candidate that is a deterministic
NON-LINEAR function of X (the degree entropy is one: the degree counts n_k are
linear in the m_ij, and the entropy is a function of the n_k) carries no
information beyond X either, yet a random forest residualizes it only
approximately, which inflates rho_nl.  Null C draws such columns:

    u = M c  (c ~ N(0, I), standardised),   z = phi(u),
    phi in {u^2, |u|^(1/2), log(1+|u|), tanh(u), u^3},

plus a Null D of entropy-type functions of the degree-count vector
(n_1..n_4)/n with random weights, which is the exact class of the degree
entropy.  Both use the residualization of scripts/73 with 100 trees, like
Null B.  Empirical p-values for the 20 candidates against Null C and D.

Outputs: results/nonlinear_nullC_draws.csv, _summary.csv, _candidates.csv
"""
from __future__ import annotations

import importlib.util
import os
import sys
import time
from collections import Counter

import numpy as np
import pandas as pd

PROJECT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, PROJECT)
_s = importlib.util.spec_from_file_location("m73", os.path.join(PROJECT, "scripts", "73_nonlinear_conditional.py"))
m73 = importlib.util.module_from_spec(_s); _s.loader.exec_module(m73)
_c = importlib.util.spec_from_file_location("cert", os.path.join(PROJECT, "scripts", "70_bid_span_certificate.py"))
cert = importlib.util.module_from_spec(_c); _c.loader.exec_module(cert)
from src.load_data import load                              # noqa: E402
from src.mol_to_graph import smiles_to_graph                # noqa: E402
from src.novel_candidates import CANDIDATE_INDICES          # noqa: E402
from src.standard_indices import compute_all                # noqa: E402

SEED = 20261005
K_C, K_D = 200, 200
m73.N_TREES = 100
RES = os.path.join(PROJECT, "results")
PHIS = {"square": lambda u: u ** 2, "sqrtabs": lambda u: np.sqrt(np.abs(u)), "log1pabs": lambda u: np.log1p(np.abs(u)),
        "tanh": np.tanh, "cube": lambda u: u ** 3}


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); q = np.empty(n); prev = 1.0
    for rank, i in zip(range(n, 0, -1), o[::-1]):
        prev = min(prev, p[i] * n / rank); q[i] = prev
    return q


def main():
    committed = pd.read_csv(os.path.join(RES, "nonlinear_noise_floor_candidates.csv"))
    draws, summ, crows = [], [], []
    for name in m73.DATASETS:
        t0 = time.time()
        df = load(name); graphs, ys = [], []
        for s, t in zip(df["smiles"], df["target"]):
            G = smiles_to_graph(s)
            if G is not None:
                graphs.append(G); ys.append(t)
        y = np.asarray(ys, float); n = len(y)
        X = np.column_stack([m73.impute(c) for c in pd.DataFrame([compute_all(G) for G in graphs]).values.T])
        Xs = (X - X.mean(0)) / np.where(X.std(0) > 0, X.std(0), 1.0)
        pairs = [cert.edge_degree_pairs(G) for G in graphs]
        M = np.array([[ps.count(k) for k in cert.P4] for ps in pairs], float)
        N = np.array([[Counter(d for _, d in G.degree()).get(k, 0) for k in (1, 2, 3, 4)] for G in graphs], float)
        P = N / N.sum(1, keepdims=True)
        ry = y - m73.rf_oof(Xs, y, SEED)
        rng = np.random.default_rng(SEED)

        def stat(z):
            zs = (z - z.mean()) / z.std()
            return m73.pearson(zs - m73.rf_oof(Xs, zs, SEED), ry)

        for k in range(K_C):
            u = M @ rng.standard_normal(M.shape[1])
            while u.std() == 0:
                u = M @ rng.standard_normal(M.shape[1])
            u = (u - u.mean()) / u.std()
            phi = list(PHIS)[k % len(PHIS)]
            z = PHIS[phi](u)
            if z.std() == 0:
                continue
            draws.append(dict(dataset=name, null="C_nonlinear_inspan", draw=k, phi=phi, rho_nl=stat(z)))
        for k in range(K_D):
            w = rng.standard_normal(4); eps = 10 ** rng.uniform(-3, -1)
            z = -(P * np.log(P + eps) * w).sum(1)          # weighted entropy-type function of the degree distribution
            if z.std() == 0:
                continue
            draws.append(dict(dataset=name, null="D_entropy_type", draw=k, phi="entropy", rho_nl=stat(z)))
        D = pd.DataFrame(draws)
        for nl in ("C_nonlinear_inspan", "D_entropy_type"):
            v = D[(D.dataset == name) & (D.null == nl)].rho_nl.values
            summ.append(dict(dataset=name, null=nl, K=len(v), sd=v.std(ddof=1), p97_5=np.percentile(v, 97.5),
                             p99=np.percentile(v, 99), p99_abs=np.percentile(np.abs(v), 99), max_abs=np.abs(v).max()))
        nullC = np.abs(D[(D.dataset == name) & (D.null == "C_nonlinear_inspan")].rho_nl.values)
        nullD = np.abs(D[(D.dataset == name) & (D.null == "D_entropy_type")].rho_nl.values)
        sub = committed[committed.dataset == name]
        for _, r in sub.iterrows():
            st = r.stat_used
            crows.append(dict(dataset=name, candidate=r.candidate, stat_used=st,
                              p_emp_nullC=(np.sum(nullC >= abs(st)) + 1) / (len(nullC) + 1),
                              p_emp_nullD=(np.sum(nullD >= abs(st)) + 1) / (len(nullD) + 1)))
        print(f"[{name}] done {time.time()-t0:.0f}s; nullC p99 |rho| {np.percentile(nullC,99):.3f}; nullD p99 {np.percentile(nullD,99):.3f}", flush=True)
    C = pd.DataFrame(crows); C["q_bh_nullC"] = bh(C.p_emp_nullC.values); C["q_bh_nullD"] = bh(C.p_emp_nullD.values)
    pd.DataFrame(draws).to_csv(os.path.join(RES, "nonlinear_nullC_draws.csv"), index=False)
    pd.DataFrame(summ).to_csv(os.path.join(RES, "nonlinear_nullC_summary.csv"), index=False)
    C.to_csv(os.path.join(RES, "nonlinear_nullC_candidates.csv"), index=False)
    print(pd.DataFrame(summ).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
