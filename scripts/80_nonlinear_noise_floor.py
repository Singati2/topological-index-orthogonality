"""
Null distribution for the non-linear conditional statistic rho_nl (Paper 2, v2.1).

scripts/73 reports rho_nl = corr(res_RF(z|X), res_RF(y|X)) with cross-fitted
random forests, calibrated against only two null references per dataset.
This script draws a proper null distribution with the SAME residualization
procedure (imported from scripts/73; only the number of trees is reduced to
N_TREES_NULL to make hundreds of draws affordable):

  Null A  K_A independent N(0,1) columns (pure noise, no X-structure);
  Null B  K_B random BID indices z = M c, c ~ N(0,1) over the ten m_ij
          counts -- exactly in the linear span of X, so any nonzero rho_nl
          measures leakage of X-structure through the imperfect RF
          residualization (the relevant null for a redundant index).

For every candidate the empirical two-sided p-value against Null B and the
Benjamini-Hochberg q across the 80 candidate-dataset pairs are reported.
Candidates with |rho_nl| >= 0.08 in the committed results are re-run with
N_TREES_NULL trees so that their statistic is like-for-like with the null.

Outputs: results/nonlinear_noise_floor_draws.csv, _summary.csv,
         _candidates.csv, _report.md
"""
from __future__ import annotations

import importlib.util
import os
import sys
import time

import numpy as np
import pandas as pd

PROJECT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, PROJECT)
_spec = importlib.util.spec_from_file_location("m73", os.path.join(PROJECT, "scripts", "73_nonlinear_conditional.py"))
m73 = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(m73)
_c = importlib.util.spec_from_file_location("cert", os.path.join(PROJECT, "scripts", "70_bid_span_certificate.py"))
cert = importlib.util.module_from_spec(_c); _c.loader.exec_module(cert)

from src.load_data import load                              # noqa: E402
from src.mol_to_graph import smiles_to_graph                # noqa: E402
from src.novel_candidates import CANDIDATE_INDICES          # noqa: E402
from src.standard_indices import compute_all                # noqa: E402

SEED = 20261005
N_TREES_NULL = 100
K_A, K_B = 100, 60
RERUN_ABS = 0.08
RES = os.path.join(PROJECT, "results")
m73.N_TREES = N_TREES_NULL   # all residualizations in this script use the reduced forest


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); q = np.empty(n)
    prev = 1.0
    for rank, i in zip(range(n, 0, -1), o[::-1]):
        prev = min(prev, p[i] * n / rank); q[i] = prev
    return q


def main():
    committed = pd.read_csv(os.path.join(RES, "nonlinear_conditional_tests.csv"))
    draws, summ, crows = [], [], []
    for name in m73.DATASETS:
        t0 = time.time()
        df = load(name)
        graphs, ys = [], []
        for s, t in zip(df["smiles"], df["target"]):
            G = smiles_to_graph(s)
            if G is None:
                continue
            graphs.append(G); ys.append(t)
        y = np.asarray(ys, float); n = len(y)
        X = np.column_stack([m73.impute(c) for c in pd.DataFrame([compute_all(G) for G in graphs]).values.T])
        Xs = (X - X.mean(0)) / np.where(X.std(0) > 0, X.std(0), 1.0)
        pairs = [cert.edge_degree_pairs(G) for G in graphs]
        M = np.array([[ps.count(k) for k in cert.P4] for ps in pairs], float)
        ry = y - m73.rf_oof(Xs, y, SEED)
        rng = np.random.default_rng(SEED)

        def stat(z):
            zs = (z - z.mean()) / z.std()
            return m73.pearson(zs - m73.rf_oof(Xs, zs, SEED), ry)

        for k in range(K_A):
            draws.append(dict(dataset=name, null="A_noise", draw=k, rho_nl=stat(rng.standard_normal(n))))
        for k in range(K_B):
            z = M @ rng.standard_normal(M.shape[1])
            while z.std() == 0:
                z = M @ rng.standard_normal(M.shape[1])
            draws.append(dict(dataset=name, null="B_inspan", draw=k, rho_nl=stat(z)))
        D = pd.DataFrame(draws)
        for nl in ["A_noise", "B_inspan"]:
            v = D[(D.dataset == name) & (D.null == nl)].rho_nl.values
            summ.append(dict(dataset=name, null=nl, K=len(v), mean=v.mean(), sd=v.std(ddof=1),
                             p2_5=np.percentile(v, 2.5), p97_5=np.percentile(v, 97.5),
                             p0_5=np.percentile(v, 0.5), p99_5=np.percentile(v, 99.5),
                             max_abs=np.abs(v).max(), p99_abs=np.percentile(np.abs(v), 99)))
        nullB = np.abs(D[(D.dataset == name) & (D.null == "B_inspan")].rho_nl.values)
        sub = committed[(committed.dataset == name) & ~committed.candidate.str.startswith("REF")]
        for _, r in sub.iterrows():
            r100 = np.nan
            if abs(r.rho_nl) >= RERUN_ABS:
                vals = []
                for G in graphs:
                    try:
                        vals.append(CANDIDATE_INDICES[r.candidate](G))
                    except Exception:
                        vals.append(np.nan)
                r100 = stat(m73.impute(vals))
            ref = r100 if np.isfinite(r100) else r.rho_nl
            p_emp = (np.sum(nullB >= abs(ref)) + 1) / (len(nullB) + 1)
            crows.append(dict(dataset=name, candidate=r.candidate, rho_nl_300=r.rho_nl, rho_nl_100=r100,
                              stat_used=ref, p_emp_nullB=p_emp))
        print(f"[{name}] n={n} done in {time.time()-t0:.0f}s; nullB 99.5th |rho| = "
              f"{np.percentile(nullB, 99.5):.3f}", flush=True)
    C = pd.DataFrame(crows); C["q_bh"] = bh(C.p_emp_nullB.values)
    S = pd.DataFrame(summ)
    pd.DataFrame(draws).to_csv(os.path.join(RES, "nonlinear_noise_floor_draws.csv"), index=False)
    S.to_csv(os.path.join(RES, "nonlinear_noise_floor_summary.csv"), index=False)
    C.to_csv(os.path.join(RES, "nonlinear_noise_floor_candidates.csv"), index=False)
    with open(os.path.join(RES, "nonlinear_noise_floor_report.md"), "w") as fh:
        fh.write(f"# Null distribution of rho_nl ({N_TREES_NULL} trees; K_A={K_A}, K_B={K_B})\n\n")
        fh.write(S.to_string(index=False) + "\n\n## Candidates with q < 0.05 against Null B\n\n")
        fh.write(C[C.q_bh < 0.05].to_string(index=False) + "\n")
    print(S.to_string(index=False)); print(C[C.q_bh < 0.05].to_string(index=False))


if __name__ == "__main__":
    main()
