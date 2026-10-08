"""
Chemistry-aware screen, hygienic version (Paper 2, v2.2).

Replaces scripts/79 for the out-of-sample comparison and adds the verdicts
against RDKit alone.  Differences from 79:
  * curated data (data/curated/<name>_curated.csv: de-duplicated, graph hash);
  * folds grouped by isomorphism class of the hydrogen-suppressed graph
    (GroupKFold / StratifiedGroupKFold), 5 folds, as in the benchmark;
  * removal of non-finite, constant and collinear RDKit columns and the
    standardisation are fitted on the training fold only;
  * ridge penalty / logistic C chosen by an inner grouped CV;
  * candidate verdicts reported against three baselines: topology (30),
    RDKit alone, RDKit + topology (on the curated data, full sample, as a
    descriptive diagnostic like Section 6).
Outputs: results/screen_vs_chemistry_grouped_cv.csv,
         results/screen_vs_chemistry_grouped_candidates.csv,
         results/screen_vs_chemistry_grouped_summary.csv
"""
from __future__ import annotations

import importlib.util
import os
import sys

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import mean_squared_error, roc_auc_score
from sklearn.model_selection import GroupKFold, StratifiedGroupKFold

PROJECT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, PROJECT)
RDLogger.DisableLog("rdApp.*")
from src import orthogonality as ortho                      # noqa: E402
from src.mol_to_graph import smiles_to_graph                # noqa: E402
from src.novel_candidates import CANDIDATE_INDICES          # noqa: E402
from src.standard_indices import compute_all                # noqa: E402

_s = importlib.util.spec_from_file_location("b72", os.path.join(PROJECT, "scripts", "72_benchmark_v2.py"))
b72 = importlib.util.module_from_spec(_s); _s.loader.exec_module(b72)
_s2 = importlib.util.spec_from_file_location("s79", os.path.join(PROJECT, "scripts", "79_screen_vs_chemistry_baseline.py"))
s79 = importlib.util.module_from_spec(_s2); _s2.loader.exec_module(s79)

SEED = 20261005
DATASETS = {"esol": "regression", "freesolv": "regression", "lipophilicity": "regression", "bbbp": "classification"}
RES = os.path.join(PROJECT, "results")
ALPHAS = np.logspace(-3, 3, 25)
CS = np.logspace(-3, 2, 11)


def fit_predict(task, Xtr, ytr, gtr, Xte, seed):
    """Column cleaning + scaling + penalty selection all inside the training fold."""
    keep = [j for j in range(Xtr.shape[1]) if np.isfinite(Xtr[:, j]).all() and Xtr[:, j].std() > 0]
    Xtr = Xtr[:, keep]; Xte = Xte[:, keep]
    kk = s79.drop_collinear(Xtr, list(range(Xtr.shape[1])))
    Xtr = Xtr[:, kk]; Xte = Xte[:, kk]
    mu, sd = Xtr.mean(0), Xtr.std(0); sd[sd == 0] = 1
    Xtr = (Xtr - mu) / sd; Xte = np.clip((Xte - mu) / sd, -10, 10)
    Xte = np.where(np.isfinite(Xte), Xte, 0.0)
    inner = (StratifiedGroupKFold if task == "classification" else GroupKFold)(3)
    best, best_score = None, -np.inf
    grid = CS if task == "classification" else ALPHAS
    for h in grid:
        sc = []
        for tr, va in inner.split(Xtr, ytr if task == "classification" else None, gtr):
            if task == "classification":
                m = LogisticRegression(C=h, max_iter=5000).fit(Xtr[tr], ytr[tr])
                sc.append(roc_auc_score(ytr[va], m.predict_proba(Xtr[va])[:, 1]))
            else:
                m = Ridge(alpha=h).fit(Xtr[tr], ytr[tr])
                sc.append(-mean_squared_error(ytr[va], m.predict(Xtr[va])))
        if np.mean(sc) > best_score:
            best_score, best = np.mean(sc), h
    if task == "classification":
        return LogisticRegression(C=best, max_iter=5000).fit(Xtr, ytr).predict_proba(Xte)[:, 1]
    return Ridge(alpha=best).fit(Xtr, ytr).predict(Xte)


def main():
    cv_rows, cand_rows, summ = [], [], []
    for name, task in DATASETS.items():
        cur = pd.read_csv(os.path.join(PROJECT, "data", "curated", f"{name}_curated.csv"))
        graphs, mols, y, g = [], [], [], []
        for s, t, h in zip(cur.canonical_smiles, cur.target, cur.graph_wl_hash):
            G = smiles_to_graph(s)
            if G is None:
                continue
            graphs.append(G); mols.append(Chem.MolFromSmiles(s)); y.append(float(t)); g.append(h)
        y = np.asarray(y); g = np.asarray(g); n = len(y)
        X_topo = np.column_stack([s79.impute(c) for c in pd.DataFrame([compute_all(G) for G in graphs]).values.T])
        X_rd = pd.DataFrame([b72.rdkit_desc(m) for m in mols]).values.astype(float)
        C = np.column_stack([s79.impute([CANDIDATE_INDICES[k](G) for G in graphs]) for k in CANDIDATE_INDICES])
        sets = {"rdkit": X_rd, "rdkit+topo": np.hstack([X_rd, X_topo]), "rdkit+topo+cands": np.hstack([X_rd, X_topo, C]), "topo": X_topo}
        outer = (StratifiedGroupKFold if task == "classification" else GroupKFold)(5)
        splits = list(outer.split(X_topo, y if task == "classification" else None, g))
        for lab, X in sets.items():
            pred = np.zeros(n)
            for tr, te in splits:
                pred[te] = fit_predict(task, X[tr], y[tr], g[tr], X[te], SEED)
            metric = roc_auc_score(y, pred) if task == "classification" else float(np.sqrt(mean_squared_error(y, pred)))
            cv_rows.append(dict(dataset=name, feature_set=lab, n=n, n_groups=len(set(g)), metric=("auc" if task == "classification" else "rmse"), value=metric))
            print(f"[{name}] {lab:18s} {metric:.4f}", flush=True)
        # descriptive full-sample verdicts (clean RDKit columns on the full curated sample)
        ok = [j for j in range(X_rd.shape[1]) if np.isfinite(X_rd[:, j]).all() and X_rd[:, j].std() > 0]
        R = X_rd[:, ok]; R = R[:, s79.drop_collinear(R, list(range(R.shape[1])))]
        both = np.hstack([X_topo, R]); both = both[:, s79.drop_collinear(both, list(range(both.shape[1])))]
        for k, z in zip(CANDIDATE_INDICES, C.T):
            if np.std(z) == 0:
                continue
            row = dict(dataset=name, candidate=k)
            for lab, X in [("topo", X_topo), ("rdkit", R), ("both", both)]:
                p, lo, hi, rr, v = s79.pcor_row(z, X, y)
                row.update({f"pcor_{lab}": p, f"ci_lo_{lab}": lo, f"ci_hi_{lab}": hi, f"verdict_{lab}": v})
            cand_rows.append(row)
        cr = pd.DataFrame(cand_rows); cr = cr[cr.dataset == name]
        summ.append(dict(dataset=name, n=n, n_groups=len(set(g)), n_rdkit=R.shape[1], n_both=both.shape[1],
                         **{f"{lab}_{v}": int((cr[f"verdict_{lab}"] == v).sum()) for lab in ("rdkit", "both") for v in ("PASS", "INCONCLUSIVE", "NEGLIGIBLE", "IN_SPAN")}))
    pd.DataFrame(cv_rows).to_csv(os.path.join(RES, "screen_vs_chemistry_grouped_cv.csv"), index=False)
    pd.DataFrame(cand_rows).to_csv(os.path.join(RES, "screen_vs_chemistry_grouped_candidates.csv"), index=False)
    pd.DataFrame(summ).to_csv(os.path.join(RES, "screen_vs_chemistry_grouped_summary.csv"), index=False)
    print(pd.DataFrame(summ).to_string(index=False))


if __name__ == "__main__":
    main()
