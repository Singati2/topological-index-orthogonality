"""
Target-aware screen against a chemistry-aware baseline (Paper 2, v2.1).

The screen of Section 6 uses a baseline of 30 topological indices.  A QSPR
practitioner's baseline also contains atom- and bond-aware descriptors.  This
script repeats the linear target-aware screen of the 20 non-BID candidates
against three baselines:
    X_topo   the 30 topological indices (as in the paper),
    X_rdkit  RDKit 2D descriptors (the same 180 used by scripts/72),
    X_both   their union, with exactly collinear columns removed,
and reports, per candidate, the partial correlation with the target given
each baseline (Fisher 95% interval, in-span certificate), and the relative
residual of the candidate on each baseline.  It also asks the reverse
question -- do the 30 topological indices carry residual target signal
beyond RDKit 2D? -- and gives an out-of-sample view with ridge regression
(regression) / ridge-penalised logistic regression (BBBP) on X_rdkit versus
X_rdkit + X_topo + candidates, 5-fold CV, preprocessing inside folds.

Outputs: results/screen_vs_chemistry_candidates.csv,
         results/screen_vs_chemistry_topo30.csv,
         results/screen_vs_chemistry_summary.csv,
         results/screen_vs_chemistry_report.md
"""
from __future__ import annotations

import importlib.util
import os
import sys

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from sklearn.linear_model import LogisticRegressionCV, RidgeCV
from sklearn.metrics import mean_squared_error, roc_auc_score
from sklearn.model_selection import KFold, StratifiedKFold
from sklearn.preprocessing import StandardScaler

PROJECT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, PROJECT)
RDLogger.DisableLog("rdApp.*")

from src import orthogonality as ortho                      # noqa: E402
from src.load_data import load                              # noqa: E402
from src.mol_to_graph import smiles_to_graph                # noqa: E402
from src.novel_candidates import CANDIDATE_INDICES          # noqa: E402
from src.standard_indices import compute_all                # noqa: E402

_spec = importlib.util.spec_from_file_location("b72", os.path.join(PROJECT, "scripts", "72_benchmark_v2.py"))
b72 = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(b72)

SEED = 20261005
TAU = 0.10
DATASETS = {"esol": "regression", "freesolv": "regression",
            "lipophilicity": "regression", "bbbp": "classification"}
RES = os.path.join(PROJECT, "results")


def impute(v):
    v = np.asarray(v, float)
    if np.any(~np.isfinite(v)):
        v = np.where(np.isfinite(v), v, np.nanmean(np.where(np.isfinite(v), v, np.nan)))
    return v


def verdict(cert, lo, hi):
    if cert["in_span"]:
        return "IN_SPAN"
    abs_lo = 0.0 if lo <= 0 <= hi else min(abs(lo), abs(hi))
    abs_hi = max(abs(lo), abs(hi))
    return "PASS" if abs_lo >= TAU else "NEGLIGIBLE" if abs_hi < TAU else "INCONCLUSIVE"


def pcor_row(z, X, y):
    c = ortho.span_certified_pcor(z, X, y)
    lo, hi = ortho.fisher_ci(c["pcor"], c["df"]) if not c["in_span"] else (0.0, 0.0)
    return c["pcor"], lo, hi, c["rel_resid"], verdict(c, lo, hi)


def drop_collinear(X, names, rtol=1e-9):
    """Greedy QR-style removal of columns that are (numerically) linear
    combinations of the columns kept before them.  Returns kept indices."""
    Xc = X - X.mean(0)
    norms = np.linalg.norm(Xc, axis=0)
    keep, Q = [], np.zeros((X.shape[0], 0))
    for j in range(X.shape[1]):
        if norms[j] == 0:
            continue
        v = Xc[:, j] / norms[j]
        r = v - Q @ (Q.T @ v)
        if np.linalg.norm(r) > rtol:
            keep.append(j)
            Q = np.column_stack([Q, r / np.linalg.norm(r)])
    return keep


def cv_metric(task, X, y, seed):
    """5-fold CV of a ridge-type linear model; scaler fitted inside folds."""
    kf = (StratifiedKFold if task == "classification" else KFold)(5, shuffle=True, random_state=seed)
    preds = np.zeros(len(y))
    for tr, te in kf.split(X, y if task == "classification" else None):
        sc = StandardScaler().fit(X[tr])
        Xtr, Xte = sc.transform(X[tr]), sc.transform(X[te])
        if task == "classification":
            m = LogisticRegressionCV(Cs=10, cv=5, penalty="l2", max_iter=5000, random_state=seed).fit(Xtr, y[tr])
            preds[te] = m.predict_proba(Xte)[:, 1]
        else:
            m = RidgeCV(alphas=np.logspace(-3, 3, 25)).fit(Xtr, y[tr])
            preds[te] = m.predict(Xte)
    return roc_auc_score(y, preds) if task == "classification" else float(np.sqrt(mean_squared_error(y, preds)))


def r2(y, X):
    A = np.column_stack([np.ones(len(y)), X])
    b, *_ = np.linalg.lstsq(A, y, rcond=None)
    return float(1 - np.sum((y - A @ b) ** 2) / np.sum((y - y.mean()) ** 2))


def main():
    cand_rows, topo_rows, summ = [], [], []
    for name, task in DATASETS.items():
        df = load(name)
        graphs, mols, y = [], [], []
        for s, t in zip(df["smiles"], df["target"]):
            G = smiles_to_graph(s)
            if G is None:
                continue
            mol = Chem.MolFromSmiles(s)
            frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
            if len(frags) > 1:
                mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())
            graphs.append(G); mols.append(mol); y.append(float(t))
        y = np.asarray(y); n = len(y)
        topo = pd.DataFrame([compute_all(G) for G in graphs])
        X_topo = np.column_stack([impute(c) for c in topo.values.T]); topo_names = list(topo.columns)
        R = pd.DataFrame([b72.rdkit_desc(m) for m in mols], columns=[nm for nm, _ in b72.RDKIT_DESC])
        ok = [c for c in R.columns if np.isfinite(R[c]).all() and R[c].std() > 0]
        X_rd = R[ok].values.astype(float); rd_names = ok
        keep_rd = drop_collinear(X_rd, rd_names)
        X_rd = X_rd[:, keep_rd]; rd_names = [rd_names[j] for j in keep_rd]
        both = np.column_stack([X_topo, X_rd]); both_names = topo_names + rd_names
        keep = drop_collinear(both, both_names)
        X_both = both[:, keep]; both_kept = [both_names[j] for j in keep]
        n_topo_dropped = sum(1 for j in range(len(topo_names)) if j not in keep)
        n_rd_dropped_vs_topo = sum(1 for j in range(len(topo_names), len(both_names)) if j not in keep)
        print(f"[{name}] n={n} rdkit cols {len(ok)} -> {len(rd_names)} after internal collinearity; "
              f"both {len(both_names)} -> {len(both_kept)} (topo dropped {n_topo_dropped}, rdkit dropped {n_rd_dropped_vs_topo})", flush=True)

        cands = {}
        for cname, fn in CANDIDATE_INDICES.items():
            vals = []
            for G in graphs:
                try:
                    vals.append(fn(G))
                except Exception:
                    vals.append(np.nan)
            cands[cname] = impute(vals)
        for cname, z in cands.items():
            if np.std(z) == 0:
                continue
            r_both = max(abs(np.corrcoef(z, X_both[:, j])[0, 1]) for j in range(X_both.shape[1]) if np.std(X_both[:, j]) > 0)
            row = dict(dataset=name, candidate=cname, max_abs_r_both=r_both)
            for lab, X in [("topo", X_topo), ("rdkit", X_rd), ("both", X_both)]:
                p, lo, hi, rr, v = pcor_row(z, X, y)
                row.update({f"pcor_{lab}": p, f"ci_lo_{lab}": lo, f"ci_hi_{lab}": hi,
                            f"relres_{lab}": rr, f"verdict_{lab}": v})
            cand_rows.append(row)
        for j, tn in enumerate(topo_names):
            p, lo, hi, rr, v = pcor_row(X_topo[:, j], X_rd, y)
            topo_rows.append(dict(dataset=name, index=tn, pcor_given_rdkit=p, ci_lo=lo, ci_hi=hi,
                                  relres_on_rdkit=rr, verdict=v))
        n_rd_signal_beyond_topo = 0
        for j in range(X_rd.shape[1]):
            p, lo, hi, rr, v = pcor_row(X_rd[:, j], X_topo, y)
            n_rd_signal_beyond_topo += (v == "PASS")
        C = np.column_stack([z for z in cands.values() if np.std(z) > 0])
        X_all = np.column_stack([X_both, C])
        vc = pd.DataFrame(cand_rows)[lambda d: d.dataset == name]
        summ.append(dict(
            dataset=name, task=task, n=n, n_rdkit_cols=len(rd_names), n_both_cols=len(both_kept),
            n_topo_dropped_collinear=n_topo_dropped, n_rdkit_dropped_collinear_with_topo=n_rd_dropped_vs_topo,
            r2_rdkit=r2(y, X_rd), r2_rdkit_topo=r2(y, X_both), r2_rdkit_topo_cands=r2(y, X_all),
            r2_topo=r2(y, X_topo),
            cv_rdkit=cv_metric(task, X_rd, y, SEED), cv_rdkit_topo=cv_metric(task, X_both, y, SEED),
            cv_rdkit_topo_cands=cv_metric(task, X_all, y, SEED), cv_topo=cv_metric(task, X_topo, y, SEED),
            cand_pass_both=int((vc.verdict_both == "PASS").sum()),
            cand_inconclusive_both=int((vc.verdict_both == "INCONCLUSIVE").sum()),
            cand_negligible_both=int((vc.verdict_both == "NEGLIGIBLE").sum()),
            cand_inspan_both=int((vc.verdict_both == "IN_SPAN").sum()),
            topo_pass_given_rdkit=int(sum(1 for r in topo_rows if r["dataset"] == name and r["verdict"] == "PASS")),
            rdkit_pass_given_topo=int(n_rd_signal_beyond_topo),
        ))
        print(summ[-1], flush=True)
    pd.DataFrame(cand_rows).to_csv(os.path.join(RES, "screen_vs_chemistry_candidates.csv"), index=False)
    pd.DataFrame(topo_rows).to_csv(os.path.join(RES, "screen_vs_chemistry_topo30.csv"), index=False)
    S = pd.DataFrame(summ); S.to_csv(os.path.join(RES, "screen_vs_chemistry_summary.csv"), index=False)
    with open(os.path.join(RES, "screen_vs_chemistry_report.md"), "w") as fh:
        fh.write("# Target-aware screen against a chemistry-aware baseline\n\n" + S.T.to_string() + "\n\n")
        fh.write("## Candidates not negligible given X_both\n\n")
        c = pd.DataFrame(cand_rows)
        fh.write(c[c.verdict_both.isin(["PASS", "INCONCLUSIVE"])].to_string(index=False) + "\n\n")
        fh.write("## Topological indices with residual signal beyond RDKit 2D (PASS)\n\n")
        t = pd.DataFrame(topo_rows)
        fh.write(t[t.verdict == "PASS"].to_string(index=False) + "\n")


if __name__ == "__main__":
    main()
