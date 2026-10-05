"""
Deduplication sensitivity of the target-aware screen (Paper 2, v2).

The screens in scripts 16/17 use the MoleculeNet files as distributed.  Those
files contain repeated molecules (same canonical SMILES after keeping the
largest fragment) and, in BBBP, conflicting labels.  This script repeats the
linear target-aware screen of the 20 alternative-family candidates on a
de-duplicated copy of each dataset and reports whether any verdict changes.

Curation rule (identical to scripts/72_benchmark_v2.py):
  * parse with RDKit, keep the largest fragment, RDKit canonical SMILES;
  * regression: average the target over duplicates whose labels agree within
    DUP_TOL = 0.5 target units, drop duplicate groups that disagree by more;
  * classification: drop duplicate groups with conflicting labels.

Output: results/dedup_sensitivity.csv, results/dedup_sensitivity_summary.txt
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

PROJECT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, PROJECT)
RDLogger.DisableLog("rdApp.*")

from src import orthogonality as ortho                      # noqa: E402
from src.load_data import load                              # noqa: E402
from src.mol_to_graph import smiles_to_graph                # noqa: E402
from src.novel_candidates import CANDIDATE_INDICES          # noqa: E402
from src.standard_indices import compute_all                # noqa: E402

DATASETS = {"esol": "regression", "freesolv": "regression",
            "lipophilicity": "regression", "bbbp": "classification"}
TAU_R, TAU_P = 0.95, 0.10
DUP_TOL = 0.5
RES = os.path.join(PROJECT, "results")


def canonical_key(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    frags = Chem.GetMolFrags(m, asMols=True)
    m = max(frags, key=lambda x: x.GetNumHeavyAtoms())
    return Chem.MolToSmiles(m)


def curate(df, task):
    df = df.copy()
    df["key"] = [canonical_key(s) for s in df["smiles"]]
    n_raw = len(df)
    df = df[df["key"].notna()]
    n_valid = len(df)
    rng = df.groupby("key")["target"].agg(lambda t: t.max() - t.min())
    bad = set(rng[rng > (DUP_TOL if task == "regression" else 0)].index)
    n_conflict = len(bad)
    g = df[~df["key"].isin(bad)].groupby("key").agg(
        smiles=("smiles", "first"), target=("target", "mean"))
    return g.reset_index(drop=True), dict(n_raw=n_raw, n_valid=n_valid,
                                          n_unique=len(g) + n_conflict,
                                          n_conflicting_dropped=n_conflict,
                                          n_curated=len(g))


def screen(df):
    graphs, y = [], []
    for s, t in zip(df["smiles"], df["target"]):
        G = smiles_to_graph(s)
        if G is None:
            continue
        graphs.append(G); y.append(float(t))
    y = np.asarray(y)
    base = pd.DataFrame([compute_all(G) for G in graphs])
    X = base.values
    rows = []
    for name, fn in CANDIDATE_INDICES.items():
        z = []
        for G in graphs:
            try:
                z.append(fn(G))
            except Exception:
                z.append(np.nan)
        z = np.asarray(z, float)
        if np.isnan(z).any():
            z = np.where(np.isnan(z), np.nanmean(z), z)
        if np.std(z) == 0:
            continue
        maxr = max(abs(np.corrcoef(z, base[c])[0, 1]) for c in base.columns if np.std(base[c]) > 0)
        cert = ortho.span_certified_pcor(z, X, y)
        lo, hi = ortho.fisher_ci(cert["pcor"], cert["df"]) if not cert["in_span"] else (0.0, 0.0)
        abs_lo = 0.0 if lo <= 0 <= hi else min(abs(lo), abs(hi))
        abs_hi = max(abs(lo), abs(hi))
        verdict = ("IN_SPAN" if cert["in_span"] else "PASS" if abs_lo >= TAU_P
                   else "NEGLIGIBLE" if abs_hi < TAU_P else "INCONCLUSIVE")
        rows.append(dict(candidate=name, max_abs_r=maxr, pcor=cert["pcor"], ci_lo=lo, ci_hi=hi,
                         pairwise="PASS" if maxr < TAU_R else "FAIL", ci_verdict=verdict, n=len(y)))
    return pd.DataFrame(rows)


def main():
    out, lines = [], []
    for name, task in DATASETS.items():
        df = load(name)
        cur, info = curate(df, task)
        r = screen(cur)
        r.insert(0, "dataset", name)
        orig = pd.read_csv(os.path.join(RES, "novel_candidates_experiment",
                                        "novel_candidates_multidataset.csv"))
        orig = orig[orig.dataset == name].set_index("candidate")
        r["orig_pcor"] = r["candidate"].map(orig["partial_corr_target"])
        r["orig_pairwise"] = r["candidate"].map(orig["pairwise_verdict"])
        r["orig_ci_verdict"] = r["candidate"].map(orig["pcor_ci_verdict"])
        out.append(r)
        changed = r[(r.pairwise != r.orig_pairwise) | (r.ci_verdict != r.orig_ci_verdict)]
        lines.append(f"{name}: {info}; max |pcor change| = "
                     f"{(r.pcor - r.orig_pcor).abs().max():.4f}; verdict changes: "
                     + (", ".join(f"{c} {a}->{b}" for c, a, b in
                                  zip(changed.candidate, changed.orig_ci_verdict, changed.ci_verdict))
                        or "none"))
    pd.concat(out).to_csv(os.path.join(RES, "dedup_sensitivity.csv"), index=False)
    with open(os.path.join(RES, "dedup_sensitivity_summary.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
