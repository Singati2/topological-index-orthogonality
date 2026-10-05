"""
Dataset verification and curated-data export (Paper 2, v2).

For each MoleculeNet file in data/ this script
  1. records its SHA-256 and, if --check-remote is given, re-downloads the
     official DeepChem/MoleculeNet file and confirms byte identity;
  2. checks the schema (SMILES and target columns present), row count,
     missing values, target range, SMILES validity, multi-fragment SMILES,
     exact canonical-SMILES duplicates and conflicting duplicate labels;
  3. writes the curated dataset used by the benchmark to
     data/curated/<name>_curated.csv (canonical SMILES, target, graph hash),
     using the curation rule of scripts/74 (identical to scripts/72);
  4. cross-checks the curated counts against results/benchmark_v2_curation.csv,
     which scripts/72 produces with an independent implementation.

Outputs: data/SHA256SUMS, data/curated/*.csv, results/dataset_verification.csv,
results/dataset_verification.md
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import os
import sys
import urllib.request

import networkx as nx
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

PROJECT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, PROJECT)
RDLogger.DisableLog("rdApp.*")

from src.load_data import DATASETS                          # noqa: E402
from src.mol_to_graph import smiles_to_graph                # noqa: E402

DATA = os.path.join(PROJECT, "data")
RES = os.path.join(PROJECT, "results")
TASK = {"esol": "regression", "freesolv": "regression",
        "lipophilicity": "regression", "bbbp": "classification"}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def load_curate():
    spec = importlib.util.spec_from_file_location(
        "dedup", os.path.join(PROJECT, "scripts", "74_dedup_sensitivity.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.curate


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check-remote", action="store_true",
                    help="re-download the official files and compare bytes")
    args = ap.parse_args()
    curate = load_curate()
    os.makedirs(os.path.join(DATA, "curated"), exist_ok=True)
    bench = None
    bp = os.path.join(RES, "benchmark_v2_curation.csv")
    if os.path.exists(bp):
        bench = pd.read_csv(bp).set_index("dataset")

    rows, sums = [], []
    for name, meta in DATASETS.items():
        path = os.path.join(DATA, meta["filename"])
        digest = sha256(path)
        sums.append(f"{digest}  {meta['filename']}")
        remote = "not checked"
        if args.check_remote:
            with urllib.request.urlopen(meta["url"]) as r:
                remote = "identical" if hashlib.sha256(r.read()).hexdigest() == digest else "DIFFERENT"
        raw = pd.read_csv(path)
        smi, tgt = meta["smiles_col"], meta["target_col"]
        assert smi in raw.columns and tgt in raw.columns, (name, list(raw.columns))
        mols = [Chem.MolFromSmiles(s) for s in raw[smi]]
        n_invalid = sum(m is None for m in mols)
        n_multi = sum(len(Chem.GetMolFrags(m)) > 1 for m in mols if m is not None)
        y = raw[tgt]
        if TASK[name] == "classification":
            assert set(y.dropna().unique()) <= {0, 1}, name

        df = raw.rename(columns={smi: "smiles", tgt: "target"})[["smiles", "target"]]
        cur, info = curate(df, TASK[name])
        keep = []
        for s, t in zip(cur["smiles"], cur["target"]):
            G = smiles_to_graph(s)
            if G is None:
                continue
            m = Chem.MolFromSmiles(s)
            frags = Chem.GetMolFrags(m, asMols=True)
            can = Chem.MolToSmiles(max(frags, key=lambda x: x.GetNumHeavyAtoms()))
            keep.append((can, float(t), nx.weisfeiler_lehman_graph_hash(G, iterations=4)))
        out = pd.DataFrame(keep, columns=["canonical_smiles", "target", "graph_wl_hash"])
        out.to_csv(os.path.join(DATA, "curated", f"{name}_curated.csv"), index=False)
        n_graphs = out.graph_wl_hash.nunique()
        shared = 100.0 * out.graph_wl_hash.duplicated(keep=False).mean()

        rec = dict(dataset=name, file=meta["filename"], sha256=digest, official_copy=remote,
                   smiles_col=smi, target_col=tgt, n_rows=len(raw),
                   n_missing_smiles=int(raw[smi].isna().sum()), n_missing_target=int(y.isna().sum()),
                   target_min=float(y.min()), target_max=float(y.max()), target_mean=float(y.mean()),
                   n_invalid_smiles=n_invalid, n_multifragment=n_multi,
                   n_conflicting_groups_dropped=info["n_conflicting_dropped"],
                   n_curated=len(out), n_distinct_graphs=n_graphs, pct_graph_shared=round(shared, 2))
        if bench is not None and name in bench.index:
            b = bench.loc[name]
            rec["matches_benchmark_curation"] = bool(
                int(b.n_final) == len(out) and int(b.n_distinct_graphs) == n_graphs
                and int(b.n_invalid_smiles) == n_invalid)
        rows.append(rec)

    with open(os.path.join(DATA, "SHA256SUMS"), "w") as fh:
        fh.write("\n".join(sums) + "\n")
    ver = pd.DataFrame(rows)
    ver.to_csv(os.path.join(RES, "dataset_verification.csv"), index=False)
    with open(os.path.join(RES, "dataset_verification.md"), "w") as fh:
        fh.write("# Dataset verification\n\n" + ver.T.to_string() + "\n")
    print(ver.T.to_string())


if __name__ == "__main__":
    main()
