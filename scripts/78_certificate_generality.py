"""
Generality of the span certificate across MoleculeNet (Paper 2, v2.1).

For every SMILES-based MoleculeNet dataset available (the four of the paper
plus BACE, ClinTox, SIDER, Tox21, ToxCast, HIV, MUV, QM8, QM9) this script
  1. parses all SMILES (largest fragment, hydrogen-suppressed graph);
  2. records the maximum degree distribution and the realised degree-pair set
     P_D, INCLUDING pairs with a degree above 4 (hypervalent S, P, ...);
  3. computes the count matrix M over P_D and the 18-index BID baseline B;
  4. checks the certificate rank[1,B] == rank[1,M] on the data, and reports the
     largest relative residual of any m_ij column on the baseline;
  5. reports whether the data-free certificate (rank Phi_B = 10 on P_4) alone
     would have sufficed, i.e. whether any molecule has degree > 4.

Only the SMILES column is used; no targets.  Datasets are read from the
directory given by --dir (default data/moleculenet_extra) or data/.

Output: results/certificate_generality.csv, results/certificate_generality.md
"""
from __future__ import annotations

import argparse
import glob
import os
import sys
from collections import Counter

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

PROJECT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, PROJECT)
RDLogger.DisableLog("rdApp.*")

from src import orthogonality as ortho                      # noqa: E402
from src.mol_to_graph import smiles_to_graph                # noqa: E402
from src.standard_indices import ALL_INDICES                # noqa: E402
import hashlib                                              # noqa: E402
import urllib.request                                       # noqa: E402

BID18 = ["M1", "M2", "mM1", "mM2", "F", "R", "SCI", "H", "GA", "AG",
         "ABC", "ABS", "AZI", "SO", "SO_red", "Alb", "Sigma", "redM1"]
RES = os.path.join(PROJECT, "results")
MIRROR = "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/"
EXTRA = {  # file name on the DeepChem/MoleculeNet mirror -> local name
    "bace.csv": "bace.csv", "clintox.csv.gz": "clintox.csv", "sider.csv.gz": "sider.csv",
    "tox21.csv.gz": "tox21.csv", "toxcast_data.csv.gz": "toxcast_data.csv", "HIV.csv": "HIV.csv",
    "muv.csv.gz": "muv.csv", "qm8.csv": "qm8.csv", "qm9.csv": "qm9.csv",
}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def download_extra(d):
    import gzip, shutil
    os.makedirs(d, exist_ok=True)
    for remote, local in EXTRA.items():
        dst = os.path.join(d, local)
        if os.path.exists(dst):
            continue
        tmp = dst + (".gz" if remote.endswith(".gz") else ".tmp")
        urllib.request.urlretrieve(MIRROR + remote, tmp)
        if remote.endswith(".gz"):
            with gzip.open(tmp, "rb") as fi, open(dst, "wb") as fo:
                shutil.copyfileobj(fi, fo)
            os.remove(tmp)
        else:
            os.rename(tmp, dst)
        print("downloaded", remote)


def find_smiles_col(df):
    for c in df.columns:
        if c.lower() in ("smiles", "mol", "canonical_smiles"):
            return c
    for c in df.columns:
        if "smiles" in c.lower():
            return c
    raise KeyError(list(df.columns))


def analyse(name, smiles):
    graphs, maxdegs, n_fail = [], [], 0
    for s in smiles:
        G = smiles_to_graph(s) if isinstance(s, str) else None
        if G is None or G.number_of_edges() == 0:
            n_fail += 1
            continue
        graphs.append(G)
        maxdegs.append(max(d for _, d in G.degree()))
    pairs = []
    for G in graphs:
        d = dict(G.degree())
        pairs.append(Counter(tuple(sorted((d[u], d[v]))) for u, v in G.edges()))
    keys = sorted({k for c in pairs for k in c})
    M = np.array([[c.get(k, 0) for k in keys] for c in pairs], float)
    B = np.array([[ALL_INDICES[b](G) for b in BID18] for G in graphs])
    rM, rB, rBM = ortho.numerical_rank(M), ortho.numerical_rank(B), ortho.numerical_rank(np.hstack([B, M]))
    max_res = max(ortho.relative_residual(M[:, j], B) for j in range(M.shape[1]))
    # which m_ij columns (if any) are NOT in the span of the baseline?
    bad = [keys[j] for j in range(M.shape[1]) if ortho.relative_residual(M[:, j], B) >= ortho.SPAN_RTOL]
    md = Counter(maxdegs)
    # For datasets with hypervalent atoms, also certify the max-degree<=4 subset,
    # which is the class the data-free certificate covers.
    restricted = {}
    if any(k > 4 for k in md):
        sel = [i for i, g in enumerate(graphs) if maxdegs[i] <= 4]
        keys4 = sorted({k for i in sel for k in pairs[i]})
        M4 = np.array([[pairs[i].get(k, 0) for k in keys4] for i in sel], float)
        B4 = B[sel]
        res4 = max(ortho.relative_residual(M4[:, j], B4) for j in range(M4.shape[1]))
        restricted = dict(restricted_n=len(sel), restricted_rank_mij=ortho.numerical_rank(M4),
                          restricted_rank_BID18=ortho.numerical_rank(B4),
                          restricted_certificate_holds=bool(ortho.numerical_rank(B4) == ortho.numerical_rank(np.hstack([B4, M4])) and res4 < ortho.SPAN_RTOL),
                          restricted_max_rel_resid=res4)
    return dict(
        dataset=name, n_smiles=len(smiles), n_graphs=len(graphs), n_unparsed_or_empty=n_fail,
        n_maxdeg_gt4=sum(v for k, v in md.items() if k > 4), max_degree=max(maxdegs),
        degree_pairs_realised=len(keys), pairs_above_4=sum(1 for k in keys if k[1] > 4),
        rank_mij=rM, rank_BID18=rB, rank_BID18_plus_mij=rBM,
        certificate_holds=bool(rB == rBM and max_res < ortho.SPAN_RTOL),
        max_rel_resid_mij_on_BID18=max_res,
        uncovered_pairs=" ".join(f"({a},{b})" for a, b in bad),
        **restricted,
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=os.path.join(PROJECT, "data", "moleculenet_extra"))
    ap.add_argument("--max-rows", type=int, default=0, help="subsample large files (0 = all)")
    args = ap.parse_args()
    download_extra(args.dir)
    files = sorted(glob.glob(os.path.join(PROJECT, "data", "*.csv"))) + sorted(glob.glob(os.path.join(args.dir, "*.csv")))
    rows = []
    for f in files:
        name = os.path.splitext(os.path.basename(f))[0]
        df = pd.read_csv(f)
        col = find_smiles_col(df)
        smi = df[col].tolist()
        if args.max_rows and len(smi) > args.max_rows:
            rng = np.random.default_rng(20261005)
            smi = [smi[i] for i in rng.choice(len(smi), args.max_rows, replace=False)]
        r = analyse(name, smi)
        r["file_sha256"] = sha256(f)
        r["source"] = "data/" if f.startswith(os.path.join(PROJECT, "data") + os.sep) and "moleculenet_extra" not in f else MIRROR
        rows.append(r)
        print({k: v for k, v in r.items()})
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(RES, "certificate_generality.csv"), index=False)
    with open(os.path.join(RES, "certificate_generality.md"), "w") as fh:
        fh.write("# Span certificate across MoleculeNet datasets\n\n" + out.to_string(index=False) + "\n")


if __name__ == "__main__":
    main()
