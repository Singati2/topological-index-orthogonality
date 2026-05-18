"""Dataset loaders for ESOL, FreeSolv, Lipophilicity (MoleculeNet).

All three are public, freely downloadable, and standard QSPR benchmarks.
We download on first call and cache under data/.

  ESOL (Delaney 2004) — aqueous solubility, n ≈ 1128
  FreeSolv (Mobley 2014, via DeepChem SAMPL) — hydration free energy, n ≈ 642
  Lipophilicity (Wenlock 2015, via DeepChem ChEMBL extract) — logD, n ≈ 4200

After loading, the returned DataFrame has columns: smiles, target.
The target name is recorded in df.attrs["target_name"] and df.attrs["dataset"].
"""
from __future__ import annotations
import os
import urllib.request
import pandas as pd

DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data"))

DATASETS = {
    "esol": {
        "url":         "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/delaney-processed.csv",
        "filename":    "esol.csv",
        "smiles_col":  "smiles",
        "target_col":  "measured log solubility in mols per litre",
        "target_name": "logS",
        "description": "Aqueous solubility (log mol/L), Delaney 2004",
    },
    "freesolv": {
        "url":         "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/SAMPL.csv",
        "filename":    "freesolv.csv",
        "smiles_col":  "smiles",
        "target_col":  "expt",
        "target_name": "dG_hyd",
        "description": "Hydration free energy (kcal/mol), Mobley 2014",
    },
    "lipophilicity": {
        "url":         "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/Lipophilicity.csv",
        "filename":    "lipophilicity.csv",
        "smiles_col":  "smiles",
        "target_col":  "exp",
        "target_name": "logD",
        "description": "Octanol/water distribution logD at pH 7.4, ChEMBL extract",
    },
}


def _download_if_missing(name: str) -> str:
    spec = DATASETS[name]
    cache_path = os.path.join(DATA_DIR, spec["filename"])
    if not os.path.exists(cache_path):
        os.makedirs(DATA_DIR, exist_ok=True)
        print(f"Downloading {name} to {cache_path} ...")
        urllib.request.urlretrieve(spec["url"], cache_path)
    return cache_path


def load(name: str) -> pd.DataFrame:
    """Load a benchmark dataset by name ('esol' | 'freesolv' | 'lipophilicity').

    Returns a DataFrame with columns ['smiles', 'target'].
    The target's chemical name is in df.attrs['target_name'].
    """
    name = name.lower()
    if name not in DATASETS:
        raise ValueError(f"Unknown dataset {name!r}. Available: {list(DATASETS)}")
    spec = DATASETS[name]
    path = _download_if_missing(name)
    df = pd.read_csv(path)
    df = df[[spec["smiles_col"], spec["target_col"]]].rename(
        columns={spec["smiles_col"]: "smiles", spec["target_col"]: "target"}
    ).reset_index(drop=True)
    df.attrs["dataset"] = name
    df.attrs["target_name"] = spec["target_name"]
    df.attrs["description"] = spec["description"]
    return df


# Backward-compatible alias used by Phase-0 scripts.
def load_esol() -> pd.DataFrame:
    df = load("esol")
    return df.rename(columns={"target": "logS"})


if __name__ == "__main__":
    for name in DATASETS:
        df = load(name)
        print(f"{name:14s}  n={len(df):>5}  target={df.attrs['target_name']:8s}  "
              f"range=[{df['target'].min():.2f}, {df['target'].max():.2f}]")
