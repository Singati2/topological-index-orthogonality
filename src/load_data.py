"""ESOL (Delaney) dataset loader: 1128 small drug-like molecules + log solubility.

We use ESOL because:
  - Real diverse chemistry (not n=18 octane isomers)
  - n=1128 — appropriate for kill-test
  - Property is meaningful (aqueous solubility, log mol/L)
  - Standard MoleculeNet benchmark, so any future modeling comparisons are reproducible.
"""
from __future__ import annotations
import os
import urllib.request
import pandas as pd

ESOL_URL = "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/delaney-processed.csv"
CACHE_PATH = os.path.join(os.path.dirname(__file__), "..", "data", "esol.csv")


def load_esol() -> pd.DataFrame:
    cache = os.path.abspath(CACHE_PATH)
    if not os.path.exists(cache):
        os.makedirs(os.path.dirname(cache), exist_ok=True)
        print(f"Downloading ESOL to {cache} ...")
        urllib.request.urlretrieve(ESOL_URL, cache)
    df = pd.read_csv(cache)
    target = "measured log solubility in mols per litre"
    df = df[["smiles", target]].rename(columns={target: "logS"})
    return df.reset_index(drop=True)


if __name__ == "__main__":
    df = load_esol()
    print(f"ESOL loaded: n={len(df)}")
    print(df.head(3))
    print(f"logS range: [{df['logS'].min():.2f}, {df['logS'].max():.2f}]")
