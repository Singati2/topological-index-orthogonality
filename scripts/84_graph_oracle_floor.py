"""
Graph-oracle floors (Paper 2, v2.2) -- generating script for
results/graph_oracle_floor.csv.

For each curated dataset (data/curated/<name>_curated.csv, grouped by the
Weisfeiler-Lehman hash of the hydrogen-suppressed graph, which scripts/77
confirms to be exact isomorphism classes):
  * in-sample oracle: every molecule predicted by the mean target of its
    isomorphism class (an upper bound on the in-sample fit of any function
    of the unlabelled graph);
  * leave-one-out oracle: each molecule predicted by the mean of the OTHER
    members of its class, or by the global mean when it has none.
Regression: RMSE and R^2; BBBP: ROC-AUC.
"""
from __future__ import annotations

import os

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

PROJECT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RES = os.path.join(PROJECT, "results")


def main():
    rows = []
    for d in ["esol", "freesolv", "lipophilicity", "bbbp"]:
        c = pd.read_csv(os.path.join(PROJECT, "data", "curated", f"{d}_curated.csv"))
        y = c.target.values.astype(float); g = c.graph_wl_hash.values
        df = pd.DataFrame({"y": y, "g": g})
        gm = df.groupby("g").y.transform("mean").values
        cnt = df.groupby("g").y.transform("count").values
        gs = df.groupby("g").y.transform("sum").values
        loo = np.where(cnt > 1, (gs - y) / np.maximum(cnt - 1, 1), y.mean())
        if d == "bbbp":
            rows.append(dict(dataset=d, n=len(y), oracle_in_sample=roc_auc_score(y, gm), oracle_loo=roc_auc_score(y, loo), metric="AUC"))
        else:
            r2 = lambda p: 1 - np.sum((y - p) ** 2) / np.sum((y - y.mean()) ** 2)
            rm = lambda p: float(np.sqrt(np.mean((y - p) ** 2)))
            rows.append(dict(dataset=d, n=len(y), oracle_in_sample_rmse=rm(gm), oracle_in_sample_r2=r2(gm),
                             oracle_loo_rmse=rm(loo), oracle_loo_r2=r2(loo), sd=float(y.std()), metric="RMSE"))
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(RES, "graph_oracle_floor.csv"), index=False)
    print(out.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
