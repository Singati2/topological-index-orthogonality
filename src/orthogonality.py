"""Orthogonality / redundancy analysis on a descriptor block.

Given a DataFrame of descriptors (columns) over molecules (rows), produce:
  - Pearson correlation matrix
  - Hierarchical clustering of descriptors by |1 - r|
  - PCA variance-explained
  - VIF for each descriptor against the rest
  - Partial correlation of each descriptor with the target, controlling for the others

The kill-criterion for a *new* candidate index:
  - If max |r| with any baseline index >= 0.95, the new index is statistically
    redundant on this chemistry set. Redesign or pivot.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler


def correlation_matrix(df_desc: pd.DataFrame) -> pd.DataFrame:
    return df_desc.corr(method="pearson")


def max_abs_corr_with_baseline(corr: pd.DataFrame, new_name: str) -> tuple[str, float]:
    """Return (baseline_index_name, |r|) for the baseline most correlated with new_name."""
    if new_name not in corr.columns:
        raise KeyError(f"{new_name!r} not in correlation matrix")
    s = corr[new_name].drop(index=new_name).abs()
    best = s.idxmax()
    return best, float(s.loc[best])


def pca_variance_explained(df_desc: pd.DataFrame) -> pd.DataFrame:
    X = StandardScaler().fit_transform(df_desc.values)
    p = PCA().fit(X)
    return pd.DataFrame({
        "component": np.arange(1, len(p.explained_variance_ratio_) + 1),
        "var_explained": p.explained_variance_ratio_,
        "cumulative": np.cumsum(p.explained_variance_ratio_),
    })


def vif(df_desc: pd.DataFrame) -> pd.DataFrame:
    """Variance inflation factor for each column regressed on the others.

    VIF = 1 / (1 - R^2). Anything > 10 is heavily collinear.
    """
    cols = list(df_desc.columns)
    X = df_desc.values
    rows = []
    for j, c in enumerate(cols):
        y = X[:, j]
        Xrest = np.delete(X, j, axis=1)
        # Center / scale not strictly needed for R^2; LinearRegression handles it.
        r2 = LinearRegression().fit(Xrest, y).score(Xrest, y)
        vif_j = float("inf") if r2 >= 0.999999 else 1.0 / (1.0 - r2)
        rows.append((c, r2, vif_j))
    return pd.DataFrame(rows, columns=["descriptor", "R2_vs_rest", "VIF"]).sort_values("VIF", ascending=False).reset_index(drop=True)


def partial_corr_with_target(df_desc: pd.DataFrame, y: pd.Series) -> pd.DataFrame:
    """For each descriptor d, partial correlation with y controlling for all other descriptors.

    Computed as: residualize d on others, residualize y on others, then pearson.
    """
    cols = list(df_desc.columns)
    X = df_desc.values
    y_arr = y.values
    rows = []
    for j, c in enumerate(cols):
        d = X[:, j]
        Xrest = np.delete(X, j, axis=1)
        d_res = d - LinearRegression().fit(Xrest, d).predict(Xrest)
        y_res = y_arr - LinearRegression().fit(Xrest, y_arr).predict(Xrest)
        # Pearson of residuals
        if np.std(d_res) < 1e-12 or np.std(y_res) < 1e-12:
            pr = 0.0
        else:
            pr = float(np.corrcoef(d_res, y_res)[0, 1])
        # Also raw correlation for comparison
        raw = float(np.corrcoef(d, y_arr)[0, 1])
        rows.append((c, raw, pr))
    return pd.DataFrame(rows, columns=["descriptor", "raw_corr_y", "partial_corr_y"]).sort_values("partial_corr_y", key=lambda s: s.abs(), ascending=False).reset_index(drop=True)


def redundancy_report(corr: pd.DataFrame, threshold: float = 0.95) -> pd.DataFrame:
    """Pairs of descriptors with |r| >= threshold."""
    cols = corr.columns
    pairs = []
    for i in range(len(cols)):
        for j in range(i + 1, len(cols)):
            r = corr.iloc[i, j]
            if abs(r) >= threshold:
                pairs.append((cols[i], cols[j], float(r)))
    return pd.DataFrame(pairs, columns=["a", "b", "r"]).sort_values("r", key=lambda s: s.abs(), ascending=False).reset_index(drop=True)


def kill_test(corr: pd.DataFrame, new_name: str, threshold: float = 0.95) -> dict:
    """Return verdict for a new candidate index vs. the baseline set."""
    best, val = max_abs_corr_with_baseline(corr, new_name)
    verdict = "FAIL — statistically redundant with existing index" if val >= threshold else "PASS — carries information not in this baseline set"
    return {
        "new_index": new_name,
        "most_correlated_baseline": best,
        "max_abs_corr": val,
        "threshold": threshold,
        "verdict": verdict,
    }
