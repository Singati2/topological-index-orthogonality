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


# ---------------------------------------------------------------------------
# Span certificate and numerically safe partial correlation (v2, 2026-10).
#
# The earlier code declared a residual "zero" only when its *absolute* standard
# deviation fell below 1e-12.  For an index that lies exactly in the column
# span of the baseline (every bond-additive degree-based index does, once the
# baseline spans the edge-degree-pair counts), the least-squares residual is
# floating-point noise of relative size ~1e-14 but absolute size well above
# 1e-12, and correlating that noise with the target produced spurious
# "partial correlations" of order 0.01-0.04.  The helpers below use a
# *relative* residual and report the in-span status explicitly.
# ---------------------------------------------------------------------------
SPAN_RTOL = 1e-8


def _design(X: np.ndarray) -> np.ndarray:
    X = np.asarray(X, dtype=float)
    return np.hstack([np.ones((X.shape[0], 1)), X])


def residualize(v: np.ndarray, X: np.ndarray) -> np.ndarray:
    """OLS residual of v on [1, X] via a rank-revealing least-squares solve."""
    A = _design(X)
    beta, *_ = np.linalg.lstsq(A, np.asarray(v, float), rcond=None)
    return np.asarray(v, float) - A @ beta


def relative_residual(v: np.ndarray, X: np.ndarray) -> float:
    """||resid(v | 1, X)|| / ||v - mean(v)||  (0 for v exactly in span)."""
    v = np.asarray(v, float)
    denom = np.linalg.norm(v - v.mean())
    if denom == 0:
        return 0.0
    return float(np.linalg.norm(residualize(v, X)) / denom)


def numerical_rank(X: np.ndarray, rtol: float = 1e-9) -> int:
    """Numerical rank of the column-centred, column-normalised matrix."""
    X = np.asarray(X, float)
    Xc = X - X.mean(0)
    norms = np.linalg.norm(Xc, axis=0)
    Xc = Xc[:, norms > 0] / norms[norms > 0]
    if Xc.shape[1] == 0:
        return 0
    s = np.linalg.svd(Xc, compute_uv=False)
    return int((s > s[0] * rtol).sum())


def span_certified_pcor(z: np.ndarray, X: np.ndarray, y: np.ndarray,
                        rtol: float = SPAN_RTOL) -> dict:
    """Partial correlation of z with y given X, with an exact-span certificate.

    Returns dict(pcor, rel_resid, in_span, df) where df = n - rank([1,X])
    (= n - k - 1 for k non-redundant conditioning columns); Fisher-z intervals
    then use SE = 1/sqrt(df - 2) = 1/sqrt(n - k - 3).  If z lies in
    span([1, X]) to relative tolerance ``rtol`` the partial correlation is
    reported as exactly 0.0 and ``in_span`` is True.
    """
    z = np.asarray(z, float); y = np.asarray(y, float)
    rr = relative_residual(z, X)
    n = len(z)
    rank = (numerical_rank(X) if X.shape[1] else 0) + 1  # + intercept
    df = n - rank
    if rr < rtol:
        return {"pcor": 0.0, "rel_resid": rr, "in_span": True, "df": df}
    zr = residualize(z, X); yr = residualize(y, X)
    if np.linalg.norm(yr) == 0:
        return {"pcor": 0.0, "rel_resid": rr, "in_span": False, "df": df}
    return {"pcor": float(np.corrcoef(zr, yr)[0, 1]), "rel_resid": rr,
            "in_span": False, "df": df}


def fisher_ci(r: float, df: int, level: float = 0.95) -> tuple[float, float]:
    """Fisher-z confidence interval for a (partial) correlation.

    ``df`` = n - rank([1, X]); the standard error of atanh(r) is
    1/sqrt(df - 2) = 1/sqrt(n - k - 3).
    """
    from scipy.stats import norm
    if df <= 3 or not np.isfinite(r):
        return (float("nan"), float("nan"))
    r = float(np.clip(r, -0.999999, 0.999999))
    z = np.arctanh(r); se = 1.0 / np.sqrt(df - 2)
    q = norm.ppf(0.5 + level / 2)
    return (float(np.tanh(z - q * se)), float(np.tanh(z + q * se)))
