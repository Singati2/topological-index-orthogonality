"""
Figures for Paper 2 (v2.1), generated from results/.

fig1  (a) conventional pairwise screen (max |r| with the 30-index baseline)
          versus the certificate (relative residual on the baseline, log scale)
          for the published BID indices, the Loyola grid and the non-BID controls;
      (b) forest plot of the candidate-dataset pairs of Table nonbid: linear
          partial correlation (hollow) and non-linear rho_nl (filled) with 95%
          intervals, the null band from results/nonlinear_noise_floor_summary.csv
          (99.5th percentile of |rho_nl| under the in-span null) if present.
fig2  information ceiling: per dataset, overall target SD, in-sample
          graph-oracle RMSE, best topology-only model, RDKit-2D reference.
"""
from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

P = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
R = os.path.join(P, "results")
OUT = os.path.join(P, "figures")
os.makedirs(OUT, exist_ok=True)
DS = ["esol", "freesolv", "lipophilicity", "bbbp"]
NICE = {"esol": "ESOL", "freesolv": "FreeSolv", "lipophilicity": "Lipophilicity", "bbbp": "BBBP"}
plt.rcParams.update({"font.size": 8.5, "axes.spines.top": False, "axes.spines.right": False})


def fig1():
    v = pd.read_csv(os.path.join(R, "bid_variant_screen.csv"))
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 3.5), gridspec_kw={"width_ratios": [1, 1.25]})
    style = {"published BID": ("o", "tab:blue", "published BID index (15)"),
             "Loyola grid": (".", "tab:gray", "Loyola grid (100 points)"),
             "non-BID control": ("s", "tab:red", "non-BID control (2)")}
    for fam, (mk, col, lab) in style.items():
        s = v[v.family == fam]
        a.scatter(s.max_abs_r_baseline30, np.log10(s.rel_resid_on_baseline30.clip(lower=1e-16)),
                  marker=mk, s=14 if fam != "Loyola grid" else 8, c=col, alpha=0.75, label=lab, linewidths=0)
    a.axvline(0.95, ls=":", c="k", lw=0.8); a.axhline(-8, ls="--", c="k", lw=0.8)
    a.text(0.9515, -10.5, "pairwise threshold", fontsize=6.5, rotation=90, va="top")
    a.text(0.893, -8.5, "certificate tolerance", fontsize=6.5, va="top")
    a.set_xlabel(r"max $|r|$ with the 30-index baseline (pairwise screen)")
    a.set_ylabel(r"$\log_{10}$ relative residual on baseline (certificate)")
    a.set_xlim(0.89, 1.005); a.set_ylim(-15.5, 0)
    a.legend(frameon=False, fontsize=6.5, loc="upper left", bbox_to_anchor=(0.0, 0.72))
    a.set_title("(a) pairwise screen vs span certificate", fontsize=9, loc="left")

    c = pd.read_csv(os.path.join(R, "novel_candidates_experiment", "novel_candidates_multidataset.csv"))
    n = pd.read_csv(os.path.join(R, "nonlinear_conditional_tests.csv"))
    m = c.merge(n, on=["dataset", "candidate"])
    sel = m[(m.partial_corr_target.abs() >= 0.10) | (m.rho_nl.abs() >= 0.10)].copy()
    sel["o"] = sel.dataset.map({d: i for i, d in enumerate(DS)}); sel = sel.sort_values(["o", "candidate"]).reset_index(drop=True)
    band = 0.05
    fp = os.path.join(R, "nonlinear_noise_floor_summary.csv")
    if os.path.exists(fp):
        s = pd.read_csv(fp); band = float(s[s.null == "B_inspan"].p99_abs.max())
    y = np.arange(len(sel))[::-1]
    b.axvspan(-band, band, color="0.9", lw=0, label=f"in-span null band (±{band:.3f})")
    b.axvline(0, c="k", lw=0.6)
    b.errorbar(sel.partial_corr_target, y + 0.18, xerr=[sel.partial_corr_target - sel.pcor_ci_lo, sel.pcor_ci_hi - sel.partial_corr_target],
               fmt="o", mfc="white", mec="tab:blue", ecolor="tab:blue", ms=4, lw=1, label="linear pcor (OLS residuals)")
    b.errorbar(sel.rho_nl, y - 0.18, xerr=[sel.rho_nl - sel.ci_lo, sel.ci_hi - sel.rho_nl],
               fmt="o", c="tab:orange", ms=4, lw=1, label=r"$\rho_{\rm nl}$ (RF residuals)")
    b.set_yticks(y); b.set_yticklabels([f"{NICE[d]}: {k.replace('_', ' ')}" for d, k in zip(sel.dataset, sel.candidate)], fontsize=7)
    b.set_xlabel("residual association with target (95% interval)")
    b.legend(frameon=False, fontsize=6.5, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=3); b.set_title("(b) non-BID candidates: linear vs non-linear", fontsize=9, loc="left")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"fig1_certificate_and_screen.{ext}"), dpi=300)


def fig2():
    cur = pd.read_csv(os.path.join(R, "benchmark_v2_curation.csv")).set_index("dataset")
    sm = pd.read_csv(os.path.join(R, "benchmark_v2_summary.csv"))
    orc = pd.read_csv(os.path.join(R, "graph_oracle_floor.csv")).set_index("dataset")
    fig, axes = plt.subplots(1, 4, figsize=(7.2, 2.9))
    for ax, d in zip(axes, DS):
        s = sm[sm.dataset == d].set_index("model")["mean"]
        topo_models = [k for k in s.index if k.startswith(("a_", "b_", "c_", "d_", "e_", "g_"))]
        if d == "bbbp":
            vals = [0.5, max(s[topo_models]), s["f_rdkit2d_rf"]]
            labs = ["chance", "topology\n(best, CV)", "RDKit-2D\n(CV)"]; cols = ["0.6", "tab:blue", "tab:green"]
            ax.set_ylabel("ROC-AUC"); ax.set_ylim(0.45, 1.02)
        else:
            vals = [cur.loc[d, "overall_target_sd"], orc.loc[d, "oracle_in_sample_rmse"], min(s[topo_models]), s["f_rdkit2d_rf"]]
            labs = ["target\nSD", "graph\noracle", "topology\n(best, CV)", "RDKit-2D\n(CV)"]; cols = ["0.6", "tab:gray", "tab:blue", "tab:green"]
            ax.set_ylabel("RMSE")
        ax.bar(range(len(vals)), vals, color=cols, width=0.7)
        for i, v in enumerate(vals):
            ax.text(i, v, f"{v:.2f}", ha="center", va="bottom", fontsize=7)
        ax.set_xticks(range(len(vals))); ax.set_xticklabels([l.replace("\n", " ") for l in labs], fontsize=6, rotation=30, ha="right")
        ax.set_title(f"{NICE[d]}\n{cur.loc[d, 'pct_mols_graph_shared']:.0f}% of molecules share a graph", fontsize=7.5)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"fig2_information_ceiling.{ext}"), dpi=300)


if __name__ == "__main__":
    fig1(); fig2(); print("figures written to", OUT)
