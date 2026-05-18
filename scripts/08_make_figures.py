"""Produce publication-quality figures for Paper 2 (Wuzi family).

Reads CSVs from results/ and writes figures/{name}.png and figures/{name}.pdf.

Figures generated (numbering matches docs/paper2_wuzi_manuscript_skeleton.md):
  fig2_pca_scree                — PCA cumulative variance, 3 datasets overlaid
  fig3_redundancy_bars          — # redundant pairs at |r|>=0.90/0.95, 3 datasets
  fig4_octane_heatmap           — |r| of each index vs each property on octanes
  fig5_wuzi_param_heatmaps      — max|r|(α,β) for 4 γ slices on ESOL
  fig5b_wuzi_param_heatmaps_all — same but 3 datasets × 4 γ slices (supplementary)
  fig6_degeneracy_bars          — % degeneracy of indices on 106 trees(10)
  fig7_structure_sensitivity    — SS, Abr, SA on 75 decanes
  fig8_correlation_heatmap      — 30x30 baseline correlation matrix (ESOL)

Usage: python scripts/08_make_figures.py
"""
from __future__ import annotations
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

THIS_DIR = os.path.abspath(os.path.dirname(__file__))
PROJECT = os.path.abspath(os.path.join(THIS_DIR, ".."))
FIG_DIR = os.path.join(PROJECT, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

# Publication style
plt.rcParams.update({
    "font.family": "serif",
    "font.size": 10,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "figure.dpi": 100,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "axes.spines.top": False,
    "axes.spines.right": False,
})

DATASETS = [("esol", "ESOL  (n=1127)"), ("freesolv", "FreeSolv  (n=639)"),
            ("lipophilicity", "Lipophilicity  (n=4200)")]


def save_both(fig, name):
    fig.savefig(os.path.join(FIG_DIR, f"{name}.png"))
    fig.savefig(os.path.join(FIG_DIR, f"{name}.pdf"))
    plt.close(fig)
    print(f"  → {name}.{{png,pdf}}")


# ===== Figure 2: PCA scree =====

def fig2_pca_scree():
    fig, ax = plt.subplots(figsize=(6, 4))
    colors = ["#1f77b4", "#d62728", "#2ca02c"]
    for (ds, label), color in zip(DATASETS, colors):
        df = pd.read_csv(os.path.join(PROJECT, "results", ds, "pca_variance.csv"))
        ax.plot(df["component"], df["cumulative"], marker="o", markersize=4,
                label=label, color=color, linewidth=1.5)
    ax.axhline(0.95, color="gray", linestyle="--", linewidth=1, alpha=0.7)
    ax.text(28, 0.955, "95%", color="gray", fontsize=9)
    ax.axhline(0.99, color="gray", linestyle=":", linewidth=1, alpha=0.7)
    ax.text(28, 0.995, "99%", color="gray", fontsize=9)
    ax.set_xlabel("PCA component")
    ax.set_ylabel("Cumulative variance explained")
    ax.set_xlim(1, 30)
    ax.set_ylim(0.4, 1.02)
    ax.set_xticks([1, 3, 5, 10, 15, 20, 25, 30])
    ax.legend(loc="lower right", frameon=True)
    ax.set_title("Effective rank of the 30-index baseline\n"
                 "(3 components capture 95% of variance on all three datasets)",
                 fontsize=10)
    save_both(fig, "fig2_pca_scree")


# ===== Figure 3: Redundancy bars =====

def fig3_redundancy_bars():
    summary = pd.read_csv(os.path.join(PROJECT, "results", "cross_dataset_summary.csv"))
    fig, ax = plt.subplots(figsize=(6, 4))
    x = np.arange(len(summary))
    width = 0.35
    bars1 = ax.bar(x - width/2, summary["pairs_r_ge_090"], width,
                   label="|r| ≥ 0.90", color="#aec7e8", edgecolor="#1f77b4")
    bars2 = ax.bar(x + width/2, summary["pairs_r_ge_095"], width,
                   label="|r| ≥ 0.95", color="#d62728", edgecolor="#8c0000")
    ax.axhline(435, color="gray", linestyle=":", linewidth=1, alpha=0.6)
    ax.text(2.4, 440, "Total = 435 pairs", fontsize=8, color="gray", ha="right")
    ax.set_xticks(x)
    ax.set_xticklabels([s.replace("_", " ").title() for s in summary["dataset"]])
    ax.set_ylabel("Number of pairs")
    ax.set_ylim(0, 480)
    ax.set_title("Redundant pairs out of C(30, 2) = 435 baseline index pairs",
                 fontsize=10)
    ax.legend(loc="upper right")
    for bar in bars1:
        ax.annotate(f"{int(bar.get_height())}", xy=(bar.get_x() + bar.get_width()/2, bar.get_height()),
                    xytext=(0, 3), textcoords="offset points", ha="center", fontsize=8)
    for bar in bars2:
        ax.annotate(f"{int(bar.get_height())}", xy=(bar.get_x() + bar.get_width()/2, bar.get_height()),
                    xytext=(0, 3), textcoords="offset points", ha="center", fontsize=8)
    save_both(fig, "fig3_redundancy_bars")


# ===== Figure 4: Octane correlation heatmap =====

def fig4_octane_heatmap():
    df = pd.read_csv(os.path.join(PROJECT, "results", "octane_prediction.csv"))
    df = df.set_index("index")
    # Reorder: sort by max |r|
    df["__rank"] = df[["r_T_B", "r_dHf", "r_dHvap", "r_S", "r_omega"]].abs().max(axis=1)
    df = df.sort_values("__rank", ascending=False).drop(columns="__rank")
    fig, ax = plt.subplots(figsize=(4, 9))
    im = ax.imshow(df.abs().values, cmap="Reds", aspect="auto", vmin=0, vmax=1)
    ax.set_xticks(range(len(df.columns)))
    ax.set_xticklabels(["T_B", "ΔH_f", "ΔH_vap", "S", "ω"], rotation=0)
    ax.set_yticks(range(len(df.index)))
    ax.set_yticklabels(df.index, fontsize=7)
    for i in range(df.shape[0]):
        for j in range(df.shape[1]):
            v = df.iloc[i, j]
            color = "white" if abs(v) > 0.6 else "black"
            ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=6, color=color)
    cbar = plt.colorbar(im, ax=ax, fraction=0.04, pad=0.04)
    cbar.set_label("|r|")
    ax.set_title("Octane prediction:\n|r| of each index with each property (n=18)", fontsize=10)
    save_both(fig, "fig4_octane_heatmap")


# ===== Figure 5: Wuzi parameter heatmaps (4 γ slices, ESOL) =====

def _wuzi_heatmap_panel(ax, ds, gamma, vmin=0.9, vmax=1.0):
    df = pd.read_csv(os.path.join(PROJECT, "results", ds, "wuzi_grid.csv"))
    slice_df = (df[df["gamma"] == gamma]
                .pivot(index="beta", columns="alpha", values="max_abs_r_baseline"))
    slice_df = slice_df.sort_index(ascending=False)  # high β at top
    cmap = LinearSegmentedColormap.from_list(
        "rg", ["#2ca02c", "#fdae61", "#d62728", "#8c0000"])
    im = ax.imshow(slice_df.values, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(slice_df.columns)))
    ax.set_xticklabels([f"{v:g}" for v in slice_df.columns])
    ax.set_yticks(range(len(slice_df.index)))
    ax.set_yticklabels([f"{v:g}" for v in slice_df.index])
    ax.set_xlabel("α")
    ax.set_ylabel("β")
    ax.set_title(f"γ = {gamma:g}")
    for i in range(slice_df.shape[0]):
        for j in range(slice_df.shape[1]):
            v = slice_df.iloc[i, j]
            color = "white" if v > 0.97 else "black"
            ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=7, color=color)
    return im


def fig5_wuzi_param_heatmaps():
    fig, axes = plt.subplots(2, 2, figsize=(9, 7))
    gammas = [0.0, 0.5, 1.0, 2.0]
    for ax, gamma in zip(axes.flat, gammas):
        im = _wuzi_heatmap_panel(ax, "esol", gamma)
    fig.suptitle("Wuzi family kill-test on ESOL (n=1127):\n"
                 "max |r| with classical 30-index baseline at each (α, β, γ)",
                 fontsize=11, y=1.00)
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label("max |r|")
    cbar.ax.axhline(0.95, color="black", linewidth=1)
    fig.text(0.94, 0.50, "kill threshold",
             rotation=270, va="center", fontsize=8)
    plt.tight_layout(rect=[0, 0, 0.9, 0.96])
    save_both(fig, "fig5_wuzi_param_heatmaps")


def fig5b_wuzi_param_heatmaps_all():
    fig, axes = plt.subplots(3, 4, figsize=(13, 9))
    gammas = [0.0, 0.5, 1.0, 2.0]
    for row, (ds, label) in enumerate(DATASETS):
        for col, gamma in enumerate(gammas):
            im = _wuzi_heatmap_panel(axes[row, col], ds, gamma)
        axes[row, 0].set_ylabel(f"{label}\n\nβ", fontsize=9)
    fig.suptitle("Wuzi family kill-test across three QSPR datasets — "
                 "no parameter point passes |r| < 0.95 on any dataset",
                 fontsize=11, y=0.995)
    cbar_ax = fig.add_axes([0.93, 0.15, 0.012, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label("max |r| with baseline")
    plt.tight_layout(rect=[0, 0, 0.92, 0.97])
    save_both(fig, "fig5b_wuzi_param_heatmaps_all")


# ===== Figure 6: Degeneracy bars =====

def fig6_degeneracy_bars():
    df = pd.read_csv(os.path.join(PROJECT, "results", "wuzi_degeneracy.csv"))
    # Show all entries, sorted by degeneracy
    df = df.sort_values("degeneracy_pct").reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(8, 0.25 * len(df) + 1))
    colors = ["#d62728" if f == "wuzi" else "#1f77b4" for f in df["family"]]
    ax.barh(range(len(df)), df["degeneracy_pct"], color=colors,
            edgecolor="black", linewidth=0.3)
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df["index"], fontsize=7)
    ax.set_xlabel("% degeneracy on 106 non-isomorphic trees of order 10")
    ax.set_xlim(0, 100)
    ax.invert_yaxis()
    # Color legend
    from matplotlib.patches import Patch
    leg = [Patch(facecolor="#1f77b4", label="baseline (classical)"),
           Patch(facecolor="#d62728", label="Wuzi parameter triple")]
    ax.legend(handles=leg, loc="lower right", frameon=True)
    ax.set_title("Index degeneracy on trees(10): lower is more discriminating",
                 fontsize=10)
    save_both(fig, "fig6_degeneracy_bars")


# ===== Figure 7: Structure sensitivity =====

def fig7_structure_sensitivity():
    df = pd.read_csv(os.path.join(PROJECT, "results", "structure_sensitivity.csv"))
    df = df.sort_values("SS", ascending=False).head(20).reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(11, 5.5))
    x = np.arange(len(df))
    width = 0.27
    ax.bar(x - width, df["SS"],  width, label="SS",  color="#1f77b4")
    ax.bar(x,         df["Abr"], width, label="Abr", color="#ff7f0e")
    ax.bar(x + width, df["SA"],  width, label="SA",  color="#2ca02c")
    ax.set_xticks(x)
    ax.set_xticklabels(df["index"], rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Value")
    ax.set_title("Structure sensitivity metrics on 75 decane isomers "
                 "(top 20 by SS)", fontsize=10)
    ax.legend(loc="upper right")
    save_both(fig, "fig7_structure_sensitivity")


# ===== Figure 8: Full correlation heatmap (ESOL baseline) =====

def fig8_correlation_heatmap():
    corr = pd.read_csv(os.path.join(PROJECT, "results", "esol", "correlation_matrix.csv"),
                       index_col=0)
    fig, ax = plt.subplots(figsize=(9, 8))
    cmap = plt.get_cmap("RdBu_r")
    im = ax.imshow(corr.values, cmap=cmap, vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(len(corr.columns)))
    ax.set_xticklabels(corr.columns, rotation=90, fontsize=7)
    ax.set_yticks(range(len(corr.index)))
    ax.set_yticklabels(corr.index, fontsize=7)
    cbar = plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_label("Pearson r")
    ax.set_title("30-index baseline pairwise correlation matrix on ESOL "
                 "(n = 1127)\nRed = highly correlated; 159 of 435 pairs have |r| ≥ 0.95",
                 fontsize=10)
    save_both(fig, "fig8_correlation_heatmap")


def main():
    print("Generating publication figures →", FIG_DIR)
    fig2_pca_scree()
    fig3_redundancy_bars()
    fig4_octane_heatmap()
    fig5_wuzi_param_heatmaps()
    fig5b_wuzi_param_heatmaps_all()
    fig6_degeneracy_bars()
    fig7_structure_sensitivity()
    fig8_correlation_heatmap()
    print("Done.")


if __name__ == "__main__":
    main()
