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

DATASET_DISPLAY = {"esol": "ESOL", "freesolv": "FreeSolv",
                   "lipophilicity": "Lipophilicity"}


def _wuzi_label(name: str) -> str:
    """Convert programmatic names like 'Wuzi_a1.0_b0.0_g0.0' or
    'Wuzi(α=1.0, β=0.0, γ=0.0)' to a clean publication-ready form
    W(1, 0, 0)."""
    if name.startswith("Wuzi_a"):
        # parse "Wuzi_a<α>_b<β>_g<γ>"
        rest = name[len("Wuzi_a"):]
        a_str, rest = rest.split("_b", 1)
        b_str, g_str = rest.split("_g", 1)
        a, b, g = float(a_str), float(b_str), float(g_str)
        return f"$W({a:g},\\, {b:g},\\, {g:g})$"
    if name.startswith("Wuzi(α="):
        s = name.replace("Wuzi(α=", "").replace(", β=", ",").replace(", γ=", ",").rstrip(")")
        a, b, g = (float(x) for x in s.split(","))
        return f"$W({a:g},\\, {b:g},\\, {g:g})$"
    return name


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
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    x = np.arange(len(summary))
    width = 0.35
    bars1 = ax.bar(x - width/2, summary["pairs_r_ge_090"], width,
                   label=r"$|r| \geq 0.90$", color="#aec7e8", edgecolor="#1f77b4")
    bars2 = ax.bar(x + width/2, summary["pairs_r_ge_095"], width,
                   label=r"$|r| \geq 0.95$", color="#d62728", edgecolor="#8c0000")
    ax.axhline(435, color="gray", linestyle=":", linewidth=1, alpha=0.7)
    ax.text(len(summary) - 0.55, 437, "Total = 435 pairs",
            fontsize=8, color="gray", ha="right", va="bottom")
    ax.set_xticks(x)
    ax.set_xticklabels([DATASET_DISPLAY[s] for s in summary["dataset"]])
    ax.set_ylabel("Number of pairs")
    ax.set_ylim(0, 480)
    ax.set_title(r"Redundant index pairs out of $\binom{30}{2}=435$ "
                 "baseline pairs", fontsize=10)
    ax.legend(loc="upper left", bbox_to_anchor=(0.01, 0.98), frameon=True)
    for bar in list(bars1) + list(bars2):
        ax.annotate(f"{int(bar.get_height())}",
                    xy=(bar.get_x() + bar.get_width()/2, bar.get_height()),
                    xytext=(0, 3), textcoords="offset points",
                    ha="center", fontsize=8)
    save_both(fig, "fig3_redundancy_bars")


# ===== Figure 4: Octane correlation heatmap =====

def fig4_octane_heatmap():
    df = pd.read_csv(os.path.join(PROJECT, "results", "octane_prediction.csv"))
    df = df.set_index("index")
    cols = ["r_T_B", "r_dHf", "r_dHvap", "r_S", "r_omega"]
    # Drop indices that are degenerate (constant) on the 18 octanes — their
    # correlation is mathematically undefined (NaN). Currently this is PI.
    finite = df[cols].notna().all(axis=1)
    degenerate = df.index[~finite].tolist()
    df = df.loc[finite].copy()
    # Reorder by max |r|
    df["__rank"] = df[cols].abs().max(axis=1)
    df = df.sort_values("__rank", ascending=False).drop(columns="__rank")
    # Relabel Wuzi entries
    df.index = [_wuzi_label(n) for n in df.index]

    fig, ax = plt.subplots(figsize=(5.2, 9))
    im = ax.imshow(df[cols].abs().values, cmap="Reds", aspect="auto",
                   vmin=0, vmax=1)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([r"$T_B$", r"$\Delta H_f$", r"$\Delta H_{\mathrm{vap}}$",
                        r"$S$", r"$\omega$"], rotation=0, fontsize=10)
    ax.set_yticks(range(len(df.index)))
    ax.set_yticklabels(df.index, fontsize=7)
    for i in range(df.shape[0]):
        for j, c in enumerate(cols):
            v = df.iloc[i][c]
            color = "white" if abs(v) > 0.6 else "black"
            ax.text(j, i, f"{v:+.2f}", ha="center", va="center",
                    fontsize=6, color=color)
    cbar = plt.colorbar(im, ax=ax, fraction=0.04, pad=0.04)
    cbar.set_label(r"$|r|$")
    note = ""
    if degenerate:
        note = (f"\nDegenerate (constant on the 18 octanes, |r| undefined): "
                f"{', '.join(degenerate)}.")
    ax.set_title("Octane prediction: signed Pearson $r$ "
                 "(color encodes $|r|$) — $n=18$" + note, fontsize=9)
    save_both(fig, "fig4_octane_heatmap")


# ===== Figure 5: Wuzi parameter heatmaps (4 γ slices, ESOL) =====

def _wuzi_heatmap_panel(ax, ds, gamma, vmin=0.9, vmax=1.0):
    df = pd.read_csv(os.path.join(PROJECT, "results", ds, "wuzi_grid.csv"))
    slice_df = (df[df["gamma"] == gamma]
                .pivot(index="beta", columns="alpha", values="max_abs_r_baseline"))
    slice_df = slice_df.sort_index(ascending=False)  # high β at top
    # Diverging at the |r|=0.95 screening threshold: green = below, red = above
    cmap = LinearSegmentedColormap.from_list(
        "screen", ["#2ca02c", "#a1d99b", "#fdae61", "#d62728", "#8c0000"])
    im = ax.imshow(slice_df.values, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(slice_df.columns)))
    ax.set_xticklabels([f"{v:g}" for v in slice_df.columns])
    ax.set_yticks(range(len(slice_df.index)))
    ax.set_yticklabels([f"{v:g}" for v in slice_df.index])
    ax.set_xlabel(r"$\alpha$")
    ax.set_ylabel(r"$\beta$")
    ax.set_title(rf"$\gamma = {gamma:g}$")
    for i in range(slice_df.shape[0]):
        for j in range(slice_df.shape[1]):
            v = slice_df.iloc[i, j]
            color = "white" if v > 0.97 else "black"
            # Mark the one passing point in bold/italic
            if v < 0.95:
                ax.text(j, i, f"{v:.2f}*", ha="center", va="center",
                        fontsize=7, color="white", fontweight="bold")
            else:
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        fontsize=7, color=color)
    return im


def fig5_wuzi_param_heatmaps():
    fig, axes = plt.subplots(2, 2, figsize=(9, 7))
    gammas = [0.0, 0.5, 1.0, 2.0]
    for ax, gamma in zip(axes.flat, gammas):
        im = _wuzi_heatmap_panel(ax, "esol", gamma)
    fig.suptitle(r"Wuzi family — redundancy screen on ESOL ($n=1127$): "
                 r"$\max|r|$ with the 30-index"
                 "\n"
                 r"baseline at each $(\alpha,\beta,\gamma)$. "
                 r"All 100 grid points are highly correlated with at least one baseline ($|r|\geq 0.95$).",
                 fontsize=10.5, y=1.00)
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label(r"$\max|r|$")
    # Reference the |r|=0.95 screening threshold inside the colorbar
    cbar.ax.plot([0, 1], [0.95, 0.95], color="black", linewidth=1.2,
                 transform=cbar.ax.transAxes if False else cbar.ax.transData)
    fig.text(0.945, 0.50, "screening threshold (0.95)",
             rotation=270, va="center", fontsize=8)
    plt.tight_layout(rect=[0, 0, 0.9, 0.95])
    save_both(fig, "fig5_wuzi_param_heatmaps")


def fig5b_wuzi_param_heatmaps_all():
    # Honest count of passing points across all three datasets
    total_pass = 0
    total_pts = 0
    pass_locations = []
    for ds, _ in DATASETS:
        d = pd.read_csv(os.path.join(PROJECT, "results", ds, "wuzi_grid.csv"))
        p = d[d["screen_verdict"] == "PASS"]
        total_pass += len(p)
        total_pts += len(d)
        for _, r in p.iterrows():
            pass_locations.append(
                f"{DATASET_DISPLAY[ds]} at $(\\alpha,\\beta,\\gamma)=({r['alpha']:g},{r['beta']:g},{r['gamma']:g})$"
                f": $\\max|r|={r['max_abs_r_baseline']:.3f}$, "
                f"partial $r$ w/ target $={r['partial_corr_target']:.3f}$"
            )

    fig, axes = plt.subplots(3, 4, figsize=(13, 9))
    gammas = [0.0, 0.5, 1.0, 2.0]
    for row, (ds, label) in enumerate(DATASETS):
        for col, gamma in enumerate(gammas):
            im = _wuzi_heatmap_panel(axes[row, col], ds, gamma)
        axes[row, 0].set_ylabel(f"{label}\n\n" + r"$\beta$", fontsize=9)
    # Honest title: report exact pass count and any marginal cases
    if total_pass == 0:
        sup = (rf"Wuzi family — redundancy screen across three QSPR datasets: "
               rf"0 of {total_pts} grid points fall below the $|r|=0.95$ threshold.")
    elif total_pass == 1:
        sup = (rf"Wuzi family — redundancy screen across three QSPR datasets: "
               rf"only 1 of {total_pts} grid points falls below the "
               rf"$|r|=0.95$ threshold (starred cell, FreeSolv).")
    else:
        sup = (rf"Wuzi family — redundancy screen across three QSPR datasets: "
               rf"{total_pass} of {total_pts} grid points fall below the $|r|=0.95$ threshold (starred cells).")
    fig.suptitle(sup, fontsize=10.5, y=0.995)
    cbar_ax = fig.add_axes([0.93, 0.15, 0.012, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label(r"$\max|r|$ with baseline")
    plt.tight_layout(rect=[0, 0, 0.92, 0.97])
    save_both(fig, "fig5b_wuzi_param_heatmaps_all")
    # Print honest pass details so they can also go in the manuscript caption
    if pass_locations:
        print("    Marginal passing points (starred in figure):")
        for loc in pass_locations:
            print(f"      - {loc}")


# ===== Figure 6: Degeneracy bars =====

def fig6_degeneracy_bars():
    df = pd.read_csv(os.path.join(PROJECT, "results", "wuzi_degeneracy.csv"))
    df["display"] = [_wuzi_label(n) for n in df["index"]]
    df = df.sort_values("degeneracy_pct").reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(8, 0.25 * len(df) + 1))
    colors = ["#d62728" if f == "wuzi" else "#1f77b4" for f in df["family"]]
    ax.barh(range(len(df)), df["degeneracy_pct"], color=colors,
            edgecolor="black", linewidth=0.3)
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df["display"], fontsize=7)
    ax.set_xlabel("% degeneracy on 106 non-isomorphic trees of order 10")
    ax.set_xlim(0, 100)
    ax.invert_yaxis()
    from matplotlib.patches import Patch
    leg = [Patch(facecolor="#1f77b4", label="baseline (classical)"),
           Patch(facecolor="#d62728", label=r"Wuzi $W(\alpha,\beta,\gamma)$")]
    ax.legend(handles=leg, loc="lower right", frameon=True)
    ax.set_title("Index degeneracy on trees of order 10\n"
                 "(lower = better discrimination among 106 distinct structures)",
                 fontsize=10)
    save_both(fig, "fig6_degeneracy_bars")


# ===== Figure 7: Structure sensitivity =====

def fig7_structure_sensitivity():
    df = pd.read_csv(os.path.join(PROJECT, "results", "structure_sensitivity.csv"))
    df["display"] = [_wuzi_label(n) for n in df["index"]]
    df = df.sort_values("SS", ascending=False).head(20).reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(11, 5.5))
    x = np.arange(len(df))
    width = 0.27
    ax.bar(x - width, df["SS"],  width, label=r"$SS$",  color="#1f77b4")
    ax.bar(x,         df["Abr"], width, label=r"$Abr$", color="#ff7f0e")
    ax.bar(x + width, df["SA"],  width, label=r"$SA$",  color="#2ca02c")
    ax.set_xticks(x)
    ax.set_xticklabels(df["display"], rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Value")
    ax.set_title(r"Structure sensitivity on 75 decane isomers — "
                 r"top 20 by $SS$ "
                 r"($SS=\sigma/\mu$, $Abr=(\max-\min)/\mu$, $SA=SS/Abr$)",
                 fontsize=9.5)
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


def fig9_ml_benchmark():
    """ML benchmark figure: per-dataset CV performance vs feature
    configuration, plus the corresponding feature count.

    Two-row layout: top row plots the primary task-appropriate metric
    (RMSE for regression, ROC-AUC for classification) as mean across
    5 CV folds with error bars; bottom row plots the mean feature
    count. One column per dataset.
    """
    df = pd.read_csv(os.path.join(PROJECT, "results", "ml_benchmark.csv"))
    config_order = ["full", "pca_95", "pairwise_pruned", "combined_pruned"]
    config_label = {
        "full": "full\n(30)",
        "pca_95": "PCA-95",
        "pairwise_pruned": "pair\n|r|<0.95",
        "combined_pruned": "combined\n+|pcor|>=0.10",
    }
    config_color = {
        "full":            "#4c72b0",
        "pca_95":          "#dd8452",
        "pairwise_pruned": "#55a868",
        "combined_pruned": "#c44e52",
    }
    datasets = ["esol", "freesolv", "lipophilicity", "bbbp"]
    metric_for = {
        "regression": ("rmse_mean", "rmse_std", "RMSE (lower is better)"),
        "classification": ("roc_auc_mean", "roc_auc_std",
                           "ROC-AUC (higher is better)"),
    }
    display_name = {
        "esol": "ESOL", "freesolv": "FreeSolv",
        "lipophilicity": "Lipophilicity", "bbbp": "BBBP",
    }

    fig, axes = plt.subplots(2, 4, figsize=(14, 6.8),
                             gridspec_kw={"hspace": 0.55, "wspace": 0.40})
    for col, ds in enumerate(datasets):
        sub = df[df["dataset"] == ds].set_index("config")
        task = sub["task_type"].iloc[0]
        mean_col, std_col, metric_label = metric_for[task]
        xs = np.arange(len(config_order))
        means = [sub.loc[c, mean_col] if c in sub.index else np.nan
                 for c in config_order]
        stds = [sub.loc[c, std_col] if c in sub.index else np.nan
                for c in config_order]
        colors = [config_color[c] for c in config_order]
        ax_top = axes[0, col]
        # Bars + error bars. NaN bars get skipped automatically.
        ax_top.bar(xs, means, color=colors, edgecolor="black", linewidth=0.5,
                   yerr=[0 if np.isnan(s) else s for s in stds],
                   capsize=3, ecolor="black")
        for x, m, s in zip(xs, means, stds):
            if np.isnan(m):
                ax_top.text(x, 0.5, "killed\n(see\nbelow)",
                            ha="center", va="center",
                            fontsize=7, color="#c44e52",
                            transform=ax_top.get_xaxis_transform())
            else:
                offset = s if not np.isnan(s) else 0
                ax_top.text(x, m + offset, f"{m:.3f}",
                            ha="center", va="bottom", fontsize=8)
        ax_top.set_xticks(xs)
        ax_top.set_xticklabels([config_label[c] for c in config_order],
                               fontsize=7)
        ax_top.set_title(f"{display_name[ds]}\n({task})", fontsize=10)
        ax_top.set_ylabel(metric_label, fontsize=9)
        if task == "classification":
            ax_top.set_ylim(0.5, 1.0)

        # Bottom row: mean feature count with min-max range as error bar.
        ax_bot = axes[1, col]
        nf_mean = [sub.loc[c, "n_features_mean"] if c in sub.index else 0
                   for c in config_order]
        nf_min = [sub.loc[c, "n_features_min"] if c in sub.index else 0
                  for c in config_order]
        nf_max = [sub.loc[c, "n_features_max"] if c in sub.index else 0
                  for c in config_order]
        # Error bar: from min to max
        lower = [m - mn for m, mn in zip(nf_mean, nf_min)]
        upper = [mx - m for m, mx in zip(nf_mean, nf_max)]
        ax_bot.bar(xs, nf_mean, color=colors, edgecolor="black",
                   linewidth=0.5,
                   yerr=[lower, upper], capsize=3, ecolor="black")
        for x, nf, mn, mx in zip(xs, nf_mean, nf_min, nf_max):
            if mn == mx:
                label = f"{int(nf)}"
            else:
                label = f"{nf:.1f}\n({mn}-{mx})"
            ax_bot.text(x, mx + 0.5, label, ha="center", va="bottom",
                        fontsize=7)
        ax_bot.set_xticks(xs)
        ax_bot.set_xticklabels([config_label[c] for c in config_order],
                               fontsize=7)
        ax_bot.set_ylim(0, 35)
        if col == 0:
            ax_bot.set_ylabel("Mean feature count\n(error bar: min-max across folds)",
                              fontsize=8)
    fig.suptitle(
        "ML benchmark across 5-fold cross-validation: redundancy screening "
        "preserves RandomForest performance while reducing input dimension\n"
        "Top: primary metric mean +/- std across folds.  "
        "Bottom: mean feature count, error bar = min-max range across folds.",
        fontsize=10.5, y=1.00,
    )
    save_both(fig, "fig9_ml_benchmark")


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
    fig9_ml_benchmark()
    print("Done.")


if __name__ == "__main__":
    main()
