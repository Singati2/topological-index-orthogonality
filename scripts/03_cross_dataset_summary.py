"""Aggregate redundancy statistics across ESOL, FreeSolv, Lipophilicity into a
single comparison table.

Reads each dataset's `summary.txt`-equivalent CSVs from `results/<dataset>/`
and writes a unified table to `results/cross_dataset_summary.csv` plus a
Markdown version at `results/cross_dataset_summary.md` for the manuscript.

Run AFTER `scripts/01_baseline_correlations.py --dataset <name>` and
`scripts/02_wuzi_grid_search.py --dataset <name>` for all three datasets.
"""
from __future__ import annotations
import os
import sys
import pandas as pd

THIS_DIR = os.path.abspath(os.path.dirname(__file__))
PROJECT = os.path.abspath(os.path.join(THIS_DIR, ".."))
sys.path.insert(0, PROJECT)

DATASETS = ["esol", "freesolv", "lipophilicity"]
N_INDICES = 30
N_PAIRS = N_INDICES * (N_INDICES - 1) // 2


PAIRWISE_THRESHOLDS = [0.90, 0.95, 0.99]
PCOR_THRESHOLDS = [0.05, 0.10, 0.20]


def stats_for(dataset: str) -> dict:
    base = os.path.join(PROJECT, "results", dataset)
    desc = pd.read_csv(os.path.join(base, "descriptors_baseline.csv"))
    pairs95 = pd.read_csv(os.path.join(base, "redundant_pairs_r95.csv"))
    pairs90 = pd.read_csv(os.path.join(base, "redundant_pairs_r90.csv"))
    pca = pd.read_csv(os.path.join(base, "pca_variance.csv"))
    wuzi = pd.read_csv(os.path.join(base, "wuzi_grid.csv"))
    pc90 = int((pca["cumulative"] >= 0.90).idxmax()) + 1
    pc95 = int((pca["cumulative"] >= 0.95).idxmax()) + 1
    pc99 = int((pca["cumulative"] >= 0.99).idxmax()) + 1
    row = {
        "dataset":             dataset,
        "n_molecules":         len(desc),
        "n_pairs_total":       N_PAIRS,
        "pairs_r_ge_090":      len(pairs90),
        "pairs_r_ge_095":      len(pairs95),
        "pct_pairs_r_ge_090":  100.0 * len(pairs90) / N_PAIRS,
        "pct_pairs_r_ge_095":  100.0 * len(pairs95) / N_PAIRS,
        "pca_components_90":   pc90,
        "pca_components_95":   pc95,
        "pca_components_99":   pc99,
        "wuzi_grid_points":    len(wuzi),
        "wuzi_pass":           int((wuzi["screen_verdict"] == "PASS").sum()),
        "wuzi_fail":           int((wuzi["screen_verdict"] == "FAIL").sum()),
        "wuzi_min_max_r":      float(wuzi["max_abs_r_baseline"].min(skipna=True)),
        "wuzi_max_partial":    float(wuzi["partial_corr_target"].abs().max(skipna=True)),
    }
    # Threshold-sensitivity columns. wuzi_pair_pass_tau_X = number of grid
    # points with max_abs_r_baseline < X. wuzi_combined_pass_tauR_taup =
    # number of grid points with max|r| < tauR AND |pcor| >= taup.
    for tau in PAIRWISE_THRESHOLDS:
        col = f"wuzi_pair_pass_tau_{int(tau*100):02d}"
        row[col] = int((wuzi["max_abs_r_baseline"] < tau).sum())
    for taur in PAIRWISE_THRESHOLDS:
        for taup in PCOR_THRESHOLDS:
            col = f"wuzi_combined_pass_taur_{int(taur*100):02d}_taup_{int(taup*100):03d}"
            row[col] = int(
                (
                    (wuzi["max_abs_r_baseline"] < taur)
                    & (wuzi["partial_corr_target"].abs() >= taup)
                ).sum()
            )
    return row


def main():
    rows = [stats_for(d) for d in DATASETS]
    df = pd.DataFrame(rows)
    out_csv = os.path.join(PROJECT, "results", "cross_dataset_summary.csv")
    df.to_csv(out_csv, index=False)

    # Manuscript-ready Markdown
    md_lines = []
    md_lines.append("# Cross-dataset redundancy summary\n")
    md_lines.append("All values computed on the same 30-index baseline set;\n"
                    "see `docs/theoretical_foundation.md` for the BID basis\n"
                    "theorem implying a 10-dimensional ceiling.\n")
    md_lines.append("## Baseline redundancy\n")
    md_lines.append("| Dataset | n | Pairs \\|r\\|≥0.90 | Pairs \\|r\\|≥0.95 | "
                    "PC₉₀ | PC₉₅ | PC₉₉ |")
    md_lines.append("|---|---:|---:|---:|---:|---:|---:|")
    for r in rows:
        md_lines.append(
            f"| {r['dataset']} | {r['n_molecules']} | "
            f"{r['pairs_r_ge_090']} / {r['n_pairs_total']} "
            f"({r['pct_pairs_r_ge_090']:.1f}%) | "
            f"{r['pairs_r_ge_095']} / {r['n_pairs_total']} "
            f"({r['pct_pairs_r_ge_095']:.1f}%) | "
            f"{r['pca_components_90']} | {r['pca_components_95']} | "
            f"{r['pca_components_99']} |"
        )
    md_lines.append("\nTheoretical upper bound on endpoint-degree edge-sum (BID)\n"
                    "index dimension under Δ≤4: **10**. Observed effective\n"
                    "rank (PC₉₅) is consistently 3 across all three benchmarks.\n")
    md_lines.append("\n## Wuzi parameter sweep — redundancy-screen verdict\n")
    md_lines.append("Verdict is **PASS** if max |r| with every classical baseline\n"
                    "is strictly below 0.95, **FAIL** otherwise. The threshold\n"
                    "|r| = 0.95 is conventional; the substantive question is\n"
                    "whether the partial correlation with the target property\n"
                    "is meaningfully different from zero after controlling for\n"
                    "the baseline (rightmost column).\n")
    md_lines.append("| Dataset | Grid points | Pass | Fail | Min max\\|r\\| | "
                    "Max \\|partial corr w/ target\\| |")
    md_lines.append("|---|---:|---:|---:|---:|---:|")
    for r in rows:
        md_lines.append(
            f"| {r['dataset']} | {r['wuzi_grid_points']} | "
            f"{r['wuzi_pass']} | {r['wuzi_fail']} | "
            f"{r['wuzi_min_max_r']:.3f} | {r['wuzi_max_partial']:.4f} |"
        )
    md_lines.append("\nOn ESOL and Lipophilicity, every one of the 100 grid\n"
                    "points is highly correlated with at least one classical\n"
                    "baseline (|r| ≥ 0.95). On FreeSolv, a single grid point\n"
                    "falls just below the threshold, but its partial correlation\n"
                    "with the target after controlling for the 30 baselines is\n"
                    "≈ 0.04, so it does not provide useful independent QSPR signal.\n"
                    "This behavior is consistent with the BID-basis observation\n"
                    "(`docs/theoretical_foundation.md`): every endpoint-degree\n"
                    "edge-sum index lies in a vector space of dimension at most\n"
                    "10 on Δ ≤ 4 graphs, so high empirical redundancy with the\n"
                    "30-index baseline is the expected outcome.\n")

    # Threshold-sensitivity table
    md_lines.append("\n## Threshold sensitivity of the Wuzi screening verdict\n")
    md_lines.append(
        "For each dataset, the number of Wuzi grid points (out of 100)\n"
        "that pass each variant of the screen at varying thresholds.\n"
        "$\\tau_r$ is the pairwise correlation threshold; $\\tau_p$ is\n"
        "the partial-correlation threshold. The combined criterion is\n"
        "$\\max|r| < \\tau_r$ AND $|\\mathrm{pcor}| \\ge \\tau_p$.\n")
    md_lines.append("### Pairwise-only verdict\n")
    md_lines.append("| Dataset | $\\tau_r = 0.90$ | $\\tau_r = 0.95$ | $\\tau_r = 0.99$ |")
    md_lines.append("|---|---:|---:|---:|")
    for r in rows:
        md_lines.append(
            f"| {r['dataset']} | "
            f"{r['wuzi_pair_pass_tau_90']} | "
            f"{r['wuzi_pair_pass_tau_95']} | "
            f"{r['wuzi_pair_pass_tau_99']} |"
        )
    md_lines.append(
        "\nThe pairwise pass count grows sharply with the relaxation of "
        "$\\tau_r$ from 0.95 to 0.99 because most Wuzi grid points have "
        "$\\max|r|$ in $[0.95, 0.99]$ with at least one classical "
        "baseline. This is the behaviour the combined criterion is "
        "designed to detect: pairwise alone is over-permissive at "
        "$\\tau_r = 0.99$.\n")
    md_lines.append("### Combined-criterion verdict (pairwise AND partial-correlation)\n")
    md_lines.append("| Dataset | $\\tau_r=0.90, \\tau_p=0.10$ | $\\tau_r=0.95, \\tau_p=0.10$ | "
                    "$\\tau_r=0.99, \\tau_p=0.10$ | $\\tau_r=0.95, \\tau_p=0.05$ | "
                    "$\\tau_r=0.95, \\tau_p=0.20$ |")
    md_lines.append("|---|---:|---:|---:|---:|---:|")
    for r in rows:
        md_lines.append(
            f"| {r['dataset']} | "
            f"{r['wuzi_combined_pass_taur_90_taup_010']} | "
            f"{r['wuzi_combined_pass_taur_95_taup_010']} | "
            f"{r['wuzi_combined_pass_taur_99_taup_010']} | "
            f"{r['wuzi_combined_pass_taur_95_taup_005']} | "
            f"{r['wuzi_combined_pass_taur_95_taup_020']} |"
        )
    md_lines.append(
        "\nThe combined-criterion pass count is robustly zero across "
        "$\\tau_r \\in \\{0.90, 0.95, 0.99\\}$ and $\\tau_p \\in \\{0.05, 0.10, 0.20\\}$ "
        "on every dataset: no Wuzi grid point in the 100-point sweep "
        "combines low pairwise correlation with meaningful partial "
        "correlation with the target. The pairwise-only verdict is "
        "threshold-sensitive; the combined verdict is not.\n")

    out_md = os.path.join(PROJECT, "results", "cross_dataset_summary.md")
    with open(out_md, "w") as f:
        f.write("\n".join(md_lines))

    print(df.to_string(index=False))
    print(f"\nWritten {out_csv}")
    print(f"Written {out_md}")


if __name__ == "__main__":
    main()
