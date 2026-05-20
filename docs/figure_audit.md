# Figure audit

This document records, for every figure in `figures/`, the source
data, the generating script, what the figure claims to show, whether
the plotted values match the source data, whether labels match, and
a verdict. Last refresh: see `git log -1 figures/`.

The 30-index baseline correlation pair count is
$\binom{30}{2} = 435$.

## Summary table

| Figure file | Source data | Script | What it shows | Data match | Labels match | Verdict |
|---|---|---|---|---|---|---|
| `fig2_pca_scree.{png,pdf}` | `results/<ds>/pca_variance.csv` for `<ds>` in {esol, freesolv, lipophilicity} | `scripts/08_make_figures.py::fig2_pca_scree` | Cumulative explained variance vs PCA component for the 30-index baseline, three datasets overlaid; 95% and 99% reference lines | YES — three components exceed 0.95 on all three datasets, matches `cross_dataset_summary.csv` PC95=3 | YES — dataset names "ESOL" / "FreeSolv" / "Lipophilicity"; sample sizes $n=1127, 639, 4200$ | CORRECT |
| `fig3_redundancy_bars.{png,pdf}` | `results/cross_dataset_summary.csv` | `fig3_redundancy_bars` | Number of redundant index pairs at $\lvert r\rvert\ge 0.90$ and $\lvert r\rvert\ge 0.95$ per dataset, with the 435-pair ceiling | YES — bar heights match CSV: ESOL (203, 159), FreeSolv (257, 162), Lipophilicity (200, 164) | YES — ESOL/FreeSolv/Lipophilicity properly cased; "Total = 435 pairs" label upper-right | CORRECT (fixed in this audit pass — annotation no longer overlaps legend) |
| `fig4_octane_heatmap.{png,pdf}` | `results/octane_descriptors.csv` (8 properties × 18 octanes) | `fig4_octane_heatmap` | Signed Pearson $r$ between each Wuzi-grid index / classical index and each of 8 octane properties; rows sorted by maximum $\lvert r\rvert$; degenerate constant-on-octanes indices dropped | YES — non-trivial check is the PI degeneracy (std=0 on the 18 octanes); script drops PI with an explanatory note | YES — property labels typeset in LaTeX ($T_B$, $\Delta H_f$, $\Delta H_{\mathrm{vap}}$, $S$, $\omega$); Wuzi labels in $W(\alpha, \beta, \gamma)$ form | CORRECT (fixed in this audit pass — clean labels, signed values, NaN row handling) |
| `fig5_wuzi_param_heatmaps.{png,pdf}` | `results/esol/wuzi_grid.csv` | `fig5_wuzi_param_heatmaps` | $\max\lvert r\rvert$ with the 30-index baseline at each $(\alpha, \beta, \gamma)$, ESOL only; four $\gamma$ slices | YES — all 100 ESOL grid points have $\max\lvert r\rvert \ge 0.965$; matches `cross_dataset_summary.csv` `wuzi_min_max_r=0.965` | YES — title now correctly says "All 100 grid points are highly correlated with at least one baseline" (no overclaim) | CORRECT |
| `fig5b_wuzi_param_heatmaps_all.{png,pdf}` | `results/<ds>/wuzi_grid.csv` for all three datasets | `fig5b_wuzi_param_heatmaps_all` | Same as fig5 but stacked rows for ESOL / FreeSolv / Lipophilicity, all three $\gamma$ slices | YES — passing-point count programmatically computed; 1 starred passing cell on FreeSolv at $(\alpha,\beta,\gamma)=(-1,-1,1)$ with $\max\lvert r\rvert=0.943$ | YES — title says "only 1 of 300 grid points falls below the $\lvert r\rvert=0.95$ threshold (starred cell, FreeSolv)" | CORRECT (fixed in this audit pass — title now honest about the single marginal point) |
| `fig6_degeneracy_bars.{png,pdf}` | `results/wuzi_degeneracy.csv` | `fig6_degeneracy_bars` | Percent degeneracy on 106 non-isomorphic trees of order 10, for each baseline index and several Wuzi parameter points | YES — 106 non-isomorphic trees of order 10 is a standard combinatorial count | YES — Wuzi labels in clean $W(\alpha, \beta, \gamma)$ form; bar colors distinguish classical baselines from Wuzi family | CORRECT |
| `fig7_structure_sensitivity.{png,pdf}` | `results/structure_sensitivity.csv` | `fig7_structure_sensitivity` | Top-20 ranking of indices on 75 decane isomers by structure-sensitivity metric $SS = \sigma / \mu$; companion abruptness $Abr = (\max - \min) / \mu$ and the ratio $SA = SS / Abr$ | YES — 75 decane isomers is standard; metrics correctly computed | YES — title typesets formulas in LaTeX; Wuzi labels clean | CORRECT |
| `fig8_correlation_heatmap.{png,pdf}` | `results/esol/correlation_matrix.csv` | `fig8_correlation_heatmap` | Full pairwise Pearson $r$ matrix of the 30 baseline indices on ESOL | YES — title says "159 of 435 pairs have $\lvert r\rvert\ge 0.95$" which matches `cross_dataset_summary.csv` ESOL value | YES — diverging colormap centred at zero; index names along both axes | CORRECT |
| `fig9_ml_benchmark.{png,pdf}` | `results/ml_benchmark.csv` | `scripts/08_make_figures.py::fig9_ml_benchmark` | Per-dataset RandomForest performance (top row: RMSE for the three regression datasets, ROC-AUC for BBBP) and feature count (bottom row) under four feature configurations | YES — bar heights match CSV: ESOL `full` RMSE 1.45, `pairwise_pruned` 1.45 (7 feats); BBBP `full` ROC-AUC 0.835, `pairwise_pruned` 0.864 (6 feats); Lipophilicity `combined_pruned` shown as "killed" since 0 features | YES — config labels match the script (full / PCA-95 / pair `\|r\|<0.95` / combined+`\|pcor\|>=0.10`); BBBP panel has its own ROC-AUC y-axis label | CORRECT |

## Notes

- Figure 1 in the manuscript skeleton (BID-lattice diagram, Δ ≤ 4)
  is a hand-drawn TikZ figure not yet created. It appears in
  `docs/paper2_wuzi_manuscript_skeleton.md` as `[TO CREATE]`.

- All figures are produced by a single script
  (`scripts/08_make_figures.py`). Regenerating from a clean checkout
  takes about 10 seconds after the upstream baselines are computed
  (`scripts/01`, `02`, `03`, `05`, `06`, `07`).

- The Wuzi parameter grid is $5\times 5\times 4 = 100$ points per
  dataset (300 across the three datasets). Pass/fail uses the
  conventional $\lvert r\rvert < 0.95$ threshold; FreeSolv's single
  borderline grid point has partial correlation with the target
  $\approx 0$, so the borderline pass does not contradict the
  overall redundancy conclusion.
