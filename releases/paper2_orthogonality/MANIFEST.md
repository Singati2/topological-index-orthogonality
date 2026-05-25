# Paper 2 Reproducibility Manifest

**Paper.** *Orthogonality Screening of Topological Indices for QSPR Modeling: A Structural and Empirical Redundancy Analysis*
**Authors.** G. Shiwakoti, A. Natarajan, M. Arockiaraj
**Repository.** <https://github.com/Singati2/topological-index-orthogonality>
**Tag.** `v1.4-paper2-orthogonality`

This file lists every script, dataset, results CSV, and figure
that Paper 2 depends on. A reviewer who clones the repository at
the tag `v1.4-paper2-orthogonality` and follows the steps below
should be able to reproduce every numerical claim and every figure
in the manuscript. The full reproduction (including BBBP
descriptor regeneration and the alt-family multi-dataset screen)
requires RDKit; the LASSO / ElasticNet / mRMR rows of
Table tab:fs_comparison can be re-run without RDKit using the
cached descriptor matrices.

The repository also hosts the companion graph-theory paper
(Paper 1, tagged `v1.4-paper1-wuzi`); Paper 2 cites Paper 1 for
the Wuzi-specific mathematical content. Paper 2 restates the
Edge-Degree-Pair Basis theorem with its own proof because the
theorem is used operationally in the screening argument.

## Source manuscript

- `docs/paper2_orthogonality_screening.tex` — full LaTeX source
  (bibliography embedded).

## Reference implementation

| File | Purpose |
|---|---|
| `src/standard_indices.py` | 30 baseline classical indices ($M_1, M_2, R, \chi, H, GA, ABC, \ldots$). |
| `src/orthogonality.py` | Correlation matrix, PCA, VIF, partial correlation. |
| `src/novel_candidates.py` | 20 alternative-family candidate indices. |
| `src/mol_to_graph.py` | SMILES → hydrogen-suppressed `networkx.Graph` (RDKit-backed). |
| `src/load_data.py` | ESOL / FreeSolv / Lipophilicity / BBBP loaders (MoleculeNet). |
| `src/wuzi_index.py` | Wuzi parametric family (shared with Paper 1, used here as a screened candidate). |

## Datasets

| Name | Target | Task | $n$ (post-RDKit-parse) | Source |
|---|---|---|---|---|
| ESOL          | $\log S$ | regression | 1,127 | Delaney 2004 via MoleculeNet |
| FreeSolv      | $\Delta G_{\mathrm{hyd}}$ | regression | 639   | Mobley 2014 via MoleculeNet |
| Lipophilicity | $\log D_{7.4}$ | regression | 4,200 | ChEMBL via MoleculeNet |
| BBBP          | blood-brain barrier penetration | binary classification | 2,039 | Martins 2012 via MoleculeNet |

All datasets download automatically via `src.load_data.load(name)`
on first use; cached CSVs are written to `data/` (gitignored). The
30-index baseline descriptors for ESOL / FreeSolv / Lipophilicity
are also cached at `results/<ds>/descriptors_baseline.csv` so the
LASSO comparator (script 15) can be re-run without RDKit.

## Reproducibility scripts

Run from the repository root. All scripts complete in seconds-to-minutes on a single core.

| Script | Output | Manuscript reference |
|---|---|---|
| `scripts/01_baseline_correlations.py --dataset <name>` | `results/<ds>/{correlation_matrix,pca_variance,vif,partial_corr_*,descriptors_baseline}.csv` | Table tab:baseline_redundancy + Section 6.1 baseline-redundancy block |
| `scripts/02_wuzi_grid_search.py --dataset <name>` | `results/<ds>/wuzi_grid.csv` | Table tab:wuzi_verdict_paper1 + §6.1 |
| `scripts/03_cross_dataset_summary.py` | `results/cross_dataset_summary.{csv,md}` | Table tab:baseline_redundancy + Table tab:threshold_sensitivity |
| `scripts/04_novel_candidates_test.py` | `results/novel_candidates_experiment/{novel_candidates.csv,novel_candidates_by_family.csv,novel_candidates_summary.txt}` | §6.2 (ESOL case study) |
| `scripts/08_make_figures.py` | `figures/{fig2,fig3,fig5b,fig8,fig9}.{png,pdf}` | All Paper 2 figures |
| `scripts/11_ml_benchmark.py` | `results/ml_benchmark.{csv,md}` | Table tab:ml_benchmark + Figure fig9 |
| `scripts/15_feature_selection_comparison.py` | `results/paper2_feature_selection_comparison.csv` | Table tab:fs_comparison (LASSO/ElasticNet/mRMR/Boruta on regression datasets; L1-Logistic/mRMR/Boruta on BBBP) |
| `scripts/16_novel_candidates_multidataset.py` | `results/novel_candidates_experiment/novel_candidates_multidataset.csv` | Cross-dataset alt-family replication on ESOL/FreeSolv/Lipophilicity (rdkit required) |
| `scripts/17_bbbp_descriptors_and_altfamily.py` | `results/bbbp/descriptors_baseline.csv`, `results/novel_candidates_experiment/novel_candidates_bbbp.csv` | BBBP baseline descriptors + BBBP row of Table tab:novel_combined |
| `scripts/18_pipeline_then_lasso.py` | `results/paper2_pipeline_then_lasso.csv` | Table tab:matched_class (pipeline-then-LASSO same-model-class benchmark + paired t-tests, §6.3.2) |

The threshold-sensitivity table (Table tab:threshold_sensitivity)
is extracted from `results/cross_dataset_summary.csv` (columns 16–27);
no separate script is required. The extracted clean form is
saved at `results/paper2_threshold_sensitivity.csv` for convenience.

## Figures referenced by the manuscript

| File | LaTeX label | Section |
|---|---|---|
| `figures/fig2_pca_scree.pdf` | `fig:pca` | §6.1 |
| `figures/fig3_redundancy_bars.pdf` | `fig:redundancy_bars` | §6.1 |
| `figures/fig5b_wuzi_param_heatmaps_all.pdf` | `fig:wuzi_grid_all` | §6.1 |
| `figures/fig8_correlation_heatmap.pdf` | `fig:corr_heatmap` | §6.1 |
| `figures/fig9_ml_benchmark.pdf` | `fig:ml_benchmark` | §6.3 |

## Test suite

The full pytest suite (`pytest -q` from the repo root) executes
130 tests: 124 covering Paper 1 mechanics (BID dimension bound,
30 baseline indices' identity tests, Wuzi reductions, ratio-bound
sanity, caterpillar predicate, closed-form values) plus 6 new
Paper-2-specific reproducibility tests in
`tests/test_paper2_reproducibility.py` that verify the schema and
sanity of `results/paper2_threshold_sensitivity.csv`,
`results/paper2_feature_selection_comparison.csv`, and
`results/novel_candidates_experiment/novel_candidates_multidataset_summary.csv`,
and check that LASSO RMSE is within $1.3\times$ of the
pairwise_pruned RandomForest RMSE on every regression dataset.
All 130 tests pass at this tag.

## Software environment

- Python 3.11+
- `networkx >= 3.0`, `numpy`, `pandas`, `scipy`, `scikit-learn`, `matplotlib`
- `rdkit` (required to regenerate `descriptors_baseline.csv` from
  SMILES, and to run scripts 01, 02, 04, 11, 16, 17; **not required** for
  scripts 03, 15, or for any pure-Python computation on the cached
  descriptor matrices)
- `mrmr_selection`, `Boruta` (required for the mRMR / Boruta rows of
  `scripts/15_feature_selection_comparison.py`; not required for the
  LASSO / ElasticNet / L1-Logistic rows)
- `pytest` for the test suite

See `requirements.txt` for exact versions.

## Build the PDF

```bash
cd docs
pdflatex paper2_orthogonality_screening.tex
pdflatex paper2_orthogonality_screening.tex   # second pass for cross-references
```

Bibliography is embedded as `\begin{thebibliography}`; no `bibtex`
pass needed. Compilation requires no extra LaTeX packages beyond
standard CTAN distributions (`amsmath`, `amssymb`, `amsthm`,
`mathtools`, `booktabs`, `array`, `graphicx`, `hyperref`, `xcolor`,
`enumitem`, `listings`).
