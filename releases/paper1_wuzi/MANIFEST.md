# Paper 1 Reproducibility Manifest

**Paper.** *The Wuzi Index Family: Graph-Theoretic Properties, Bounds, Extremal Graphs, and Sensitivity Analysis*
**Authors.** A. Natarajan, G. Shiwakoti, M. Arockiaraj
**Repository.** <https://github.com/Singati2/topological-index-orthogonality>
**Tag.** `v1.4-paper1-wuzi`

This file lists every script, dataset, results CSV, and figure that
this paper depends on. A reviewer who clones the repository at the
tag `v1.4-paper1-wuzi` and follows the steps below should be able
to reproduce every numerical claim and every figure in the
manuscript.

The repository is named `topological-index-orthogonality` because
it also hosts the companion methodology paper (Paper 2). Paper 1
uses only the subset listed here.

## Source manuscript

- `docs/paper1_wuzi.tex` — the full LaTeX source (bibliography embedded).

## Reference implementation

| File | Purpose |
|---|---|
| `src/wuzi_index.py` | The Wuzi index function $W(G;\alpha,\beta,\gamma)$. |
| `src/standard_indices.py` | The 30 classical baseline indices ($M_1, M_2, R, \chi, H, GA, ABC, \ldots$). |
| `src/wuzi_analytical.py` | Closed-form values on $K_n, C_n, P_n, S_n, K_{p,q}, Q_k, W_n, F_p, r$-regular (cross-validated against the definition). |
| `src/mol_to_graph.py` | SMILES → hydrogen-suppressed `networkx.Graph` (RDKit-backed). |

## Datasets

| File | Source | Used for |
|---|---|---|
| `results/octane_descriptors.csv` | NIST WebBook + standard physical-chemistry references; 18 octane isomers, 5 physicochemical properties. | Section 6.1 (Table 3, Figure 1). |
| in-script: 106 non-isomorphic trees of order 10 | Generated via `networkx.generators.nonisomorphic_trees(10)`. | Section 6.2 (Figure 2 degeneracy). |
| in-script: 75 decane isomers | Built from canonical SMILES of all C${}_{10}$H${}_{22}$ alkane isomers. | Section 6.3 (Figure 3 structure-sensitivity). |

## Reproducibility scripts

Run from the repository root. All scripts complete in seconds on a single core.

| Script | Output | Manuscript reference |
|---|---|---|
| `scripts/05_octane_prediction.py` | `results/octane_competitor_comparison.{csv,md}`, `results/octane_descriptors.csv` | Table 3 (Pearson correlations), Figure 1 |
| `scripts/06_wuzi_degeneracy.py` | `results/wuzi_degeneracy.{csv,md}` | Figure 2 |
| `scripts/07_structure_sensitivity.py` | `results/structure_sensitivity.{csv,md}` | Figure 3 |
| `scripts/09_wuzi_bounds_tables.py` | `results/wuzi_bounds_ratio_tables.{csv,md}` | Remark 4.11 (ratio-bound constants) |
| `scripts/10_wuzi_extremal_search.py` | `results/wuzi_extremal_trees.{csv,md}` | Section 5 narrative (8 parameter triples at $n \le 12$) |
| `scripts/13_caterpillar_max_n5_to_20.py` | `results/wuzi_caterpillar_max_n5_to_20.csv` | Table 2 (extended n=5..20 enumeration at $(-1,-1,1)$) |
| `scripts/14_octane_bootstrap.py` | `results/octane_competitor_bootstrap.csv` | Table 3 caption (95% bootstrap CIs, B = 10,000) |

The two trailing scripts (13 and 14) are committed in this
release so a reviewer pulling the `v1.4-paper1-wuzi` tag can
reproduce the extended n=5..20 enumeration and the bootstrap CIs
end-to-end. Both scripts run in under one minute on a single
core.

## Figures referenced by the manuscript

| File | LaTeX label | Section |
|---|---|---|
| `figures/fig4_octane_heatmap.pdf` | `fig:octanes` | §6.1 |
| `figures/fig6_degeneracy_bars.pdf` | `fig:degeneracy` | §6.2 |
| `figures/fig7_structure_sensitivity.pdf` | `fig:ss` | §6.3 |
| `figures/fig10_caterpillar_maximizers.pdf` | `fig:caterpillar_maximizers` | §5 (computational extremal-graph search) |

## Test suite

The full pytest suite (`pytest -q` from the repo root) executes
124 tests covering the Wuzi index definition, the closed-form
identities, the BID-basis dimension count, the ratio-bound
edge-positivity hypothesis, and the per-script result-file
consistency. All 124 tests pass at this tag.

## Software environment

- Python 3.11+
- `networkx >= 3.0`
- `numpy`, `pandas`, `matplotlib`
- `rdkit` (only required to regenerate `results/octane_descriptors.csv` from SMILES; not required for the headline reproduction since the descriptor CSV is committed)
- `pytest` for the test suite

See `requirements.txt` for exact versions.

## Build the PDF

```
cd docs
pdflatex paper1_wuzi.tex
pdflatex paper1_wuzi.tex   # second pass for cross-references
```

The bibliography is embedded in `paper1_wuzi.tex` as a
`\begin{thebibliography}` block; no `bibtex` pass is required.

The `\graphicspath{{../figures/}}` directive points to the
top-level `figures/` directory. Compilation requires no extra
LaTeX packages beyond standard CTAN distributions (`amsmath`,
`amssymb`, `amsthm`, `mathtools`, `booktabs`, `array`, `graphicx`,
`hyperref`, `xcolor`, `enumitem`, `tikz`).
