# Orthogonality Screening of Topological Indices for QSPR Modeling

A reproducible pipeline and a parametric case study for assessing
whether a proposed scalar topological index contributes
**non-redundant** information beyond the standard classical baseline
on real chemistry datasets.

The repository combines (i) a fast open-source screening pipeline
(correlation, PCA, partial correlation, VIF) replicated across three
QSPR benchmarks, with (ii) a worked parametric case study — the
**Wuzi index family** — used as the running example.

---

## What this repository is

- A reproducible orthogonality- / redundancy-screening pipeline that
  anyone proposing a new scalar topological index can run, in
  seconds, before claiming QSPR utility.
- A parametric case study — the Wuzi index family — used as a worked
  example of how a multi-parameter bond-incident-degree (BID)
  generalization can be screened against a 30-index classical
  baseline.
- An explicit *Edge-Degree-Pair Basis* observation giving a
  structural reason why endpoint-degree edge-sum indices on hydrogen-
  suppressed molecular graphs with maximum degree $\Delta\le 4$ lie
  in a vector space of dimension at most $10$.

## What this repository is *not*

- It is **not** a claim that the Wuzi family is a superior
  topological index. Across the three benchmarks tested, the family
  is statistically redundant with the classical 30-index baseline
  at every grid point examined except one borderline case on
  FreeSolv (which does not carry independent QSPR signal).
- It is **not** a claim that classical topological indices are
  useless. Many of them are widely used, well-studied, and
  appropriate for specific problems.
- It is **not** a claim that orthogonality alone proves QSPR
  usefulness. Orthogonality with the baseline is **necessary but
  not sufficient** — an index can clear the redundancy threshold
  and still add effectively zero predictive signal once partial
  correlations are inspected.
- The BID dimension-bound observation is not itself new to the
  mathematical-chemistry community; the contribution is the
  **explicit framing** plus the **redundancy-screening consequence**
  implemented end-to-end here.

---

## Two-paper strategy

This repository supports two related papers, developed in sequence.

**Paper 2 — prioritized, drafted first.**
*"The Wuzi Index Family: Graph-Theoretic Properties, Bounds,
Extremal Graphs, and Redundancy Analysis."*

A mathematical-chemistry paper in the tradition of Movahedi, Gutman,
Redžepović & Furtula, *MATCH Commun. Math. Comput. Chem.* 95(1),
2026, 141–162. Includes the definition, closed-form values on
standard graph families, bounds in graph parameters and in
classical indices, extremal graph characterizations, and a
numerical section reporting the orthogonality-screening results.
Target venues: *J. Math. Chem.*, *MATCH*, *AKCE Int. J. Graphs
Comb.*, *SAR & QSAR Env. Res.*

**Paper 1 — drafted after Paper 2.**
*"Orthogonality Screening of Topological Indices for QSPR Modeling:
A Structural and Empirical Redundancy Analysis."*

The broader methodology and benchmark paper. Cites Paper 2 as the
worked case study. Target venues: *Molecular Informatics*,
*J. Cheminformatics*.

See `docs/phase2_plan.md` and `docs/paper2_wuzi_manuscript_skeleton.md`
for current section drafts.

---

## Current status

| Component                                                    | Status |
|---|---|
| Wuzi definition + identity tests vs $M_2$, $R$, SCI, $H$     | ✅ done |
| Closed-form Wuzi values on 9 standard graph families         | ✅ done (numerically verified) |
| Edge-contribution function analysis                          | ✅ done |
| Edge-Degree-Pair Basis observation (dim $\le 10$ on $\Delta\le 4$) | ✅ done |
| 30-index baseline on ESOL / FreeSolv / Lipophilicity         | ✅ done |
| Wuzi 100-point grid sweep on each dataset                    | ✅ done |
| Cross-dataset redundancy comparison                          | ✅ done |
| Octane prediction (analog of MATCH §5.1)                     | ✅ done |
| Degeneracy on 106 trees of order 10                          | ✅ done |
| Structure sensitivity on 75 decane isomers                   | ✅ done |
| Reproducible end-to-end pipeline                             | ✅ done |
| Publication-quality figures (PNG + PDF, 8 figures)           | ✅ done |
| Paper 2 §3 bounds in $(n, m, \Delta, \delta)$                | ⏳ math pending |
| Paper 2 §4 bounds via classical indices                      | ⏳ math pending |
| Paper 2 §5 extremal-graph characterization                   | ⏳ math pending |
| Paper 2 §1 introduction + §7 discussion                      | ⏳ writing |
| Paper 1 manuscript                                           | ⏳ scheduled after Paper 2 |

---

## Mathematical formulation

### Wuzi index family

For a simple connected graph $G = (V, E)$ with vertex degrees
$d_u, d_v$ on each edge $uv \in E(G)$, the **Wuzi index family**
is the three-parameter scalar functional

$$W(G; \alpha, \beta, \gamma) = \sum_{uv \in E(G)} (d_u d_v)^{\alpha} (d_u + d_v)^{\beta} \exp\!\left( \gamma \cdot \frac{|d_u - d_v|}{d_u + d_v} \right),$$

with $\alpha, \beta, \gamma \in \mathbb{R}$. Several classical
bond-incident-degree indices are recovered as special cases:

| Parameters                                  | Recovered index                |
|---|---|
| $\alpha = 1$, $\beta = 0$, $\gamma = 0$     | Second Zagreb $M_2$            |
| $\alpha = -1/2$, $\beta = 0$, $\gamma = 0$  | Randić $R$                     |
| $\alpha = 0$, $\beta = -1/2$, $\gamma = 0$  | Sum-connectivity $\mathrm{SCI}$|
| $\alpha = 0$, $\beta = -1$, $\gamma = 0$    | (Half) Harmonic $H/2$          |

When $\gamma = 0$ the family reduces to the standard two-parameter
$(\alpha, \beta)$ degree-based family; the $\gamma$ term introduces
explicit dependence on the normalized degree imbalance
$|d_u - d_v| / (d_u + d_v) \in [0, 1)$ along each edge.

### Edge-Degree-Pair Basis observation

For hydrogen-suppressed molecular graphs with maximum degree
$\Delta \le 4$, the only unordered degree pairs on an edge are
$(i, j)$ with $1 \le i \le j \le 4$, giving $10$ admissible
edge-degree pairs. Let

$$m_{ij}(G) = \bigl| \{\, uv \in E(G) \,:\, \{d_u, d_v\} = \{i, j\} \,\} \bigr|, \qquad 1 \le i \le j \le 4.$$

Then every *endpoint-degree edge-sum* (BID) index
$I_f(G) = \sum_{uv \in E(G)} f(d_u, d_v)$ with $f$ symmetric is a
linear functional of these $10$ counts:

$$I_f(G) = \sum_{1 \le i \le j \le 4} f(i, j) \cdot m_{ij}(G).$$

Hence the vector space spanned by all BID indices has dimension at
most $10$ on this graph class; any $11$ such indices are linearly
dependent. See `docs/theoretical_foundation.md` for the full
statement and proof. The observation that BID indices are linear in
$m_{ij}$ counts is standard in the mathematical-chemistry
literature; the contribution here is the explicit dimension bound
and its application to redundancy screening.

---

## Key results

Baseline pairwise redundancy on the 30-index set
($\binom{30}{2} = 435$ pairs):

| Dataset       | $n$    | Pairs $\lvert r\rvert\ge 0.90$ | Pairs $\lvert r\rvert\ge 0.95$ | PC$_{95}$ |
|---|---:|---:|---:|---:|
| ESOL          | 1127   | 203 / 435 (46.7%)              | 159 / 435 (36.6%)              | 3         |
| FreeSolv      | 639    | 257 / 435 (59.1%)              | 162 / 435 (37.2%)              | 3         |
| Lipophilicity | 4200   | 200 / 435 (46.0%)              | 164 / 435 (37.7%)              | 3         |

PC$_{95}$ denotes the number of principal components needed to
capture 95% of the variance of the 30 standardized indices.
Three components suffice on all three benchmarks; five capture 99%.

Wuzi-family 100-point grid sweep against the same 30-index baseline:

| Dataset       | Pass (max $\lvert r\rvert < 0.95$) | Fail | Min max $\lvert r\rvert$ | Max $\lvert$partial $r$ w/ target$\rvert$ |
|---|---:|---:|---:|---:|
| ESOL          | 0    | 100  | 0.965 | 0.0282 |
| FreeSolv      | 1    | 99   | 0.943 | 0.0362 |
| Lipophilicity | 0    | 100  | 0.965 | 0.0093 |

The single borderline grid point on FreeSolv occurs at
$(\alpha, \beta, \gamma) = (-1, -1, 1)$ with $\max \lvert r\rvert = 0.943$;
its partial correlation with $\Delta G_{\mathrm{hyd}}$ after
controlling for the 30 baseline indices is $\approx 0$, so it does
not carry useful independent QSPR signal. The empirical pattern is
consistent with the dimension-bound observation above: high
redundancy with the classical baseline is the expected behavior of
BID-family indices on this graph class.

Full per-dataset tables are in `results/<dataset>/` and the
cross-dataset summary in `results/cross_dataset_summary.md`.

---

## How to reproduce results

```bash
git clone https://github.com/Singati2/topological-index-orthogonality.git
cd topological-index-orthogonality
pip install -r requirements.txt
```

### Baseline + Wuzi grid on each dataset

```bash
python scripts/01_baseline_correlations.py --dataset esol            # ~ 1 s
python scripts/01_baseline_correlations.py --dataset freesolv        # ~ 1 s
python scripts/01_baseline_correlations.py --dataset lipophilicity   # ~ 10 s

python scripts/02_wuzi_grid_search.py --dataset esol                 # ~ 2 s
python scripts/02_wuzi_grid_search.py --dataset freesolv             # ~ 1 s
python scripts/02_wuzi_grid_search.py --dataset lipophilicity        # ~ 15 s
```

Outputs land in `results/<dataset>/`.

### Cross-dataset comparison and numerical sections

```bash
python scripts/03_cross_dataset_summary.py     # results/cross_dataset_summary.{csv,md}
python scripts/05_octane_prediction.py         # results/octane_prediction.{csv,md}
python scripts/06_wuzi_degeneracy.py           # results/wuzi_degeneracy.{csv,md}
python scripts/07_structure_sensitivity.py     # results/structure_sensitivity.{csv,md}
python scripts/08_make_figures.py              # figures/*.{png,pdf}
```

### Adding a new candidate index

1. Implement a function `f(G: nx.Graph) -> float` on the
   H-suppressed molecular graph produced by
   `src.mol_to_graph.smiles_to_graph`.
2. Register it in `src/standard_indices.py` (`ALL_INDICES`) or in
   `src/novel_candidates.py` (`CANDIDATE_INDICES`).
3. Rerun `scripts/01_baseline_correlations.py` per dataset and
   inspect `results/<dataset>/correlation_matrix.csv`,
   `partial_corr_<target>.csv`, and `vif.csv`.

If the candidate is a BID index, the dimension-bound observation
implies it lies in a vector space of dimension at most $10$; the
screening pipeline quantifies how much of that subspace the
30-index baseline already covers.

---

## Datasets

| Name          | Property                                  | $n$    | Source                       |
|---|---|---:|---|
| ESOL          | log aqueous solubility (mol / L)          | 1 128  | Delaney 2004 / MoleculeNet   |
| FreeSolv      | hydration free energy (kcal / mol)        | 642    | Mobley 2014 / MoleculeNet    |
| Lipophilicity | octanol-water $\log D$ at pH 7.4          | 4 200  | ChEMBL / MoleculeNet         |

All datasets download automatically via `src.load_data.load(name)`
on first use; cached CSVs are written to `data/` (gitignored).

---

## Implemented baseline indices (30)

**Degree-based (18):** First and Second Zagreb ($M_1$, $M_2$),
Modified Zagreb ($mM_1$, $mM_2$), Forgotten ($F$), Randić ($R$),
Sum-connectivity (SCI), Harmonic ($H$), Geometric-arithmetic ($GA$),
Arithmetic-geometric ($AG$), Atom-bond connectivity (ABC),
Atom-bond-sum connectivity (ABS), Augmented Zagreb (AZI),
Sombor ($SO$), Reduced Sombor ($SO_{\mathrm{red}}$),
Albertson irregularity (Alb), Sigma ($\sigma$),
Reduced First Zagreb ($\mathrm{redM}_1$).

**Distance-based (8):** Wiener ($W$), Hyper-Wiener (WW),
Harary (HR), Schultz (MTI), Gutman, Mostar, Szeged (Sz),
Padmakar–Ivan (PI).

**Spectral (3):** Estrada (EE), Graph Energy, Spectral Radius.

**Other (1):** Balaban $J$.

---

## Repository layout

```
src/
  mol_to_graph.py          SMILES -> NetworkX (H-suppressed molecular graph)
  standard_indices.py      30 baseline classical indices
  wuzi_index.py            Wuzi parametric family with special-case identity tests
  wuzi_analytical.py       Closed-form Wuzi values on standard graph families
  novel_candidates.py      20 alternative-family candidate indices (info-theoretic, centrality, ...)
  load_data.py             ESOL / FreeSolv / Lipophilicity loaders
  orthogonality.py         Correlation matrix, PCA, VIF, partial correlation

scripts/
  01_baseline_correlations.py    Per-dataset baseline orthogonality
  02_wuzi_grid_search.py         Per-dataset 100-point Wuzi sweep
  03_cross_dataset_summary.py    Cross-dataset redundancy table
  04_novel_candidates_test.py    Screening of 20 alternative-family candidates on ESOL
  05_octane_prediction.py        MATCH §5.1 analog on 18 octanes
  06_wuzi_degeneracy.py          MATCH §5.3 analog on 106 trees of order 10
  07_structure_sensitivity.py    MATCH §5.4 analog on 75 decane isomers
  08_make_figures.py             Generates the 8 manuscript figures

docs/
  theoretical_foundation.md      Edge-Degree-Pair Basis observation
  edge_contribution_analysis.md  Edge-contribution function analysis
  phase2_plan.md                 Section-by-section paper plan
  paper2_wuzi_manuscript_skeleton.md   Paper 2 skeleton with math gaps marked
  literature_notes.md            Citation discipline document
  figure_audit.md                Figure-by-figure consistency record
  figure_captions.md             Manuscript-ready figure captions
  advisor_update.md              Update document for Arockiaraj sir
  collaboration.md               Contributor workflow

figures/                          Publication-quality figures (PNG + PDF)
  fig2_pca_scree                  Cross-dataset PCA scree (3 components -> 95%)
  fig3_redundancy_bars            Per-dataset redundant-pair counts at three thresholds
  fig4_octane_heatmap             Octane prediction: signed correlation with 8 properties
  fig5_wuzi_param_heatmaps        Wuzi redundancy screen on ESOL
  fig5b_wuzi_param_heatmaps_all   Wuzi redundancy screen across all three datasets
  fig6_degeneracy_bars            Index degeneracy on 106 trees of order 10
  fig7_structure_sensitivity      Structure sensitivity on 75 decane isomers
  fig8_correlation_heatmap        Full 30-index pairwise correlation matrix (ESOL)

results/                          Reproducible outputs by dataset
data/                             Cached download CSVs (gitignored)
```

---

## What remains unfinished

1. Paper 2 §3, §4, §5 — bounds in graph parameters, bounds in
   classical indices, and extremal-graph characterizations.
   Graph-theoretic derivations following the MATCH 95:141–162
   template.
2. Paper 2 manuscript prose — Introduction, Discussion, Conclusion.
3. Resolution of every `[REFERENCE NEEDED]` placeholder in
   `docs/literature_notes.md` by reading the source.
4. Paper 1 manuscript, scheduled after Paper 2 lands.

---

## Citation

```
[CITATION PLACEHOLDER — manuscript in preparation]

A. Natarajan, G. Shiwakoti, M. Arockiaraj.
"The Wuzi Index Family: Graph-Theoretic Properties, Bounds,
 Extremal Graphs, and Redundancy Analysis."
In preparation, 2026.
```

## License

MIT License. See `LICENSE`. Code is openly available; please cite
the manuscript above (once published) for academic use.
