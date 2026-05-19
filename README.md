# Orthogonality Screening of Topological Indices for QSPR Modeling

An open-source pipeline and a parametric case-study aimed at testing
whether a proposed scalar topological index contributes
**non-redundant information** beyond the standard classical 30-index
baseline on real chemistry datasets.

---

## What this repository is *not*

- This is **not** a "Wuzi index is superior" paper.
- We do **not** claim that classical topological indices are useless.
- We do **not** claim that information-theoretic, motif, centrality,
  or other "alternative-family" indices are newly invented — they
  appear in this repository only as candidates for screening, with
  full attribution to the existing literature.
- A single dataset is **not** evidence of general redundancy; we
  replicate across three independent QSPR benchmarks.

## What this repository *is*

- A reproducible, fast (~ seconds per benchmark) open-source
  **orthogonality-screening pipeline** that anyone proposing a new
  scalar topological index can run before claiming QSPR utility.
- A **parametric case study** — the **Wuzi index family** — used as a
  worked example demonstrating how a multi-parameter degree-based
  generalization can be screened.
- An explicit **Edge-Degree-Pair Basis theorem** giving a structural
  reason why endpoint-degree edge-sum indices on Δ ≤ 4 molecular
  graphs lie in a ≤ 10-dimensional subspace.

---

## Mathematical formulation

### Wuzi index family

For a simple connected graph $G = (V, E)$ with vertex degrees
$d_u, d_v$ on each edge $uv \in E(G)$, the **Wuzi index family** is
the three-parameter scalar functional

$$
W(G;\ \alpha, \beta, \gamma) \;=\; \sum_{uv \in E(G)}
  (d_u \, d_v)^{\alpha}
  \;(d_u + d_v)^{\beta}
  \;\exp\!\left(\gamma \cdot \frac{|d_u - d_v|}{d_u + d_v}\right),
$$

with $\alpha, \beta, \gamma \in \mathbb{R}$. The family is a strict
generalization of several classical bond-incident-degree (BID)
indices recovered as special cases:

| Parameters | Recovered index |
|---|---|
| $\alpha = 1,\ \beta = 0,\ \gamma = 0$        | Second Zagreb $M_2$ |
| $\alpha = -\tfrac{1}{2},\ \beta = 0,\ \gamma = 0$ | Randić $R$ |
| $\alpha = 0,\ \beta = -\tfrac{1}{2},\ \gamma = 0$ | Sum-connectivity $\mathrm{SCI}$ |
| $\alpha = 0,\ \beta = -1,\ \gamma = 0$        | (Half) Harmonic $H/2$ |

When $\gamma = 0$ the family reduces to the standard two-parameter
$(\alpha, \beta)$ degree-based family; the $\gamma$ term introduces
explicit dependence on the *normalized degree imbalance*
$|d_u - d_v| / (d_u + d_v) \in [0, 1)$ along each edge.

### Edge-Degree-Pair Basis theorem (informal)

Let $\mathcal{F}_{\mathrm{BID}}$ denote the class of
**bond-incident-degree** scalar indices, i.e. real-valued functionals
on graphs of the form

$$
I(G) \;=\; \sum_{uv \in E(G)} f(d_u, d_v)
\qquad
\text{for some symmetric } f : \mathbb{N}^2 \to \mathbb{R}.
$$

On the subclass of **H-suppressed molecular graphs** with maximum
degree $\Delta \le 4$, the only possible unordered degree pairs on an
edge are $(i, j)$ with $1 \le i \le j \le 4$, giving 10 admissible
edge-degree pairs. Hence every $I \in \mathcal{F}_{\mathrm{BID}}$ is
a linear functional of the 10 edge-degree-pair counts
$m_{ij}(G) = \big|\{uv \in E(G) : \{d_u, d_v\} = \{i, j\}\}\big|$,
i.e.

$$
I(G) \;=\; \sum_{1 \le i \le j \le 4}
  f(i, j) \cdot m_{ij}(G).
$$

**Consequence.** On Δ ≤ 4 graphs the vector space spanned by all BID
indices has dimension at most 10; any 11 BID indices are necessarily
linearly dependent. This is the structural reason why high empirical
redundancy is *unsurprising*, and motivates the orthogonality
screening protocol.

See `docs/theoretical_foundation.md` for the full statement and proof.

---

## Two-paper strategy

This repository supports two related papers being developed in
sequence.

**Paper 2 (prioritized — written first):**
*"The Wuzi Index Family: Graph-Theoretic Properties, Bounds, Extremal
Graphs, and Redundancy Analysis"*

Follows the mathematical-chemistry tradition (analog of Movahedi,
Gutman, Redžepović & Furtula, *MATCH* 95(1), 2026, 141–162). Includes:

- Definition + closed-form values on K_n / C_n / P_n / S_n / K_{p,q} /
  Q_k / W_n / F_p / regular graphs (✅ done)
- Edge-contribution function analysis (✅ done)
- Bounds in terms of (n, m, Δ, δ) (⏳ pending — math)
- Bounds in terms of classical indices (M_1, M_2, R, SO, GA, H, …)
  (⏳ pending — math)
- Extremal-graph characterization for trees, unicyclic, bicyclic
  (⏳ pending — math)
- Numerical Section 6: octane prediction, intercorrelations across
  three QSPR datasets, PCA effective rank, degeneracy, structure
  sensitivity (✅ done)

Target venues (realistic): J. Math. Chem., MATCH, AKCE Int. J. Graphs
Comb., SAR & QSAR Env. Res. See `docs/phase2_plan.md` and
`docs/paper2_wuzi_manuscript_skeleton.md` for details.

**Paper 1 (drafted later):**
*"Orthogonality Screening of Topological Indices for QSPR Modeling:
A Structural and Empirical Redundancy Analysis"*

The broader methodology / pipeline / cross-dataset benchmark paper.
Cites Paper 2 as the worked example. More ambitious target
(J. Cheminformatics or J. Chem. Inf. Model.); depends on Paper 2
landing first.

---

## Current status

| Component                                       | Status      |
|---|---|
| Wuzi definition + identity tests vs M_2, R, SCI, H | ✅ done    |
| Closed-form Wuzi values on 9 standard graph classes | ✅ done (numerically verified) |
| Edge-contribution function ψ analysis              | ✅ done   |
| Edge-Degree-Pair Basis theorem                     | ✅ done   |
| 30-index baseline on ESOL / FreeSolv / Lipophilicity | ✅ done |
| Wuzi 100-point grid sweep on each dataset          | ✅ done   |
| Cross-dataset redundancy comparison                | ✅ done   |
| Octane prediction (MATCH §5.1 analog)              | ✅ done   |
| Degeneracy on trees of order 10                    | ✅ done   |
| Structure sensitivity on 75 decanes                | ✅ done   |
| Open-source pipeline + reproducibility             | ✅ done   |
| Bibliographic discipline document                  | ✅ done   |
| Paper 2 § 3 bounds in (n, m, Δ, δ)                 | ⏳ awaiting math |
| Paper 2 § 4 bounds via classical indices           | ⏳ awaiting math |
| Paper 2 § 5 extremal graphs                        | ⏳ awaiting math |
| Paper 2 § 1 introduction + § 7 discussion          | ⏳ writing |
| Paper 1 manuscript                                  | ⏳ scheduled after Paper 2 |
| Figures (production quality, PNG + PDF, 8 figures)  | ✅ done    |

---

## What remains unfinished

1. The three math sections of Paper 2 (bounds in graph parameters,
   bounds in classical indices, extremal graphs). These require
   graph-theoretic derivations following the MATCH 95:141-162 template.
2. Manuscript writing (intro, discussion, conclusion) for Paper 2.
3. Full reference bibliography (every `[REFERENCE NEEDED]` placeholder
   in `docs/literature_notes.md` must be resolved by reading the source).
4. Paper 1 manuscript (deliberately deferred until Paper 2 lands).

---

## How to reproduce results

```bash
git clone https://github.com/Singati2/topological-index-orthogonality.git
cd topological-index-orthogonality
pip install -r requirements.txt
```

### Baseline + Wuzi grid on each dataset

```bash
python scripts/01_baseline_correlations.py --dataset esol            # ≈ 1 s
python scripts/01_baseline_correlations.py --dataset freesolv        # ≈ 1 s
python scripts/01_baseline_correlations.py --dataset lipophilicity   # ≈ 10 s

python scripts/02_wuzi_grid_search.py --dataset esol            # ≈ 2 s
python scripts/02_wuzi_grid_search.py --dataset freesolv        # ≈ 1 s
python scripts/02_wuzi_grid_search.py --dataset lipophilicity   # ≈ 15 s
```

Outputs land in `results/<dataset>/`.

### Cross-dataset comparison

```bash
python scripts/03_cross_dataset_summary.py
```
Outputs `results/cross_dataset_summary.{csv,md}`.

### Numerical sections of Paper 2

```bash
python scripts/05_octane_prediction.py       # §6.1 analog (MATCH §5.1)
python scripts/06_wuzi_degeneracy.py         # §6.4 analog (MATCH §5.3)
python scripts/07_structure_sensitivity.py   # §6.5 analog (MATCH §5.4)
```

### Optional: 20-candidate "alternative families" experiment

```bash
python scripts/04_novel_candidates_test.py   # output in results/novel_candidates_experiment/
```

This was the experiment that confirmed alternative-family candidates
(Shannon degree entropy, mean clustering coefficient, etc.) also have
prior chemistry literature — kept in the repository as supporting
empirical evidence; the candidates themselves are not proposed as new
indices.

### Verifying closed-form Wuzi values

```bash
python -m src.wuzi_analytical    # checks closed-form vs direct compute
python -m src.standard_indices   # checks 30-index implementations on P4, C6
python -m src.mol_to_graph       # checks SMILES → NetworkX
```

---

## How to plug in a new candidate index

1. Implement a function `f(G: nx.Graph) -> float` taking the
   H-suppressed molecular graph produced by `src.mol_to_graph.smiles_to_graph`.
2. Register it in the `ALL_INDICES` dict of `src/standard_indices.py`
   (or add to `CANDIDATE_INDICES` in `src/novel_candidates.py`).
3. Rerun `scripts/01_baseline_correlations.py --dataset <name>` on
   each dataset.
4. Inspect:
   - `results/<dataset>/correlation_matrix.csv` — max |r| with any baseline ≥ 0.95 ⇒ statistically redundant.
   - `results/<dataset>/partial_corr_<target>.csv` — rank of new index ⇒ unique signal after controlling for the baseline.
   - `results/<dataset>/vif.csv` — VIF > 10 indicates heavy
     collinearity; VIF = ∞ indicates exact linear dependence.

If your candidate is an **endpoint-degree edge-sum** (BID) index, the
Edge-Degree-Pair Basis theorem implies it lies in a ≤ 10-dimensional
subspace. The screening quantifies how much of that subspace existing
indices already cover.

---

## Datasets

| Name           | Property                                  | n      | Source                          |
|----------------|-------------------------------------------|--------|---------------------------------|
| ESOL           | log aqueous solubility (mol / L)          | 1 128  | Delaney 2004 / MoleculeNet      |
| FreeSolv       | hydration free energy (kcal / mol)        | 642    | Mobley 2014 / MoleculeNet       |
| Lipophilicity  | octanol / water log D at pH 7.4           | 4 200  | ChEMBL extract / MoleculeNet    |

All download automatically via `src/load_data.load(name)` on first use.

---

## Implemented baseline indices (30)

**Degree-based (18):** First & Second Zagreb (M_1, M_2), Modified Zagreb
(mM_1, mM_2), Forgotten (F), Randić (R), Sum-connectivity (SCI),
Harmonic (H), Geometric-arithmetic (GA), Arithmetic-geometric (AG),
Atom-bond-connectivity (ABC), Atom-bond-sum-connectivity (ABS),
Augmented Zagreb (AZI), Sombor (SO), Reduced Sombor (SO_red),
Albertson irregularity (Alb), Sigma (σ), Reduced First Zagreb (redM_1)

**Distance-based (8):** Wiener (W), Hyper-Wiener (WW), Harary (HR),
Schultz (MTI), Gutman, Mostar, Szeged (Sz), Padmakar–Ivan (PI)

**Spectral (3):** Estrada (EE), Graph Energy, Spectral Radius

**Other (1):** Balaban J

---

## Repository layout

```
src/
  mol_to_graph.py          SMILES → NetworkX (H-suppressed molecular graph)
  standard_indices.py      30 baseline classical indices
  wuzi_index.py            Wuzi parametric family with special-case identity tests
  wuzi_analytical.py       Closed-form Wuzi values on standard graph families
  novel_candidates.py      20 "alternative-family" candidates (info-theoretic, centrality, etc.)
  load_data.py             ESOL / FreeSolv / Lipophilicity loaders
  orthogonality.py         Correlation matrix, PCA, VIF, partial correlation

scripts/
  01_baseline_correlations.py    Per-dataset baseline orthogonality
  02_wuzi_grid_search.py         Per-dataset 100-point Wuzi sweep
  03_cross_dataset_summary.py    Aggregate redundancy stats across datasets
  04_novel_candidates_test.py    20-candidate kill-test experiment
  05_octane_prediction.py        MATCH §5.1 analog on 18 octanes
  06_wuzi_degeneracy.py          MATCH §5.3 analog on 106 trees of order 10
  07_structure_sensitivity.py    MATCH §5.4 analog on 75 decanes

docs/
  theoretical_foundation.md      Edge-Degree-Pair Basis theorem
  edge_contribution_analysis.md  ψ(x, y; α, β, γ) analysis
  phase2_plan.md                 Section-by-section paper plan
  paper2_wuzi_manuscript_skeleton.md   Full Paper-2 skeleton with math gaps marked
  literature_notes.md            Bibliographic discipline document
  advisor_update.md              Two-page update for Arockiaraj sir
  email_to_adviK_update.md       Email draft for Advik
  collaboration.md               Contributor workflow

results/                          Reproducible outputs by dataset
figures/                          Publication-quality figures (PNG + PDF)
  fig2_pca_scree                  Cross-dataset PCA scree (3 components ⇒ 95%)
  fig3_redundancy_bars            Per-dataset redundant-pair counts at three |r| thresholds
  fig4_octane_heatmap             Octane prediction: signed correlation with 8 properties
  fig5_wuzi_param_heatmaps        Wuzi kill-test on ESOL (γ ∈ {0, 0.5, 1, 2})
  fig5b_wuzi_param_heatmaps_all   Same kill-test across all three datasets
  fig6_degeneracy_bars            Index degeneracy on 106 trees of order 10
  fig7_structure_sensitivity      Structure sensitivity on 75 decane isomers
  fig8_correlation_heatmap        30-index pairwise correlation heatmap (ESOL)
data/                             Cached download CSVs (gitignored)
```

---

## Citation

```
[CITATION PLACEHOLDER — manuscript in preparation]

A. Natarajan, G. Shiwakoti, M. Arockiaraj.
"The Wuzi Index Family: Graph-Theoretic Properties, Bounds, Extremal
Graphs, and Redundancy Analysis."
In preparation, 2026.
```

## License

See `LICENSE`. Code is openly available; please cite the manuscript
above (once published) for academic use.
