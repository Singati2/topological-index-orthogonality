# Orthogonality Screening of Topological Indices for QSPR Modeling

An open-source pipeline for testing whether a proposed topological index
contributes **non-redundant information** beyond standard indices on real
chemistry datasets, prior to any QSPR validity claim.

This is **not** a "new superior topological index" project. It is a
**quality-control framework** intended to sit between formula proposal and
QSPR application — a pre-flight check authors and reviewers can run on any
candidate index in seconds.

A parametric **Wuzi index family** is used throughout as a worked case
study: we define it, derive its mathematical properties, run it through
the screening pipeline, and find that — at every parameter setting tested
— it is statistically redundant with classical degree-based indices. This
is a feature of the paper, not a defect: it illustrates *why* such
screening is needed.

---

## What this repository contains

```
src/
  mol_to_graph.py        SMILES → NetworkX (H-suppressed molecular graph)
  standard_indices.py    30 baseline indices (degree-, distance-, spectral-based, Balaban J)
  wuzi_index.py          Wuzi parametric family W(G; α, β, γ) + special-case identity tests
  novel_candidates.py    20 candidate indices from 5 alternative design families
                         (info-theoretic, centrality, spectral-hybrid, motif, eccentricity).
                         All are previously-published quantities under their own
                         literature; included to demonstrate the screening pipeline
                         on indices outside the BID subspace.
  load_data.py           ESOL / FreeSolv / Lipophilicity loaders (MoleculeNet)
  orthogonality.py       Correlation matrix, PCA / effective rank, VIF, partial correlation

scripts/
  01_baseline_correlations.py    30-index baseline orthogonality analysis
  02_wuzi_grid_search.py         100-point sweep of (α, β, γ) under the kill-test
  04_novel_candidates_test.py    Screening across the 20 alternative-family candidates

docs/
  theoretical_foundation.md      Edge-Degree-Pair Basis theorem
                                 (10-D ceiling on endpoint-degree edge-sum indices)
  [literature_notes.md, manuscript_draft.md, etc. — Phase 3]

results/                          Empirical outputs (CSVs + summary.txt) per dataset
LICENSE, requirements.txt, .gitignore
```

---

## The core claim

For hydrogen-suppressed molecular graphs with maximum degree `Δ ≤ 4` —
essentially all of organic chemistry — every *endpoint-degree edge-sum
index* of the form

```
I_f(G) = Σ_{uv ∈ E(G)} f(d_u, d_v)
```

(for symmetric `f`) lies in a real vector space of dimension **at most
10**, spanned by the ten edge-degree-pair counts `m_ij`. This is a
restatement of the well-known **Bond-Incident-Degree (BID)** framework
([REFERENCE NEEDED: Borovićanin et al. on BID bounds; Gutman & Tošović
2013]). The new emphasis here is the **explicit dimension bound** and the
**redundancy consequence**: any collection of ≥ 11 BID indices is
necessarily linearly dependent on this graph family, so proposing a new
BID index without a redundancy check is structurally limited.

The screening pipeline operationalizes this by testing whether a candidate
index has `|r| ≥ 0.95` (a "kill-test failure") with any of 30 standard
classical indices on a real chemistry benchmark, and also reports
effective rank, VIF, and partial-correlation diagnostics.

Full statement, proof, and discussion: [`docs/theoretical_foundation.md`](docs/theoretical_foundation.md).

---

## Quick start

```bash
git clone https://github.com/Singati2/topological-index-orthogonality.git
cd topological-index-orthogonality
pip install -r requirements.txt

# Phase-1 baseline + Wuzi screening on ESOL (default)
python scripts/01_baseline_correlations.py --dataset esol
python scripts/02_wuzi_grid_search.py     --dataset esol

# Replicate on FreeSolv and Lipophilicity (downloads on first run)
python scripts/01_baseline_correlations.py --dataset freesolv
python scripts/02_wuzi_grid_search.py     --dataset freesolv
python scripts/01_baseline_correlations.py --dataset lipophilicity
python scripts/02_wuzi_grid_search.py     --dataset lipophilicity
```

Outputs land in `results/<dataset>/`. ESOL takes ≈ 2 s; FreeSolv ≈ 1 s;
Lipophilicity ≈ 8 s (4 200 molecules).

### Datasets

| Name           | Property                                | n      | Source                          |
|----------------|-----------------------------------------|--------|---------------------------------|
| ESOL           | log aqueous solubility (mol/L)          | 1 128  | Delaney 2004 / MoleculeNet      |
| FreeSolv       | hydration free energy (kcal/mol)        | 642    | Mobley 2014 / MoleculeNet       |
| Lipophilicity  | octanol/water logD at pH 7.4            | 4 200  | ChEMBL extract / MoleculeNet    |

All are downloaded automatically on first call to `src/load_data.load(name)`.

---

## How to test your own candidate index

1. Implement a function `f(G: nx.Graph) -> float` taking a NetworkX graph
   (the H-suppressed molecular graph produced by
   `src.mol_to_graph.smiles_to_graph`).
2. Register it in the `ALL_INDICES` dict in `src/standard_indices.py`
   (or add to `CANDIDATE_INDICES` in `src/novel_candidates.py`).
3. Rerun `scripts/01_baseline_correlations.py` on each dataset.

Interpretation of outputs:

- `correlation_matrix.csv`: look at your index's row. If `max |r| ≥ 0.95`
  with **any** baseline, your index is statistically redundant with that
  baseline on this chemistry — the "kill-test failure" verdict.
- `partial_corr_target.csv`: rank of your index. A high rank means it
  adds unique signal *after* controlling for all baselines — the
  strongest case for a new index.
- `vif.csv`: VIF > 10 indicates heavy collinearity; VIF = ∞ indicates
  exact linear dependence (i.e., the index is a linear combination of
  the others within numerical precision).
- `pca_variance.csv`: how many PCA components your block needs for 95 %
  variance is the effective rank of the descriptor set.

If your candidate is an **endpoint-degree edge-sum index** (BID family),
the basis theorem implies it lies in a ≤ 10-dimensional subspace; the
screening will quantify how much of that subspace existing indices
already cover.

---

## What this framework does *not* prove

- Orthogonality is **necessary but not sufficient** for QSPR usefulness.
  An index can pass the screening and still add no predictive value.
- A single dataset is **not** evidence of general redundancy; we
  replicate across three independent QSPR benchmarks (ESOL, FreeSolv,
  Lipophilicity) precisely because cross-dataset consistency is the
  load-bearing claim.
- The pipeline targets **scalar topological indices**. Modern descriptor
  blocks (Mordred ≈ 1 800 descriptors, RDKit ≈ 200, fingerprints, graph
  neural networks) are out of scope here — the screening is intended for
  the math-chem / cheminformatics tradition of proposing named scalar
  invariants.
- Bond order is ignored (standard convention for classical topological
  indices). Weighting by bond order is an open extension.

---

## Status

| Phase | Scope                                                                  | Status         |
|-------|------------------------------------------------------------------------|----------------|
| 1     | BID basis theorem; multi-dataset replication; reframed README          | **in progress**|
| 2     | Wuzi family bounds, extremal-graph characterization, relations to known indices | pending  |
| 3     | Manuscript drafting, figures, journal targeting                        | pending        |

See [`docs/project_status.md`](docs/project_status.md) once Phase 3 begins.

---

## Citation

```
[CITATION PLACEHOLDER — manuscript in preparation]

Singati2 (GitHub), Natarajan, A., Arockiaraj, M., and collaborators.
"Orthogonality Screening of Topological Indices for QSPR Modeling:
An Open-Source Redundancy Benchmark with a Parametric Index Case Study."
In preparation, 2026.
```

## License

See [`LICENSE`](LICENSE).
