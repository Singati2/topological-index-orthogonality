# Exact redundancy of degree-based topological indices

[![tests](https://github.com/Singati2/topological-index-orthogonality/actions/workflows/test.yml/badge.svg)](https://github.com/Singati2/topological-index-orthogonality/actions/workflows/test.yml)

Code, data and results for the manuscript

> G. Shiwakoti, A. Natarajan, P. Chalise, M. Arockiaraj,
> *Degree-Based Topological Indices Are Exactly Redundant on Molecular Data:
> A Span Certificate and a Target-Aware Screen for New Descriptors* (v2, in preparation).

Source: [`docs/paper2_orthogonality_screening.tex`](docs/paper2_orthogonality_screening.tex)
(tables are generated into `docs/tables_v2/` by `scripts/75_make_v2_tables.py`).
The companion mathematical paper on the Loyola index family lives in
[Singati2/LOYOLA-PAPER](https://github.com/Singati2/LOYOLA-PAPER).

## What the paper shows

1. **Span certificate.** Every bond-incident-degree (BID) index
   `I_f = Σ_{uv∈E} f(d_u, d_v)` is a linear combination of the edge-degree-pair
   counts `m_ij` (classical; formalised by Rada, MATCH 96 (2026) 841–863).
   We turn this into a rank check: if the BID block of a baseline has the same
   rank as the realised `m_ij` matrix, every BID index — any `f`, any parameters,
   any future variant — has zero residual on the baseline and partial correlation
   0 with any target. The 18 degree-based indices of the standard 30-index
   baseline satisfy it for **all** graphs of maximum degree ≤ 4.
2. **Census.** 14 published degree-based indices (15 instances, 8 Sombor-type)
   and 400 points of the Loyola family are exactly in span on ESOL, FreeSolv,
   Lipophilicity and BBBP (relative residual ≤ 2×10⁻¹²); non-degree-based Sombor
   controls are not.
3. **Numerical pitfall.** Partial correlations of exactly redundant indices
   computed with absolute tolerances are floating-point noise (an earlier draft
   reported 0.028/0.036/0.009; the true value is 0).
4. **Target-aware screen for other indices.** Of 20 non-BID candidates, 15–17 per
   dataset pass a pairwise |r| < 0.95 screen, but only two candidate–dataset
   pairs keep residual target association under both linear (with intervals) and
   cross-fitted non-linear residualization.
5. **Information ceiling.** Many molecules share their hydrogen-suppressed graph
   (FreeSolv: 639 molecules, 236 graphs), which bounds every topology-only
   descriptor; a chemistry-aware RDKit-2D reference outperforms all topological
   models in a graph-grouped, corrected benchmark.

## Repository layout

| Path | Contents |
|---|---|
| `data/` | The four MoleculeNet files (verified byte-identical to the official mirror), `SHA256SUMS`, curated copies with graph hashes; see `data/README.md` |
| `src/` | Graph construction, the 30 baseline indices, 20 candidate indices, the Loyola index, redundancy diagnostics incl. `span_certified_pcor` |
| `scripts/` | Analyses (below) |
| `results/` | All outputs used by the manuscript |
| `docs/` | Manuscript and generated tables |
| `tests/` | Unit and consistency tests (`pytest`) |

## Reproduce

```bash
pip install -r requirements.txt
python scripts/76_verify_datasets.py --check-remote   # data integrity + curated sets   (~2 min)
python scripts/70_bid_span_certificate.py             # certificate, census, identities (~2 min)
for d in esol freesolv lipophilicity; do python scripts/02_wuzi_grid_search.py --dataset $d; done
python scripts/16_novel_candidates_multidataset.py    # candidate screen, 3 regression sets
python scripts/17_bbbp_descriptors_and_altfamily.py   # candidate screen, BBBP
python scripts/73_nonlinear_conditional.py            # cross-fitted non-linear check  (~50 min)
python scripts/74_dedup_sensitivity.py                # de-duplication sensitivity
python scripts/72_benchmark_v2.py                     # curation, collisions, benchmark (~80 min)
python scripts/77_isomorphism_and_r2_bound.py         # exact isomorphism of graph groups, R^2 bound
python scripts/75_make_v2_tables.py                   # LaTeX tables from results/
pytest -q
cd docs && pdflatex paper2_orthogonality_screening.tex   # or tectonic
```

Results were produced with Python 3.11, RDKit 2022.09, scikit-learn 1.8 and
NumPy 2.4. Under that RDKit/NumPy combination 28 RDKit EState/VSA-EState/
Balaban/Bertz/Ipc descriptors cannot be computed and are excluded from the
RDKit-2D reference (documented in `scripts/72_benchmark_v2.py`).

## History

This repository contains only Paper 2 (v2). Version 1 of the paper and the Paper 1
material that previously lived here have been removed; Paper 1 is maintained in
[Singati2/LOYOLA-PAPER](https://github.com/Singati2/LOYOLA-PAPER).

## License

MIT (code). Datasets are redistributed from MoleculeNet; see `data/README.md`
for their original sources.
