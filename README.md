# Topological Index Kill-Test Pilot

Pre-flight infrastructure for evaluating whether a *new* topological index
carries information not already captured by the existing ~30 standard indices,
before committing to either a graph-theoretic bounds paper or a full QSPR
benchmark.

**Status:** Baseline established 2026-05-14 on ESOL (1127 drug-like molecules).
Awaiting a new candidate index from collaborator (Advik Natarajan, Loyola
College Chennai) for kill-test.

---

## Why this exists

Proposing a new topological index by tweaking an existing formula is easy.
Showing it carries *novel* information vs. the 3000+ existing indices is the
hard part. This pipeline runs that test in ~2 seconds on real drug-like
chemistry. If the candidate has |r| ≥ 0.95 with any classical baseline, it is
statistically redundant on this set and any downstream paper (bounds *or*
applied) rests on a dead claim.

## Headline baseline finding (n=1127 ESOL molecules)

The 30 standard indices implemented here carry the **effective rank of ~3
dimensions** on real chemistry:

- 95% of total variance captured by **3 PCA components**
- 99% captured by **5 components**
- **159 of 435 pairs** of indices (37%) have |r| ≥ 0.95
- Multiple indices have VIF = ∞ (perfect linear dependencies)
- Partial correlations with target property (logS) max out at |r| ≈ 0.12

In other words: the 30-index block is doing the work of ~3-5. Most "new"
indices are downstream of that ceiling. This is also a publishable result on
its own — see "Pivot option" below.

## Quick start

```bash
pip install -r requirements.txt
python scripts/01_baseline_correlations.py
```

Outputs land in `results/`. Runs in 2 seconds.

## Adding a new candidate index

1. Implement a function `f(G: nx.Graph) -> float` for your candidate.
2. Register it: add to `ALL_INDICES` dict in `src/standard_indices.py`
   (or create `src/new_index.py` and merge in the script).
3. Rerun `scripts/01_baseline_correlations.py`.
4. Inspect:
   - `results/correlation_matrix.csv` — look at the new index's row.
     Max |r| with any baseline ≥ 0.95 → **redundant, redesign**.
   - `results/partial_corr_logS.csv` — find your index. High rank means it
     adds unique signal *after* controlling for all baselines. This is the
     strongest possible case for a new index.

## Three "good index" definitions (pick one before writing)

A new index can be "better" in three non-equivalent senses. The paper should
commit to one explicitly:

1. **Orthogonal information** — adds predictive signal in a regression
   already containing the standard descriptors. Measured by partial
   correlation; pipeline reports this directly.
2. **Higher standalone correlation** with target property, consistent across
   multiple datasets (not one cherry-picked n=20 set).
3. **Interpretable** — has a clean chemical meaning (branching, electronic
   structure, surface area) that existing indices muddle.

Definition (1) is the strongest, hardest, and most defensible at JCIM / J.
Cheminformatics tier.

## Pivot option (if the new-index path looks unpromising)

The baseline results above (3 effective dimensions out of 30 indices,
n=1127, real chemistry) already constitute a publishable empirical critique
of the index-proliferation literature. Reframing the paper as:

> *"Empirical orthogonality analysis of 30 commonly-used topological indices
> on real drug-like chemistry shows the family carries the predictive
> information of ~3 dimensions and contains 159 statistically redundant
> pairs (|r| ≥ 0.95). We argue that new indices should clear an orthogonality
> bar before publication."*

...lands cleanly in J. Chem. Inf. Model. or J. Cheminformatics, uses
Ganesh's biostat methodology training, and is more useful to the field than
yet another tweaked-formula paper.

## Files

```
src/mol_to_graph.py        SMILES -> NetworkX (H-suppressed molecular graph)
src/standard_indices.py    30 baseline indices (degree, distance, spectral, other)
src/load_data.py           ESOL dataset loader (Delaney, 1128 molecules)
src/orthogonality.py       Correlation, PCA, VIF, partial correlation tools
scripts/01_baseline_correlations.py   End-to-end runner
results/                    Outputs (CSVs + summary.txt)
data/                       Cached ESOL CSV
```

## Implemented indices (30)

**Degree-based (18):** First & Second Zagreb (M1, M2), Modified Zagreb (mM1,
mM2), Forgotten (F), Randić (R), Sum-connectivity (SCI), Harmonic (H),
Geometric-arithmetic (GA), Arithmetic-geometric (AG), Atom-bond-connectivity
(ABC), Atom-bond-sum-connectivity (ABS), Augmented Zagreb (AZI), Sombor
(SO), Reduced Sombor (SO_red), Albertson irregularity (Alb), Sigma, Reduced
First Zagreb (redM1)

**Distance-based (8):** Wiener (W), Hyper-Wiener (WW), Harary (HR), Schultz
(MTI), Gutman, Mostar, Szeged (Sz), Padmakar–Ivan (PI)

**Spectral (3):** Estrada (EE), Graph Energy, Spectral Radius

**Other (1):** Balaban J

## Methodology notes

- H-suppressed molecular graph (heavy atoms only, bonds as unweighted edges)
  — standard convention for classical topological indices.
- Multi-component SMILES (salts) reduced to largest fragment.
- Single-atom / disconnected graphs dropped.
- Floyd–Warshall for all-pairs shortest paths (fast enough for ≤50-atom
  molecules; ESOL run completes in 2s on a laptop).
- Mostar / Szeged / PI use strict-inequality vertex partitioning (standard).
- VIF = ∞ values indicate exact linear dependence within the descriptor
  block; this is a real signal of severe redundancy, not a numerical error.

## What is *not* yet built (next phases)

- Phase 2: scaffold-split k-fold CV on ESOL + FreeSolv + Lipophilicity
- Phase 3: Modern baselines (Chemprop, RDKit2D+RF, Mordred+XGB)
- Phase 4: Graph-theoretic bounds for new index (parallel track)
- Phase 5: OECD QSAR applicability domain (leverage / Williams plot)

Each is conditional on the new index passing the kill-test in phase 1.
