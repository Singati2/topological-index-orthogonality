# Implementation priority plan

Compiled after the second external review of the project
(round 2 of Woz + Jobs reviews; harsh outside-review).
Ranked by impact-per-effort. The first three are scoped to land
inside the current manuscript drafts; the remaining items are
larger and require either further empirical work or senior-author
contributions.

## Tier 1 — high impact, in-scope this iteration

| # | Action | Status | Notes |
|---|---|---|---|
| 1 | Add an ML benchmark with one new dataset | **DONE (this commit)** | `scripts/11_ml_benchmark.py`, `results/ml_benchmark.{csv,md}`. BBBP added as the additional classification dataset. RandomForest on four feature configurations on four datasets; pairwise-pruned matches `full` performance with 4--5x fewer features. On BBBP the pairwise-pruned model marginally improves ROC-AUC. On Lipophilicity the strict combined screen kills every feature, which is reported as an honest failure mode. |
| 2 | Position the pipeline against existing feature-selection methods | **pending** | Paper 2 should cite and briefly compare to mRMR, Boruta, mutual information, LASSO. Paper 2 §2 ("Related work") already mentions intercorrelation analyses by Movahedi et al. but the standard ML feature-selection toolbox is not addressed. One paragraph + 4-5 citations is sufficient. |
| 3 | Fill `[REFERENCE NEEDED]` placeholders | **pending** | Todeschini & Consonni 2009; Mauri / Dragon; Vukičević--Furtula 2009 (GA); Movahedi 2025 arXiv ID; Gutman--Borovićanin BID survey citation. Resolve before any submission. |
| 4 | Frame the Wuzi family explicitly as a stress test | **partially done** | Already framed as "parametric case study"; Paper 1 abstract sharpened from "novel" to "distinguishing feature" of the imbalance factor. One more sentence in Paper 2 introduction explicitly identifying the Wuzi family as the negative-control test of the screening pipeline would close this. |
| 5 | Add `tests/` with reproducibility checks | **pending** | At minimum: test that the closed-form Wuzi values match `src/wuzi_analytical.py` on a fixed seed; test that the 30 baseline indices reproduce frozen reference values on a known graph (e.g. P_5 or C_6). 3-4 tests are sufficient for a methodology paper. |

## Tier 2 — important, scoped for the next iteration

| # | Action | Status | Notes |
|---|---|---|---|
| 6 | Threshold-sensitivity table | **DONE in round 2** | `results/cross_dataset_summary.{csv,md}` now contains the 3x3x3 grid of pairwise-and-combined pass counts. Conclusion: the pairwise-only screen is over-permissive at $\tau_r = 0.99$; the combined screen is robustly zero across all 27 threshold combinations on the three regression datasets. The new BBBP dataset has not yet been added to the threshold-sensitivity grid; doing so is a half-day's work. |
| 7 | Resolve the Paper 2 / Paper 1 circular cite | **partially done** | Paper 2 now has its own self-contained 3-row Wuzi verdict table (`tab:wuzi_verdict_paper1`); the cross-reference to Paper 1 remains as a pointer. Paper 2 still cannot be submitted before Paper 1 has at least an arXiv identifier. Enforce: do not submit Paper 2 until Paper 1 has an arXiv ID. |
| 8 | Authorship / handle clarification | **pending** | The GitHub account `Singati2` does not match the cited authors. One-line note in `README.md` "Repository maintainer: Ganesh Shiwakoti (GitHub: Singati2)" is sufficient. |
| 9 | Update `Conjecture 6.4` after the falsification | **DONE in round 2** | The original caterpillar conjecture was falsified by the `is_caterpillar` predicate added in `scripts/10_wuzi_extremal_search.py`. The conjecture has been weakened to: at $(-1, -1, 1)$ and $n \ge 7$, the maximizer is neither path, star, nor caterpillar. Empirically supported across all 8 sampled tree orders. |

## Tier 3 — larger asks, gated on senior-author engagement

| # | Action | Status | Notes |
|---|---|---|---|
| 10 | Complete Paper 1 §4 (bounds in $(n, m, \Delta, \delta)$) | **pending; needs Arockiaraj sir** | Currently 4 rigorous theorems + 4 [PROOF SKETCH] markers (suppressed in PDF). MATCH-style submission requires the sharp Nordhaus--Gaddum bound and the sharp $\Delta$-dependent upper bound. |
| 11 | Complete Paper 1 §5 (bounds via classical indices) | **pending; needs Arockiaraj sir** | The ratio-bound theorem (rigorous) gives a uniform template; sharper per-index bounds remain [PROOF SKETCH]. |
| 12 | Complete Paper 1 §6 (extremal characterization) | **pending; needs Arockiaraj sir** | All three graph classes (trees, unicyclic, bicyclic) currently have conjectures only. The most interesting open problem is the non-classical extremal at $(-1, -1, 1)$ for trees of order $n \ge 7$ (Conjecture 6.4). |
| 13 | Expand BBBP threshold sensitivity | **pending** | `scripts/03_cross_dataset_summary.py` operates on the three regression datasets only. Extending to BBBP requires re-running `scripts/02_wuzi_grid_search.py` on BBBP (the target column needs to be the binary `p_np` for the partial-correlation computation, which uses linear regression — re-evaluation needed against logistic-regression residuals for the classification setting). Roughly a half-day. |
| 14 | One additional regression dataset for the combined-screen failure mode | **pending** | Lipophilicity is the only one where the combined screen kills all features; replicating on a second wide regression target (e.g. QM9 internal energy on a 10k subset) would test whether this is a Lipophilicity quirk or a general behaviour of QSPR with weak topology-target correlation. Not in this iteration. |

## What was explicitly NOT adopted from the external review

The external review recommended several aggressive changes that are
declined here:

- **Drop the math paper / Paper 1 entirely.** Declined: this would
  remove the senior author's contribution and the realistic path to
  a mid-tier specialty publication (J. Math. Chem., MATCH).
- **Expand to 6--8 datasets, 50--60 indices, 500-point grid.** Declined
  as scope expansion. One additional dataset (BBBP) was added per the
  in-scope directive.
- **Reframe MATCH and J. Math. Chem. as low-impact venues.** Declined:
  these are the realistic venues for the math program of Paper 1.
  The DSO paper template (Movahedi et al. 2026, MATCH 95:141--162)
  is taken seriously by the math-chem community.

## Project publishability assessment (updated after the BBBP + ML benchmark addition)

- **Paper 2 (methodology) before this iteration:** conditionally
  publishable at Mol. Informatics / SAR & QSAR Env. Res. tier.
  Main referee concerns: the 15/20 -> 4/20 finding was on a single
  dataset (ESOL); no downstream ML evidence; no positioning against
  existing feature-selection methods.
- **Paper 2 (methodology) after this iteration:** the ML-benchmark
  result strengthens the paper to the point where J. Cheminformatics
  becomes a defensible submission target. Specifically: the
  pairwise-pruned screen now demonstrably matches `full`-feature
  RandomForest performance on three regression datasets, marginally
  improves ROC-AUC on a classification dataset, and the strict
  combined screen has one honest failure mode (Lipophilicity).
  The BBBP addition diversifies the empirical base across task type.
  Items 2--3 of the Tier 1 list above remain.
- **Paper 1 (math) assessment unchanged:** publishability still gated
  on Sections 4--6 (Tier 3) being completed by Arockiaraj sir.
