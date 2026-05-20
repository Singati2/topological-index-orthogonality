# Advisor-ready summary

Concise summary of where the manuscripts stand after Arockiaraj
sir's restructure directive. Intended as a 5-minute orientation
note before reading the drafts.

## 0. What we are asking from you (one screen)

If you have time for only three actions in this round:

1. **Look at Conjecture 6.4 in `docs/paper1_wuzi.tex` §6.4** (the
   non-classical extremal at $(\alpha, \beta, \gamma) = (-1, -1, 1)$
   for trees of order $n \ge 7$). The computational evidence
   is in `results/wuzi_extremal_trees.md`. If this looks like a
   characterisable structure to you, it would be the most
   plausibly novel mathematical contribution of Paper 1; if not,
   we will drop it and present the result only as a remark.

2. **Decide whether you want to co-author Paper 1 §4 (Nordhaus-Gaddum
   and sharp $\Delta$-dependent bounds), §5 (sharp per-index bounds
   for $R$, $SO$, $GA$, $H$, $ABC$), and §6 (formal extremal
   characterizations) with Advik**, or whether to ship Paper 1
   with the current `[PROOF SKETCH]` placeholders (now made visible
   in the rendered PDF) converted to "open problem" remarks or
   numbered Conjectures. Either choice unblocks submission.

3. **Confirm the affiliation line** "Department of Mathematics and
   Statistics, Florida Atlantic University" (Ganesh) and "Department
   of Mathematics, Loyola College, Chennai" (Advik, yourself) are
   correct as they appear in both `.tex` files and in `CITATION.cff`.

The remainder of this document is the longer-form context.

## 1. What changed after your feedback

Per your direction, the two manuscripts have been restructured so
they have minimal overlap. Previously, both papers carried the
Wuzi family math AND the cross-dataset screening / ML benchmark
content; this duplication has been removed. The Wuzi family is
now the central object of Paper 1 only; the screening pipeline is
the central object of Paper 2 only.

The two `.tex` files were also renamed to match your ordering:

- `docs/paper1_wuzi.tex` -- Wuzi index paper (graph theory)
- `docs/paper2_orthogonality_screening.tex` -- methodology paper

The previous filenames (`paper2_wuzi.tex`, `paper1_orthogonality_screening.tex`)
have been retired via `git mv` (history preserved).

## 2. How Paper 1 and Paper 2 are now separated

Both papers cite each other but do not duplicate work. Paper 1
references Paper 2 in exactly three places (abstract,
introduction, conclusion) with a single sentence: *"A companion
methodology paper studies multi-dataset redundancy screening of
topological indices."* Paper 2 references Paper 1 in its
introduction and Section 6 ("Case study 1: the Wuzi parametric
family") to cite the math, then proceeds to use Wuzi only as one
screened candidate among others.

## 3. What Paper 1 contains

A mathematical-chemistry / graph-theory paper in the tradition of
MATCH-style index papers:

- Section 1: Introduction
- Section 2: Preliminaries (notation, Wuzi definition, special-case
  reductions to $M_2$, Randi\'c, sum-connectivity, half-harmonic,
  edge-contribution function $\psi$, standard lemmas)
- Section 3: The Edge-Degree-Pair Basis observation
  (Theorem~3.1: any BID index on $\Delta \le 4$ molecular graphs
  lies in a 10-dimensional subspace)
- Section 4: Bounds in $(n, m, \Delta, \delta)$
  -- bracket bound (rigorous, tight on regular graphs)
  -- Cauchy-Schwarz geometric-mean inequality (rigorous)
  -- Jensen-type bounds via $M_1$ (rigorous, both $\beta \ge 1$
  lower and $0 < \beta < 1$ upper)
  -- Upper bound in $(\Delta, m)$ in the non-negative parameter
  region (rigorous, Proposition 4.7)
  -- Lower bound in $(\delta, m)$ in the non-negative parameter
  region (rigorous, Proposition 4.8)
  -- Nordhaus-Gaddum bound (CONJECTURE 4.9, open in mixed-sign
  $\gamma \ne 0$ regime)
- Section 5: Bounds via classical indices
  -- ratio-bound theorem (rigorous, applies uniformly to any
  classical BID index; Theorem 5.2)
  -- bound via generalized $R$ and $\chi$ (rigorous,
  Proposition 5.4)
  -- indirect bound via Sombor $SO$ and $M_1$ for $0 < \beta < 1$
  (rigorous, Proposition 5.5)
  -- sharp two-sided bounds via $GA$, $H$, $ABC$ (CONJECTURES 5.6,
  5.7, open)
- Section 6: Extremal graphs
  -- Conjectures for trees, unicyclic, bicyclic in positive
  parameter region
  -- Conjecture 6.4 (mixed-sign non-classical extremal at
  $(\alpha, \beta, \gamma) = (-1, -1, 1)$): computationally
  verified across $n = 5, \ldots, 12$ via
  `scripts/10_wuzi_extremal_search.py`; for $n \ge 7$ the
  observed maximum is neither $P_n$, nor $S_n$, nor a caterpillar
- Section 7: Sensitivity analysis
  -- octane competitor comparison (Wuzi at 8 parameter triples vs
  $M_1, M_2, R, SO, GA, H, ABC, AZI$ on the 18-octane benchmark,
  Pearson correlations with 5 physicochemical properties)
  -- degeneracy on the 106 non-isomorphic trees of order 10
  -- structure sensitivity on the 75 decane isomers
- Section 8: Discussion
- Section 9: Conclusion

## 4. What Paper 2 contains

A software/methodology paper in the tradition of *J. Cheminformatics*
/ *Molecular Informatics* methods papers:

- Section 1: Introduction (leads with the 15-of-20-vs-4-of-20
  over-permissiveness finding on ESOL)
- Section 2: Related work
- Section 3: Datasets and implementation (ESOL, FreeSolv,
  Lipophilicity, BBBP)
- Section 4: The Edge-Degree-Pair Basis observation
  (cites Paper 1 for the proof)
- Section 5: The screening protocol (correlation, PCA, partial
  correlation, VIF, combined criterion)
- Section 6: Case study 1 -- the Wuzi parametric family
  (300-grid screening across three datasets)
- Section 7: Case study 2 -- 20 alternative-family candidates on
  ESOL
- Section 8: Downstream ML benchmark (5-fold CV across 4 datasets,
  Figure 9 ml_benchmark)
- Section 9: Discussion
- Section 10: Conclusion

## 5. What was removed to avoid overlap

From Paper 1 (Wuzi):
- ESOL / FreeSolv / Lipophilicity cross-dataset screening
- BBBP dataset
- Downstream RandomForest ML benchmark, fig9
- PCA / pairwise-pruned / combined-pruned feature configurations
- All other methodology / cheminformatics-pipeline content
  (zero hits in a strict forbidden-term grep)

From Paper 2 (Screening):
- Wuzi closed-form derivations
- Wuzi bounds (bracket, ratio, Cauchy-Schwarz, Jensen)
- Wuzi extremal-graph analysis
- Octane / decane / standalone-tree sensitivity sections

## 6. What still needs your mathematical input in Paper 1

These are the items currently labelled as Conjectures in the LaTeX
source. None present a fabricated theorem; the round-10 work
upgraded every "Theorem -- skeleton" to either a rigorous
Proposition or an honest Conjecture.

Section 4 (Bounds in $(n, m, \Delta, \delta)$):
- Nordhaus-Gaddum bound $W(G) + W(\bar G) \ge f_n(\alpha, \beta, \gamma)$
  (Conjecture~4.9, `conj:NG`) -- open in the mixed-sign $\gamma \ne 0$
  regime; the explicit $f_n$ is the missing object.
- The upper bound in $(\Delta, m)$ (Proposition~4.7) and the lower
  bound in $(\delta, m)$ (Proposition~4.8) are now *fully proved*
  in the non-negative parameter region; the corresponding sharp
  versions in the mixed-sign region remain natural follow-ups.

Section 5 (Bounds via classical indices):
- Sharp two-sided bound via $GA$ (Conjecture~5.6, `conj:GA`).
- Sharp two-sided bounds via $H$ and $ABC$ (Conjecture~5.7,
  `conj:H_ABC`).
- The ratio-bound theorem (Theorem~5.2, `thm:ratio_bound`), the
  Cauchy--Schwarz bound via generalized $R$ and $\chi$
  (Proposition~5.4, `thm:R_chi`), and the indirect Sombor / $M_1$
  bound for $0 < \beta < 1$ (Proposition~5.5, `thm:sombor`) are
  now *fully proved*.
- (Optional) Closed-form bounds paired with the Albertson
  irregularity $\mathrm{Alb}(G) = \sum_e |d_u - d_v|$ for the
  $\gamma$-axis (not yet stated; listed in §9 future work).

Section 6 (Extremal graphs):
- Formal extremal-tree characterization in the positive parameter
  region (Conjecture~6.1, `conj:trees_beta_pos`).
- Formal extremal characterization for unicyclic graphs
  (Conjecture~6.2, `conj:unicyclic`).
- Bicyclic extremals (Conjecture~6.3, `conj:bicyclic` --
  candidate-family-level only).
- Mixed-sign non-classical extremal at $(\alpha, \beta, \gamma) = (-1, -1, 1)$
  (Conjecture~6.4, `conj:caterpillar_max` -- computationally
  verified for $n = 5, \ldots, 12$; the formal structural
  characterization is the central open problem).

The companion strategy document `docs/wuzi_bounds_strategy.md`
contains rigorous theorem templates, proof techniques, and
explicit ratio-bound computations that you and Advik can use as
scaffolding.

## 7. What is already complete in Paper 2

- Full screening protocol with combined criterion
- Cross-dataset baseline-redundancy table (3 regression datasets)
- 300-grid Wuzi screening verdict + commentary on the FreeSolv
  borderline grid point
- 20-candidate alternative-family case study with combined-pass
  breakdown
- 5-fold CV ML benchmark across 4 datasets (regression and
  classification), reported as mean $\pm$ std
- Software-reproducibility scaffolding: `tests/` (124 passing
  assertions), `.github/workflows/test.yml` (pytest-only CI),
  `CITATION.cff`

Paper 2 still needs final prose polishing in §9 (Discussion)
and §10 (Conclusion), and final reference resolution (a few
`[REFERENCE NEEDED]` slots in the bibliography).

## 8. Exact files to read first

In this order:

1. `docs/paper1_wuzi.tex` -- the Wuzi math paper
2. `docs/paper2_orthogonality_screening.tex` -- the methodology paper
3. `results/octane_competitor_comparison.md` -- Paper 1's main
   sensitivity-analysis table, in standalone Markdown form
4. `results/wuzi_extremal_trees.md` -- the computational evidence
   for Conjecture~6.4
5. `results/ml_benchmark.md` -- Paper 2's main empirical evidence,
   in standalone Markdown form

If only one file: `docs/paper1_wuzi.tex` (this is where your
expertise is most needed).

---

Repository: <https://github.com/Singati2/topological-index-orthogonality>

Latest commit on `main`: see `git log -1` at the time of reading.

Both papers compile only via Overleaf in the local environment
used here; no `pdflatex` is installed locally, so static checks
have been used (env-balance, ref/cite resolution, no `\\end{xxx>`
typos). Overleaf compilation is the recommended next step before
substantive reading.
