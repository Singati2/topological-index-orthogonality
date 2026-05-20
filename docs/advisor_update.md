# Project Update for Dr. Michael Arockiaraj

**Project:** Topological-index redundancy analysis with parametric case study
**Repository:** <https://github.com/Singati2/topological-index-orthogonality>
**Collaborators:** Advik Natarajan (Loyola College, Chennai), Ganesh Shiwakoti
(FAU, Department of Mathematics and Statistics)
**Date:** 2026-05-20

## TL;DR (three sentences)

We have a complete numerical pipeline plus LaTeX-ready Paper 2
(methodology) and Paper 1 (Wuzi mathematical case study with closed
forms, bracket bound, ratio bound, and an empirical observation
about a non-classical extremal at $(\alpha,\beta,\gamma)=(-1,-1,1)$).
What we need from you, sir, is the formal proof / refinement of
Sections 4, 5, and 6 of Paper 1 (bounds in $(n, m, \Delta, \delta)$,
bounds via classical indices, and the extremal-graph
characterization, including the non-classical extremal at
$(-1,-1,1)$). With those sections completed Paper 1 is submittable
to a math-chem venue within a few weeks; Paper 2 follows once
Paper 1 has at minimum an arXiv ID.

## Most interesting finding to-date (one paragraph)

An exhaustive enumeration of non-isomorphic trees of order
$5 \le n \le 12$ shows that the path / star extremal pattern (path
= minimizer, star = maximizer) holds in the positive-parameter
region $\alpha,\beta,\gamma \ge 0$; the extremals flip in the
Randić-style sign region; and at the mixed-sign triple
$(\alpha,\beta,\gamma) = (-1,-1,1)$ the maximizer is *neither*
the path *nor* the star for $n \ge 6$, and for $n \ge 7$ it is
also *not a caterpillar* (verified by the `is_caterpillar`
predicate in `scripts/10_wuzi_extremal_search.py`). The structural
characterization of this non-classical extremal regime
(`docs/paper1_wuzi.tex`, Conjecture 6.4 / label
`conj:caterpillar_max`) is the most plausibly novel mathematical
observation in the project and is the natural anchor for
Section 6 of Paper 1.

---

## 1. What the project is about (brief)

We began with the question Advik raised — *"how do we show that a new
topological index is better than existing ones?"* — and developed an
**open-source orthogonality-screening pipeline** that tests whether any
proposed scalar topological index contributes information not already
captured by the standard 30-index baseline on real chemistry.

Along the way, Advik proposed the **Wuzi parametric index family**:

```
W(G; α, β, γ) = Σ_{uv ∈ E(G)} (d_u d_v)^α · (d_u + d_v)^β · exp(γ |d_u − d_v|/(d_u + d_v))
```

a three-parameter generalization of degree-based edge-sum indices that
contains M_1, M_2, Randić, sum-connectivity, harmonic, and Sombor-like
behavior as special cases.

## 2. Why we reframed the project (and what it is *not*)

After implementing the pipeline and running it across three benchmark
datasets (ESOL, FreeSolv, Lipophilicity), our 100-point grid search on
(α, β, γ) showed that the Wuzi family is **statistically redundant** with
the classical 30-index baseline (max |r| ≥ 0.96 with mM_1 / mM_2 / SO_red
at every grid point on every dataset).

A literature check then confirmed that several alternative "novel index"
candidates we considered (Shannon degree entropy, mean clustering
coefficient, network-science centralities, etc.) are also already
published as named chemical-graph invariants in the
information-theoretic, centrality, and motif-based topological-index
literature (Bonchev and collaborators on information indices;
Cao--Dehmer--Shi on entropy indices; Ghorbani and collaborators
on degree-entropy variants; the MOLTOP graph-classification
descriptor). A curated reference list is maintained in
`docs/literature_notes.md`.

We therefore **do not claim** that:

- Wuzi is a superior topological index.
- Entropy, motif, centrality, or other "alternative-family" indices are
  newly invented.
- All topological indices are useless.

We **do claim**, and have evidence for:

- Endpoint-degree edge-sum (BID) indices on Δ ≤ 4 molecular graphs lie
  in a vector space of dimension ≤ 10, spanned by the 10 edge-degree-pair
  counts m_ij (well-known BID structure; we make the dimension bound
  explicit and draw the redundancy consequence).
- The empirical effective rank of a 30-index baseline is 3 PCA components
  for 95 % variance, consistently across ESOL (n=1127), FreeSolv (n=639),
  and Lipophilicity (n=4200).
- The Wuzi family lies inside this 10-D subspace at every (α, β, γ); its
  100-point sweep fails the |r|≥0.95 orthogonality criterion against the
  baseline on all three datasets.
- An open-source pipeline can mechanize this redundancy check in seconds
  for any future scalar topological index proposal.

## 3. Why two papers make sense

The work naturally splits into two contributions for two different
audiences:

**Paper 1 (we recommend writing first):**
*"The Wuzi Index Family: Graph-Theoretic Properties, Bounds, Extremal
Graphs, and Sensitivity Analysis"*

Follows the mathematical-chemistry tradition (analog of Movahedi, Gutman,
Redžepović & Furtula, *MATCH Commun. Math. Comput. Chem.* 95(1), 2026,
141–162):

- Definition + special-graph closed-forms + edge-contribution analysis
- BID-basis observation (10-dim ceiling on Δ ≤ 4)
- Bounds in (n, m, Δ, δ)
- Bounds via classical indices (ratio-bound theorem)
- Extremal graph characterization (trees, unicyclic, bicyclic; the
  headline non-classical extremal at (−1, −1, 1))
- Sensitivity Section: octane physicochemical-property correlations,
  degeneracy on 106 trees of order 10, structure sensitivity on
  the 75 decane isomers

The cross-dataset redundancy analysis and the orthogonality-screening
pipeline are NOT in Paper 1 — they are the central content of the
methodology paper (Paper 2).

Target venues: J. Math. Chem., MATCH, AKCE, Iranian J. Math. Chem.,
SAR & QSAR Env. Res. Realistic.

**Paper 2 (drafted in parallel, submitted second):**
*"Orthogonality Screening of Topological Indices for QSPR Modeling: A
Structural and Empirical Redundancy Analysis"*

The broader methodology / software / benchmark paper. References Paper 1
as the worked case study. Target venues: J. Cheminformatics, J. Chem. Inf.
Model., Molecular Informatics. More ambitious; depends on Paper 1 going
through smoothly.

## 4. What is already complete

All numerical infrastructure is in place: the $30$-index baseline,
the $100$-point Wuzi grid sweep, the cross-dataset comparison, the
octane prediction (MATCH §5.1 analog on the $18$ octane isomers),
the degeneracy analysis on $106$ non-isomorphic trees of order
$10$, the structure-sensitivity analysis on the $75$ decane
isomers, and the reproducible end-to-end pipeline (runs in well
under a minute on a single laptop). The Wuzi family closed-form
values on nine standard graph classes have been derived
analytically and verified numerically by the test code in
`src/wuzi_analytical.py`. The edge-contribution function $\psi$
and the Edge-Degree-Pair Basis observation (dimension $\le 10$
on $\Delta \le 4$) are documented in
`docs/edge_contribution_analysis.md` and
`docs/theoretical_foundation.md` respectively. A scaffolding doc
(`docs/wuzi_bounds_strategy.md`) gives rigorous templates for
the bracket bound, the ratio-bound theorem in terms of any
classical BID index, and parameter-monotonicity statements;
computational extremal searches are in
`scripts/10_wuzi_extremal_search.py`. The bibliographic
discipline is maintained in `docs/literature_notes.md`.

## 5. What remains — the math sections of Paper 1

This is where your input would be most valuable, sir. Three sections need
graph-theoretic derivations following the template of MATCH 95:141–162
and Movahedi (2025) arXiv preprint:

| Section | Content needed                                                     |
|---|---|
| **4** | Sharp $\Delta$-dependent and $\delta$-dependent bounds, Nordhaus–Gaddum, Polya–Szego (currently `[PROOF SKETCH]` placeholders) |
| **5** | Sharp per-index bounds for $R$, $SO$, $GA$, $H$, $ABC$ (ratio-bound theorem is rigorous; per-index sharp forms are `[PROOF SKETCH]` placeholders) |
| **6** | Extremal-graph characterization among trees, unicyclic, bicyclic graphs, and the non-classical extremal at $(-1,-1,1)$ (Conjecture 6.4) |

Each is comparable in scope to a single section of the MATCH DSO paper.
The sensitivity-analysis section (§7) and all numerical infrastructure
are already in place so the math sections can drop directly into a
prepared manuscript skeleton.

## 5a. Mathematical scaffolding added (2026-05-19)

To make it as easy as possible for Advik and yourself, sir, to fill
in Sections 3, 4, and 5, we have added rigorous templates,
ratio-bound calculations, and computational extremal searches.
These are **not** final theorem replacements for the bounds work;
they are starting material.

| Artifact | What it provides |
|---|---|
| `docs/wuzi_bounds_strategy.md` | Bracket bound (rigorous), ratio-bound theorem (rigorous), parameter monotonicity, extremal-graph conjectures with heuristic motivation, list of proof techniques. |
| `scripts/09_wuzi_bounds_tables.py` | Brute-forces the ratio constants $c_{\min}^h, c_{\max}^h$ for $h \in \{M_1, M_2, R, SO, GA, H, ABC\}$ on $\Delta \le 4$ at eight parameter triples; sanity-checks against exact Wuzi-classical reductions. |
| `scripts/10_wuzi_extremal_search.py` | Enumerates all non-isomorphic trees of order $5 \le n \le 12$ and identifies the observed minimizer / maximizer of $W$ at each parameter triple. |
| `results/wuzi_bounds_ratio_tables.{csv,md}` | Output of script 09. |
| `results/wuzi_extremal_trees.{csv,md}` | Output of script 10; includes whether the observed extremals match $P_n / S_n$. |

**Empirical findings of the extremal search** (trees only, $5 \le n \le 12$):
- For $(\alpha, \beta, \gamma) \in \{(1,0,0),\ (0,0,1),\ (0,0,2),\ (1,1,1)\}$:
  observed minimizer = $P_n$ and observed maximizer = $S_n$ at every $n$.
- For $(\alpha, \beta, \gamma) \in \{(-1/2,0,0),\ (0,-1/2,0),\ (0,-1,0)\}$
  (the Randić / sum-connectivity / harmonic sign region):
  observed minimizer = $S_n$ and observed maximizer = $P_n$ at every $n$
  (the extremals flip with the sign of the parameters).
- For the mixed-sign triple $(-1, -1, 1)$: observed extremals are
  *neither* $P_n$ nor $S_n$. This is the most interesting case
  computationally and the natural target for a "genuinely new"
  result in Section 6 of Paper 1 (Conjecture 6.4).

**What this still does not provide.** Formal proofs of any of the
extremal conjectures above; equality-condition analysis for the
Jensen / Cauchy-Schwarz bounds; bounds in the full $(\alpha, \beta, \gamma)$
parameter region rather than just at the eight sampled triples.
These are exactly the contributions reserved for the senior author.

## 6. Honest publishability status

- Paper 1 **is not currently publishable**. Sections 3–5 (the math) are
  unwritten; without them the manuscript would be desk-rejected by a
  mathematical-chemistry venue. Numerical Section 6 alone does not make
  a math-chem paper.
- With Sections 3–5 derived to the standard of the MATCH DSO paper, the
  manuscript is realistically targetable at *J. Math. Chem.*, *MATCH*,
  *AKCE Int. J. Graphs Comb.*, or *SAR & QSAR Env. Res.* within 6–8 weeks.
- Paper 2 is more ambitious and depends on Paper 1 landing first.

## 7. What we request

If you would be open to it, we would like to:

1. Discuss the Wuzi family definition and confirm the family is worth
   the bounds derivation (versus folding into the screening paper only).
2. Co-author Sections 3, 4, and 5 of Paper 1 — Advik as lead, you as
   senior author providing the math-chem derivations and connections to
   your group's prior work, Ganesh on the numerical / methodology side
   (or in acknowledgments only, at your discretion).
3. Advise on venue choice given the BID-family parametric framing.

The full repository is public and all code is reproducible end-to-end
(scripts run in well under a minute except Lipophilicity, which takes ~15 s).

We look forward to your thoughts.

— Advik Natarajan & Ganesh Shiwakoti
