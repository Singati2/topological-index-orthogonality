# Project Update for Dr. Michael Arockiaraj

**Project:** Topological-index redundancy analysis with parametric case study
**Repository:** <https://github.com/Singati2/topological-index-orthogonality>
**Collaborators:** Advik Natarajan (Loyola College, Chennai), Ganesh Shiwakoti
(FAU, Biostatistics PhD)
**Date:** 2026-05-18

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
published as named chemical-graph invariants ([REFERENCE NEEDED:
Ghorbani et al. 2018, Cao-Dehmer-Shi 2014, Aslam et al. 2023,
Bonchev-Mekenyan-Trinajstić 1981 "superindex," MOLTOP ECAI 2024]).

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

**Paper 2 (we recommend writing first):**
*"The Wuzi Index Family: Graph-Theoretic Properties, Bounds, Extremal
Graphs, and Redundancy Analysis"*

Follows the mathematical-chemistry tradition (analog of Movahedi, Gutman,
Redžepović & Furtula, *MATCH Commun. Math. Comput. Chem.* 95(1), 2026,
141–162):

- Definition + special-graph closed-forms + edge-contribution analysis
- Bounds in (n, m, Δ, δ)
- Bounds in classical indices
- Extremal graph characterization (trees, unicyclic, bicyclic)
- Numerical Section: octane prediction, intercorrelations across three
  datasets, PCA effective rank, degeneracy, structure sensitivity

The redundancy analysis is included **inside Paper 2** as Section 6
("Numerical work"), framed constructively, and positioned the same way
the MATCH DSO paper Section 5.2 positions its own intercorrelation
observation.

Target venues: J. Math. Chem., MATCH, AKCE, Iranian J. Math. Chem.,
SAR & QSAR Env. Res. Realistic.

**Paper 1 (drafted in parallel, submitted second):**
*"Orthogonality Screening of Topological Indices for QSPR Modeling: A
Structural and Empirical Redundancy Analysis"*

The broader methodology / software / benchmark paper. References Paper 2
as the worked case study. Target venues: J. Cheminformatics, J. Chem. Inf.
Model., Molecular Informatics. More ambitious; depends on Paper 2 going
through smoothly.

## 4. What is already complete

| Component                                       | Status      | Location |
|---|---|---|
| Wuzi definition + identity tests vs M_2/R/SCI/H | ✅ done     | `src/wuzi_index.py` |
| Closed-form Wuzi values: K_n, C_n, P_n, S_n, K_{p,q}, Q_k, W_n, F_p, regular | ✅ done (numerically verified) | `src/wuzi_analytical.py` |
| Edge-contribution function ψ(x,y; α,β,γ) analysis | ✅ done   | `docs/edge_contribution_analysis.md` |
| Edge-Degree-Pair Basis theorem (Δ ≤ 4: dim ≤ 10) | ✅ done   | `docs/theoretical_foundation.md` |
| 30-index baseline on ESOL (n=1127)              | ✅ done     | `results/esol/` |
| 30-index baseline on FreeSolv (n=639)           | ✅ done     | `results/freesolv/` |
| 30-index baseline on Lipophilicity (n=4200)     | ✅ done     | `results/lipophilicity/` |
| Wuzi 100-point grid sweep on each dataset       | ✅ done     | `results/<dataset>/wuzi_grid.csv` |
| Cross-dataset redundancy comparison             | ✅ done     | `results/cross_dataset_summary.md` |
| Octane prediction (n=18, MATCH §5.1 analog)     | ✅ done     | `results/octane_prediction.md` |
| Degeneracy on trees of order 10 (§5.3 analog)   | ✅ done     | `results/wuzi_degeneracy.md` |
| Structure sensitivity on 75 decanes (§5.4 analog) | ✅ done   | `results/structure_sensitivity.md` |
| Open-source pipeline (Python; runs in seconds)  | ✅ done     | `scripts/*.py` |
| Bibliographic discipline document               | ✅ done     | `docs/literature_notes.md` |

## 5. What remains — the math sections of Paper 2

This is where your input would be most valuable, sir. Three sections need
graph-theoretic derivations following the template of MATCH 95:141–162
and Movahedi (2025) arXiv preprint:

| Section | Content needed                                                     |
|---|---|
| **3** | Bounds on W(G; α, β, γ) in terms of (n, m, Δ, δ), with sharp regular-graph extremals |
| **4** | Bounds on W in terms of M_1, M_2, R, SO, GA, H, ABC (Movahedi 2025 template) |
| **5** | Extremal-graph characterization among trees, unicyclic graphs, bicyclic graphs |

Each is comparable in scope to a single section of the MATCH DSO paper.
The numerical sections (§6) and all infrastructure are already in place
so the math sections can drop directly into a prepared manuscript skeleton.

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
  result in Section 5 of Paper 2.

**What this still does not provide.** Formal proofs of any of the
extremal conjectures above; equality-condition analysis for the
Jensen / Cauchy-Schwarz bounds; bounds in the full $(\alpha, \beta, \gamma)$
parameter region rather than just at the eight sampled triples.
These are exactly the contributions reserved for the senior author.

## 6. Honest publishability status

- Paper 2 **is not currently publishable**. Sections 3–5 (the math) are
  unwritten; without them the manuscript would be desk-rejected by a
  mathematical-chemistry venue. Numerical Section 6 alone does not make
  a math-chem paper.
- With Sections 3–5 derived to the standard of the MATCH DSO paper, the
  manuscript is realistically targetable at *J. Math. Chem.*, *MATCH*,
  *AKCE Int. J. Graphs Comb.*, or *SAR & QSAR Env. Res.* within 6–8 weeks.
- Paper 1 is more ambitious and depends on Paper 2 landing first.

## 7. What we request

If you would be open to it, we would like to:

1. Discuss the Wuzi family definition and confirm the family is worth
   the bounds derivation (versus folding into the screening paper only).
2. Co-author Sections 3, 4, and 5 of Paper 2 — Advik as lead, you as
   senior author providing the math-chem derivations and connections to
   your group's prior work, Ganesh on the numerical / methodology side
   (or in acknowledgments only, at your discretion).
3. Advise on venue choice given the BID-family parametric framing.

The full repository is public and all code is reproducible end-to-end
(scripts run in well under a minute except Lipophilicity, which takes ~15 s).

We look forward to your thoughts.

— Advik Natarajan & Ganesh Shiwakoti
