# Paper 2 — Wuzi family manuscript skeleton

Working template for the parametric-index / math-chem paper. Section
numbering mirrors **Movahedi, Gutman, Redžepović & Furtula, MATCH 95(1),
2026, 141-162** ("Diminished Sombor Index") for direct reviewer
comparability. Math sections marked `[NEEDS DERIVATION ...]` are honest
gaps awaiting Advik + Arockiaraj sir's input.

---

## Title — three options

1. *"The Wuzi Parametric Index Family: Definition, Bounds, Extremal
   Graphs, and Redundancy Analysis on Molecular Chemistry Benchmarks"*
2. *"A Three-Parameter Bond-Incident-Degree Index Family: Mathematical
   Properties and Empirical Redundancy"*
3. *"Wuzi Index: Closed Forms, Sharp Bounds, and the Limits of
   Parametric Extensions in the Bond-Incident-Degree Family"*

(Working title for now: **option 1**.)

## Authors

Advik Natarajan¹, Ganesh Shiwakoti², Michael Arockiaraj¹
¹ Department of Mathematics, Loyola College, Chennai
² Department of Biostatistics, Florida Atlantic University

(Author order to be confirmed with Arockiaraj sir; placeholder above.)

---

## Abstract (working draft)

> A new parametric degree-based topological index, the **Wuzi index
> family**, is introduced for a simple graph G:
>
> `W(G; α, β, γ) = Σ_{uv ∈ E(G)} (d_u d_v)^α (d_u + d_v)^β
>                 · exp(γ · |d_u − d_v|/(d_u + d_v))`
>
> for parameters (α, β, γ) ∈ ℝ³. Several established indices —
> first/second Zagreb, Randić, sum-connectivity, harmonic — are
> recovered as special cases. Closed-form values are derived for the
> Wuzi family on standard graph classes (complete, cycle, path, star,
> complete bipartite, hypercube, wheel, friendship, and regular
> graphs). Sharp upper and lower bounds are established in terms of
> the graph order n, size m, maximum and minimum degrees Δ and δ, and
> in terms of classical degree-based topological indices. Extremal
> graphs achieving equality are characterized among trees, unicyclic
> graphs, and bicyclic graphs. We additionally apply a transparent
> orthogonality-screening protocol to a 100-point grid in (α, β, γ)
> against a 30-index classical baseline on three QSPR benchmark
> datasets (ESOL, FreeSolv, Lipophilicity) totaling 6 000 molecules
> and report the empirical redundancy structure of the family. The
> open-source code and all numerical results are publicly available.

(~ 220 words. To be tightened for the venue limit.)

---

## 1. Introduction

`[NEEDS WRITING — ~1.5 pages]`

- The role of topological indices in QSPR; their proliferation; the
  recurring question of "what makes a new index a meaningful one?"
- Brief survey of the parametric-index tradition: general Randić R_α
  (Bollobás–Erdős 1998), general sum-connectivity χ_α (Zhou–Trinajstić
  2010), generalized Sombor SO_α, reduced/co-Sombor variants, Banhatti
  Sombor, Diminished Sombor [Movahedi/Gutman/Redžepović/Furtula 2026].
- Motivation for the Wuzi family: combine the Randić-style product
  axis with the sum-connectivity axis and the (novel here) exponential
  irregularity axis into one parametric family that subsumes
  multiple existing indices and introduces a new degree-of-freedom.
- Outline of contributions and paper structure.

## 2. Preparations

### 2.1 Notation

- Graph-theoretic notation (V, E, n, m, d_u, Δ, δ, K_n, C_n, P_n, S_n,
  K_{p,q}, Q_k, W_n, F_p, etc.). Identical to MATCH 95:141-162 p. 142.

### 2.2 Definition and special cases

`[CONTENT EXISTS — embed `src/wuzi_index.py` definition and the special
cases table from `docs/edge_contribution_analysis.md` §1.]`

The Wuzi family contains as special cases:

| (α, β, γ)         | Reduces to                                       |
|---|---|
| (0, 0, 0)          | `|E|` (number of edges)                         |
| (1, 0, 0)          | Second Zagreb M_2                                |
| (-1/2, 0, 0)       | Randić R                                         |
| (0, -1/2, 0)       | Sum-connectivity SCI                             |
| (0, -1, 0)         | Harmonic / 2                                     |
| (1/2, 0, 0)        | Σ_e √(d_u d_v)  (geometric-mean index variant)   |

### 2.3 Closed-form values on standard graphs

`[CONTENT EXISTS in src/wuzi_analytical.py, numerically verified.]`

Present as **Proposition 1** (analog of MATCH 95:141-162 Prop 1) with
items 1)–9) covering K_n, C_n, P_n, S_n, K_{p,q}, Q_k, W_n, D_p^(q)
(Dutch windmill — mechanical to add), F_p (friendship).

### 2.4 Edge-contribution function

`[CONTENT EXISTS in docs/edge_contribution_analysis.md.]` Migrate
inline. Includes the ψ(x, y; α, β, γ) decomposition, the 10-row BID
lattice table, partial-derivative monotonicity discussion, and the
Lemma 2.3.1 / Corollary 2.3.2 m·ψ_min ≤ W ≤ m·ψ_max bracket.

### 2.5 Lemmas used in proofs

`[NEEDS WRITING ~half page]` — Cauchy-Schwarz, Radon's inequality,
Polya-Szegő, Jensen on convex functions, Maclaurin sums-of-products,
power-mean. Same set used in Movahedi 2025.

## 3. Bounds in terms of n, m, Δ, δ

`[NEEDS DERIVATION BY ADVIK / AROCKIARAJ SIR — ~3 pages]`

Template: Section 3 of MATCH 95:141-162 has 13 theorems of this kind.
Suggested goals:

- **Theorem 3.1.** Sharp upper and lower bounds of the form
  m · ψ_min(δ, Δ) ≤ W ≤ m · ψ_max(δ, Δ) with regular-graph equality.
- **Theorem 3.2.** Bounds when α, β, γ are restricted to specific
  sign regions (positive, negative).
- **Theorem 3.3 (Cauchy-Schwarz type).** W² ≤ (something) ·
  (something), using one of the standard lemmas.
- **Theorem 3.4 (Nordhaus-Gaddum type).** W(G) + W(Ḡ) ≥ f(n).
- **Theorem 3.5–3.6.** Sharper bounds for trees / specific structures.

## 4. Bounds in terms of classical topological indices

`[NEEDS DERIVATION BY ADVIK / AROCKIARAJ SIR — ~3-4 pages]`

Template: Movahedi 2025 arXiv preprint (DSO bounds paper). Suggested
goals:

- **Theorem 4.1.** W in terms of M_1 (first Zagreb).
- **Theorem 4.2.** W in terms of M_2 (second Zagreb).
- **Theorem 4.3.** W in terms of Randić R.
- **Theorem 4.4.** W in terms of geometric-arithmetic GA.
- **Theorem 4.5.** W in terms of Sombor SO.
- **Theorem 4.6.** W in terms of harmonic H and ABC.
- **Theorem 4.7.** W in terms of Albertson irregularity Alb (the
  γ-axis pairing).

For special parameter values where Wuzi reduces to a classical
index, these bounds collapse to exact identities — flag those.

## 5. Extremal graphs

`[NEEDS DERIVATION BY ADVIK / AROCKIARAJ SIR — ~2-3 pages]`

Template: Section 4 of MATCH 95:141-162. Suggested goals:

- **Theorem 5.1.** Extremal graphs in the class of all connected
  graphs of order n.
- **Theorem 5.2.** Extremal trees of order n.
- **Theorem 5.3.** Extremal unicyclic graphs of order n.
- **Theorem 5.4.** Extremal bicyclic graphs of order n.

Identify which structures (path P_n, star S_n, cycle C_n, specific
"caterpillar" or "broom" trees) attain the min/max in each parameter
region. Many of these mirror the Albertson / Randić extremal
literature when γ dominates / α dominates respectively.

## 6. Numerical work

`[CONTENT EXISTS — full results in results/*.md]`

### 6.1 Prediction ability on octanes

`[CONTENT EXISTS in results/octane_prediction.md]`

18 octane isomers; properties bp, ΔH_f, ΔH_vap, S, ω from NIST
WebBook. Wuzi at canonical parameters and the 30-index baseline.
Threshold |r| > 0.8 ("predicts well") matches MATCH §5.1 convention.

### 6.2 Cross-dataset orthogonality screening

`[CONTENT EXISTS in results/cross_dataset_summary.md]`

100-point Wuzi grid against the 30-index baseline on **three**
QSPR benchmark datasets (ESOL n=1127, FreeSolv n=639, Lipophilicity
n=4200, total 5 966 molecules). Cross-dataset agreement on
effective rank (3 PCA components for 95 %) is a key finding. Frame
constructively: we extend the analytical move of Movahedi et al.
2026 §5.2 to a much larger, replicated setting with a defined
threshold (|r| ≥ 0.95) and PCA / partial-correlation diagnostics.

### 6.3 PCA effective rank

`[CONTENT EXISTS in results/<dataset>/pca_variance.csv]` — three
components capture 95 % of variance on all three datasets; five
components capture 99 %.

### 6.4 Degeneracy

`[CONTENT EXISTS in results/wuzi_degeneracy.md]` — 106
non-isomorphic trees of order 10; Wuzi at γ = 2 attains ~25 %
degeneracy, comparable to the Sombor-family values (17–22 %) in
MATCH 95:141-162 Figure 5.

### 6.5 Structure sensitivity

`[CONTENT EXISTS in results/structure_sensitivity.md]` — 75 decane
isomers; SS, Abr, SA metrics. Wuzi values in the same regime as
Sombor variants in MATCH 95:141-162 Table 2.

## 7. Discussion

`[NEEDS WRITING ~1.5 pages]`

Key points to make:

- The Wuzi family demonstrates the general phenomenon (formal in
  §2 BID basis observations, empirical in §6) that adding parameters
  to an endpoint-degree edge-sum index does not escape the BID
  family's structural redundancy ceiling.
- This is a *constructive* observation — the family is mathematically
  rich (the bounds, extremal graphs, and special-case identities of
  §3-§5 are non-trivial), but it is *complementary* to the existing
  zoo, not orthogonal to it.
- A formal screening protocol (correlation-based + PCA + partial
  correlation across multiple datasets) would help the field
  distinguish parametric extensions that add complementary
  mathematical structure from those that add genuine descriptor
  information.

## 8. Conclusion

`[NEEDS WRITING ~0.5 page]` Restate contributions: definition,
closed-forms, bounds (§3–4), extremal graphs (§5), honest
redundancy analysis (§6).

## 9. Data and code availability

> All code, results, and a reproducible end-to-end pipeline are
> publicly available at
> <https://github.com/Singati2/topological-index-orthogonality>.

## Acknowledgements

`[TO FILL — funding, institutional support if applicable]`

## References

`[FILL FROM docs/literature_notes.md — every REFERENCE NEEDED slot
must be resolved before submission. Do NOT submit with placeholders.]`

---

## Figures and tables (proposed)

| Figure / Table         | Source                                                  | Status |
|---|---|---|
| Fig. 1 — BID lattice diagram (Δ ≤ 4)              | `[TO CREATE — TikZ]`                 | pending |
| Fig. 2 — PCA scree plot, three datasets            | `results/<dataset>/pca_variance.csv` | `figures/fig2_pca_scree.{png,pdf}` [VERIFIED AGAINST CSV] |
| Fig. 3 — Cross-dataset redundancy bars             | `results/cross_dataset_summary.csv`  | `figures/fig3_redundancy_bars.{png,pdf}` [VERIFIED AGAINST CSV] |
| Fig. 4 — Octane correlation heatmap                | `results/octane_descriptors.csv`     | `figures/fig4_octane_heatmap.{png,pdf}` [VERIFIED AGAINST CSV] |
| Fig. 5 — Wuzi parameter sweep heatmaps (γ slices, ESOL) | `results/esol/wuzi_grid.csv`    | `figures/fig5_wuzi_param_heatmaps.{png,pdf}` [VERIFIED AGAINST CSV] |
| Fig. 5b — Wuzi parameter sweep across three datasets | `results/<dataset>/wuzi_grid.csv`  | `figures/fig5b_wuzi_param_heatmaps_all.{png,pdf}` [VERIFIED AGAINST CSV] |
| Fig. 6 — Degeneracy bar chart (trees of order 10)  | `results/wuzi_degeneracy.csv`        | `figures/fig6_degeneracy_bars.{png,pdf}` [VERIFIED AGAINST CSV] |
| Fig. 7 — Structure sensitivity bar chart (decanes) | `results/structure_sensitivity.csv`  | `figures/fig7_structure_sensitivity.{png,pdf}` [VERIFIED AGAINST CSV] |
| Fig. 8 — Full 30-index correlation matrix (ESOL)   | `results/esol/correlation_matrix.csv`| `figures/fig8_correlation_heatmap.{png,pdf}` [VERIFIED AGAINST CSV] |
| Table 1 — Special-case identities table            | `docs/edge_contribution_analysis.md` | done |
| Table 2 — Closed-form values on standard graphs    | `src/wuzi_analytical.py`             | done |
| Table 3 — Bounds in (n, m, Δ, δ) summary           | [NEEDS DERIVATION BY ADVIK / AROCKIARAJ SIR] | pending |
| Table 4 — Octane prediction R values               | `results/octane_prediction.md`       | done |

Figure captions are kept manuscript-ready in
`docs/figure_captions.md`; the figure ↔ data consistency record is
in `docs/figure_audit.md`.

---

## Notes to collaborators

- **Math style:** match MATCH 95:141-162 conventions exactly (theorem
  format, equality-condition discussion at end of each statement,
  proof environment with ∎ at end). Don't reinvent the wheel.
- **Variable naming:** use d_u, d_v not d(u), d(v) — Movahedi/Gutman
  use the subscript form.
- **No "we propose" without checking:** any phrase like "we introduce
  the Wuzi family" must be paired with a sentence citing the
  parametric-index tradition (general Randić, general Sombor) so we
  do not look isolated.
- **Section 6 framing:** constructive, not critical. We extend the
  intercorrelation analysis of Movahedi et al. 2026 §5.2 to three
  datasets with a formal criterion. We are not attacking that paper.

When math sections 3, 4, 5 are drafted, drop them into the
correspondingly-numbered placeholder blocks above. All citation slots
must be resolved from `docs/literature_notes.md` before submission.
