# Paper 1 — Wuzi family manuscript skeleton

Working template for the parametric-index / math-chem paper. Section
numbering mirrors **Movahedi, Gutman, Redžepović & Furtula, MATCH 95(1),
2026, 141-162** ("Diminished Sombor Index") for direct reviewer
comparability. Math sections marked `[NEEDS DERIVATION ...]` are honest
gaps awaiting Advik + Arockiaraj sir's input.

---

## Title

*"The Wuzi Index Family: Graph-Theoretic Properties, Bounds, Extremal
Graphs, and Sensitivity Analysis"*

(Matches `docs/paper1_wuzi.tex` title verbatim.)

## Authors

Advik Natarajan¹, Ganesh Shiwakoti², Michael Arockiaraj¹
¹ Department of Mathematics, Loyola College, Chennai
² Department of Mathematics and Statistics, Florida Atlantic University

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
> graphs). A rigorous ratio-bound theorem and a Cauchy–Schwarz
> geometric-mean inequality are established for the family in terms
> of any classical BID index with strictly positive edge contribution;
> further sharp bounds in (n, m, Δ, δ) and per-classical-index
> bounds are stated with proof strategies. Extremal-graph behaviour
> among trees, unicyclic, and bicyclic graphs is conjectured by
> parameter-sign region; a non-classical extremal regime at
> (α, β, γ) = (−1, −1, 1) is documented computationally and stated
> as Conjecture 6.4. A sensitivity-analysis section reports octane
> property correlations, degeneracy on the 106 trees of order 10,
> and structure sensitivity on the 75 decane isomers. The
> open-source code and numerical results are publicly available.
> (Cross-dataset orthogonality screening — ESOL, FreeSolv,
> Lipophilicity, BBBP — is the subject of the companion methodology
> paper, not this paper.)

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

Companion strategy document: `docs/wuzi_bounds_strategy.md` (sections A, C).

**Rigorous results available to drop in:**

- **Theorem 3.1 (Bracket bound, rigorous).** m · ψ_min(δ, Δ) ≤ W ≤
  m · ψ_max(δ, Δ), with equality iff every edge of G has the same
  unordered degree pair. Proof: trivial sum bracket. Already
  stated as Theorem A.1 in `docs/wuzi_bounds_strategy.md` and as
  Theorem~\ref{thm:bracket} in `docs/paper1_wuzi.tex`.
- **Theorem 3.2 (Jensen on $M_1$, rigorous for $\beta \ge 1$).**
  For $\alpha = \gamma = 0$, $\beta \ge 1$: $W \ge m^{1-\beta} M_1^\beta$.
  Proof in `docs/paper1_wuzi.tex` Theorem~\ref{thm:m1_bound}.
- **Theorem 3.3 (Geometric mean inequality, rigorous).** For
  $\alpha_1+\alpha_2 = 2\alpha$, etc., $W(G; \alpha, \beta, \gamma)
  \le \sqrt{W(G; \alpha_1, \beta_1, \gamma_1) \cdot W(G; \alpha_2, \beta_2, \gamma_2)}$.
  Proof: Cauchy–Schwarz, in `docs/paper1_wuzi.tex` Theorem~\ref{thm:cs}.

**Now complete (round-10 update):**

- **Proposition 4.7 / 4.8 (upper bound in $(\Delta, m)$ / lower
  bound in $(\delta, m)$ in the non-negative parameter region,
  rigorous).** Termwise factor-wise bounds, summed; equality
  conditions stated. In `docs/paper1_wuzi.tex` as
  `thm:upper_Delta` and `thm:lower_delta`.

**Still needed (NEEDS DERIVATION BY ADVIK / AROCKIARAJ SIR):**

- **Conjecture 4.9 (Nordhaus–Gaddum type).** $W(G) + W(\bar G) \ge f_n(\alpha,\beta,\gamma)$ with an explicit $f_n$; only the
  classical $M_1, M_2$ case has a known template.
- **Sharper $\Delta$-dependent bounds in the mixed-sign region**
  (the current Proposition 4.7 covers only $\alpha, \beta, \gamma \ge 0$).
- **Polya–Szegő-type inverse-Cauchy-Schwarz bounds** using
  $\delta, \Delta$ as Polya–Szegő input.

## 4. Bounds in terms of classical topological indices

Companion strategy document: `docs/wuzi_bounds_strategy.md` (section B).

**Verified ratio-bound theorem (rigorous).**

For every classical BID index $I_h$ with strictly positive edge
contribution $h(i, j)$ on the admissible degree pairs of $G$, the
following holds:

> **Theorem (Ratio bound, rigorous).** $c_{\min}^h \cdot I_h(G) \le
> W(G; \alpha, \beta, \gamma) \le c_{\max}^h \cdot I_h(G)$, where
> $c_{\min}^h, c_{\max}^h$ are the extrema of $\psi(i, j)/h(i, j)$
> over the admissible degree pairs.

This is rigorously stated and proved as Theorem B.1 in
`docs/wuzi_bounds_strategy.md`. The constants $c_{\min}^h$ and
$c_{\max}^h$ are computed by `scripts/09_wuzi_bounds_tables.py`
for $h \in \{M_1, M_2, R, SO, GA, H, ABC\}$ and a list of parameter
triples; output in `results/wuzi_bounds_ratio_tables.{csv,md}`.

The ratio bound is rigorous and now included in Paper 1 (§5,
Theorem~5.2 `thm:ratio_bound`): it covers all seven classical indices in one
unified statement (one theorem + one table per index), and at
exact Wuzi-classical reductions the constants collapse to a single
identity (sanity-checked in the output: e.g.\ $(\alpha,\beta,\gamma) = (1,0,0)$
gives $c_{\min}^{M_2} = c_{\max}^{M_2} = 1$).

**Still needed (NEEDS DERIVATION BY ADVIK / AROCKIARAJ SIR):**

- **Theorem 4.x (sharper than ratio bound).** Bounds that interpolate
  over the distribution of $m_{ij}$ counts rather than using a
  single $c_{\min}, c_{\max}$ pair. These would be tight on a wider
  family of graphs than just the regular ones.
- **Theorem 4.y (Albertson irregularity Alb pairing).** The
  $\gamma$-axis of Wuzi pairs naturally with the irregularity index
  $\mathrm{Alb}(G) = \sum_e |d_u - d_v|$; a clean bound here is the
  natural opening for the $\gamma$ component of the family.

## 5. Extremal graphs

Companion strategy document: `docs/wuzi_bounds_strategy.md` (section D).
Computational support: `scripts/10_wuzi_extremal_search.py`
exhaustively enumerates trees of order $5 \le n \le 12$ and tests
the path-vs-star conjecture at each parameter triple. Output:
`results/wuzi_extremal_trees.{csv,md}`.

**Empirically verified across $5 \le n \le 12$ (Conjecture D.1).**

For trees of order $n \ge 3$ in the sign region
$\alpha \ge 0, \beta \ge 0, \gamma \ge 0$:

> **Conjecture 5.1 (Extremal trees, positive parameter region).**
> The star $S_n$ uniquely attains $\max W$; the path $P_n$ uniquely
> attains $\min W$.

Empirical verification at $(\alpha, \beta, \gamma) \in \{(1,0,0),
(0,0,1), (0,0,2), (1,1,1)\}$ over all $5 \le n \le 12$ shows
$P_n$ as observed minimizer and $S_n$ as observed maximizer in
$100\%$ of cases. See `results/wuzi_extremal_trees.md` for the full
table. This is empirical evidence, *not* proof; the formal proof
requires an edge-transformation argument.

**Empirically observed sign-flip behaviour.**

For the Randić-region triples $(-0.5, 0, 0)$, $(0, -0.5, 0)$,
$(0, -1, 0)$, the path $P_n$ becomes the observed *maximizer* and
the star $S_n$ the observed *minimizer*. This is consistent with
the classical Bollobás–Erdős extremal-trees result for general
Randić $R_\alpha$ with $\alpha < 0$.

**Mixed-sign region $(-1, -1, 1)$.**

In this region (the FreeSolv "borderline pass" triple) the
observed extremals are *neither* $P_n$ nor $S_n$ for $n \ge 5$.
This is the most interesting case for genuine new mathematics in
Paper 1 §6: the extremal structures here are not the classical
ones and require characterization. The script output identifies
specific extremal trees by structure for each $n$.

**Still needed (NEEDS DERIVATION BY ADVIK / AROCKIARAJ SIR):**

- **Theorem 5.1 (formal proof of Conjecture 5.1).** Edge-shifting
  argument moving a pendant from a low-degree vertex to a
  high-degree vertex; equality conditions.
- **Theorem 5.2 (extremal unicyclic graphs).** Compare $C_n$ to
  hub-concentrated unicyclic graphs in each sign region. Enumeration
  for unicyclic / bicyclic is *not* done by script 10 (networkx has
  no built-in non-isomorphic unicyclic enumerator); a SageMath
  or nauty implementation is left as a follow-up.
- **Theorem 5.3 (extremal bicyclic graphs).** Three sub-classes
  (θ-graph, bowtie-like, two-cycles-joined-by-path) and their
  pendant-augmented variants must be considered separately.
- **Characterization in the mixed-sign region.** What is the
  extremal tree of order $n$ for $(\alpha, \beta, \gamma) = (-1, -1, 1)$?
  The script identifies it computationally; the formal
  characterization is open.

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
closed-forms, BID-basis observation (§3), bounds (§4–5),
extremal graphs (§6), sensitivity analysis (§7).

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
| Fig. 5b — Wuzi parameter sweep across three datasets, 4 γ slices | `results/<dataset>/wuzi_grid.csv` | `figures/fig5b_wuzi_param_heatmaps_all.{png,pdf}` [VERIFIED AGAINST CSV] |
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
