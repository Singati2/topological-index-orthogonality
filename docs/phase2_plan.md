# Phase 2 Plan — Wuzi Family as a Parametric Topological Index

> **SUPERSEDED — kept for history only.** This document describes the
> pre-restructure single-paper layout (Wuzi math + screening + ML
> benchmark in one 10-section paper). Arockiaraj sir's directive split
> the project into two non-overlapping manuscripts: Paper 1 (Wuzi
> mathematical chemistry, `docs/paper1_wuzi.tex`) and Paper 2
> (orthogonality-screening methodology, `docs/paper2_orthogonality_screening.tex`).
> The authoritative manuscript plan is now `docs/paper1_wuzi_manuscript_skeleton.md`,
> `docs/advisor_ready_summary.md`, and `docs/final_10_of_10_upgrade_plan.md`.
> Do not use this file's section layout as a guide — it no longer matches the project.

Structural template: **Movahedi, Gutman, Redžepović, Furtula (2026),
"Diminished Sombor Index,"** MATCH Commun. Math. Comput. Chem. 95(1),
141–162. We mirror that paper's section structure for the Wuzi family
because (a) the math-chem community accepts that template, (b) it lets
reviewers verify completeness section-by-section, and (c) the DSO origin
paper's Section 5.2 intercorrelation observation provides the
load-bearing motivation for our orthogonality-screening framework.

The companion **Movahedi (2025), "Diminished Sombor index and its
relationship with topological indices"** (arXiv) is the template for our
Section 4 (bounds in terms of classical indices) — it shows the genre
of `DSO ≤ f(Alb, GA, H, R, ...)` inequalities.

---

## Status legend

| Marker | Meaning |
|---|---|
| ✅ DONE | Already implemented and verified in this repo |
| 🔵 PARTIAL | Some content exists, needs extension |
| 🟡 MECHANICAL | I can derive this; mostly algebra/computation |
| 🔴 NEEDS MATH | Real graph-theoretic work; best done by Advik/Arockiaraj sir |
| 🟣 EMPIRICAL | Already covered in our Phase-1 data |

---

## Proposed section structure

### Section 1 — Introduction (3–4 paragraphs)

Goal: motivate the Wuzi family as a parametric extension of the Sombor /
Diminished-Sombor / Randić / Zagreb tradition that absorbs all three
endpoint-degree-based design axes (product, sum, irregularity) into a
single tunable formula.

- 🔴 NEEDS MATH/WRITING — opening framing, citation walk-through.
  Cite [Gutman 2021 Sombor], [Movahedi/Gutman/Redžepović/Furtula 2026
  DSO], [Movahedi 2025 DSO bounds], [Randić 1975], [Gutman & Trinajstić
  1972 Zagreb], [Albertson 1997].

### Section 2 — Preparations

#### 2.1 — Notation and the Wuzi definition

- ✅ DONE — `src/wuzi_index.py` and our README give the formal definition.

#### 2.2 — Closed-form Wuzi values on simple graphs

(Direct analog of Proposition 1 in MATCH 95:141–162.)

- ✅ DONE — `src/wuzi_analytical.py` provides closed-form expressions
  with full numerical verification for:

  - W(K_n) and general k-regular graphs of order n
  - W(C_n) — cycle
  - W(P_n) — path
  - W(S_n) — star
  - W(K_{p,q}) — complete bipartite
  - W(Q_k) — hypercube
  - W(W_n) — wheel
  - W(F_p) — friendship graph

  All identities are verified against the direct `wuzi(G, α, β, γ)`
  computation on a NetworkX instance for six representative
  parameter triples.

- 🟡 MECHANICAL — add Dutch windmill graph D_p^(q), Petersen-type small
  graphs, dendrimer-like specific structures if needed for the
  manuscript. Mechanical — ~30 min each.

#### 2.3 — Analysis of the edge-contribution function

The MATCH DSO paper studies γ(x, y) = √(x²+y²)/(x+y) and its
monotonicity. Our analog is the Wuzi edge-contribution function

> `ψ(x, y; α, β, γ) := (xy)^α · (x+y)^β · exp(γ · |x-y| / (x+y))`

- 🟡 MECHANICAL — derive partial derivatives ∂ψ/∂x, ∂ψ/∂y; identify
  parameter regions where ψ is monotone in each argument; locate the
  maximum and minimum on the degree-pair lattice {1, ..., n-1}². This
  is the closest formal analog of Proposition 1 / Eqs. (3)–(4) in
  MATCH 95:141–162.

#### 2.4 — Inequalities used in the proofs

- 🟡 MECHANICAL — list lemmas (Cauchy–Schwarz, Radon's inequality,
  Polya–Szegő, Jensen on convex functions, Maclaurin sums of products,
  power-mean inequalities). Same lemma set used in Movahedi 2025.

### Section 3 — Bounds in terms of (n, m, Δ, δ)

Analog of Section 3 in MATCH DSO and Movahedi 2025.

- 🔴 NEEDS MATH — typical statements:
  - Lower bound: `W(G; α, β, γ) ≥ m · g_min(Δ, δ; α, β, γ)`
    with equality on regular graphs.
  - Upper bound: `W(G; α, β, γ) ≤ m · g_max(Δ, δ; α, β, γ)`
    with equality on regular graphs.
  - Refinements when sign of α, β, γ is restricted.
  - Nordhaus–Gaddum-type relations W(G) + W(Ḡ) ≥ ... for the complement.

  Estimated 5–10 theorems with proofs, ~3–6 pages.

### Section 4 — Bounds via classical indices

Analog of Movahedi 2025 (the DSO bounds-paper).

- 🔴 NEEDS MATH — Wuzi bounds in terms of Zagreb M_1, M_2, Randić R,
  Sombor SO, Albertson Alb, GA, harmonic H, ABC, sum-connectivity SCI.
  For specific (α, β, γ) values, exact identities reduce Wuzi to a
  known index — those become corollaries.

  Estimated 6–8 theorems with proofs, ~4–6 pages.

### Section 5 — Extremal graphs

Analog of Section 4 of MATCH DSO.

- 🔴 NEEDS MATH — characterize among connected graphs of given order
  which graphs maximize/minimize W in each parameter region. Cover:
  general connected graphs, trees, unicyclic, bicyclic. Cite the
  "shift toward pendant" arguments common in degree-based extremal
  results.

  Estimated 5–7 theorems with proofs, ~3–5 pages.

### Section 6 — Numerical work (the orthogonality-screening contribution)

Analog of Section 5 of MATCH DSO, but **substantially extended**. This
is the load-bearing section of our paper because it differentiates us
from a routine new-index proposal.

#### 6.1 — Wuzi prediction ability on physicochemical properties

🟡 MECHANICAL — analog of MATCH Section 5.1.
- Apply the standard octane benchmark (n = 18) used in the MATCH paper:
  measure correlation of Wuzi(best (α,β,γ)) against boiling point (T_B),
  enthalpy of formation (ΔH_f), enthalpy of vaporization (ΔH_vap),
  entropy (S), acentric factor (ω).
- Threshold "predicts well" at r > 0.8 (same convention as DSO paper).
- **Why the small-n octane benchmark is now defensible**: we're not
  *claiming* Wuzi is superior; we're including this set for direct
  comparability with the DSO Section 5.1 numbers.

#### 6.2 — Intercorrelations among Wuzi parameter points

🟣 EMPIRICAL — already done for the full Wuzi 100-point grid on three
benchmark datasets (ESOL n=1127, FreeSolv n=639, Lipophilicity n=4200).
See `results/<dataset>/wuzi_grid.csv` and `results/cross_dataset_summary.md`.

**Key finding to feature:** Even at parameter points where Wuzi is
"least correlated" with the baseline, max |r| stays ≥ 0.94 across all
three datasets. This *directly replicates and extends* the MATCH DSO
Section 5.2 observation (DSO vs other Sombor variants: r ∈ [0.92, 0.97])
on a 50× larger benchmark and with proper partial-correlation analysis.

#### 6.3 — Effective rank / PCA

🟣 EMPIRICAL — done. Cross-dataset PCA finds effective rank 3 for the
30-index baseline on all three datasets. This is our novel contribution
beyond what the DSO paper reports.

#### 6.4 — Degeneracy

🟡 MECHANICAL — analog of MATCH 5.3. Compute % of distinct Wuzi values
among trees of size 10 (or similar small classes) for representative
parameter triples. Easy to add (~30 min).

#### 6.5 — Structure sensitivity (SS, Abr, SA)

🟡 MECHANICAL — analog of MATCH 5.4. Use the [12] decanes set
(reference in the DSO paper). Easy to add (~1 hour).

### Section 7 — The orthogonality screening framework (the methodology contribution)

🟣 EMPIRICAL — `scripts/01_baseline_correlations.py`,
`scripts/02_wuzi_grid_search.py`, `scripts/03_cross_dataset_summary.py`
are all written, tested, and replicated across three benchmark datasets.

Required additions:
- 🟡 MECHANICAL — `docs/usage_guide.md`: walk a hypothetical user through
  "I want to propose index X, how do I run the screening?" — concrete
  CLI recipe, interpretation of outputs, common pitfalls.
- 🟡 MECHANICAL — minor refactor exposing `screen_new_index(G_list, f)`
  as a one-call public API.

### Section 8 — Discussion

- 🔴 NEEDS WRITING — what redundancy means; why formula-tweaking inside
  the BID family can't escape the 10-D ceiling; why orthogonality should
  precede QSPR claims; reframing the DSO paper's Section 5.2 observation
  formally.

### Section 9 — Limitations

- 🟡 MECHANICAL — drug-like chemistry only; bond order ignored; scalar
  indices; orthogonality is necessary but not sufficient; etc.

### Section 10 — Conclusion

- 🔴 NEEDS WRITING — short. Restate contribution.

---

## Realistic effort estimate

| Section | Effort | Owner |
|---|---|---|
| 1 (intro) + 2 (preparations) | ~2 days | Advik + me (writing) |
| 3 (n,m,Δ,δ bounds) | ~1 week | Advik + Arockiaraj sir (math) |
| 4 (classical-index bounds) | ~1 week | Advik + Arockiaraj sir (math) |
| 5 (extremal graphs) | ~1 week | Advik + Arockiaraj sir (math) |
| 6 (numerical) — additions only | ~1 day | me |
| 7 (screening framework) | ~1 day | me |
| 8–10 (discussion, limitations, conclusion) | ~3 days | all of us |

**Total wall-clock if everyone is reasonably available:** 5–7 weeks.
That hits the "8–10 weeks to JCIM-tier" estimate I gave earlier.

---

## Order of operations

1. **First (this week, by me):** finish Section 2 (preparations) — done
   for special-graph values, todo: 2.3 (edge-contribution analysis,
   mechanical) and 2.4 (inequality list, mechanical). Also Section 6.4
   (degeneracy) and 6.5 (structure sensitivity).
2. **Next (Advik + sir):** Sections 3, 4, 5. The math-chem heavy
   lifting. These three sections are independent so they can be
   parallelized among collaborators.
3. **Then (collaborative):** Section 6.1 (octane prediction), Section 7
   (usage guide), Section 8 (discussion writeup), Sections 9–10.

I'll prepare the LaTeX skeleton with all section headings and the
"already-derived" results pre-filled so the math sections can be
written into ready-made slots without setup friction.
