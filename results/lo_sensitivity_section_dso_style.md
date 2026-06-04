# Sensitivity analysis (DSO-style) — Loyola index — honest rewrite

Modeled on Movahedi–Gutman–Redžepović–Furtula, *MATCH* 95 (2026) 141–162 (DSO).
Compared indices (fixed order): the reductions **M1, M2, HM, ᵐM2, R, χ, H/2, ISI,
GA, AG** plus the two **pre-specified** pure-γ points **LO(0,0,1), LO(0,0,2)**.
Edge count m = LO(0,0,0) omitted (constant on octanes). No decimal-tuned triples.

## 5.1 Prediction ability (18 octanes, NIST properties)
Best |r| per property:

| property | best index | r |
|---|---|---|
| T_B | ᵐM₂ | +0.861 |
| ΔH_f | ᵐM₂ | +0.888 |
| ΔH_vap | **LO(0,0,2)** | −0.967 |
| S | HM | −0.975 |
| ω | M₂ | −0.987 |

→ Loyola is competitive and, at a pre-specified setting (LO(0,0,2)), improves the
correlation on ΔH_vap. (NIST source; DSO used a different source — ΔH_vap/S may
differ by source, not by index computation; see property-source note.)
Table: `lo_sensitivity_prediction_table_latex.tex`; figure:
`fig_lo_prediction_correlation_heatmap.{png,pdf}`.

## 5.2 Intercorrelation (18 octanes)
Max |r| of each pure-γ point with a classical reduction: **LO(0,0,1) = 0.983,
LO(0,0,2) = 0.993** → the pure-γ points are highly collinear with the classical
reductions, so the **γ-imbalance term does not add independent linear information
on this small, degree-limited benchmark.** Figure:
`fig_lo_intercorrelation_heatmap.{png,pdf}`.

## 5.3 Degeneracy (106 order-10 trees)
Degeneracy(%) = 100(1 − #distinct/N). LO(0,0,1)=LO(0,0,2)=**25.5%** (among the
lowest = best discrimination, ties R/GA/AG); M1 = 83.0%, M2 = 67.0% (most
degenerate). Formula to be verified against DSO's source [8]. Figure:
`fig_lo_degeneracy_order10_trees.{png,pdf}`.

## 5.4 Structure sensitivity (75 decanes)
SS = σ/μ, Abr = (max−min)/μ, SA = SS/Abr. **LO(0,0,2) has the highest SS (0.165)
and Abr** among the compared indices; LO(0,0,1) also high → the γ term's
contribution shows up in structure sensitivity. ⚠ Standard formulas; they do **not**
reproduce DSO Table 2 (our SO = 0.098/0.446/0.219 vs DSO 0.193/0.396/0.488) — to be
reconciled with Rakić–Furtula [12]. Figure:
`fig_lo_structure_sensitivity_decanes.{png,pdf}`.

## Safe summary sentence (manuscript-ready)
> "The Loyola family unifies several classical degree-based indices as exact
> reductions. In the octane benchmark, selected pre-specified Loyola values were
> competitive with classical reductions (best on ΔH_vap), while the pure-γ values
> showed high structure sensitivity but limited independent correlation information
> on this small, degree-limited set."

## Claims to AVOID
- "Loyola dominates / generally outperforms the classical indices."
- "Small decimal tweaks improve prediction."
- "We tuned Loyola to beat the best index."
- "The DSO property table was recovered" / "aligned to DSO" / "SO row matches DSO."
- "The γ term adds independent information on octanes" (the data say the opposite).
- "Following DSO exactly" for structure sensitivity (formulas not verified).
