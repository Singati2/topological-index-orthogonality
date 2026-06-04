# Paper 1 — math core audit (closed-forms & bounds)

Audited `docs/paper1_wuzi_match.tex` §3 (closed-forms) and §4 (bounds), plus the
ceiling proposition in §6. **Verdict: all 12 formal results are mathematically
correct. No errors found.** (Audited the repo copy; the canonical Loyola manuscript
is on Overleaf — the math should be identical under the Wuzi→LO rename.)

## Results checked (each verified line by line)
| # | Result | Status |
|---|---|---|
| Prop 3.x | Closed forms: K_n, C_n, P_n, S_n, K_{p,q}, Q_k, wheel, friendship, r-regular | ✅ all 9 correct |
| Prop (bracket) | m·ψ_min ≤ W ≤ m·ψ_max + equality conditions | ✅ correct |
| Thm (CS split) | W(α,β,γ) ≤ √(W₁·W₂) for additive parameter splits | ✅ correct (imbalance term is parameter-independent, exponents add) |
| Thm (M₁, β≥1) | W ≥ m^{1−β}M₁^β via convex Jensen | ✅ correct |
| Thm (M₁, 0<β<1) | W ≤ m^{1−β}M₁^β via concave Jensen | ✅ correct |
| Thm (M₂, α≥1) | R_α ≥ m^{1−α}M₂^α | ✅ correct |
| Prop (Δ upper) | W ≤ m·Δ^{2α}(2Δ)^β·exp(γ(Δ−1)/(Δ+1)) | ✅ valid (see note 1) |
| Prop (δ lower) | W ≥ m·δ^{2α}(2δ)^β, equality iff δ-regular | ✅ correct |
| Thm (ratio bound) | c_min·I_h ≤ W ≤ c_max·I_h | ✅ correct |
| Prop (R, χ) | W(α,β,0) ≤ √(R_{2α}·χ_{2β}) | ✅ correct |
| Prop (Sombor) | χ_β ≤ m^{1−β}(√2·SO)^β via M₁≤√2·SO | ✅ correct (x²+y²≥(x+y)²/2) |
| Prop (γ-axis) | γ→±∞ asymptotics to extremal-imbalance edge sets | ✅ correct |
| Prop (ceiling) | Δ≤4 ⇒ 10 degree-pair classes ⇒ span dim ≤ min(10,rank) | ✅ correct |

## Minor polish notes (NOT errors — optional tightening for reviewers)
1. **Δ-upper bound is valid but not sharp.** It multiplies three *separately*
   maximised factors whose maxima conflict: d_ud_v=Δ² and d_u+d_v=2Δ require the
   pair (Δ,Δ) (imbalance 0), while the imbalance max (Δ−1)/(Δ+1) needs (1,Δ).
   So for γ>0 the bound is never attained (e.g. on an r-regular graph it overshoots
   the true value by exp(γ(r−1)/(r+1))). Suggest one sentence acknowledging it is a
   non-sharp upper bound, or stating no equality case (unlike the δ-lower bound,
   which correctly gives equality iff δ-regular).
2. **Bracket bound admissible set.** P(δ,Δ)={(i,j):δ≤i≤j≤Δ} reads as a real box;
   degrees are integers, so specifying (i,j)∈ℤ² makes ψ_min/ψ_max well-defined and
   the equality conditions exactly meaningful. Cosmetic.
3. **Star notation.** S_n is used with n vertices / (n−1) leaves; confirm this matches
   the intro's S_n convention so the reduction table is consistent. Cosmetic.

## Bottom line
The mathematical contribution — the part a MATCH reviewer weighs most — is sound.
The only substantive caveats in the paper live in the (reproducibility-limited)
sensitivity section, now reframed honestly; the theorems themselves are solid.
