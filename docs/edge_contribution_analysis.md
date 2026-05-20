# Edge-contribution function analysis (Paper 1 §2)

This document is the Wuzi analog of the γ(x, y) analysis in
**Movahedi, Gutman, Redžepović, Furtula 2026** (MATCH 95:141-162),
Eqs. (3)-(4), pp. 143-144. It characterizes the monotonicity and
extremal behavior of the Wuzi *edge contribution function*

> `ψ(x, y; α, β, γ) := (x · y)^α · (x + y)^β · exp(γ · |x − y| / (x + y))`,

so that bounds in Sections 3 and 4 can be derived by direct substitution.

Throughout, `x, y ∈ ℝ_{>0}` are the degrees of adjacent vertices. For
molecular graphs under Assumption (A) (`Δ ≤ 4`), we have `(x, y) ∈
{1, 2, 3, 4}²`. By symmetry of ψ in its first two arguments, we may
restrict attention to `1 ≤ x ≤ y ≤ 4`.

---

## 1. Decomposition

Write `ψ = ψ_P · ψ_S · ψ_I` as a product of three factors:

| Factor             | Formula                              | Encodes              |
|---|---|---|
| Product factor `ψ_P` | `(xy)^α`                             | Randić-family axis   |
| Sum factor `ψ_S`     | `(x+y)^β`                            | Sum-connectivity axis|
| Irregularity factor `ψ_I` | `exp(γ · |x-y| / (x+y))`         | Albertson/σ axis     |

For the special parameter slices:

| (α, β, γ)       | ψ reduces to              | I_f reduces to         |
|---|---|---|
| (0, 0, 0)        | 1                         | `|E|` (edge count)     |
| (1, 0, 0)        | xy                        | Second Zagreb M₂       |
| (-1/2, 0, 0)     | `1/√(xy)`                 | Randić R               |
| (0, -1/2, 0)     | `1/√(x+y)`                | Sum-connectivity SCI   |
| (0, -1, 0)       | `1/(x+y)`                 | Harmonic / 2           |
| (1/2, 0, 0)      | √(xy)                     | (geometric mean)       |

---

## 2. Closed-form values on the (1, 1)…(4, 4) lattice

For γ = 0 (no irregularity weighting), `ψ` is symmetric and factors as
`ψ(x, y; α, β, 0) = (xy)^α · (x+y)^β`. The 10 ordered pairs and their
contributions for any fixed (α, β):

| (x, y) | xy | x + y | `\|x − y\|/(x+y)` |
|---|---:|---:|---:|
| (1, 1) | 1  | 2 | 0     |
| (1, 2) | 2  | 3 | 1/3 ≈ 0.333 |
| (1, 3) | 3  | 4 | 2/4 = 0.5   |
| (1, 4) | 4  | 5 | 3/5 = 0.6   |
| (2, 2) | 4  | 4 | 0           |
| (2, 3) | 6  | 5 | 1/5 = 0.2   |
| (2, 4) | 8  | 6 | 2/6 ≈ 0.333 |
| (3, 3) | 9  | 6 | 0           |
| (3, 4) | 12 | 7 | 1/7 ≈ 0.143 |
| (4, 4) | 16 | 8 | 0           |

The right-most column is the input to the exponential `exp(γ · ·)`.
The maximum irregularity ratio on the BID lattice is `3/5 = 0.6`,
attained at the (1, 4) pair (pendant attached to a quaternary atom).

---

## 3. Partial derivatives

Let `f(x, y) := ln ψ(x, y; α, β, γ) = α ln(xy) + β ln(x+y) + γ · |x-y|/(x+y)`.
For `x ≠ y`, taking partial derivatives:

> `∂f/∂x = α/x + β/(x+y) + γ · sign(x-y) · 2y / (x+y)²`

> `∂f/∂y = α/y + β/(x+y) − γ · sign(x-y) · 2x / (x+y)²`

For `y > x` (sign(x-y) = -1):
> `∂f/∂x = α/x + β/(x+y) − 2γy / (x+y)²`

> `∂f/∂y = α/y + β/(x+y) + 2γx / (x+y)²`

**Consequence (γ ≥ 0):** For `y > x`, `∂ψ/∂y > 0` whenever
`α/y + β/(x+y) + 2γx/(x+y)² > 0` — in particular, when α, β, γ are
all non-negative. So `ψ` is monotone non-decreasing in the larger
argument when all three parameters are non-negative.

---

## 4. Extremal pairs on `[1, 4]²` for canonical parameter regions

By inspection of the closed-form values above, ψ attains its
extremes on the BID lattice at the following pairs in each canonical
parameter region:

| Region                  | min over (x, y) ∈ {1..4}²  | max over (x, y) ∈ {1..4}² |
|---|---|---|
| α = β = γ = 0            | (any) → 1                  | (any) → 1                 |
| α > 0, β = γ = 0         | (1, 1) → 1                 | (4, 4) → 16^α             |
| α < 0, β = γ = 0         | (4, 4) → 16^α              | (1, 1) → 1                |
| α = 0, β > 0, γ = 0      | (1, 1) → 2^β               | (4, 4) → 8^β              |
| α = 0, β < 0, γ = 0      | (4, 4) → 8^β               | (1, 1) → 2^β              |
| α = β = 0, γ > 0         | regular (γ = 0)            | (1, 4) → exp(0.6 γ)       |

This table is what enables the *m · ψ_min(Δ, δ) ≤ W ≤ m · ψ_max(Δ, δ)*
template of bounds in Section 3.

**Lemma 2.3.1.** For any (α, β, γ) ∈ ℝ³ and any G with Δ(G) ≤ 4,

> `m · min_{(x,y) ∈ T} ψ(x, y; α, β, γ) ≤ W(G; α, β, γ) ≤ m · max_{(x,y) ∈ T} ψ(x, y; α, β, γ)`

where `T = {(x,y) : 1 ≤ x ≤ y ≤ 4, m_xy(G) > 0}` is the set of edge-type
pairs actually realized in G. Equality on either side holds if and only
if all edges of G have endpoint degrees attaining that extremum.

**Proof.** Direct from `W = Σ_e ψ(d_u, d_v)` and the fact that the
average of values is between their min and max. ∎

**Corollary 2.3.2 (general Δ, δ bounds).** For any connected G with
maximum degree Δ and minimum degree δ,

> `m · ψ_min(δ, Δ; α, β, γ) ≤ W(G; α, β, γ) ≤ m · ψ_max(δ, Δ; α, β, γ)`,

where ψ_min and ψ_max are taken over the rectangle `[δ, Δ] × [δ, Δ]`.
The Δ ≤ 4 case reduces to the discrete maximum above.

---

## 5. Regular-graph value (sharp lower bound for many parameter regions)

For k-regular G with `n` vertices and `m = nk/2`:

> `W(G; α, β, γ) = m · k^{2α + β} · 2^β = (nk/2) · k^{2α + β} · 2^β`

This is the special case where every edge has `|x − y| = 0`, so the
irregularity exponential is `exp(0) = 1` and the parameter γ drops
out. *On regular graphs, the Wuzi family is degenerate in γ.*

**Implication for §6 (extremal graphs):** the regular-graph
"all γ collapse" degeneracy means that any γ ≠ 0 sensitivity has to
manifest on irregular graphs — pendants on trees, hub-rim distinction
on wheels, etc. The γ axis is what distinguishes Wuzi from a pure
(α, β)-parametrized member of the general Randić / Sombor family.
