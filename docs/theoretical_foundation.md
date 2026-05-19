# Theoretical Foundation: The Edge-Degree-Pair Basis and Redundancy

This document establishes the mathematical reason that many proposed
topological indices end up statistically redundant on real chemistry, and
why the Wuzi family in particular cannot escape that redundancy at any
parameter setting.

The core observation is **not new**: the **Bond-Incident-Degree (BID)** index
family has been studied for decades by Gutman, Borovićanin, Furtula and
collaborators (see [REFERENCE NEEDED: Borovićanin, Das, Furtula, Gutman,
"Bounds for Zagreb indices," MATCH Commun. Math. Comput. Chem. 78, 2017]
and [REFERENCE NEEDED: Gutman & Tošović, "Testing the quality of molecular
structure descriptors," J. Serb. Chem. Soc. 78, 2013] for the canonical
edge-partition framing). Gutman's own survey [I. Gutman, "Degree-based
topological indices," *Croatica Chem. Acta* 86(4), 2013, 351-361] catalogs
the family.

The redundancy implication is foreshadowed by the explicit intercorrelation
observation in [F. Movahedi, I. Gutman, I. Redžepović, B. Furtula,
"Diminished Sombor Index," *MATCH Commun. Math. Comput. Chem.* 95(1),
2026, 141-162], Section 5.2, where the authors report correlations
0.92–0.99 between the Diminished Sombor index and other Sombor variants
on octanes (n = 18) and infer informally that the lower correlation
"indicates that the diminished Sombor index captures some additional
information." Our pipeline replicates and formalizes this kind of
analysis on three independent QSPR benchmark datasets (ESOL n = 1127,
FreeSolv n = 639, Lipophilicity n = 4200) with a defined |r| ≥ 0.95
threshold, partial-correlation diagnostics, and PCA-based effective rank.

What this document does is therefore: **state the dimension bound
explicitly** in the form of Theorem 1 below, and **draw the redundancy
consequence** that motivates the orthogonality-screening pipeline.

---

## 1. Setup: hydrogen-suppressed molecular graphs

Let `G = (V, E)` be a connected, simple, undirected hydrogen-suppressed
molecular graph: vertices are heavy atoms and edges are covalent bonds
(bond order ignored). Under standard organic valence rules, every heavy
atom has at most four heavy-atom neighbors, so

> **Assumption (A).** `Δ(G) ≤ 4` for every graph `G` in the family
> considered. (This is essentially all of organic chemistry; the bound
> is reached only by quaternary carbons and similar.)

For an edge `uv ∈ E(G)`, the *endpoint-degree pair* is the **unordered**
pair `{d_u, d_v}` where `d_v := deg_G(v)`. Under Assumption (A), the
endpoint-degree pair belongs to the 10-element set

> `P := { (i, j) : 1 ≤ i ≤ j ≤ 4 } = {(1,1), (1,2), (1,3), (1,4),
>          (2,2), (2,3), (2,4), (3,3), (3,4), (4,4)}`.

For each `(i, j) ∈ P`, define the **edge-degree-pair count**

> `m_ij(G) := | { uv ∈ E(G) : {d_u, d_v} = {i, j} } |`.

The `m_ij` are integer-valued graph invariants whose sum equals `|E(G)|`.

---

## 2. The basis theorem

**Definition (endpoint-degree edge-sum index).** Let `f : {1,2,3,4}² → ℝ`
be a function that is **symmetric** in its arguments
(`f(i,j) = f(j,i)`). The *endpoint-degree edge-sum index associated to `f`*
is

> `I_f(G) := Σ_{uv ∈ E(G)} f(d_u, d_v)`.

A graph invariant is an *endpoint-degree edge-sum index* (equivalently, a
*BID index*) if it can be written in this form for some symmetric `f`.

**Theorem 1 (Edge-Degree-Pair Basis).** *Under Assumption (A), for every
symmetric `f : {1,…,4}² → ℝ`,*

> `I_f(G) = Σ_{(i,j) ∈ P} f(i,j) · m_ij(G)`. (★)

*Consequently, the real vector space*

> `V := { I_f : f symmetric } ⊆ ℝ^{Graphs}`

*spanned by all endpoint-degree edge-sum indices (viewed as real-valued
functions on the family of `Δ ≤ 4` molecular graphs) has dimension at
most 10.*

**Proof.** For each edge `uv ∈ E(G)`, the contribution `f(d_u, d_v)`
depends only on the unordered pair `{d_u, d_v}` because `f` is symmetric.
Grouping edges by their endpoint-degree pair gives

> `Σ_{uv∈E(G)} f(d_u, d_v) = Σ_{(i,j)∈P} f(i,j) · m_ij(G)`,

which is (★). The right-hand side is a linear combination (with
coefficients `f(i,j)` depending only on `f`, not on `G`) of the 10
quantities `{ m_ij : (i,j) ∈ P }`. Hence every `I_f` lies in the span of
these 10 functions, so `dim V ≤ 10`. ∎

---

## 3. Corollaries

**Corollary 2 (linear dependence).** *Any collection of 11 or more
endpoint-degree edge-sum indices is linearly dependent on the family of
`Δ ≤ 4` molecular graphs. Equivalently, given any 10 linearly independent
such indices, every other endpoint-degree edge-sum index is an exact
linear combination of those 10.*

**Corollary 3 (the Wuzi family is in `V`).** *The Wuzi parametric family*

> `W(G; α, β, γ) = Σ_{uv∈E(G)} (d_u d_v)^α · (d_u + d_v)^β
>                            · exp( γ · |d_u − d_v| / (d_u + d_v) )`

*lies in `V` for every choice of `(α, β, γ) ∈ ℝ³`. In particular, no
parameter setting produces a `W(·; α, β, γ)` linearly independent from any
10 linearly independent classical BID indices.*

**Proof.** The summand depends only on `d_u, d_v` and is symmetric in
them. Apply Theorem 1. ∎

**Corollary 4 (effective rank ceiling).** *On any subfamily `F` of `Δ ≤ 4`
molecular graphs, the empirical rank of the descriptor matrix formed from
any set of endpoint-degree edge-sum indices is at most `min(10,
rank(M_F))`, where `M_F` is the matrix of edge-degree-pair counts `m_ij`
over `F`.*

In practice the rank of `M_F` is **strictly less than 10** on drug-like
chemistry because several `m_ij` are sparse or near-constant:

- `m_11`: a degree-1–degree-1 edge requires an isolated bond (e.g., HCN
  with hydrogens suppressed reduces to `(C ≡ N)` — a single edge `m_11`).
  Vanishingly rare in drug-like chemistry.
- `m_14`, `m_44`: require quaternary carbons connected to degree-1
  pendants or to other quaternary atoms — also rare.

We observe this directly: on ESOL (n = 1127), the empirical effective
rank of the 30-index baseline is **≈ 3** (PCA: 95 % of variance in 3
components). Three is consistent with — and bounded above by — the
theoretical limit of 10, with the gap explained by sparse `m_ij`.

---

## 4. Which classical indices are in `V` (and which are not)

**Inside `V`** (endpoint-degree edge-sum form):

| Index           | `f(i, j)`                              |
|-----------------|----------------------------------------|
| `|E|`           | `1`                                    |
| First Zagreb M₁ | `i + j`  (uses identity `M₁ = Σ_{uv∈E}(d_u + d_v)`) |
| Second Zagreb M₂| `i · j`                                |
| Randić R        | `(i · j)^{−1/2}`                       |
| Sum-conn. SCI   | `(i + j)^{−1/2}`                       |
| Harmonic H      | `2 / (i + j)`                          |
| GA              | `2√(i·j) / (i + j)`                    |
| AG              | `(i + j) / (2√(i·j))`                  |
| ABC             | `√((i + j − 2) / (i · j))`             |
| ABS             | `√(1 − 2/(i + j))`                     |
| AZI             | `(i j / (i + j − 2))³` for `i+j>2`     |
| Sombor SO       | `√(i² + j²)`                           |
| Reduced Sombor  | `√((i−1)² + (j−1)²)`                   |
| Forgotten F     | `i² + j²`  (uses `F = Σ_{uv∈E}(d_u² + d_v²)`) |
| Albertson       | `\|i − j\|`                            |
| Sigma σ         | `(i − j)²`                             |
| Wuzi `W(α,β,γ)` | `(ij)^α (i+j)^β exp(γ\|i−j\|/(i+j))`   |

**Outside `V`** (use information beyond endpoint degrees):

- **Distance-based:** Wiener `W`, hyper-Wiener `WW`, Harary, Schultz,
  Gutman, Mostar, Szeged, Padmakar-Ivan (PI) — depend on shortest-path
  structure beyond the endpoints' degrees.
- **Spectral:** Estrada `EE`, graph energy, spectral radius, algebraic
  connectivity (Fiedler value) — depend on the adjacency or Laplacian
  spectrum.
- **Other:** Balaban `J` (uses vertex distance sums), eccentricity-based
  indices (use longest shortest paths), centrality-based indices
  (betweenness, closeness, etc.).

These indices live *outside* the 10-D subspace and can therefore — in
principle — supply information not reducible to edge-degree-pair counts.
This is consistent with their generally lower pairwise correlation with
the BID block in our empirical baseline.

---

## 5. Implications for new-index proposals

1. **Any new BID index is trivially redundant given enough existing BID
   indices on the same graph family.** Proposing a new endpoint-degree
   edge-sum index without checking its linear dependence on existing ones
   is not a meaningful contribution: the dimension is capped at 10.

2. **Parametric BID families (general Randić, general Sombor, Wuzi, …)
   cannot escape the 10-D subspace by adding knobs.** The knobs change
   the coefficients `f(i, j)` but not the basis. The empirical screening
   result we observed for Wuzi — every one of the 100 grid points has
   max |r| with at least one of the 30 classical baselines above 0.95
   on ESOL and Lipophilicity, and 99 of 100 on FreeSolv — is consistent
   with this structural ceiling.

3. **Genuinely novel indices must use structure beyond endpoint degrees**:
   distance, spectrum, centrality, motif, eccentricity, information
   content, or some hybrid. Each of these design philosophies has prior
   literature in chemoinformatics (see `docs/literature_notes.md`), and
   the screening pipeline tests whether a candidate from these families
   adds non-redundant information.

4. **The screening pipeline is therefore the right pre-flight check** for
   any proposed topological index. It is a low-cost test (≈ seconds on a
   1000-molecule benchmark) that filters out redundant proposals before
   they consume reviewer attention or research time.

---

## 6. What we do *not* claim

- We do not claim the BID structure observation is novel. The fact that
  endpoint-degree edge-sum indices are determined by `m_ij` counts is
  well-known to the math-chem community (Gutman, Borovićanin, Furtula
  and collaborators have used this framing extensively in BID-bound
  papers). Our contribution is the **explicit dimension-10 statement
  with redundancy consequences** and the **screening pipeline** that
  operationalizes it.

- We do not claim that being in `V` is *bad*. Many useful classical
  indices are in `V`. The point is that **any new BID index proposal
  must be tested for linear dependence on existing BID indices**, since
  adding parameters or non-linear transformations cannot escape the
  10-D ceiling on the molecular-graph family considered.

- We do not claim the screening pipeline guarantees QSPR usefulness.
  Orthogonality is **necessary but not sufficient** for predictive value.
  An index can clear the |r|<0.95 redundancy threshold and still add
  effectively zero predictive signal (cf. our partial-correlation
  analyses).

---

## 7. Open questions

- For `Δ ≤ 3` (acyclic, no quaternary atoms) the dimension drops to 6
  (`P_3 = {(i,j) : 1 ≤ i ≤ j ≤ 3}`). For `Δ ≤ 5` (rare; phosphorus,
  sulfur in some oxidation states) the dimension grows to 15. Worth
  noting for materials chemistry / inorganic extensions.
- Replacing the simple-graph model with bond-order-weighted graphs
  expands the basis dimension. Whether this expansion buys meaningful
  orthogonal information is an open empirical question.
- Whether the dimension argument generalizes to "edge-sum over
  *k*-vertex neighborhood degrees" (Mondal-style neighborhood degree
  indices) is an interesting extension; the basis would be larger but
  finite.

[Phase 2 of this project will derive sharp bounds for the Wuzi family
over trees, unicyclic, and regular graph classes — the standard
parametric-index treatment.]
