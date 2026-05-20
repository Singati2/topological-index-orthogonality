# Bounds-strategy document for the Wuzi index family

This document records the **rigorous theorem templates**, the
**ratio-bound principle**, the **monotonicity discussion**, and the
**candidate extremal graph list** that should help Sections 4, 5,
and 6 of Paper 1 (`docs/paper1_wuzi.tex`) be completed by Advik
Natarajan and Dr. Michael Arockiaraj. (The .tex section numbering
is the authoritative one: §3 is the BID-basis observation, §4 is
bounds in $(n, m, \Delta, \delta)$, §5 is bounds via classical
indices, §6 is extremal graphs.)

> **Position.** This document is scaffolding. It contains theorem
> statements that are *true and rigorously proved* (ratio bounds and
> the bracket bound), and **conjectures** marked clearly as such for
> the extremal-graph characterization. It does **not** replace the
> mathematical work that needs to be done by the senior author. In
> particular: nothing here should be treated as a final theorem of
> Paper 1 without a careful re-derivation and equality-condition
> analysis.

---

## A. Degree-pair contribution framework

### A.1 Definition of the edge-contribution function

For positive reals $x, y$ and parameters
$\alpha, \beta, \gamma \in \mathbb{R}$, define

$$\psi(x, y; \alpha, \beta, \gamma) = (x\,y)^{\alpha}\,(x + y)^{\beta}\,\exp\!\left(\gamma\cdot\frac{|x - y|}{x + y}\right).$$

The Wuzi index of a graph $G$ is

$$W(G; \alpha, \beta, \gamma) = \sum_{uv \in E(G)} \psi\bigl(d_u, d_v; \alpha, \beta, \gamma\bigr).$$

### A.2 Admissible degree pairs on $\mathcal{G}_4$

Hydrogen-suppressed molecular graphs have $\Delta \le 4$. The
admissible edge degree pairs $(i, j)$ with $1 \le i \le j \le 4$
form a set of size $\binom{4+1}{2} = 10$:

| $(i,j)$ | Note (typical chemistry) |
|---|---|
| $(1, 1)$ | An isolated bond. Empty in any connected molecular graph of $\ge 3$ heavy atoms. |
| $(1, 2)$ | Pendant on a degree-2 carbon (terminal methyl-like). Very common. |
| $(1, 3)$ | Pendant on a tertiary carbon. Common. |
| $(1, 4)$ | Pendant on a quaternary carbon. Less common. |
| $(2, 2)$ | Interior of a chain. Very common in alkanes. |
| $(2, 3)$ | Branch-point junction. Common. |
| $(2, 4)$ | Chain meets a quaternary carbon. Less common. |
| $(3, 3)$ | Two tertiary carbons adjacent (heavily branched alkanes). |
| $(3, 4)$ | Branched--quaternary junction. |
| $(4, 4)$ | Two quaternary carbons adjacent. Rare; appears in cage/cluster molecules. |

### A.3 Universal bracket bound (rigorous)

**Theorem A.1 (Bracket bound).** *Let $G$ be a simple connected graph
of size $m \ge 1$ whose degrees lie in $[\delta, \Delta]$. Define*

$$\psi_{\min}(\delta, \Delta) = \min\{\psi(i, j; \alpha, \beta, \gamma) : \delta \le i \le j \le \Delta\},$$

$$\psi_{\max}(\delta, \Delta) = \max\{\psi(i, j; \alpha, \beta, \gamma) : \delta \le i \le j \le \Delta\}.$$

*Then*

$$m \cdot \psi_{\min}(\delta, \Delta) \le W(G; \alpha, \beta, \gamma) \le m \cdot \psi_{\max}(\delta, \Delta).$$

*The upper inequality is an equality if and only if every edge of
$G$ has a degree pair* $\{d_u, d_v\}$ *attaining $\psi_{\max}$ over
the admissible set; the lower inequality is dual. If the
extremum is attained at a unique admissible pair, this collapses to
"every edge of $G$ has the same unordered degree pair", which holds
for regular graphs and some biregular graphs.*

**Proof.** Per-edge, $\psi_{\min} \le \psi(d_u, d_v; \cdot) \le \psi_{\max}$.
Summing over the $m$ edges yields the inequalities. Equality on
the upper side forces every term in the sum to equal $\psi_{\max}$,
i.e. every edge to have a degree pair achieving the maximum; lower
side is dual.  $\blacksquare$

> The bracket is tight on $K_n$, $C_n$, $Q_k$, and any $r$-regular
> graph. It is usually loose on chemistry datasets because most
> molecular graphs are non-regular.

---

## B. Ratio bounds in terms of classical indices

This section gives a **rigorous template** for bounding $W$ above
and below in terms of a chosen classical BID index $I_h$, provided
$h(i, j) > 0$ on every admissible degree pair appearing in $G$.

### B.1 Edge contributions of seven classical indices

For each classical index $I_h(G) = \sum_{uv \in E(G)} h(d_u, d_v)$,
the corresponding edge contribution $h(i, j)$ is recorded below.

| Index | $h(i, j)$ | Positive on $\mathcal{G}_4$? |
|---|---|---|
| First Zagreb $M_1$        | $i + j$                              | yes ($\ge 2$) |
| Second Zagreb $M_2$       | $i\,j$                               | yes ($\ge 1$) |
| Randić $R$                | $1 / \sqrt{i\,j}$                    | yes |
| Sombor $SO$               | $\sqrt{i^2 + j^2}$                   | yes |
| Geometric-arithmetic $GA$ | $2\sqrt{i\,j}\,/(i + j)$             | yes |
| Harmonic $H$              | $2\,/(i + j)$                        | yes |
| Atom-bond conn.\ $ABC$    | $\sqrt{(i + j - 2)/(i\,j)}$          | **0 at $(1,1)$, positive elsewhere** |

The $ABC$ ratio bound therefore applies on graphs that do not
contain a $(1,1)$ edge, i.e. on any connected graph of order $\ge 3$
(since a connected graph with a $(1,1)$ edge has only that single
edge: $P_2$).

### B.2 The ratio principle (rigorous)

**Theorem B.1 (Ratio bound).** *Let $G$ be a simple graph with
maximum degree $\Delta$, and let $I_h(G) = \sum_{uv \in E(G)} h(d_u, d_v)$
be a BID index whose edge contribution $h$ is strictly positive on
every degree pair $(i, j)$ appearing in $G$. Define*

$$c_{\min}^{h} = \min_{(i,j)\,\text{admissible}}\,\frac{\psi(i, j; \alpha, \beta, \gamma)}{h(i, j)},
\qquad
c_{\max}^{h} = \max_{(i,j)\,\text{admissible}}\,\frac{\psi(i, j; \alpha, \beta, \gamma)}{h(i, j)}.$$

*Then*

$$c_{\min}^{h} \cdot I_h(G) \;\le\; W(G; \alpha, \beta, \gamma) \;\le\; c_{\max}^{h} \cdot I_h(G).$$

*Equality on either side holds if and only if every edge of $G$
attains the corresponding extremum of the ratio
$\psi/h$.*

**Proof.** Fix an edge $uv \in E(G)$ with degree pair $(d_u, d_v)$.
By the definition of the minimum,
$\psi(d_u, d_v) \ge c_{\min}^h \cdot h(d_u, d_v)$. Summing over
edges gives $W(G) \ge c_{\min}^h \cdot I_h(G)$. The upper bound and
the equality conditions are analogous.  $\blacksquare$

> **Practical computation.** $c_{\min}^{h}$ and $c_{\max}^{h}$ are
> finite-set extrema over the at most $10$ admissible degree pairs
> on $\mathcal{G}_4$, so they can be computed by brute force given
> any $(\alpha, \beta, \gamma)$. The script
> `scripts/09_wuzi_bounds_tables.py` produces these constants in
> tabular form for a selection of parameter triples and the
> seven classical indices above.

### B.3 Why these bounds are useful

The ratio bound is the natural template for "bounds via classical
indices" because:

1. It is mathematically rigorous --- the only hypothesis is
   positivity of $h$ on the admissible degree pairs.
2. The constants $c_{\min}^{h}, c_{\max}^{h}$ depend only on the
   parameter triple $(\alpha, \beta, \gamma)$ and the maximum
   degree $\Delta$, *not* on the graph $G$ itself.
3. At the special parameter values where $\Wuzi$ reduces to a
   classical index (e.g. $\alpha = 0, \beta = -1, \gamma = 0$
   gives $W = H/2$), the corresponding ratio constants collapse to
   the relevant constant of the identity, providing a check.
4. The bound is sharp when $\psi/h$ is constant on $G$'s edge degree
   pairs, in particular on regular graphs.

### B.4 Limitation: the bound is not always tight

The ratio bound is tight on regular graphs but can be loose on
graphs with many distinct degree pairs. Tighter bounds may exist
that interpolate over the distribution of $m_{ij}$ counts; this is
left for the formal derivation in Paper 1.

---

## C. Parameter-region monotonicity

This section discusses how $\psi(x, y; \alpha, \beta, \gamma)$
varies as a function of its three parameters and as a function of
the degree pair $(x, y)$.

### C.1 Monotonicity in parameters (rigorous)

Fix a positive degree pair $(x, y) \in \mathbb{R}_{>0}^2$. Then:

- $\psi$ is **strictly increasing** in $\alpha$ if $xy > 1$, strictly
  decreasing if $xy < 1$, and constant in $\alpha$ if $xy = 1$
  (only the trivial $(1,1)$ pair on $\mathcal{G}_4$).
- $\psi$ is **strictly increasing** in $\beta$ if $x + y > 1$
  (always true on $\mathcal{G}_4$ since $x + y \ge 2$), so $\psi$
  is strictly increasing in $\beta$ on any non-trivial graph.
- $\psi$ is **strictly increasing** in $\gamma$ if $x \ne y$
  (imbalance $> 0$), and constant in $\gamma$ if $x = y$ (regular
  edge).

**Proof.** Direct from the partial derivatives:
$\partial_\alpha \log\psi = \log(xy)$,
$\partial_\beta \log\psi = \log(x+y)$,
$\partial_\gamma \log\psi = |x - y|/(x + y)$,
all evaluated at $(x, y; \alpha, \beta, \gamma)$.  $\blacksquare$

### C.2 Monotonicity in degree pair (rigorous in part)

In the **product axis** $u := xy$:

- $\psi$ is strictly increasing in $u$ if $\alpha > 0$, strictly
  decreasing if $\alpha < 0$, and constant in $u$ if $\alpha = 0$
  (with all else fixed).

In the **sum axis** $s := x + y$:

- $\psi$ is strictly increasing in $s$ if $\beta > 0$, strictly
  decreasing if $\beta < 0$, and constant in $s$ if $\beta = 0$
  (with $u$ and the imbalance fixed).

In the **imbalance axis** $r := |x - y|/(x + y) \in [0, 1)$:

- $\psi$ is strictly increasing in $r$ if $\gamma > 0$, strictly
  decreasing if $\gamma < 0$, and constant in $r$ if $\gamma = 0$.

Note however that the three axes $(u, s, r)$ are **not algebraically
independent** on $\mathbb{R}_{>0}^2$ (given $u, s$, the pair $\{x, y\}$
is uniquely determined as the roots of $t^2 - s\,t + u = 0$ and the
imbalance is then $r = \sqrt{s^2 - 4u}\,/\,s$). So the
"monotonicity-in-one-axis-with-the-others-fixed" statements above
are useful descriptive language but **do not** combine into a
joint monotonicity result on the unconstrained $(x, y)$ plane.

**Requires checking.** The full joint monotonicity behaviour of
$\psi(x, y; \alpha, \beta, \gamma)$ on the unconstrained $(x, y)$
plane, in different sign regions of $(\alpha, \beta, \gamma)$, is
**non-trivial** because of the implicit dependence between $u, s, r$.
A careful case analysis is needed for the extremal-graph proofs of
Section~D.

### C.3 Boundary behaviour of $\gamma$ (rigorous)

- At $\gamma = 0$: $\psi$ reduces to the two-parameter degree-based
  function $(xy)^\alpha (x+y)^\beta$. The Wuzi family then reduces to
  the standard two-parameter degree-based index family.
- As $\gamma \to +\infty$ (with other parameters fixed): $\psi$ is
  asymptotically dominated by the edges with maximum imbalance
  $|x - y|/(x + y)$. On $\mathcal{G}_4$ the maximum imbalance on
  the admissible pairs is $3/5$ at $(1, 4)$, giving
  $\psi \sim (xy)^\alpha (x+y)^\beta e^{3\gamma/5}$.
- As $\gamma \to -\infty$: $\psi$ is dominated by the regular pairs
  (imbalance $= 0$), and the family is asymptotically supported on
  the diagonal of the admissible-pair lattice.

---

## D. Extremal graphs: candidates and strategy

This section lists **candidate extremal graphs** by graph class
and parameter sign region. None of the statements below is a
theorem; they are stated as conjectures with the standard
heuristic motivation. Formal proofs are left for Paper 1 and
should be derived by adapting the classical edge-transformation
arguments from the Randić / Sombor extremal literature.

### D.1 Trees

For a tree $T$ of order $n \ge 3$, two natural extremal candidates
are the **path** $P_n$ and the **star** $S_n$. They differ
structurally in degree distribution:

- $P_n$: two pendants of degree $1$, $n - 2$ interior vertices of
  degree $2$. Edge degree pairs: 2 pendants $(1, 2)$ and $n - 3$
  interiors $(2, 2)$ (for $n \ge 4$).
- $S_n$: one centre of degree $n - 1$, $n - 1$ pendants of degree
  $1$. All $n - 1$ edges have degree pair $(1, n - 1)$.

**Conjecture D.1 (Trees, $\alpha \ge 0$, $\beta \ge 0$, $\gamma \ge 0$).**
*Among all trees of order $n \ge 3$:*
1. *the star $S_n$ attains $\max W(\cdot; \alpha, \beta, \gamma)$;*
2. *the path $P_n$ attains $\min W(\cdot; \alpha, \beta, \gamma)$.*

*Equality (uniqueness) is plausible but unverified.*

**Heuristic.** The star concentrates degree on one centre, producing
$n - 1$ edges of pair $(1, n - 1)$ which are large in
$(xy)^\alpha$ (high product) and in $(x + y)^\beta$ (high sum) for
positive parameters, and large in the imbalance factor for
positive $\gamma$ (high $|x - y|/(x + y)$). The path distributes
degrees uniformly at $2$, minimizing all three factors. The
conjecture should be verified by a careful edge-transformation
argument (move a pendant from a low-degree vertex to a high-degree
vertex; show $W$ strictly increases when the parameters are in the
stated region).

The conjecture is computationally tested by
`scripts/10_wuzi_extremal_search.py` on all non-isomorphic trees
of order $5 \le n \le 12$. See the output file
`results/wuzi_extremal_trees.md`.

**Other parameter sign regions** are expected to give different
extremals; in particular:
- For $\alpha < 0$ (e.g. Randić $\alpha = -1/2$): the path $P_n$
  typically becomes the maximizer and the star the minimizer.
- For $\beta < 0$ and $\alpha = \gamma = 0$ (sum-connectivity
  region): similarly reversed.
- The mixed region $\alpha > 0, \beta < 0$ is non-trivial.

### D.2 Unicyclic graphs

A unicyclic graph of order $n$ is a connected graph with exactly
one cycle, equivalently it has $n$ edges and one cycle. Two
natural extremal candidates:

- $C_n$: the cycle of order $n$; $2$-regular, all edge degree
  pairs $(2, 2)$.
- $U_n^{(1)}$: the unicyclic graph obtained by attaching $n - 3$
  pendant vertices to one vertex (the "hub") of $C_3$. The hub
  then has degree $2 + (n - 3) = n - 1$, while the two non-hub
  triangle vertices have degree $2$. The triangle edges are
  therefore $(2, n - 1)$ (hub-to-non-hub, two such edges) and
  $(2, 2)$ (the single non-hub triangle edge); each of the $n - 3$
  pendant edges has degree pair $(1, n - 1)$.

**Conjecture D.2 (Unicyclic, $\alpha \ge 0$, $\beta \ge 0$, $\gamma \ge 0$).**
*Among unicyclic graphs of order $n \ge 4$:*
1. *the maximum is attained by a "hub-concentrated" unicyclic graph*
   *such as $U_n^{(1)}$ (or a closely related $C_3$-based hub);*
2. *the minimum is attained by the cycle $C_n$.*

**Heuristic.** Same intuition as trees: concentrating degree
increases $\psi$ for positive parameters; the regular cycle
minimizes.

### D.3 Bicyclic graphs

A bicyclic graph of order $n$ has $n + 1$ edges and exactly two
cycles. The class is structurally richer than trees and unicyclic:

- **Type (a):** two cycles sharing one edge (a $\theta$-graph).
- **Type (b):** two cycles sharing one vertex (a "bowtie"
  generalization).
- **Type (c):** two cycles joined by a path.

Each of (a), (b), (c) may contain additional pendant vertices.
Candidate extremals:

- The "uniform" bicyclic graph $C_3 + C_3$ sharing a vertex (the
  bowtie graph $B_n$ on $n$ vertices with $n - 5$ pendants on a
  central vertex), for the maximum in the positive-parameter region.
- A near-regular bicyclic graph for the minimum.

**Conjecture D.3 (Bicyclic, $\alpha \ge 0$, $\beta \ge 0$, $\gamma \ge 0$).**
*Among bicyclic graphs of order $n \ge 5$, the maximum and minimum
of $\Wuzi$ are attained at structures yet to be characterized,
plausibly the "hub-concentrated" bicyclic and a near-regular
bicyclic respectively.*

**To be characterized.** The bicyclic case requires the most
careful work and is the natural target for a multi-page Section 5
in Paper 1.

### D.4 Computational support for the conjectures

`scripts/10_wuzi_extremal_search.py` enumerates all non-isomorphic
trees of order $5 \le n \le 12$ (1, 1, 1, 2, 3, 6, 11, 23, 47,
106, 235, 551 non-isomorphic trees for $n = 1, \ldots, 12$
respectively) and records the minimizer and maximizer of $\Wuzi$
at each of eight parameter triples. The output file
`results/wuzi_extremal_trees.md` reports the canonical (Weisfeiler--
Lehman) label and degree sequence of the observed extremals.
This empirical evidence is **conjecture-generating**, not proof.

**Three sign regions observed.** Aggregating across $n = 5, \ldots, 12$:

| Parameter region | Observed minimizer | Observed maximizer |
|---|---|---|
| $\alpha, \beta, \gamma \ge 0$ (four triples sampled) | $P_n$ (every $n$) | $S_n$ (every $n$) |
| $\alpha \le 0$ or $\beta \le 0$ with $\gamma = 0$ (three triples) | $S_n$ (every $n$) | $P_n$ (every $n$) |
| $(\alpha, \beta, \gamma) = (-1, -1, 1)$ (mixed-sign) | $S_n$ for $n \ge 6$ | non-classical (see §D.5) |

### D.5 The non-classical regime at $(-1, -1, 1)$

The most interesting empirical finding of the computational search.
For $n \ge 6$ the observed maximizer of $W(T; -1, -1, 1)$ is
*neither* $P_n$ *nor* $S_n$. An is-caterpillar check on the
observed maximizer (verified by `is_caterpillar` in
`scripts/10_wuzi_extremal_search.py`) further distinguishes two
sub-cases:

| $n$ | Argmax degree sequence | Is caterpillar? |
|---|---|---|
| $6$  | $(3, 2, 2, 1, 1, 1)$               | yes |
| $7$  | $(3, 2, 2, 2, 1, 1, 1)$            | no  |
| $8$  | $(3, 2, 2, 2, 2, 1, 1, 1)$         | no  |
| $9$  | $(4, 2, 2, 2, 2, 1, 1, 1, 1)$      | no  |
| $10$ | $(3, 3, 2, 2, 2, 2, 1, 1, 1, 1)$   | no  |
| $11$ | $(5, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1)$| no  |
| $12$ | $(4, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1)$ | no  |

For $n = 6$ the observed argmax is a caterpillar (specifically a
broom-like tree with a single branch vertex of degree 3 on a chain
of length 2). For every $n \ge 7$ tested the argmax is NOT a
caterpillar: removing the leaves leaves a non-path subgraph.
Inspection of the search-script output (canonical labels in the
CSV) shows that the argmax for $n \ge 7$ is a *spider*-like or
*branched* tree where multiple short paths meet at one or two
internal vertices. The full structural description is left for
the formal derivation.

**Conjecture D.5 (Mixed-sign non-classical extremal).**
For the parameter triple $(\alpha, \beta, \gamma) = (-1, -1, 1)$
and tree order $n \ge 7$, the maximum of $W$ over all trees of
order $n$ is attained at a tree that is neither the path $P_n$,
the star $S_n$, nor a caterpillar. The precise structural
characterization is an open problem.

This conjecture is the most plausibly *new* mathematical
observation that the parametric Wuzi family produces: a parameter
triple in which the standard star / path extremals fail and the
maximizer escapes the caterpillar class entirely. A formal
characterization should adapt the Kelmans transformation and the
pendant-shift argument of the classical extremal-trees literature
to the joint weight
$\psi(d_u, d_v; -1, -1, 1) = (d_u d_v)^{-1}(d_u + d_v)^{-1} \exp(|d_u - d_v|/(d_u + d_v))$
at $(\alpha, \beta, \gamma) = (-1, -1, 1)$.

---

## E. Proof techniques (cheat-sheet)

The following standard techniques are expected to be needed in the
formal derivations of Paper 1 Sections 3, 4, and 5.

1. **Edge-degree pair counts $m_{ij}$.** Express any BID index as
   $I_f(G) = \sum_{1 \le i \le j \le \Delta} f(i, j)\,m_{ij}(G)$
   (Theorem 3.1 of Paper 1). This converts edge sums into finite
   linear combinations.

2. **Edge-moving transformations.** The Kelmans transformation, the
   "shift a pendant from $u$ to $v$" move, branch contraction, and
   edge subdivision are the standard tools for tree- and
   unicyclic-extremal arguments. Each preserves the order $n$ while
   changing the degree sequence in a controlled way.

3. **Majorization and Schur convexity.** When $\psi$ can be written
   as a Schur-convex (resp. Schur-concave) symmetric function of
   the degree sequence, majorization comparisons between graph
   degree sequences give immediate ordering of $\Wuzi$.

4. **Jensen's inequality.** For convex / concave functions of the
   edge degree pair, Jensen gives the kind of bound used in
   the $M_1$-via-Jensen bound of Paper 1 §4 (label `thm:m1_bound`).

5. **Chebyshev's sum inequality.** For two similarly-ordered
   sequences $(a_i), (b_i)$,
   $N \sum a_i b_i \ge (\sum a_i)(\sum b_i)$. Useful for products
   of edge contributions.

6. **Cauchy-Schwarz inequality.** Used in the Cauchy–Schwarz
   theorem of Paper 1 §4 (label `thm:cs`) to derive the
   geometric-mean inequality
   $W(G; \bar\alpha, \bar\beta, \bar\gamma) \le \sqrt{W(G; \alpha_1, \beta_1, \gamma_1) \cdot W(G; \alpha_2, \beta_2, \gamma_2)}$
   where the bars denote arithmetic averages.

7. **Power-mean inequality.** For different exponents $p, q$ with
   $p > q$, the $L^p$ and $L^q$ norms of the edge-contribution
   sequence are comparable.

8. **Discrete verification over the degree-pair lattice for $\Delta \le 4$.**
   On $\mathcal{G}_4$ the degree-pair lattice has only $10$ admissible
   pairs. Any inequality involving $\psi$ on a single edge can be
   verified by brute-force enumeration over the lattice; this is
   how $c_{\min}^h$ and $c_{\max}^h$ in Theorem B.1 are computed.

9. **Exhaustive enumeration for small $n$.** All non-isomorphic
   trees up to about $n = 18$ are enumerable in seconds with
   `networkx.algorithms.tree.nonisomorphic_trees`. This makes
   small-$n$ computational tests of Conjectures D.1--D.3 cheap.

10. **The Polya--Szegő inequality.** Useful for inverting
    Cauchy-Schwarz-type bounds when degree pairs are bounded
    above and below.

---

## F. What is *not* in this document

- A final theorem statement for any bound of Paper 1 §4 (bounds in
  $n, m, \Delta, \delta$ other than the bracket).
- A final extremal-graph characterization for Paper 1 §6.
- Specific numerical constants in upper / lower bounds beyond what
  $c_{\min}^h, c_{\max}^h$ provide.
- Equality-condition analysis for any bound other than the bracket
  and the ratio bound.

These items are reserved for Advik Natarajan and
Dr. Michael Arockiaraj.

---

## G. Pointers

- `scripts/09_wuzi_bounds_tables.py` --- generates the ratio
  constants $c_{\min}^h, c_{\max}^h$ for $(\alpha, \beta, \gamma)$
  in a selected list and $h \in \{M_1, M_2, R, SO, GA, H, ABC\}$.
- `scripts/10_wuzi_extremal_search.py` --- enumerates trees of
  order $5 \le n \le 12$, computes $\Wuzi$, records observed
  minimizers and maximizers, and checks whether they agree with
  $P_n / S_n$ at each parameter setting.
- `results/wuzi_bounds_ratio_tables.{csv,md}` --- outputs of
  script 09.
- `results/wuzi_extremal_trees.{csv,md}` --- outputs of script 10.
