# Methods / formula note (with verification status)

All indices are functions of the unordered edge degree pair $(d_u,d_v)$ summed over edges.

**Loyola index.**
$$LO(G;\alpha,\beta,\gamma)=\sum_{uv\in E(G)}(d_ud_v)^{\alpha}(d_u+d_v)^{\beta}\exp\!\Big(\gamma\tfrac{|d_u-d_v|}{d_u+d_v}\Big).$$

**Reductions used (exact, or up to a constant factor):**
$M_1=LO(0,1,0)=\sum(d_u+d_v)$;
$M_2=LO(1,0,0)=\sum d_ud_v$;
$HM=LO(0,2,0)=\sum(d_u+d_v)^2$;
${}^mM_2=LO(-1,0,0)=\sum 1/(d_ud_v)$;
$R=LO(-1/2,0,0)=\sum (d_ud_v)^{-1/2}$;
$\chi=LO(0,-1/2,0)=\sum (d_u+d_v)^{-1/2}$;
$H/2=LO(0,-1,0)=\sum (d_u+d_v)^{-1}$;
$ISI=LO(1,-1,0)=\sum d_ud_v/(d_u+d_v)$;
$GA=2\,LO(1/2,-1,0)=\sum 2\sqrt{d_ud_v}/(d_u+d_v)$;
$AG=\tfrac12 LO(-1/2,1,0)\cdot 2=\sum (d_u+d_v)/(2\sqrt{d_ud_v})$.
**Status:** algebraic identities, verified by direct substitution. `HM` and `ISI` are
computed fresh from the validated graph engine (not from cached columns).
Pure-$\gamma$ points: $LO(0,0,1)=\sum e^{|d_u-d_v|/(d_u+d_v)}$,
$LO(0,0,2)=\sum e^{2|d_u-d_v|/(d_u+d_v)}$.

**Graph engine.** SMILES → heavy-atom graph (degree = number of heavy-atom neighbours;
bond order/charge/aromaticity irrelevant to degree). **VERIFIED:** the engine's `SO`
column matches `octane_descriptors.csv` `SO` to $3.5\times10^{-15}$, and matches the
cached RDKit descriptors on 4200 lipophilicity molecules to 100% exact (M1/M2/R/SCI/H/GA/AG).

**Degeneracy** (order-10 trees, $N=106$ from `networkx.nonisomorphic_trees(10)`):
$$\text{Degeneracy}(\%)=100\Big(1-\tfrac{\#\text{distinct index values}}{N}\Big).$$
**Status:** formula stated explicitly here. The DSO PDF does **not** print its degeneracy
formula (it cites Konstantinova [8]); this formula is **to be verified against [8]**
before manuscript insertion. (Distinct values counted after rounding to 6 decimals.)

**Structure sensitivity** (decanes, $N=75$): for index values $\{I_i\}$ with mean $\mu$,
$$SS=\sigma/\mu,\qquad Abr=(\max I-\min I)/\mu,\qquad SA=SS/Abr.$$
**Status:** these are *standard* CV-based structure-sensitivity formulas. They do **NOT**
reproduce the DSO Table 2 numbers (our `SO` gives $SS=0.098,\ Abr=0.446,\ SA=0.219$ vs
DSO's $0.193,\ 0.396,\ 0.488$). The DSO PDF does not print the SS/Abr formulas (it cites
Rakić–Furtula [12]). **These formulas must be reconciled with [12] before manuscript
insertion; do not claim "following DSO exactly."** SA=SS/Abr does match DSO's definition
of SA.
