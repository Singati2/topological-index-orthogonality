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
**Status: CONFIRMED.** This formula reproduces the DSO paper's Figure-5 percentages
*exactly* as integer distinct-counts out of $N=106$: DSO $19.81\%\to85$ distinct,
ESO (least degenerate) $16.98\%\to88$, RSO (most degenerate) $22.64\%\to82$
(each $=100(1-k/106)$ for integer $k$). So our convention is the Konstantinova [8]
convention DSO used; no further reconciliation needed. (Distinct values counted after
rounding to 6 decimals.)

**Structure sensitivity** (decanes, $N=75$): for index values $\{I_i\}$ with mean $\mu$,
we use two simple **dispersion** measures defined entirely from the descriptor values,
$$SS=\sigma/\mu\ \text{(coefficient of variation)},\qquad Abr=(\max I-\min I)/\mu\ \text{(relative range)},\qquad SA=SS/Abr.$$
**Status: these are OUR OWN measures, NOT the DSO / Rakić–Furtula procedure.** DSO's
method [12] (Rakić–Furtula 2019, *J. Chemometrics*) is **fingerprint-based** — it builds
Morgan circular fingerprints for each molecule and weights descriptor differences by
pairwise **Tanimoto** structural similarity. That requires cheminformatics fingerprinting
(RDKit-class tooling) and is **not reproducible in this no-install environment, and not
validatable against DSO Table 2** (our `SO` gives $SS=0.098$ vs DSO's $0.193$ — different
method, as expected). **Therefore the manuscript must NOT claim to follow DSO or
Rakić–Furtula for §5.4, and the values are NOT directly comparable to DSO Table 2.** We
present $SS,Abr,SA$ only as a self-contained, reproducible within-study comparison; the
relative ranking (LO(0,0,2) most sensitive) is the only claim made. SA$=SS/Abr$ does
coincide with DSO's *definition* of SA, but the underlying SS/Abr differ.
