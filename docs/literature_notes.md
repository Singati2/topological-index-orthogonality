# Literature notes — references and positioning

This document enforces citation discipline for the project. **Every
claim of novelty must point to a `[REFERENCE NEEDED]` slot that someone
fills in by reading the cited paper, not by guessing.** Several
candidate "novel" framings considered during early scoping turned out
to overlap substantially with prior literature; the entries below
record those checks so the same ground does not get re-walked.

Citations below are grouped by topic. Where a citation has been verified
from a PDF in this repo's reference materials, it is marked **[verified]**.
All others are `[REFERENCE NEEDED]` and must be filled before submission.

---

## 1. Classical topological indices

- **Wiener index.** [REFERENCE NEEDED — H. Wiener, "Structural Determination
  of Paraffin Boiling Points," *J. Am. Chem. Soc.* 69(1), 1947, 17-20.]
- **Randić index.** [verified — see refs of MATCH 95:141-162 / arXiv DSO
  paper]: M. Randić, "Characterization of molecular branching," *J. Am.
  Chem. Soc.* 97(23), 1975, 6609-6615.
- **First Zagreb M_1.** [verified]: I. Gutman, N. Trinajstić, "Graph
  theory and molecular orbitals. Total φ-electron energy of alternant
  hydrocarbons," *Chem. Phys. Lett.* 17(4), 1972, 535-538.
- **Albertson (third Zagreb) Alb.** [verified]: M. O. Albertson, "The
  irregularity of a graph," *Ars Comb.* 46, 1997, 219-225.
- **Geometric-arithmetic GA.** [verified]: D. Vukičević, B. Furtula,
  "Topological index based on the ratios of geometrical and arithmetic
  means of end-vertex degrees of edges," *J. Math. Chem.* 46, 2009,
  1369-1376.
- **Harmonic.** [verified]: S. Fajtlowicz, "On conjectures of
  Graffiti-II," *Congr. Numer.* 60, 1987, 187-197.
- **Atom-bond connectivity ABC.** [REFERENCE NEEDED — E. Estrada et al.,
  "An atom-bond connectivity index: Modelling the enthalpy of formation
  of alkanes," *Indian J. Chem.* 37A, 1998, 849-855.]
- **Sum-connectivity χ.** [verified]: B. Zhou, N. Trinajstić, "On a novel
  connectivity index," *J. Math. Chem.* 46(4), 2009, 1252-1270.
- **Forgotten index F.** [verified]: B. Furtula, I. Gutman, "A forgotten
  topological index," *J. Math. Chem.* 53, 2015, 1184-1190.
- **Balaban J.** [REFERENCE NEEDED — A. T. Balaban, "Highly discriminating
  distance-based topological index," *Chem. Phys. Lett.* 89(5), 1982,
  399-404.]
- **Estrada index EE.** [REFERENCE NEEDED — E. Estrada, "Characterization
  of 3D molecular structure," *Chem. Phys. Lett.* 319, 2000, 713-718.]
- **Mostar index.** [REFERENCE NEEDED — T. Došlić et al., "Mostar
  index," *J. Math. Chem.* 56, 2018, 2995-3013.]
- **Szeged index.** [REFERENCE NEEDED — I. Gutman, "A formula for the
  Wiener number of trees and its extension to graphs containing cycles,"
  *Graph Theory Notes NY* 27, 1994, 9-15.]
- **General Randić R_α.** [REFERENCE NEEDED — B. Bollobás, P. Erdős,
  "Graphs of extremal weights," *Ars Comb.* 50, 1998, 225-233.]
- **General sum-connectivity χ_α.** [verified]: B. Zhou, N. Trinajstić,
  "On general sum-connectivity index," *J. Math. Chem.* 47(1), 2010,
  210-218.

## 2. BID (Bond-Incident-Degree) structure

This is the antecedent literature for our Edge-Degree-Pair Basis
theorem. The BID framework explicitly notes that endpoint-degree
edge-sum indices are determined by edge-type counts m_ij — what we
state as the dimension bound is a corollary that the BID literature
has not emphasized.

- **Borovićanin, Das, Furtula, Gutman.** [REFERENCE NEEDED — "Bounds for
  Zagreb indices," *MATCH Commun. Math. Comput. Chem.* 78, 2017,
  17-100.] Comprehensive treatment of BID-style bounds.
- **Gutman & Tošović.** [REFERENCE NEEDED — I. Gutman, J. Tošović,
  "Testing the quality of molecular structure descriptors,
  vertex-degree-based topological indices," *J. Serb. Chem. Soc.* 78,
  2013, 805-810.] Edge-partition framing of degree-based indices.
- **Gutman degree-based topological indices survey.** [verified]: I.
  Gutman, "Degree-based topological indices," *Croatica Chemica Acta*
  86(4), 2013, 351-361.

## 3. Sombor variants (Wuzi's design antecedents)

- **Sombor index SO.** [verified]: I. Gutman, "Geometric approach to
  degree-based topological indices: Sombor indices," *MATCH Commun.
  Math. Comput. Chem.* 86, 2021, 11-16.
- **Sombor extremal review.** [verified]: H. Liu, I. Gutman, L. You, Y.
  Huang, "Sombor index: review of extremal results and bounds," *J.
  Math. Chem.* 60, 2022, 771-798.
- **Reduced / Co-Sombor.** [verified]: D. T. Rajathagiri, "Enhanced
  mathematical models for the Sombor index: Reduced and co-Sombor index
  perspectives," *Data Anal. Artif. Intell.* 1(2), 2021, 215-228.
- **Diminished Sombor index (DSO).** [verified]: F. Movahedi, I. Gutman,
  I. Redžepović, B. Furtula, "Diminished Sombor Index," *MATCH Commun.
  Math. Comput. Chem.* 95(1), 2026, 141-162. doi:10.46793/match95-1.14125.

  **This is the structural template for our Phase 2 paper.** See
  `docs/phase2_plan.md`.

- **DSO bounds via classical indices.** [verified]: F. Movahedi,
  "Diminished Sombor index and its relationship with topological
  indices," August 2025, arXiv preprint. This is the template for our
  Section 4 (Wuzi bounds in classical indices).
- **Banhatti-Sombor.** [verified]: V. R. Kulli, "On Banhatti-Sombor
  indices," *Int. J. Appl. Chem.* 8(2021), 21-25.

## 4. Intercorrelation analyses (motivating our screening framework)

This is the load-bearing citation set. Each of these papers does
*informal* redundancy analysis on Sombor / Zagreb variants. Our
contribution is to formalize the methodology.

- **The DSO Section 5.2 observation.** [verified] — Movahedi, Gutman,
  Redžepović, Furtula 2026 (above) Section 5.2, page 157: report
  correlations 0.92-0.99 between DSO and other Sombor variants on
  octanes, conclude DSO "captures some additional information." No
  threshold, no partial correlation, no PCA, n=18.

  **This is the example we cite to motivate the need for orthogonality
  screening. We do *not* frame it as a critique of the paper — its
  authors are field-leading. We cite it to show the field is asking the
  right question informally, and our framework provides the formal
  answer.**

- **Aslam et al. degree-entropy QSPR.** [REFERENCE NEEDED — verify
  citation details]: M. Aslam et al., "Degree-Based Graph Entropy in
  Structure-Property Modeling," *Entropy* 25(7), 2023, 1092. Combines
  graph entropy with degree-based indices for QSPR; informal
  correlation comparison.
- **Sombor extremal review intercorrelation.** [verified — see Liu et al.
  2022 above]. Includes informal correlation discussion across Sombor
  variants.

## 5. Graph entropy literature (foreclosed our "InfoH_deg as new index" plan)

These pre-empt our 20-candidate experiment's Shannon-degree-entropy
result. We cite them in the paper to make clear we are *not* claiming
novelty for degree entropy:

- **Rashevsky 1955.** [REFERENCE NEEDED — N. Rashevsky, "Life,
  information theory, and topology," *Bull. Math. Biophys.* 17, 1955,
  229.]
- **Mowshowitz 1968.** [REFERENCE NEEDED — A. Mowshowitz, "Entropy and
  the complexity of graphs, I-IV," *Bull. Math. Biophys.* 30, 1968.]
- **Bonchev information indices.** [REFERENCE NEEDED — D. Bonchev,
  *Information Theoretic Indices for Characterization of Chemical
  Structures,* Research Studies Press, 1983.]
- **Cao, Dehmer, Shi 2014.** [REFERENCE NEEDED — S. Cao, M. Dehmer, Y.
  Shi, "Extremality of degree-based graph entropies," *Information
  Sciences* 278, 2014, 22-33.]
- **Ghorbani et al. 2018.** [REFERENCE NEEDED — M. Ghorbani, M. Dehmer
  et al., "First degree-based entropy of graphs," *J. Appl. Math.
  Comput.* 59, 2018-2019, 303-310.] Paper title is literally our
  candidate index name.
- **Dehmer-Mowshowitz 2011.** [REFERENCE NEEDED — M. Dehmer, A.
  Mowshowitz, "A history of graph entropy measures," *Information
  Sciences* 181(1), 2011, 57-78.]
- **Sabirov-Shepelevich 2021.** [REFERENCE NEEDED — D. Sh. Sabirov, I.
  S. Shepelevich, "Information entropy in chemistry: An overview,"
  *Entropy* 23(10), 2021, 1240.]

## 6. Superindex / composite-index literature (foreclosed our composite-family escape hatch)

- **Bonchev, Mekenyan, Trinajstić 1981.** [REFERENCE NEEDED — D.
  Bonchev, O. Mekenyan, N. Trinajstić, "Isomer discrimination by
  topological information approach," *J. Comput. Chem.* 2(2), 1981,
  127-148.] Coined the term "superindex."
- **Dehmer-Mowshowitz 2013.** [REFERENCE NEEDED — M. Dehmer, A.
  Mowshowitz, "The discrimination power of structural superindices,"
  *PLOS ONE* 8(8), 2013, e70551.] Explicit cross-family composite
  construction.
- **Pogliani LCXCI.** [REFERENCE NEEDED — L. Pogliani, "From molecular
  connectivity indices to semiempirical connectivity terms," *Chem.
  Rev.* 100(10), 2000, 3827-3858.]

## 7. Cheminformatics descriptor tools (defines the descriptor universe)

- **Mordred.** [REFERENCE NEEDED — H. Moriwaki, Y. S. Tian, N.
  Kawashita, T. Takagi, "Mordred: a molecular descriptor calculator,"
  *J. Cheminformatics* 10, 2018, Article 4.]
- **RDKit descriptors.** [REFERENCE NEEDED — G. Landrum et al., "RDKit:
  Open-source cheminformatics," https://www.rdkit.org/]
- **alvaDesc / Dragon.** [REFERENCE NEEDED — A. Mauri, "alvaDesc:
  A Tool to Calculate and Analyze Molecular Descriptors and
  Fingerprints," in *Ecotoxicological QSARs*, Springer 2020.]
- **ECFP / Morgan fingerprints.** [REFERENCE NEEDED — D. Rogers, M.
  Hahn, "Extended-connectivity fingerprints," *J. Chem. Inf. Model.*
  50, 2010, 742-754.]
- **MOLTOP (network-science descriptors on ESOL).** [REFERENCE NEEDED —
  J. Adamczyk et al., "MOLTOP — Molecular topological descriptors,"
  ECAI 2024.] Direct prior art for using network-science descriptors
  on the same datasets we use; cite as concurrent and complementary
  rather than competing.

## 8. QSPR validation methodology

- **MoleculeNet.** [REFERENCE NEEDED — Z. Wu, B. Ramsundar, E. N. Feinberg
  et al., "MoleculeNet: a benchmark for molecular machine learning,"
  *Chem. Sci.* 9, 2018, 513-530.] Provides ESOL / FreeSolv /
  Lipophilicity benchmarks.
- **ESOL / Delaney.** [REFERENCE NEEDED — J. S. Delaney, "ESOL:
  estimating aqueous solubility directly from molecular structure,"
  *J. Chem. Inf. Comput. Sci.* 44(3), 2004, 1000-1005.]
- **FreeSolv.** [REFERENCE NEEDED — D. L. Mobley, J. P. Guthrie, "FreeSolv:
  a database of experimental and calculated hydration free energies,"
  *J. Comput. Aided Mol. Des.* 28, 2014, 711-720.]
- **Lipophilicity.** [REFERENCE NEEDED — A. M. Wenlock, "Defining you
  ADME territory: the Wenlock dataset," ChEMBL.]
- **OECD QSAR validation principles.** [REFERENCE NEEDED — OECD,
  "Guidance Document on the Validation of (Q)SAR Models," 2014.]
- **Scaffold split / applicability domain.** [REFERENCE NEEDED — for
  applicability domain we can cite [F. Sahigara et al., "Comparison of
  different approaches to define the applicability domain of QSAR
  models," *Molecules* 17, 2012, 4791-4810].]

---

## Citation discipline

1. No claim of "we propose [thing]" without an explicit check against
   the relevant sections above. If the section contains a paper that
   already proposed `[thing]` (even under a different name), reframe.
2. No PDF reading is "evidence" — if a citation is `[REFERENCE NEEDED]`,
   it stays that way until somebody pastes the bib entry from the actual
   reference.
3. Where we cite a paper *constructively* (DSO MATCH 95:141-162 as
   template; Aslam 2023 as related QSPR-entropy work; MOLTOP 2024 as
   complementary network-science approach), make that explicit in the
   manuscript — never frame respected work as the target of a critique.
