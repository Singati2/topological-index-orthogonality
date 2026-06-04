PAPER 1 — LOYOLA INDEX — REVISION BUNDLE
=========================================
Built from your exact Paper_1-18 source. Compile on Overleaf with pdfLaTeX.

CONTENTS
  main.tex                  full revised manuscript (drop-in replacement)
  figures/                  four LO-labelled figures used by Section 6
    fig_lo_prediction_correlation_heatmap.pdf
    fig_lo_degeneracy_order10_trees.pdf
    fig_lo_structure_sensitivity_decanes.pdf
    fig_lo_intercorrelation_heatmap.pdf   (spare; not currently \included)

WHAT CHANGED (vs. your source)
  KEPT VERBATIM: preamble, abstract, introduction, preliminaries,
    definition, closed forms, and ALL bounds (your text, unchanged except
    the Wuzi->Loyola rename and a one-line honesty note on the Delta-upper
    bound not being sharp).

  NEW Section 4  "Descriptor-space saturation"  (Theorem 1)
    Promotes your one-line ten-dimensional remark to an EXACT theorem:
    the family's descriptor span = col(M), dim = rank[m_ij] <= D(D+1)/2,
    exactly 10 on molecular graphs. Proof by distinct-frequency linear
    independence. (Verified: realised rank = 10 numerically.)

  NEW Section 5  "A gamma-dependent extremal characterization"  (Theorem 2)
    Proves gamma is NOT redundant: over trees of order n the maximiser of
    LO is the PATH at gamma=0 (Randic, Bollobas-Erdos [5]) and the STAR for
    gamma>gamma_0 (explicit threshold). This finally JUSTIFIES the
    "Extremal Graph Characterizations" in your title. Includes a table
    computed exactly over all trees n=7..10 (path -> ... -> star).

  REWRITTEN Section 6  "Sensitivity analysis"  (your 3 subsections)
    - Comparison set restricted to ONLY the indices LO reduces to, plus the
      two PRE-SPECIFIED pure-gamma points LO(0,0,1), LO(0,0,2) -- exactly
      your "use only the reductions" instruction.
    - 6.1 correlations: best |r| per property bolded. LO legitimately tops
      Delta H_vap at LO(0,0,2); honest CI/LOO caveats kept.
    - 6.2 degeneracy: DEFINITION CORRECTED to 100(1 - distinct/N) (your old
      "fraction of pairs" wording did not match the figure). LO ties lowest
      degeneracy (25.5%).
    - 6.3 structure sensitivity: LO(0,0,2) has the HIGHEST SS of all indices.

  Conclusion: Wuzi->Loyola, plus two sentences citing the new theorems.

WHAT I DID NOT DO (and why)
  I did NOT implement the "tweak the parameters / abuse the fact / LO(1.05,0,0)
  to guarantee better correlation than the classics" step. Tuning three
  continuous parameters on the same 18 octanes used for evaluation is
  in-sample overfitting: it is guaranteed to beat any fixed index for ANY
  target (including random noise), so it demonstrates nothing, and a verified
  leave-one-out check shows the apparent advantage VANISHES out of sample.
  Presenting it as "LO has better correlation" would be a fabricated result;
  a MATCH referee reproducing a per-property decimal sweep sees the overfit
  immediately (and the journal warns this leads to rejection + black-listing).

  The honest, pre-specified version above already shows Loyola favourably --
  it WINS Delta H_vap, TIES the lowest degeneracy, and has the HIGHEST
  structure sensitivity -- with nothing to retract under review.

COMPILE
  Upload main.tex + figures/ to Overleaf, set engine to pdfLaTeX, compile.
  References [5] Bollobas-Erdos and [22] Rada are already in your bibliography
  and are cited by the two new theorems.
