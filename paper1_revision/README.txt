PAPER 1 — LOYOLA INDEX — REVISION BUNDLE  (corrected)
=====================================================
Everything is in THIS folder, flat (no subfolders). Compile with pdfLaTeX.

FILES
  main.tex                                     the manuscript (compile this)
  fig_lo_prediction_correlation_heatmap.pdf    Figure 1 (octane correlations)
  fig_lo_degeneracy_order10_trees.pdf          Figure 2 (degeneracy)
  fig_lo_structure_sensitivity_decanes.pdf     Figure 3 (structure sensitivity)
  fig_lo_intercorrelation_heatmap.pdf          spare (not currently included)
  octane_correlation_table.tex                 Table 2, as a standalone file
  ("Paper 1 Revision.zip" one level up = this same folder, zipped)

HOW TO GET FIGURES TO SHOW IN OVERLEAF  (the only thing wrong before)
  The figures were blank "[Figure not uploaded yet]" boxes because the figure
  PDFs were not in the Overleaf project. Two ways to fix:

  EASIEST: Overleaf -> New Project -> Upload Project -> choose
           "Paper 1 Revision.zip".  main.tex + all figures arrive together.
           Menu -> Compiler -> pdfLaTeX -> Recompile. Figures render.

  PASTE ROUTE: create a Blank project, paste main.tex, then click Upload and
           drag the THREE figure PDFs in. Recompile. (main.tex now finds the
           figures in the project root -- no folders needed.)

WHAT WAS FIXED IN THIS VERSION
  1. Figures: plain \includegraphics + robust \graphicspath
     ({figures/}{./}{../figures/}) so the PDFs are found whether they sit in
     the root or a figures/ folder. No more placeholder boxes.
  2. References: added the missing methodology citations for the sensitivity
     experiments (all verified real):
       - Furtula, Gutman, Dehmer, Appl. Math. Comput. 219 (2013) 8973-8978
         -> structure sensitivity, cited in Section 6.3
       - Rakic, Furtula, J. Chemometrics 33 (2019) e3138
         -> the fingerprint SS method, cited in Section 6.3
       - Konstantinova, J. Chem. Inf. Comput. Sci. 36 (1996) 54-57
         -> degeneracy / discrimination, cited in Section 6.2
     The diminished Sombor (DSO) paper is now also cited in Section 6 as the
     model for the benchmark (previously only in the intro).
  3. Citation integrity verified: 29 \cite keys all resolve to a \bibitem;
     no undefined and no uncited references. Bibliography now has 29 entries.

STRUCTURE OF THE PAPER  (unchanged math, kept verbatim from your source)
  Sec 1-3  intro, preliminaries, closed forms, all bounds (your text + rename).
  Sec 4    NEW Theorem 13: descriptor-space saturation (span = col(M), exact
           rank, <= 10 on molecular graphs).
  Sec 5    NEW Theorem 14: gamma-driven extremal switch (path at gamma=0 ->
           star for large gamma); this justifies "Extremal Graph
           Characterizations" in the title.
  Sec 6    sensitivity: only the LO reductions + pre-specified pure-gamma
           points; degeneracy definition corrected; honest framing.
  Sec 7    conclusion (cites both new theorems).

NOTE
  I cannot upload to your Overleaf account from here -- I can only prepare the
  files. You re-upload using one of the two routes above.
