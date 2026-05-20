# Final 10/10 Upgrade Plan

**Compiled:** 2026-05-20 after the round-10 final-audit sweep.
**State of repo:** branch `main`; baseline commit before this round was `58a57ac`; this plan lives in commit `a757e27` together with the round-10 fixes; round-12 meta-doc updates appended in the next commit.

**Important note on Proposition / Conjecture numbering.** Paper 1 uses a per-section shared counter (`[theorem]`), so the propositions added in round 10 receive numbers different from what an earlier draft of this plan stated. The actual compiled numbers are:

- Proposition **4.7** (upper bound in $\Delta$, `thm:upper_Delta`).
- Proposition **4.8** (lower bound in $\delta$, `thm:lower_delta`).
- Conjecture **4.9** (Nordhaus--Gaddum, `conj:NG`).
- Proposition **5.4** (bound via generalized $R$ and $\chi$, `thm:R_chi`).
- Proposition **5.5** (indirect bound via Sombor, `thm:sombor`).
- Conjecture **5.6** (sharp bound via $GA$, `thm:GA`).
- Conjecture **5.7** (sharp bounds via $H$, $ABC$, `thm:H_ABC`).
- Conjecture **6.4** (non-classical extremal at $(-1,-1,1)$, `conj:caterpillar_max`).

This document tells (a) Arockiaraj sir what to act on, (b) Ganesh what to do before arXiv, and (c) anyone reading the repo cold what is still missing for a clean journal submission.

---

## 1. What is now fixed

Round-10 work, driven by the final Woz + Jobs + figure + cross-number audits:

**Paper 1 (`docs/paper1_wuzi.tex`).**
1. The previous "Theorem (upper bound in $\Delta$)" with a `\todo{PROOF SKETCH}` marker is now **Proposition 4.7** (`thm:upper_Delta`) with an explicit termwise proof in the non-negative parameter region.
2. The previous "Theorem (lower bound in $\delta$)" is now **Proposition 4.8** (`thm:lower_delta`) with explicit proof and equality condition.
3. The previous "Theorem (R/χ bound) -- skeleton" with a `[PROOF SKETCH -- TO BE COMPLETED]` marker is now **Proposition 5.4** (`thm:R_chi`) with a proper `\begin{proof}` block, the explicit $(2\alpha,0,0)/(0,2\beta,0)$ Cauchy--Schwarz split, and the equality condition inherited from `thm:cs`. (Round-6 silently fixed a spurious "$+1$" exponent in this bound; that fix is preserved.)
4. The previous "Theorem (Sombor) -- skeleton" had a wrong-direction Jensen citation in the round-8 pass; corrected to **Proposition 5.5** (`thm:sombor`) valid for $0 < \beta < 1$ using the concave Jensen `thm:m1_concave`, giving $W(G; 0, \beta, 0) \le m^{1-\beta} (\sqrt{2}\, SO(G))^{\beta}$.
5. The previous "Theorem (GA) -- skeleton" and "Theorem (H/ABC) -- skeleton" stated no actual bound; converted to **Conjectures 5.6 and 5.7** with honest framing ("are open"). (The legacy label names `thm:GA` and `thm:H_ABC` are retained for backward-compatibility of cross-references; the environments are correctly typeset as Conjectures.)
6. Conjecture 6.3 (bicyclic) carried an in-body `\todo{TO BE CHARACTERIZED}` marker; rewritten to name two candidate families (theta-graph with pendants; bowtie-with-pendants) and to honestly state that a small-$n$ computational verification (e.g. via `nauty`/`geng`) is still open.
7. Remark 5.3 reworded to make the per-$G$ sharp / $\mathcal{G}_4$-wide uniform distinction explicit, matching what `results/wuzi_bounds_ratio_tables.csv` actually computes.
8. §8 octane wording at lines ~1138-43 corrected: "two of five properties" via $W(0, 0, 2)$ ($\gamma \ne 0$) and $W(1, 0, 0) = M_2$ ($\gamma = 0$), not "at parameter triples that include a non-trivial $\gamma$-axis contribution" (which is true of only one of the two).
9. Contribution list (4) rewritten to enumerate the now-rigorous Propositions explicitly; GA and H/ABC honestly named as Conjectures.
10. §9 Future-work paragraph rewritten — no more "[PROOF SKETCH]" reference.

**Paper 2 (`docs/paper2_orthogonality_screening.tex`).**
1. Abstract listed three feature configurations but Table 7.1 has four; corrected to "four (full, pca_95, pairwise_pruned, combined_pruned)".
2. Abstract gained an explicit "single dataset (ESOL, $\log S$); we do not generalize" caveat on the alt-family screen.
3. §6 candidate-disclaimer `\todo{[REFERENCES TO BE FILLED]}` replaced with a prose pointer to `docs/literature_notes.md`.
4. §7.1 "0.04--0.05 higher" softened to "up to 0.04 higher" (matches the actual range across regression datasets).
5. Missing `\bibitem{Martins2012BBBP}` added (Martins et al. 2012, *J. Chem. Inf. Model.*) since BBBP was textually cited but had no bibliography entry.

**Figure documentation.**
1. `docs/figure_captions.md` Figure 4 caption corrected: "34 indices / 5 properties / less the constant $PI$" (was wrongly "27 / 8 / less three").
2. `docs/figure_audit.md` fig4 source corrected to `results/octane_prediction.csv` (was `octane_descriptors.csv`; the script actually reads the correlation table).

**Results.**
1. `results/wuzi_extremal_trees.md` line 102 "TODO" reframed to "future work" with explicit acknowledgement that Paper 1's bicyclic conjecture does not depend on this computation.

**Cross-doc.**
1. Paper 1 keyword line dropped "QSPR" (Paper 2's lexicon); Paper 1 keywords now: topological index, BID index, Randić, Sombor, Wuzi family, sensitivity analysis, extremal graph.
2. README citation block has both manuscripts as proper `@unpublished` BibTeX entries.
3. README "What remains unfinished" rewritten to match current Paper 1 §4–§6 state and the arXiv-first submission plan.
4. `CITATION.cff` is well-formed YAML with project-level title, 8 keywords, two-paragraph abstract covering both papers.
5. Repo description and 6 topics set on GitHub.

After this commit: **Paper 1 has zero `[PROOF SKETCH -- TO BE COMPLETED]` markers inside theorem environments.** Only 4 `\todo{PROOF STRATEGY}` markers remain, all inside Conjecture bodies or between conjectures and the next subsection — these are visible-in-PDF outlines for the senior-author proof work, not hidden gaps.

---

## 2. What still prevents Paper 1 from being 10/10

In rough priority order:

1. **Seven Conjectures remain open:** 4.9 (Nordhaus--Gaddum, `conj:NG`), 5.6 (GA, `thm:GA`), 5.7 (H/ABC, `thm:H_ABC`), 6.1 (trees positive region, `conj:trees_beta_pos`), 6.2 (unicyclic, `conj:unicyclic`), 6.3 (bicyclic, `conj:bicyclic`), 6.4 (caterpillar_max at $(-1,-1,1)$, `conj:caterpillar_max`). This is the structural ceiling for a MATCH or J. Math. Chem. submission. Each can be closed by Arockiaraj sir + Advik via the transformation-argument template of Movahedi--Gutman--Redžepović--Furtula 2026 (cited).
2. **5 bibliography entries carry "to be added at submission" venue placeholders** (not missing entries — they compile and resolve cleanly via cite/bibitem, but the venue/page/identifier fields are still placeholder text): `Movahedi2025arXivDSO` (arXiv id), `GutmanBorovicanin2018BID` (full survey citation), `GutmanFurtulaBook2017` (publisher), `HansenMelot2003unicyclicRandic` and `HansenMelot2005bicyclicBounds` (venues). All resolvable in 30 minutes of literature lookup.
3. **No comparison to closely related parametric BID families** (general Randić $R_\alpha$, sum-connectivity $\chi_\beta$, parametric Sombor variants, "Misc-Sombor" 2024-2025). The ratio-bound theorem makes a *uniform* statement against any classical BID index but does not explicitly position Wuzi vs other parametric extensions.
4. **`\todo{PROOF STRATEGY}` markers in Conjecture bodies will render as red bold text.** This is intentional for advisor draft; remove the `\renewcommand{\todo}{...}` in the preamble (swap to `\newcommand{\todo}[1]{}`) for the submission-ready compile.
5. **Closed-form bounds involving the Albertson irregularity $\mathrm{Alb}(G) = \sum_e |d_u - d_v|$** are listed as future work (§9 item 4) but not even stated as conjectures. If Arockiaraj sir likes the $\gamma$-axis framing, this is a natural mini-section to add.

---

## 3. What still prevents Paper 2 from being 10/10

In rough priority order:

1. **One citation entry with placeholder venue text: `GutmanBorovicanin2018BID`** — same as in Paper 1; fill the venue/page fields once.
2. **No reproducibility command block.** For *J. Cheminformatics*, the methodology section should name the exact script (`scripts/11_ml_benchmark.py`) and the expected output. Add a short "Reproducibility" block at the end of §4.
3. **VIF defined but never reported in body.** §4.3 defines VIF; `results/{esol,freesolv,lipo}/vif.csv` contain the numbers. Add a one-sentence summary in §5.1 (e.g., "Of the 30 baseline indices, $k$ have $\mathrm{VIF} = \infty$ on ESOL").
4. **The alternative-family screen is ESOL-only.** The abstract now flags this, but extending to one additional dataset (e.g., FreeSolv) would close the most likely reviewer attack at J. Cheminformatics.
5. **No formal independence proof of the BID-basis 10-dim claim.** The repository contains a sharp test (`tests/test_bid_basis.py::test_bid_dimension_bound_on_delta_4`) but Paper 2 §3 (Theorem 3.1 in Paper 2's numbering) restates the observation without a separate proof. Either link explicitly to Paper 1's Theorem 3.1 or move the one-line linear-algebra proof into Paper 2's body (one paragraph).

---

## 4. Exact actions Arockiaraj sir must take

1. **Look at Conjecture 6.4 (caterpillar_max) in `docs/paper1_wuzi.tex` §6.4** and decide if the non-classical extremal at $(-1,-1,1)$ for $n \ge 7$ is characterizable. If yes, Paper 1 has a publishable new structural theorem. If no, the conjecture stays as the headline empirical observation. The computational evidence is in `results/wuzi_extremal_trees.md` (degree sequences for $n = 5..12$).
2. **Decide whether to co-author Paper 1 §4, §5, §6 conjectures** (Nordhaus--Gaddum 4.11, GA 5.8, H/ABC 5.9, trees-positive 6.1, unicyclic 6.2, bicyclic 6.3, caterpillar 6.4), or to ship Paper 1 with them as labelled conjectures + companion future-work paper. Either choice unblocks submission.
3. **Confirm the affiliation line** "Department of Mathematics, Loyola College, Chennai" (Advik + sir) and "Department of Mathematics and Statistics, Florida Atlantic University" (Ganesh) are correct as they appear in both `.tex` files and `CITATION.cff`.

---

## 5. Exact actions Ganesh must take before arXiv

1. **Fill the 5 Paper 1 bibliography stubs** (~30 minutes literature lookup).
2. **Fill the 1 Paper 2 bibliography stub** (`GutmanBorovicanin2018BID`).
3. **Add a Reproducibility block at the end of Paper 2 §4** naming `scripts/01_baseline_correlations.py`, `scripts/02_wuzi_grid_search.py`, `scripts/11_ml_benchmark.py`, `scripts/08_make_figures.py`.
4. **Add one sentence to Paper 2 §5.1 reporting VIF results.**
5. **Compile both papers on Overleaf, visually inspect the rendered PDF**, and confirm the `\todo{PROOF STRATEGY}` red markers in Paper 1 conjectures appear where expected (for advisor read). For the arXiv-submission compile, swap `\newcommand{\todo}[1]{...visible...}` to `\newcommand{\todo}[1]{}` in the preamble.
6. **Tag a `v0.1.0` release** on GitHub so `CITATION.cff`'s `version: 0.1.0` resolves to a real release.

---

## 6. Target venue recommendation

**Paper 1.** Realistic primary target: *MATCH Commun. Math. Comput. Chem.* — the paper sits squarely in the MATCH tradition (parametric BID family + closed forms + bounds + conjectures + sensitivity analysis on octanes/decanes), explicitly templates from Movahedi--Gutman--Redžepović--Furtula 2026 (MATCH 95:141-162), and is structurally complete with rigorous Theorems 4.1/4.4/4.5/4.6/4.7/4.8/5.5, Propositions 4.4/4.5/5.6/5.7, and clearly-labelled Conjectures 5.8/5.9/conj:NG/6.1/6.2/6.3/6.4. *J. Math. Chem.* is plausible. *AKCE Int. J. Graphs Comb.* and *SAR & QSAR Env. Res.* are safer fall-throughs.

**Paper 2.** Realistic primary target: *Molecular Informatics* (Wiley) — methodology paper with honest empirical findings, four-dataset RandomForest benchmark, explicit failure-mode reporting, accepting MIT-licensed open-source companion. *J. Cheminformatics* (BMC) is reachable once the reproducibility block + cross-dataset alt-family extension land. *J. Chem. Inf. Model.* is a stretch without a stronger predictive payoff.

---

## 7. Honest one-paragraph pitch for Paper 1

*A three-parameter parametric topological index family, the Wuzi index family, is introduced for a simple connected graph and parametrically interpolates between the second Zagreb $M_2$, the Randić $R$, the sum-connectivity $\chi$, the harmonic $H$, and a new exponential imbalance factor controlled by $\gamma$. Closed-form values are derived on nine standard graph classes. The central rigorous result is a uniform ratio-bound theorem: for every classical bond-incident-degree index $I_h$ with strictly positive edge contribution, the Wuzi family is sandwiched between $c_{\min}^h \cdot I_h$ and $c_{\max}^h \cdot I_h$ with explicit equality conditions. Companion results give a Cauchy--Schwarz / geometric-mean inequality, a Jensen-type bound via the first Zagreb $M_1$, and explicit upper and lower bounds in $(\Delta, m)$ and $(\delta, m)$ for the non-negative parameter region. Computational enumeration of all non-isomorphic trees of order $5 \le n \le 12$ identifies a non-classical extremal regime at $(\alpha, \beta, \gamma) = (-1, -1, 1)$: the maximum is attained for $n \ge 7$ by a tree that is neither a path, a star, nor a caterpillar (Conjecture 6.4, computationally supported), whose formal structural characterization is the central open problem. A sensitivity analysis on three classical chemical compound classes ($18$ octanes, $106$ trees of order $10$, $75$ decane isomers) shows the Wuzi family is in the same correlation regime as established BID indices and is competitive on two of the five octane properties.*

## 8. Honest one-paragraph pitch for Paper 2

*Across twenty alternative-family candidate topological indices screened on the ESOL solubility benchmark, fifteen clear the conventional pairwise $|r| < 0.95$ orthogonality threshold but only four also clear a partial-correlation threshold $|\mathrm{pcor}(z, y \mid X)| \ge 0.10$ once the standard $30$-index baseline is regressed out: eleven candidates that "look orthogonal" in the pairwise sense carry no measurable QSPR signal. We present an open-source, reproducible orthogonality-screening pipeline combining four diagnostics — pairwise correlation, principal-component effective rank, partial correlation with the target, and variance inflation factor — into a single combined criterion. A structural observation, the Edge-Degree-Pair Basis, gives a $10$-dimensional ceiling for the bond-incident-degree (BID) family on hydrogen-suppressed molecular graphs with maximum degree $\Delta \le 4$. Replication across three independent QSPR benchmarks (ESOL, FreeSolv, Lipophilicity; total $5{,}966$ molecules) gives a stable baseline-redundancy pattern: approximately $37\%$ of the $\binom{30}{2} = 435$ pairwise correlations exceed $|r| \ge 0.95$ on every dataset; three principal components capture $95\%$ of variance. A five-fold cross-validated RandomForest benchmark across all four datasets (including BBBP classification) confirms that the pairwise-pruned screen matches the full $30$-index baseline within CV noise with roughly four-fold fewer features, while exposing an honest failure mode of the strict combined screen on Lipophilicity (zero features retained in 4 of 5 folds). Orthogonality with the classical baseline is positioned as a necessary but not sufficient condition for QSPR usefulness.*

---

## 9. Do not claim

- **Do not** claim Paper 1 proves the non-classical $(-1, -1, 1)$ maximum is some specific tree class. It is a computational observation at $n \le 12$.
- **Do not** claim Paper 1 has sharp bounds via $GA$, $H$, or $ABC$. These are Conjectures.
- **Do not** claim Paper 1 has formal extremal-graph characterizations for trees / unicyclic / bicyclic graphs. These are Conjectures with proof strategies.
- **Do not** claim Paper 1 has Nordhaus--Gaddum bounds. This is a Conjecture.
- **Do not** claim Paper 2 demonstrates that orthogonality screening improves QSPR performance. The combined screen *matches* the full baseline within CV noise on three of four datasets and *underperforms* on Lipophilicity (the honest failure mode). The contribution is feature reduction, not predictive improvement.
- **Do not** claim the BBBP ROC-AUC delta of $0.014$ is a real performance gain — it is "within CV noise", with overlapping std bands. Paper 2 hedges this correctly; do not undo the hedge.
- **Do not** claim Paper 2's alternative-family findings (15/20 vs 4/20) generalize to other QSPR targets. They are ESOL-only and the manuscript flags this.
- **Do not** describe the Wuzi family as a "superior" or "novel" topological index. It is a parametric BID extension with explicit reductions to classical indices; the contribution is the family + bounds + extremal study, not the claim of superiority.

---

## 10. Ready to send?

**To Arockiaraj sir as an advisor draft: YES.** Both manuscripts are internally consistent, honestly hedge proven vs conjectured vs computational results, are LaTeX-clean (every \ref/\cite/\bibitem resolves; env balance 84/84 in Paper 1, 34/34 in Paper 2), have 124/124 local tests passing, and have 6 consecutive green CI runs. The advisor-ready summary `docs/advisor_ready_summary.md` lists three concrete actions for sir to take. Send.

**To a journal referee: NO, not yet.** Action items in §5 above plus the §2 / §3 Paper-1 / Paper-2 gap closures are needed first. After those, *MATCH* (Paper 1) and *Molecular Informatics* (Paper 2) are realistic primary targets. Estimated time from advisor sign-off to journal-ready: 1-2 weeks of focused author work, longer if Arockiaraj sir chooses to close the §4-§6 conjectures rigorously.
