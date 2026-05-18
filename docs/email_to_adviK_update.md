# Email draft to Advik — project update + two-paper split

**Subject:** Project update — Phase 2 scaffolding pushed, two-paper plan, what we need

---

Hey Advik,

Substantial progress to share since my last note, plus a strategic
recommendation that I'd like your read on before we move further.

**Where the work currently stands** (all on GitHub at
<https://github.com/Singati2/topological-index-orthogonality>):

- Full Wuzi parametric implementation, plus closed-form values verified
  numerically on K_n, C_n, P_n, S_n, K_{p,q}, Q_k, W_n, F_p, and general
  k-regular graphs (`src/wuzi_analytical.py`). These are the
  Proposition-1 analog of the Movahedi/Gutman/Redžepović/Furtula 2026
  MATCH DSO paper you sent.
- Edge-contribution function ψ(x, y; α, β, γ) analyzed in
  `docs/edge_contribution_analysis.md` (monotonicity, partial
  derivatives, extremal pairs on the {1, ..., 4}^2 lattice, regular-graph
  degeneracy of γ). This is our Eqs. (3)-(4) analog from MATCH 95:141-162.
- Edge-Degree-Pair Basis theorem (`docs/theoretical_foundation.md`): for
  Δ ≤ 4 molecular graphs, every endpoint-degree edge-sum index lies in
  a 10-dimensional space spanned by edge-degree-pair counts m_ij. We
  attribute this to the BID-index framework that Gutman, Borovićanin,
  Furtula and others developed, and only claim the explicit dimension
  bound + redundancy consequence + screening operationalization as
  novel emphasis.
- Wuzi 100-point grid sweep on ESOL (n=1127), FreeSolv (n=639), and
  Lipophilicity (n=4200). Cross-dataset summary in
  `results/cross_dataset_summary.md`. Effective rank of the 30-index
  baseline is consistently 3 across all three. Wuzi fails the
  orthogonality kill-test at all 100 grid points on each dataset (one
  marginal pass on FreeSolv at (-1, -1, 1) with max |r| = 0.943 and
  partial correlation = 0).
- Octane prediction (Section 5.1 of MATCH DSO analog) on the 18 octane
  isomers using NIST WebBook physicochemical properties: in
  `results/octane_prediction.md`.
- Degeneracy on the 106 non-isomorphic trees of order 10 (Section 5.3
  analog): in `results/wuzi_degeneracy.md`. Wuzi at γ = 2 reaches ~25 %
  degeneracy, comparable to the Sombor variants reported in MATCH
  95:141-162 Figure 5.
- Structure sensitivity on the 75 decane isomers (Section 5.4 analog):
  in `results/structure_sensitivity.md`. Wuzi at (0, 0, 2) gives
  SS ≈ 0.165, comparable to the Sombor-family values reported in
  MATCH 95:141-162 Table 2.

**Strategic recommendation — two papers, Paper 2 first**

I now think the work should ship as two papers rather than one:

**Paper 2 (write first):** *"The Wuzi Index Family: Graph-Theoretic
Properties, Bounds, Extremal Graphs, and Redundancy Analysis"* —
follows the standard math-chem template (your MATCH PDF, exactly).
Definition + special-graph values + edge-contribution analysis +
bounds + extremal graphs, then a Section 6 (numerical) that includes
the redundancy / orthogonality result on real chemistry as part of
the family's structural honest evaluation. Realistic venues:
J. Math. Chem., MATCH, AKCE, SAR & QSAR Env. Res.

**Paper 1 (write second):** *"Orthogonality Screening of Topological
Indices for QSPR Modeling: A Structural and Empirical Redundancy
Analysis"* — the broader methodology / open-source pipeline paper that
cites Paper 2 as the worked example. More ambitious target (J.
Cheminformatics or J. Chem. Inf. Model.); depends on Paper 2 going
through cleanly first.

The reason I want Paper 2 first is exactly what you said in your last
note — that's the math-chem norm: structural / theoretical paper before
the methodology paper. It also lets us settle the Wuzi family's
publishability question on its own terms before claiming the broader
screening contribution.

**One condition I'd like us to commit to up front:** Paper 2 has to
include the redundancy result inside its Section 6. If we ship a pure
"new index" paper without addressing what our own data shows about
orthogonality, a careful reviewer will catch us and the result is worse
than including it constructively. The MATCH 95:141-162 Section 5.2
intercorrelation observation (r = 0.92 – 0.99 on octanes) is our
template here — we extend it formally to three datasets with PCA and
partial correlation. We are not attacking that paper or any other
existing work — we are building on their analytical move.

**Three things I need from you to proceed**

1. **Your GitHub handle.** I'll add you as a collaborator with push
   access. If you don't have a GitHub account, take ten minutes to
   create one; the username `advik-natarajan` or similar is fine.

2. **Co-development of Sections 3, 4, 5 of Paper 2 (the math) with
   Arockiaraj sir.** Specifically:

   - **Section 3** — bounds W(G; α, β, γ) in terms of n, m, Δ, δ with
     equality conditions. Template: Section 3 of MATCH 95:141-162.
   - **Section 4** — bounds W in terms of classical indices M_1, M_2,
     R, SO, GA, H, ABC. Template: Movahedi 2025 arXiv preprint.
   - **Section 5** — extremal graphs among trees, unicyclic, bicyclic.
     Template: Section 4 of MATCH 95:141-162.

   Each is approximately 5-10 theorems with proofs. I've prepared a
   manuscript skeleton (`docs/paper2_wuzi_manuscript_skeleton.md`)
   where the math sections drop in directly; all numerical content is
   already filled in.

3. **Send the advisor update to Arockiaraj sir.** I've prepared a
   2-page advisor-facing version at `docs/advisor_update.md`. Feel
   free to forward it, modify it, or use it as talking points for an
   in-person conversation. The tone is respectful and asks for his
   collaboration; nothing is presented as a fait accompli.

**Timeline if everyone is reasonably available**

Paper 2 sections 3-5 derived in 3-4 weeks → manuscript polished in
2 weeks → first submission in ~6 weeks. Paper 1 drafted in parallel,
first submission in ~10 weeks.

Take your time on any of this. I'd rather we land it cleanly than fast.

Best,
Ganesh
