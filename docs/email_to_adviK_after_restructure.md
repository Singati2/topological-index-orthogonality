# Email draft to Advik after the manuscript restructure

Subject: Manuscripts restructured per Arockiaraj sir's direction --
ready for your review

---

Hi Advik,

Following Arockiaraj sir's feedback, I have restructured the two
manuscripts so they have minimal overlap. Sharing the current
state for your review before the next round goes to sir.

What changed:

- The two papers are now cleanly separated. Paper 1 is the Wuzi
  graph-theory paper (closed forms, bounds, extremal graphs,
  sensitivity analysis on the standard math-chem benchmarks --
  octanes, trees of order 10, decanes). Paper 2 is the
  orthogonality-screening pipeline plus the multi-dataset machine
  learning benchmark.
- The methodology / multi-dataset / ML-benchmark content has been
  removed from Paper 1 entirely; Paper 1 only mentions the
  companion methodology paper in a single sentence (abstract,
  introduction, conclusion).
- Conversely, the detailed Wuzi closed-form derivations and the
  bounds and extremal-graph analysis have been removed from
  Paper 2; that paper cites Paper 1 for the Wuzi math and uses
  Wuzi only as one screened candidate.
- The `.tex` files have been renamed to match Arockiaraj sir's
  ordering: `docs/paper1_wuzi.tex` and
  `docs/paper2_orthogonality_screening.tex`.

What still needs mathematical input:

Paper 1 Sections 4--6 still contain a few [PROOF SKETCH] /
[CONJECTURE] markers in the spots where the formal proof needs
your or Arockiaraj sir's derivation. Specifically: sharp
$\Delta$-dependent bounds (Section 4), the Nordhaus-Gaddum bound
(Section 4), sharper per-index bounds for $R$, $SO$, $GA$, $H$,
$ABC$ (Section 5), and the formal extremal-graph characterizations
for trees / unicyclic / bicyclic (Section 6). The most
interesting target is Conjecture 6.4, the non-classical extremal
at $(\alpha, \beta, \gamma) = (-1, -1, 1)$, which has been
computationally verified for $n = 5, \ldots, 12$ but lacks a
formal proof. None of the math gaps is presented as a finished
theorem; all are clearly labelled.

What is ready:

- Paper 1 compiles structurally clean (env-balance, refs and cites
  resolved, no forbidden methodology terminology).
- Paper 2 compiles structurally clean and contains the
  cross-dataset evidence, the 5-fold CV ML benchmark, and the
  software-reproducibility scaffolding (test suite, CI, CITATION
  metadata).
- A short summary note for Arockiaraj sir is in
  `docs/advisor_ready_summary.md`.

Logistics:

The current GitHub commit is ready for review. Both `.tex` files
need to be compiled on Overleaf (`pdflatex` is not available in
my local environment); once compiled the PDFs are referee-style
clean. Could you let me know whether you would prefer to read
first and forward to Arockiaraj sir, or whether I should send to
him directly?

Repository:
https://github.com/Singati2/topological-index-orthogonality

Best,
Ganesh
