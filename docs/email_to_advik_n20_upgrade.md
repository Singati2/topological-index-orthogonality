# Email draft to Advik — recent Paper 1 update (n=20 enumeration + structural conjecture)

Subject: Paper 1 update — extended enumeration to n=20 reveals a
sharp structural pattern at (−1, −1, 1)

---

Hi Advik,

Quick update on Paper 1 ("The Wuzi Index Family"). I've pushed
several rounds of revisions over the last few days; the repo is
now at tag `v1.2-paper1-wuzi` on `main`. Sharing the substantive
changes before the next round goes to sir.

**The headline mathematical upgrade.**
I extended the exhaustive enumeration of the (−1, −1, 1) extremal
regime from $n \le 12$ to $n \le 20$ (823,065 non-isomorphic trees
at $n = 20$; total ~1.3 million trees enumerated; runs end-to-end
in about 30 seconds on a single core). The pattern is much sharper
than the n ≤ 12 data alone suggested:

- For odd $n \in \{9, 11, 13, 15, 17, 19\}$: the unique maximizer
  is the *spider* $S^{(2)}_{(n-1)/2}$ — single hub of degree
  $(n-1)/2$, every spoke a length-2 path.
- For even $n \in \{10, 12, 14, 16, 18, 20\}$: the unique
  maximizer is the *balanced double-spider* $DS(a, b)$ with
  $a + b = n/2 + 1$ and $|a - b| \le 1$; equivalently
  $a = \lceil (n+2)/4 \rceil$, $b = \lfloor (n+2)/4 \rfloor$.

Every one of the 12 enumerated cases at $n \ge 9$ matches the
closed form. The earlier draft framed the conjecture *negatively*
("maximizer is neither path nor star nor caterpillar"), which
several harsh reviewers correctly attacked as too weak — the
upgraded Conjecture in §7.3 is *positive structural*: it specifies
the exact graph family at each $n$ and gives a concrete
transformation-lemma proof target. The case $n = 8$ is genuinely
exceptional (spider with one length-3 spoke); $n = 7$ fits the
spider formula at the boundary ($h = 3$). Table 2 now spans
$n = 5..20$ with both degree sequences and canonical `graph6`
strings for unique tree identification, and there is a new
Figure 10 visualising the maximizers at $n = 10, 11, 12$ (drawn
as rooted trees with the highest-degree vertex at the top).

**Other substantive changes this round.**

- *Bootstrap CIs on Table 3.* $B = 10{,}000$ percentile bootstrap
  on the octane correlations (seed fixed for reproducibility).
  The caption now reports headline CIs — e.g.\
  $r(\Delta H_{\mathrm{vap}}, GA) = +0.959$ with 95% CI
  $[+0.917, +0.986]$, and the 0.008 margin between $W(0, 0, 2)$
  and $GA / H$ on dHvap is well within the 0.034 CI half-width.
  The "competitive on $\Delta H_{\mathrm{vap}}$ and $\omega$"
  framing has been softened to "comparable in magnitude"; the
  $\omega$ "win" at $|r| = 0.987$ is the classical
  $W(1, 0, 0) = M_2$ reduction, so we no longer claim it as a
  Wuzi-specific result.

- *Section 4 bounds, language tightening.* The bracket "Theorem"
  was demoted to a "Proposition" — it really is a one-line
  consequence of min/max applied edgewise, so calling it a
  Theorem was overclaim. The upper bound in $\Delta$ is now
  explicitly labelled as a coarse product bound (per-factor
  maxima are not simultaneously attained), and the earlier
  "tight precisely on regular graphs" wording was relaxed to
  acknowledge that certain biregular and semiregular graphs can
  also attain the bracket equality depending on the parameter
  region. The math is unchanged; only the labelling.

- *Bibliography corrections (please double-check three of these).*
  A verification pass on the references found that three earlier
  citations did not actually point to existing works under the
  titles we'd given them. They have been replaced with verified
  alternatives that serve the same textual role:
  - "Gutman & Borovićanin 2018 BID survey" was retired in favour
    of **B. Borovićanin, K. C. Das, B. Furtula, I. Gutman,
    "Bounds for Zagreb indices", MATCH Commun. Math. Comput. Chem.
    78 (2017) 17–100** (also reprinted as a chapter in
    MCM-20, Univ. Kragujevac, 2017, pp.\ 67–153).
  - "Gutman & Furtula 2017 Randić book" was retired in favour of
    **I. Gutman & B. Furtula (eds.), "Recent Results in the
    Theory of Randić Index", Mathematical Chemistry Monographs
    MCM-6, Univ. Kragujevac, 2008**.
  - "Hansen & Mélot 2005 bicyclic Randić" was retired in favour
    of **H. Liu, M. Lu, F. Tian, "On the Randić index",
    J. Math. Chem. 44 (2008) 301–310** (DOI
    10.1007/s10910-007-9311-1).
  These are *my* identifications based on the role each citation
  plays in the text; please skim the three replacements when you
  read the manuscript and let me know if any of them is not the
  reference you would have cited. The two other corrected entries
  (Movahedi 2025 → arXiv:2508.06532; Todeschini–Consonni 2009
  full details filled in) are easier to verify and I'm confident
  those are right.

- *Section 7 (Open problems and conjectures).* Placement is
  unchanged. The bicyclic statement is now framed as a "Remark /
  Open question" rather than a formal Conjecture, because we
  have no computational evidence for it (the unicyclic and tree
  conjectures have analogue-from-classical-Randić support; the
  bicyclic one didn't, so it was honest to demote it).

**Reproducibility / infrastructure.**
Repo is tagged `v1.2-paper1-wuzi`. The `releases/paper1_wuzi/MANIFEST.md`
file lists every script, dataset, results CSV, and figure that
Paper 1 depends on, together with the reproduction recipe. Two new
scripts were committed so the n=5..20 enumeration and the
bootstrap CIs are reproducible end-to-end from a clean checkout:

- `scripts/13_caterpillar_max_n5_to_20.py` regenerates Table 2.
- `scripts/14_octane_bootstrap.py` regenerates the bootstrap CIs
  used in the Table 3 caption.

All 124 pytest tests pass under the new state.

**What I think comes next.**
If you have time to read the manuscript before sir's physical-print
review, the three highest-yield places to focus are:

1. The upgraded Conjecture `conj:caterpillar_max` in §7.3 — does
   the spider / double-spider closed-form statement read correctly
   to you, and is the boundary-case Remark (n = 6, 7, 8)
   precise enough?
2. The three replaced bibliography entries above — please confirm
   each is the reference you would have cited.
3. The bootstrap CI footnote on Table 3 and the rewritten
   interpretation paragraph in Section 6 — does the new tone
   ("comparable in magnitude", "statistically indistinguishable
   on n = 18") read as appropriately hedged?

If you'd prefer to wait and read after sir's review, that is also
fine — the manuscript is in what I would call defensibly
MATCH-ready shape.

Compiles from `docs/paper1_wuzi.tex` on Overleaf without extra
setup; if you want a pre-compiled zip with the figures path
already adjusted for the flat Overleaf layout, I have one ready
on my Desktop and can share.

Let me know how you'd like to proceed.

Best,
Ganesh
