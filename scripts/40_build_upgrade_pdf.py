#!/usr/bin/env python3
"""Build the upgraded Paper-1 (Loyola, probe-first) as a downloadable PDF.
reportlab for layout + matplotlib mathtext for typeset display equations.
No LaTeX engine required. Output: ~/Downloads/Paper1_Loyola_Upgraded.pdf
"""
import os, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import inch
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_JUSTIFY, TA_CENTER, TA_LEFT
from reportlab.lib import colors
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Image,
                                Table, TableStyle, HRFlowable, KeepTogether)
from reportlab.lib.utils import ImageReader
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

EQDIR = "/Users/ganeshshiwakoti/topological_index_pilot/figures/_eq"
os.makedirs(EQDIR, exist_ok=True)
OUT = "/Users/ganeshshiwakoti/Downloads/Paper1_Loyola_Upgraded.pdf"

# ---- fonts (DejaVu has Greek + math glyphs for inline unicode) ----
mpl_ttf = os.path.join(matplotlib.get_data_path(), "fonts", "ttf")
pdfmetrics.registerFont(TTFont("DJ", os.path.join(mpl_ttf, "DejaVuSans.ttf")))
pdfmetrics.registerFont(TTFont("DJ-B", os.path.join(mpl_ttf, "DejaVuSans-Bold.ttf")))
pdfmetrics.registerFont(TTFont("DJ-I", os.path.join(mpl_ttf, "DejaVuSans-Oblique.ttf")))
pdfmetrics.registerFont(TTFont("DJ-BI", os.path.join(mpl_ttf, "DejaVuSans-BoldOblique.ttf")))
pdfmetrics.registerFontFamily("DJ", normal="DJ", bold="DJ-B", italic="DJ-I", boldItalic="DJ-BI")

styles = getSampleStyleSheet()
BODY = ParagraphStyle("body", parent=styles["Normal"], fontName="DJ", fontSize=9.5,
                      leading=13.5, alignment=TA_JUSTIFY, spaceAfter=5)
H1 = ParagraphStyle("h1", parent=styles["Heading1"], fontName="DJ-B", fontSize=14,
                    leading=17, spaceBefore=10, spaceAfter=5, textColor=colors.HexColor("#1a1a2e"))
H2 = ParagraphStyle("h2", parent=styles["Heading2"], fontName="DJ-B", fontSize=11.5,
                    leading=14, spaceBefore=8, spaceAfter=3, textColor=colors.HexColor("#16213e"))
TITLE = ParagraphStyle("title", parent=styles["Title"], fontName="DJ-B", fontSize=17,
                       leading=21, alignment=TA_CENTER, spaceAfter=4)
SUB = ParagraphStyle("sub", parent=BODY, alignment=TA_CENTER, fontSize=9, textColor=colors.grey)
THM = ParagraphStyle("thm", parent=BODY, fontName="DJ-I", backColor=colors.HexColor("#eef2fb"),
                     borderPadding=6, leftIndent=4, rightIndent=4, spaceBefore=4, spaceAfter=6)
PROOF = ParagraphStyle("proof", parent=BODY, fontSize=9, leading=12.5, leftIndent=8,
                       textColor=colors.HexColor("#222222"))
NOTE = ParagraphStyle("note", parent=BODY, fontSize=9, backColor=colors.HexColor("#fff6e6"),
                      borderPadding=5, spaceBefore=3, spaceAfter=5)

import re as _re
_eqn = [0]
FAILED = []   # equations that could not render as math (must stay empty)
def eq(latex, fs=15, maxw=430):
    """Render a display equation to PNG via mathtext; guaranteed to produce a file."""
    _eqn[0] += 1
    path = os.path.join(EQDIR, f"eq_{_eqn[0]:03d}.png")
    def _san(s):
        s = _re.sub(r"\\text\{([^}]*)\}",
                    lambda m: r"\mathrm{" + m.group(1).replace(" ", r"\ ") + "}", s)
        s = _re.sub(r"\\xrightarrow(\[[^\]]*\])?(\{[^}]*\})?", r"\\to", s)
        s = _re.sub(r"\\le(?![a-zA-Z])", r"\\leq", s)   # mathtext wants \leq not \le
        s = _re.sub(r"\\ge(?![a-zA-Z])", r"\\geq", s)   # ...and \geq not \ge (guard \left)
        for tok in (r"\Bigl", r"\Bigr", r"\bigl", r"\bigr", r"\Big", r"\big"):
            s = s.replace(tok, "")
        return (s.replace(r"\tfrac", r"\frac").replace(r"\mathbb", r"\mathbf")
                 .replace(r"\!", "").replace(r"\dim", r"\mathrm{dim}")
                 .replace(r"\qquad", r"\quad").replace(r"\longrightarrow", r"\to"))
    variants = [latex, _san(latex)]
    rendered = False
    for v in variants:
        try:
            fig = plt.figure(figsize=(0.01, 0.01))
            fig.text(0, 0, f"${v}$", fontsize=fs)
            fig.savefig(path, dpi=220, bbox_inches="tight", pad_inches=0.04, transparent=True)
            plt.close(fig); rendered = True; break
        except Exception:
            plt.close("all"); continue
    if not rendered:                      # plain-text fallback: always yields a file
        FAILED.append(latex)
        plain = _re.sub(r"[\\${}^_]", " ", latex)
        fig = plt.figure(figsize=(0.01, 0.01))
        fig.text(0, 0, plain, fontsize=fs-2)
        fig.savefig(path, dpi=220, bbox_inches="tight", pad_inches=0.04, transparent=True)
        plt.close(fig)
    iw, ih = ImageReader(path).getSize()
    sc = 72.0/220.0
    w, h = iw*sc, ih*sc
    if w > maxw:
        h *= maxw/w; w = maxw
    img = Image(path, width=w, height=h)
    img.hAlign = "CENTER"
    return img

def P(t, st=BODY): return Paragraph(t, st)
def H(t): return P(t, H1)
def h(t): return P(t, H2)
def sp(x=4): return Spacer(1, x)
def rule(): return HRFlowable(width="100%", thickness=0.6, color=colors.HexColor("#bbbbbb"),
                              spaceBefore=4, spaceAfter=6)

flow = []
A = flow.append

# ===================== TITLE =====================
A(P("Descriptor-Space Saturation of a Parametric Bond-Incident-Degree Family:", TITLE))
A(P("the Loyola Index as a Probe", TITLE))
A(sp(4))
A(P("Revised &amp; upgraded manuscript — probe-first reframe", SUB))
A(P("Advik Natarajan, Ganesh Shiwakoti, Micheal Arockiaraj", SUB))
A(P("<i>Prepared as a drop-in upgrade: two new theorems (saturation; γ-driven extremal switch), "
    "completed Loyola naming, corrected sensitivity section, expanded literature review, and a "
    "referee report on this revised version.</i>", SUB))
A(rule())

# ===================== ABSTRACT =====================
A(h("Abstract"))
A(P("Degree-based topological indices proliferate, frequently by adjoining parameters to a "
    "bond-incident-degree (BID) template. We ask a structural question: can continuous parametric "
    "freedom ever yield a descriptor independent of the classical BID basis? Using a maximally "
    "flexible three-parameter family — the Loyola index "
    "LO(G;α,β,γ), which recovers a dozen classical indices as exact reductions — we answer in the "
    "negative and make the answer precise. <b>(i)</b> We prove a <b>saturation theorem</b>: on graphs "
    "of bounded maximum degree the family's descriptor span equals the column space of the "
    "edge-degree-pair count matrix, of dimension rank[m<sub>ij</sub>] ≤ Δ(Δ+1)/2 — exactly ten on "
    "molecular graphs (Δ ≤ 4) — so no choice of (α,β,γ) escapes a fixed finite-dimensional space. "
    "<b>(ii)</b> We prove the imbalance parameter γ is nonetheless <b>not redundant for extremal "
    "problems</b>: over trees of order n the maximizer of LO switches from the path (small γ) to the "
    "star (large γ), so the extremal graph genuinely depends on γ. Closed forms and a suite of sharp "
    "bounds confirm the family is analytically controlled, and a sensitivity study on octane, tree, "
    "and decane benchmarks shows — consistently with the saturation theorem — that the family occupies "
    "the same correlation regime as the classical indices rather than above it. The Loyola family is "
    "thus best understood not as a competitor index but as a controlled probe that renders the BID "
    "descriptor ceiling explicit.", BODY))
A(rule())

# ===================== 1. INTRODUCTION + LITERATURE REVIEW =====================
A(H("1.  Introduction and literature review"))
A(P("The theory of degree-based topological indices begins with the Wiener index and matures through "
    "the Randić connectivity index and the Zagreb indices of Gutman and Trinajstić. A first wave of "
    "<i>parametric</i> generalisations followed: the general Randić index R<sub>a</sub>, the general "
    "sum-connectivity index χ<sub>b</sub>, the forgotten index F, the geometric–arithmetic (GA) and "
    "atom-bond-connectivity (ABC) indices, and the Albertson irregularity measure. Each replaces a "
    "fixed edge weight by a tunable one, and a substantial extremal-graph literature characterises the "
    "trees, unicyclic and bicyclic graphs that optimise these weights.", BODY))
A(P("A second, larger wave is geometric. Gutman's Sombor index recast the edge weight as a Euclidean "
    "length √(d<sub>u</sub>²+d<sub>v</sub>²), and a now-crowded family of relatives has appeared — "
    "reduced, average, elliptic, hyperbolic and mean Sombor indices, the diminished Sombor index, and "
    "vertex–edge KG-Sombor variants — together with bound, spectral and QSPR studies of each. This "
    "rapid diversification is increasingly consolidated by structural machinery: M-polynomial and "
    "matrix formulations that compute whole index families at once, and, most relevant here, Rada's "
    "finite-dimensional linear framework, which views every vertex-degree-based index on a bounded-degree "
    "class as a linear functional of a fixed, finite set of edge-degree-pair counts.", BODY))
A(P("That linear viewpoint exposes a tension the proliferation literature rarely confronts. If, on the "
    "molecular class, <i>every</i> BID index is a linear combination of the same finite basis of "
    "edge-degree-pair counts, then introducing a new parametric index cannot enlarge the descriptor "
    "space — it can only move within it. The natural way to test this is not to add yet another index, "
    "but to take a family with the <i>most</i> parametric freedom one can reasonably write down and ask "
    "whether that freedom buys any independent descriptor. This is the programme of the present paper.", BODY))
A(P("<b>Contributions.</b> We introduce the three-parameter Loyola family as such a probe. Its product "
    "(d<sub>u</sub>d<sub>v</sub>)<sup>α</sup>, sum (d<sub>u</sub>+d<sub>v</sub>)<sup>β</sup>, and "
    "exponential-imbalance exp(γ·imb) factors recover the Zagreb, Randić, sum-connectivity, harmonic, "
    "inverse-sum-indeg, geometric–arithmetic and arithmetic–geometric indices as exact reductions, so it "
    "is a faithful chart of the BID template. We then prove (Thm 3) that its descriptor span is exactly "
    "the column space of the edge-pair-count matrix — a hard ten-dimensional ceiling on molecular graphs "
    "— and (Thm 4) that despite this, the imbalance parameter γ controls extremal structure, inverting the "
    "maximising tree from the path to the star. Closed forms (§2), bounds (§3) and a sensitivity study (§5) "
    "place and interpret the family inside that fixed space. The conclusion is deliberately a negative "
    "one, stated exactly: parametric flexibility reparametrises a saturated descriptor space rather than "
    "enlarging it.", BODY))

# ===================== 2. THE FAMILY =====================
A(H("2.  The Loyola family and its reductions"))
A(P("Let G=(V,E) be a simple connected graph; d<sub>v</sub> is the degree of v, with maximum Δ and "
    "minimum δ. For α,β,γ ∈ ℝ the <b>Loyola index family</b> is", BODY))
A(eq(r"LO(G;\alpha,\beta,\gamma)=\sum_{uv\in E(G)}(d_u d_v)^{\alpha}\,(d_u+d_v)^{\beta}\,"
     r"\exp\!\left(\gamma\,\frac{|d_u-d_v|}{d_u+d_v}\right)."))
A(P("The summand depends only on the unordered degree pair {d<sub>u</sub>,d<sub>v</sub>}; we write "
    "ψ(i,j;α,β,γ) for it and imb(i,j)=|i−j|/(i+j) ∈ [0,1) for the normalised imbalance. At specific "
    "triples LO reduces, exactly or up to a constant, to: the edge count m=LO(0,0,0), first Zagreb "
    "M₁=LO(0,1,0), second Zagreb M₂=LO(1,0,0), hyper-Zagreb HM=LO(0,2,0), modified second Zagreb "
    "ᵐM₂=LO(−1,0,0), Randić R=LO(−½,0,0), sum-connectivity χ=LO(0,−½,0), half-harmonic "
    "H/2=LO(0,−1,0), inverse-sum-indeg ISI=LO(1,−1,0), and (up to a factor) GA and AG. The whole "
    "general Randić family is the β=γ=0 axis and the general sum-connectivity family the α=γ=0 axis. "
    "The single new ingredient is the exponential imbalance factor, bounded above by e<sup>γ</sup> and "
    "equal to 1 on regular graphs.", BODY))
A(P("<b>Closed forms (summary).</b> Direct edge enumeration gives LO on the standard families: "
    "K<sub>n</sub>, C<sub>n</sub>, P<sub>n</sub>, the star S<sub>n</sub>, K<sub>p,q</sub>, the hypercube "
    "Q<sub>k</sub>, the wheel, the friendship graph, and any r-regular graph "
    "(LO = m·2<sup>β</sup>r<sup>2α+β</sup>). These are routine but exact and serve as the reference values "
    "used below; their derivations are unchanged from the standard computation and were verified "
    "symbolically and numerically.", BODY))

# ===================== 3. BOUNDS (compact) =====================
A(H("3.  Bounds (compact statement)"))
A(P("The family admits a coherent suite of sharp bounds, each a short proof from a standard inequality; "
    "we state them compactly since the techniques are classical and the proofs were verified line by line.", BODY))
A(P("• <b>Bracket bound.</b> m·ψ<sub>min</sub>(δ,Δ) ≤ LO ≤ m·ψ<sub>max</sub>(δ,Δ), with equality iff every "
    "edge attains the extremum; sharp on regular graphs.<br/>"
    "• <b>Cauchy–Schwarz splitting.</b> For additive parameter splits, "
    "LO(α,β,γ) ≤ √(LO(α₁,β₁,γ₁)·LO(α₂,β₂,γ₂)); the imbalance term is parameter-independent so exponents add.<br/>"
    "• <b>Jensen bounds.</b> LO(0,β,0) ≥ m<sup>1−β</sup>M₁<sup>β</sup> for β≥1 (≤ for 0&lt;β&lt;1); "
    "LO(α,0,0)=R<sub>α</sub> ≥ m<sup>1−α</sup>M₂<sup>α</sup> for α≥1.<br/>"
    "• <b>Degree bounds.</b> For α,β,γ≥0: LO ≤ m·Δ<sup>2α</sup>(2Δ)<sup>β</sup>exp(γ(Δ−1)/(Δ+1)); "
    "LO ≥ m·δ<sup>2α</sup>(2δ)<sup>β</sup>, equality iff δ-regular.<br/>"
    "• <b>Ratio bound.</b> For any BID index I<sub>h</sub> with positive edge weight, "
    "c<sub>min</sub>I<sub>h</sub> ≤ LO ≤ c<sub>max</sub>I<sub>h</sub>.<br/>"
    "• <b>Randić/sum-connectivity &amp; Sombor.</b> LO(α,β,0) ≤ √(R<sub>2α</sub>·χ<sub>2β</sub>); "
    "LO(0,β,0) ≤ m<sup>1−β</sup>(√2·SO)<sup>β</sup>.", BODY))
A(NOTE and P("<b>Honest note (kept from the referee feedback).</b> The Δ-upper bound multiplies three "
    "per-factor maxima that cannot be attained simultaneously (the product/sum maxima need (Δ,Δ) but the "
    "imbalance maximum needs (1,Δ)); for γ&gt;0 it is therefore valid but not sharp. The ratio bound holds "
    "for any positive edge function, so §3 establishes regularity of the family rather than anything "
    "specific to it — which is exactly why the contribution of this paper is §4, not §3.", NOTE))

# ===================== 4. NEW THEOREMS =====================
A(H("4.  Saturation and the role of γ  (new results)"))

A(h("4.1  The descriptor-space saturation theorem"))
A(P("Fix Δ₀ and let 𝒢={G₁,…,G_N} be graphs with Δ ≤ Δ₀. For 1 ≤ i ≤ j ≤ Δ₀ let "
    "m<sub>ij</sub>(G)=|{uv∈E(G):{d<sub>u</sub>,d<sub>v</sub>}={i,j}}| and collect these into the "
    "N×D matrix M, where D = Δ₀(Δ₀+1)/2 is the number of admissible degree-pair classes.", BODY))
A(P("<b>Theorem 3 (BID saturation).</b> For every (α,β,γ) the descriptor vector "
    "(LO(G<sub>k</sub>;α,β,γ))<sub>k</sub> equals M·ψ(α,β,γ), and", THM))
A(eq(r"\mathrm{span}\{(LO(G_k;\alpha,\beta,\gamma))_k:(\alpha,\beta,\gamma)\in\mathbb{R}^3\}"
     r"=\mathrm{col}(M),\quad \mathrm{dim}=\mathrm{rank}(M)\leq \frac{\Delta_0(\Delta_0+1)}{2}.", fs=14))
A(P("In particular, on molecular graphs (Δ₀ = 4) the entire three-parameter family lives in a space of "
    "dimension at most ten, and it fills that space exactly: the bound is rank(M), attained.", THM))
A(P("<b>Proof.</b> Partitioning E(G) into its degree-pair classes and using that ψ depends only on the "
    "class gives LO(G;α,β,γ)=Σ<sub>(i,j)</sub> ψ(i,j;α,β,γ)·m<sub>ij</sub>(G), i.e. the vector identity "
    "W(α,β,γ)=M·ψ(α,β,γ) with ψ(α,β,γ)=(ψ(i,j;α,β,γ))<sub>(i,j)</sub> ∈ ℝ<sup>D</sup>. Writing the edge "
    "weight in exponential form,", PROOF))
A(eq(r"\psi(i,j;\alpha,\beta,\gamma)=\exp(\alpha\ln(ij)+\beta\ln(i+j)+\gamma\,\mathrm{imb}(i,j))"
     r"=\exp(\theta\cdot v_{ij}),\quad \theta=(\alpha,\beta,\gamma)", fs=13))
A(P("with frequency vector v<sub>ij</sub>=(ln(ij), ln(i+j), imb(i,j)). The D pairs give D <i>distinct</i> "
    "frequency vectors: the sum i+j and product ij determine {i,j}, so distinct classes have distinct "
    "(ln(ij),ln(i+j)) already. Exponentials θ↦exp⟨θ,v⟩ with distinct frequencies v are linearly "
    "independent functions of θ; hence the only covector annihilating every ψ(θ) is zero, so "
    "{ψ(θ):θ∈ℝ³} spans all of ℝ<sup>D</sup>. Therefore span{M·ψ(θ)} = M·ℝ<sup>D</sup> = col(M), of "
    "dimension rank(M) ≤ D. □", PROOF))
A(P("<b>Verification.</b> For Δ₀=4 the ten frequency vectors are distinct and the 10×K matrix of ψ-values "
    "over a parameter grid has numerical rank 10, confirming the span is the full column space. "
    "On the 18 octanes the count matrix has rank 8, so the octane Loyola space is 8-dimensional and is "
    "already spanned by the classical reductions — exactly the redundancy the sensitivity study observes.", PROOF))

A(h("4.2  γ is not redundant: a tree-extremal switch"))
A(P("Saturation says γ adds no descriptor <i>dimension</i>. It does not say γ is inert: for extremal "
    "problems the parameter is decisive.", BODY))
A(P("<b>Theorem 4 (γ-driven extremal switch).</b> Let n ≥ 4 and α,β ∈ ℝ. "
    "<b>(a)</b> There is a threshold γ₀(α,β,n) such that for every γ &gt; γ₀ the star S<sub>n</sub> is the "
    "unique maximizer of LO(·;α,β,γ) over all trees of order n. "
    "<b>(b)</b> For (α,β)=(−½,0) the value at γ=0 is the ordinary Randić index, whose unique maximizer "
    "over trees of order n is the path P<sub>n</sub>; by continuity the path remains the maximizer for γ "
    "in a right-neighbourhood of 0. Hence the extremal tree depends on γ.", THM))
A(P("<b>Proof of (a).</b> Every edge of the star is the pair (1,n−1) with imbalance (n−2)/n. A non-star "
    "tree T has Δ(T) ≤ n−2, so each of its edges has imbalance at most (n−3)/(n−1), and", PROOF))
A(eq(r"\frac{n-2}{n}-\frac{n-3}{n-1}=\frac{2}{n(n-1)}>0,", fs=13))
A(P("i.e. every star edge strictly exceeds every non-star edge in imbalance. Let "
    "m<sub>n</sub>=max{(ij)<sup>α</sup>(i+j)<sup>β</sup> : 1 ≤ i ≤ j ≤ n−2} (finite, positive). Since a "
    "tree has n−1 edges and γ ≥ 0,", PROOF))
A(eq(r"LO(T)\leq (n-1)\,m_n\,e^{\gamma(n-3)/(n-1)},\quad "
     r"LO(S_n)=(n-1)(n-1)^{\alpha}n^{\beta}\,e^{\gamma(n-2)/n}", fs=12))
A(eq(r"\frac{LO(S_n)}{LO(T)}\ \geq\ \frac{(n-1)^{\alpha}n^{\beta}}{m_n}\,"
     r"\exp(\gamma\cdot 2/(n(n-1)))\ \to\ \infty\quad(\gamma\to\infty).", fs=12))
A(P("As there are finitely many non-star trees, one threshold "
    "γ₀ = max(0, (n(n−1)/2)·ln(m<sub>n</sub>/((n−1)<sup>α</sup>n<sup>β</sup>))) works for all of them, "
    "and the inequality is strict, giving uniqueness. <b>Proof of (b).</b> At γ=0, "
    "LO(·;−½,0,0)=Σ(d<sub>u</sub>d<sub>v</sub>)<sup>−½</sup>=R, the Randić index, whose maximum over "
    "trees of order n is uniquely attained by P<sub>n</sub> (Bollobás–Erdős). Each γ↦LO(T;−½,0,γ) is "
    "continuous and there are finitely many trees, so P<sub>n</sub> remains the strict maximizer on a "
    "right-neighbourhood of γ=0. □", PROOF))
A(P("<b>Numerical confirmation</b> (maximizer of LO over all trees of order n at α=−½, β=0):", BODY))

# Fix B validation table
tb = [["γ", "n=7", "n=8", "n=9", "n=10"],
      ["0",   "path","path","path","path"],
      ["0.5", "other","other","other","other"],
      ["1",   "star","other","other","other"],
      ["2",   "star","star","star","star"],
      ["3",   "star","star","star","star"],
      ["5",   "star","star","star","star"]]
t = Table(tb, colWidths=[50,70,70,70,70], hAlign="CENTER")
t.setStyle(TableStyle([
    ("FONT",(0,0),(-1,-1),"DJ",8.5),
    ("FONT",(0,0),(-1,0),"DJ-B",8.5),
    ("BACKGROUND",(0,0),(-1,0),colors.HexColor("#16213e")),
    ("TEXTCOLOR",(0,0),(-1,0),colors.white),
    ("GRID",(0,0),(-1,-1),0.4,colors.grey),
    ("ALIGN",(0,0),(-1,-1),"CENTER"),
    ("ROWBACKGROUNDS",(0,1),(-1,-1),[colors.white,colors.HexColor("#f2f4fa")]),
    ("TEXTCOLOR",(1,1),(-1,1),colors.HexColor("#0a6")),  # path row
]))
A(KeepTogether([t]))
A(P("The maximizer travels path → intermediate tree → star as γ grows: γ is load-bearing, and the "
    "intermediate regime (a family of broom-like trees) is itself an open extremal question.", PROOF))

# ===================== 5. SENSITIVITY (Fix C) =====================
A(H("5.  Sensitivity analysis  (corrected)"))
A(P("We benchmark LO on three standard testbeds — the 18 octane isomers (five NIST physicochemical "
    "properties), the 106 non-isomorphic trees of order 10, and the 75 decane isomers — using a fully "
    "self-contained protocol. Throughout, the family is written LO(·) and never abbreviated to W(·) (which "
    "would collide with the Wiener index in the comparison set).", BODY))
A(h("5.1  Property correlations on the octanes"))
A(P("Across the five properties the strongest correlations are split between classical indices "
    "(the augmented Zagreb index on T<sub>B</sub> and ΔH<sub>f</sub>; the second Zagreb index on ω) and "
    "the pure-γ point LO(0,0,2) on ΔH<sub>vap</sub> (|r|=0.967). <b>However, every Loyola margin lies "
    "within the n=18 bootstrap confidence half-width (≈0.03), so no parameter triple significantly "
    "outperforms the classical indices on any property.</b> The exact reductions reproduce their classical "
    "counterparts to floating-point precision, a useful internal check rather than a finding.", BODY))
predfig = "/Users/ganeshshiwakoti/topological_index_pilot/figures/fig_lo_prediction_correlation_heatmap.png"
if os.path.exists(predfig):
    iw, ih = ImageReader(predfig).getSize(); w = 380; A(Image(predfig, width=w, height=w*ih/iw, hAlign="CENTER"))
A(h("5.2  Degeneracy on trees of order 10  (definition corrected)"))
A(P("Degeneracy measures how often an index fails to distinguish non-isomorphic graphs. We use the "
    "standard distinct-value convention (lower is better):", BODY))
A(eq(r"\mathrm{Degeneracy}(\%)=100\left(1-\frac{k}{N}\right),\quad k=\#\{\mathrm{distinct\ values}\},\quad N=106.", fs=13))
A(NOTE and P("<b>Correction.</b> An earlier draft defined degeneracy as “the fraction of pairs of distinct "
    "trees on which the index takes the same value”, which does <i>not</i> equal the plotted bars (those are "
    "1 − #distinct/N). The definition above is the one consistent with the figure and is used throughout. "
    "On the 106 trees the pure-γ points LO(0,0,1), LO(0,0,2) attain ≈25% degeneracy — among the lowest of "
    "all indices tested, i.e. good discrimination; M₁ and M₂ are the most degenerate.", NOTE))
A(h("5.3  Structure sensitivity on decanes"))
A(P("As self-contained dispersion measures over the 75 decanes we report the coefficient of variation, the "
    "relative range, and their ratio:", BODY))
A(eq(r"SS=\sigma/\mu,\qquad Abr=(\max-\min)/\mu,\qquad SA=SS/Abr.", fs=13))
A(P("These are within-study measures (not a fingerprint-similarity procedure), so no cross-paper "
    "comparison is claimed. The one defensible positive: the pure-γ point LO(0,0,2) attains the "
    "<b>highest structure sensitivity</b> in the set — the imbalance term's contribution surfaces here, "
    "in responsiveness, rather than in added correlation. This is consistent with the saturation theorem: "
    "γ moves the descriptor within col(M) (changing its dispersion) without adding a dimension.", BODY))

# ===================== 6. CONCLUSION =====================
A(H("6.  Conclusion"))
A(P("The Loyola family is a maximally flexible three-parameter BID index. We used it not to compete with "
    "the classical indices but to probe the descriptor space they inhabit, and proved that the probe "
    "saturates: its span is exactly the column space of the edge-degree-pair count matrix, a hard "
    "ten-dimensional ceiling on molecular graphs, attained but never exceeded. The imbalance parameter γ, "
    "though it adds no descriptor dimension, is not inert — it inverts the extremal tree from the path to "
    "the star. Closed forms and sharp bounds show the family is analytically well behaved, and the "
    "sensitivity study confirms the predicted picture: on real testbeds the family sits squarely in the "
    "classical correlation regime. The take-away for the proliferation literature is a structural one: on "
    "bounded-degree graphs, continuous parametric freedom reparametrises a saturated space rather than "
    "enlarging it, and a genuinely new descriptor must come from outside the BID template, not from more "
    "parameters within it.", BODY))
A(P("<b>Open problems.</b> (1) Characterise the intermediate-γ tree maximizers between path and star. "
    "(2) Sharpen the Δ-upper bound to an attained form. (3) Identify any property whose optimal BID "
    "predictor provably requires a non-reduction triple — or prove none exists on Δ ≤ 4.", BODY))

A(rule())
# ===================== APPENDIX: REFEREE REPORT =====================
A(H("Appendix.  Referee report on this revised manuscript"))
A(P("<b>Recommendation: accept with minor revision (realistic MATCH); weak accept / minor-to-major at a "
    "harsh desk.</b> The revision repairs the structural problems of the prior version.", BODY))
A(P("<b>What is now fixed.</b> (i) The paper no longer contradicts itself: it is a coherent negative "
    "result, with the saturation statement promoted from an aside to a proved theorem (Thm 3) and "
    "upgraded from “≤10” to the exact rank — this is the contribution. (ii) The chief objection to the "
    "prior draft — that the new parameter γ does nothing — is answered: Thm 4 shows γ controls extremal "
    "structure (path↔star), so the third parameter earns its place even though saturation holds. (iii) The "
    "naming is consistent (LO throughout; the W/Wiener collision is gone). (iv) The degeneracy definition "
    "now matches its figure. (v) The prediction claims are stated at their true (non-significant) strength.", BODY))
A(P("<b>Residual weaknesses (honest).</b> (a) Thm 3 is elementary — distinct-frequency linear "
    "independence — so its value is framing and exactness, not depth; a hostile referee may still call it "
    "“an observation, well packaged.” Mitigated by the exact-rank form and the explicit tightness. "
    "(b) Thm 4(b) is anchored at one reduction (α=−½,β=0); the general (α,β) transition is only partially "
    "described. Stating at least the qualitative two-regime picture for a parameter band would strengthen "
    "it. (c) §3 remains classical technique; correctly demoted, but it adds little. (d) n=18 octanes is a "
    "small chemometric base — acceptable for a structural paper whose point is redundancy, but it cannot "
    "carry a predictive claim, and the paper now correctly makes none. (e) The index is still named after "
    "an institution; harmless, but a descriptive name would read better.", BODY))
A(P("<b>Verdict.</b> Previously this was “major revision trending reject” on novelty grounds. With a proved "
    "saturation theorem as the headline and a theorem making γ load-bearing, the paper is now a defensible, "
    "self-consistent contribution to the structural side of the degree-based-index literature — the kind a "
    "referee can respect rather than merely tolerate. The remaining work is to deepen Thm 4 beyond a single "
    "reduction and to decide whether §3 is worth its length.", BODY))

# ===================== BUILD =====================
doc = SimpleDocTemplate(OUT, pagesize=A4, topMargin=0.7*inch, bottomMargin=0.7*inch,
                        leftMargin=0.8*inch, rightMargin=0.8*inch,
                        title="Loyola Index — Upgraded Manuscript")
# strip accidental booleans from the `NOTE and P(...)` idiom
flow = [f for f in flow if not isinstance(f, bool)]
doc.build(flow)
print("WROTE", OUT)
if FAILED:
    print("!!! EQUATIONS THAT FAILED TO RENDER AS MATH (FIX THESE):")
    for f in FAILED:
        print("   ", f)
else:
    print("OK: all equations rendered as typeset math (none hit the text fallback).")
