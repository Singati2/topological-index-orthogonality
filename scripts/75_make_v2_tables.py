"""Generate the LaTeX tables of Paper 2 (v2) directly from result CSVs.

Writes docs/tables_v2/*.tex, which the manuscript \\input's, so that every
number in these tables is produced by code (no hand transcription).
"""
import os
import numpy as np
import pandas as pd

P = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
R = os.path.join(P, "results")
OUT = os.path.join(P, "docs", "tables_v2")
os.makedirs(OUT, exist_ok=True)
DS = ["esol", "freesolv", "lipophilicity", "bbbp"]
NICE = {"esol": "ESOL", "freesolv": "FreeSolv", "lipophilicity": "Lipophilicity", "bbbp": "BBBP"}


def sci(x):
    if x == 0:
        return "$0$"
    e = int(np.floor(np.log10(abs(x))))
    m = x / 10 ** e
    return f"${m:.1f}\\times10^{{{e}}}$"


def cert_table():
    c = pd.read_csv(os.path.join(R, "bid_span_certificate.csv")).set_index("dataset")
    v = pd.read_csv(os.path.join(R, "bid_variant_screen.csv"))
    rows = []
    for d in DS:
        r = c.loc[d]
        lo = v[(v.dataset == d) & (v.family == "Loyola grid")]
        pb = v[(v.dataset == d) & (v.family == "published BID")]
        rows.append(f"{NICE[d]} & ${int(r.n)}$ & ${int(r.degree_pairs_realised)}$ & ${int(r.rank_mij)}$ & "
                    f"${int(r.rank_BID18)}$ & ${int(r.rank_BID18_plus_mij)}$ & {sci(r.max_rel_resid_mij_on_BID18)} & "
                    f"{sci(pb.rel_resid_on_baseline30.max())} & {sci(lo.rel_resid_on_baseline30.max())}\\\\")
    body = "\n".join(rows)
    tex = r"""\begin{table}[h]
\centering
\small
\caption{Span certificate on the four datasets. $|P_D|$: realised degree pairs; ranks are numerical ranks of the column-centred, column-normalised matrices (relative tolerance $10^{-9}$). ``Resid.'' columns give the largest relative residual $\|r\|/\|z-\bar z\|$ on the baseline of the ten counts $m_{ij}$, of the $15$ published BID instances of Table~\ref{tab:census}, and of the $100$ Loyola grid points.}
\label{tab:cert}
\setlength{\tabcolsep}{4pt}
\begin{tabular}{lrrrrrccc}
\toprule
Dataset & $n$ & $|P_D|$ & rk $M$ & rk $B$ & rk $[B,M]$ & resid.\ $m_{ij}$ & resid.\ published & resid.\ Loyola\\
\midrule
""" + body + r"""
\bottomrule
\end{tabular}
\end{table}
"""
    open(os.path.join(OUT, "cert.tex"), "w").write(tex)


CITE = {
    "mSO (modified Sombor)": ("Modified Sombor", r"$1/\sqrt{a^2+b^2}$", "KulliGutman2021Modified"),
    "ESO (elliptic Sombor)": ("Elliptic Sombor", r"$(a+b)\sqrt{a^2+b^2}$", "GutmanFurtulaOz2024Elliptic"),
    "EU (Euler Sombor)": ("Euler Sombor", r"$\sqrt{a^2+b^2+ab}$", "Gutman2024Euler"),
    "DSO (diminished Sombor)": ("Diminished Sombor", r"$\sqrt{a^2+b^2}/(a+b)$", "MovahediGutmanRedzepovicFurtula2026"),
    "HSO (hyperbolic Sombor)": ("Hyperbolic Sombor", r"$\sqrt{a^2+b^2}/\min(a,b)$", "BarmanDas2026Hyperbolic"),
    "SO_p, p=1/2 (p-Sombor)": ("$p$-Sombor, $p=\\tfrac12$", r"$(a^{1/2}+b^{1/2})^{2}$", "RetiDoslicAli2021"),
    "SO_p, p=3 (p-Sombor)": ("$p$-Sombor, $p=3$", r"$(a^3+b^3)^{1/3}$", "RetiDoslicAli2021"),
    "KG-SO (KG-Sombor)": ("KG-Sombor", r"$\sqrt{a^2+e^2}+\sqrt{b^2+e^2}$", "KulliHarishChaluvarajuGutman2022KG"),
    "ISI (inverse sum indeg)": ("Inverse sum indeg", r"$ab/(a+b)$", "VukicevicGasperov2010Adriatic"),
    "SDD (symm. division deg)": ("Symmetric division deg", r"$a/b+b/a$", "VukicevicGasperov2010Adriatic"),
    "N (Nirmala)": ("Nirmala", r"$\sqrt{a+b}$", "Kulli2021Nirmala"),
    "IN1 (first inverse Nirmala)": ("First inverse Nirmala", r"$\sqrt{1/a+1/b}$", "KulliLokeshaNirupadi2021InvNirmala"),
    "IN2 (second inverse Nirmala)": ("Second inverse Nirmala", r"$\sqrt{ab/(a+b)}$", "KulliLokeshaNirupadi2021InvNirmala"),
    "GG1 (first Gourava)": ("First Gourava", r"$a+b+ab$", "Kulli2017Gourava"),
    "HM (hyper-Zagreb)": ("Hyper-Zagreb", r"$(a+b)^2$", "ShirdelRezapourSayadi2013HyperZagreb"),
}


def census_table():
    v = pd.read_csv(os.path.join(R, "bid_variant_screen.csv"))
    rows = []
    for lab, (name, f, key) in CITE.items():
        sub = v[v["index"] == lab].set_index("dataset")
        rs = " & ".join(f"${sub.loc[d, 'max_abs_r_baseline30']:.3f}$" for d in DS)
        rows.append(f"{name}~\\cite{{{key}}} & {f} & {rs} & {sci(sub.rel_resid_on_baseline30.max())}\\\\")
    rows.append(r"\midrule")
    for lab, name in [("SO coindex", "Sombor coindex (control)"), ("SO_ecc (eccentric Sombor)", "Eccentric Sombor (control)")]:
        sub = v[v["index"] == lab].set_index("dataset")
        rs = " & ".join(f"${sub.loc[d, 'max_abs_r_baseline30']:.3f}$" for d in DS)
        rows.append(f"{name} & non-BID & {rs} & {sci(sub.rel_resid_on_baseline30.max())}\\\\")
    tex = r"""\begin{table}[h]
\centering
\small
\caption{Census of published BID indices ($f(a,b)$ with $a,b$ the end degrees of an edge) and two non-BID Sombor-type controls. Columns 3--6: largest $|r|$ with any of the $30$ baseline indices on each dataset (the conventional pairwise screen). Last column: largest relative residual on the baseline over the four datasets. For KG-Sombor $e=a+b-2$ is the edge degree. Every BID index is in the span to machine precision; the controls are not.}
\label{tab:census}
\footnotesize
\setlength{\tabcolsep}{3pt}
\begin{tabular}{llccccc}
\toprule
Index & $f(a,b)$ & ESOL & FreeSolv & Lipo & BBBP & max resid.\\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}
\end{table}
"""
    open(os.path.join(OUT, "census.tex"), "w").write(tex)


def nonbid_table():
    c = pd.read_csv(os.path.join(R, "novel_candidates_experiment", "novel_candidates_multidataset.csv"))
    n = pd.read_csv(os.path.join(R, "nonlinear_conditional_tests.csv"))
    m = c.merge(n, on=["dataset", "candidate"], how="left")
    counts = []
    for d in DS:
        g = c[c.dataset == d]
        vc = g.pcor_ci_verdict.value_counts()
        counts.append(f"{NICE[d]} & ${(g.pairwise_verdict == 'PASS').sum()}$ & ${(g.combined_verdict == 'PASS').sum()}$ & "
                      f"${int(((g.pcor_ci_verdict == 'PASS') & (g.pairwise_verdict == 'PASS')).sum())}$ & "
                      f"${int(vc.get('INCONCLUSIVE', 0))}$ & ${int(vc.get('NEGLIGIBLE', 0))}$ & ${int(vc.get('IN_SPAN', 0))}$\\\\")
    sel = m[(m.partial_corr_target.abs() >= 0.10) | (m.rho_nl.abs() >= 0.10)].copy()
    sel["o"] = sel.dataset.map({d: i for i, d in enumerate(DS)})
    sel = sel.sort_values(["o", "candidate"])
    det = []
    for _, r in sel.iterrows():
        nm = r.candidate.replace("_", r"\_")
        det.append(f"{NICE[r.dataset]} & \\texttt{{{nm}}} & ${r.max_abs_r_baseline:.3f}$ & "
                   f"${r.partial_corr_target:+.3f}$ $[{r.pcor_ci_lo:+.3f},{r.pcor_ci_hi:+.3f}]$ & {r.pcor_ci_verdict.lower()} & "
                   f"${r.rho_nl:+.3f}$ $[{r.ci_lo:+.3f},{r.ci_hi:+.3f}]$\\\\")
    refs = n[n.candidate.str.startswith("REF")]
    rng = lambda key: (refs[refs.candidate == key].rho_nl.min(), refs[refs.candidate == key].rho_nl.max())
    a, b = rng("REF_SCI_inspan"); c1, c2 = rng("REF_noise"); l1, l2 = rng("REF_leak")
    tex = r"""\begin{table}[h]
\centering
\small
\caption{Target-aware screen of the $20$ non-BID candidates. Top: counts per dataset --- pairwise pass ($\max|r|<0.95$), point-estimate combined pass ($|\pcor|\ge0.10$ and pairwise pass), interval pass (pairwise pass and $95\%$ interval of $|\pcor|$ above $0.10$), inconclusive, negligible, and in span. Bottom: every candidate--dataset pair with $|\pcor|\ge0.10$ or $|\rho_{\mathrm{nl}}|\ge0.10$, with $95\%$ Fisher-$z$ intervals.}
\label{tab:nonbid}
\begin{tabular}{lcccccc}
\toprule
Dataset & pairwise & point comb. & interval pass & inconcl. & neglig. & in span\\
\midrule
""" + "\n".join(counts) + r"""
\bottomrule
\end{tabular}

\vspace{0.6em}
\setlength{\tabcolsep}{3pt}
\begin{tabular}{llcccc}
\toprule
Dataset & Candidate & $\max|r|$ & linear $\pcor$ [95\% CI] & verdict & $\rho_{\mathrm{nl}}$ [95\% CI]\\
\midrule
""" + "\n".join(det) + r"""
\bottomrule
\end{tabular}

\vspace{0.3em}
\noindent\footnotesize{Calibration of $\rho_{\mathrm{nl}}$ across the four datasets: exactly redundant BID reference ($\sum(d_u+d_v)^{-1/2}$, i.e.\ SCI) """ + f"${a:+.3f}$ to ${b:+.3f}$; pure noise ${c1:+.3f}$ to ${c2:+.3f}$; leaky column $y+\\varepsilon$ ${l1:.3f}$ to ${l2:.3f}$." + r"""}
\end{table}
"""
    open(os.path.join(OUT, "nonbid.tex"), "w").write(tex)


def data_table():
    c = pd.read_csv(os.path.join(R, "benchmark_v2_curation.csv")).set_index("dataset")
    target = {"esol": r"$\log S$ (mol/L)", "freesolv": r"$\Delta G_{\mathrm{hyd}}$ (kcal/mol)",
              "lipophilicity": r"$\log D_{7.4}$", "bbbp": "BBB penetration (0/1)"}
    rows = []
    for d in DS:
        r = c.loc[d]
        parsed = int(r.n_valid) - int(r["n_no_graph_dropped(<2 heavy atoms)"])
        rows.append(f"{NICE[d]} & {target[d]} & ${int(r.n_raw)}$ & ${int(r.n_invalid_smiles)}$ & "
                    f"${int(r.n_multifragment_stripped)}$ & ${parsed}$ & ${int(r.n_rows_in_dup_groups)}$ & "
                    f"${int(r.n_conflicting_dup_groups_dropped)}$ & ${int(r.n_final)}$ & "
                    f"${int(r.n_distinct_graphs)}$ & ${r.pct_mols_graph_shared:.1f}$\\\\")
    tex = r"""\begin{table}[h]
\centering
\small
\caption{Datasets and curation. Raw rows of the MoleculeNet files; invalid SMILES; multi-fragment SMILES reduced to the largest fragment; molecules parsed to a graph with at least two heavy atoms (the sets used by the screens); rows in exact canonical-SMILES duplicate groups; duplicate groups dropped for conflicting labels; molecules after curation (also dropping graphs with fewer than two heavy atoms); distinct hydrogen-suppressed unlabelled graphs (Weisfeiler--Lehman hash); and the percentage of curated molecules whose graph is shared with another molecule.}
\label{tab:data}
\setlength{\tabcolsep}{3pt}
\begin{tabular}{llrrrrrrrrr}
\toprule
Dataset & Target & raw & invalid & multi-frag. & parsed & dup.\ rows & conflicts & curated & graphs & shared (\%)\\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}
\end{table}
"""
    tex = tex.replace("\\begin{tabular}", "\\resizebox{\\textwidth}{!}{%\n\\begin{tabular}").replace("\\end{tabular}", "\\end{tabular}}")
    open(os.path.join(OUT, "data.tex"), "w").write(tex)


MODELS = [("a_full30_lasso", "full 30, LASSO / $\\ell_1$-logistic"),
          ("b_pairwise_lasso", "pairwise-screened, LASSO"),
          ("c_combined_lasso", "combined-screened, LASSO"),
          ("d_full30_rf", "full 30, random forest"),
          ("e_pairwise_rf", "pairwise-screened, random forest"),
          ("g_mij_rf", "ten counts $m_{ij}$, random forest"),
          ("g_mij_lasso", "ten counts $m_{ij}$, LASSO"),
          ("f_rdkit2d_rf", "RDKit 2D descriptors, random forest")]
COMPS = [("b_pairwise_lasso - a_full30_lasso", "pairwise vs full (LASSO)"),
         ("c_combined_lasso - a_full30_lasso", "combined vs full (LASSO)"),
         ("e_pairwise_rf - d_full30_rf", "pairwise vs full (RF)"),
         ("g_mij_rf - d_full30_rf", "$m_{ij}$ vs full 30 (RF)"),
         ("f_rdkit2d_rf - d_full30_rf", "RDKit 2D vs full 30 (RF)")]


def bench_table():
    sm = pd.read_csv(os.path.join(R, "benchmark_v2_summary.csv"))
    ts = pd.read_csv(os.path.join(R, "benchmark_v2_tests.csv"))
    rows = []
    for key, lab in MODELS:
        cells = []
        for d in DS:
            r = sm[(sm.dataset == d) & (sm.model == key)]
            if r.empty:
                cells.append("--"); continue
            r = r.iloc[0]
            star = "$^{\\ast}$" if r.n_intercept_only_folds > 0 else ""
            cells.append(f"${r['mean']:.3f}$ ({r.mean_n_features:.1f}){star}")
        rows.append(lab + " & " + " & ".join(cells) + "\\\\")
    trows = []
    for key, lab in COMPS:
        cells = []
        for d in DS:
            r = ts[(ts.dataset == d) & (ts.comparison.str.replace(" ", "") == key.replace(" ", ""))]
            if r.empty:
                cells.append("--"); continue
            r = r.iloc[0]
            eq = "eq." if r.tost_verdict == "equivalent" else ""
            ph = "<0.001" if r.p_holm < 0.001 else f"{r.p_holm:.3f}"
            cells.append(f"${r.mean_diff:+.3f}$ $[{r.ci90_lo:+.3f},{r.ci90_hi:+.3f}]$ {eq} ($p_H$={ph})")
        trows.append(lab + " & " + " & ".join(cells) + "\\\\")
    nint = int(sm[(sm.dataset == "lipophilicity") & (sm.model == "c_combined_lasso")].n_intercept_only_folds.iloc[0])
    tex = r"""\begin{table}[h]
\centering
\footnotesize
\caption{Downstream benchmark on the curated data: $5\\times5$ repeated cross-validation with folds grouped by graph, all preprocessing and selection inside the training folds. Top: mean RMSE (ESOL, FreeSolv, Lipophilicity) or ROC-AUC (BBBP), with the mean number of active features in parentheses. Bottom: mean paired difference (model minus reference; lower RMSE and higher AUC are better) with Nadeau--Bengio corrected $90\\%$ intervals; ``eq.'' marks equivalence (TOST) within $5\\%$ of the reference RMSE or $5\\%$ of the reference AUC above $0.5$; $p_H$ is the Holm-adjusted $p$-value of the corrected paired $t$-test within the dataset (equivalence and a nonzero difference can both hold).""" + f" $^{{\\ast}}$On Lipophilicity the combined screen kept no feature in {nint} of 25 folds (intercept-only model)." + r"""}
\label{tab:bench}
\setlength{\tabcolsep}{2.5pt}
\begin{tabular}{lcccc}
\toprule
Model & ESOL & FreeSolv & Lipophilicity & BBBP (AUC)\\\\
\midrule
""" + "\n".join(rows) + r"""
\midrule
Comparison & \multicolumn{4}{c}{difference [corrected 90\\% CI]}\\\\
\midrule
""" + "\n".join(trows) + r"""
\bottomrule
\end{tabular}
\end{table}
"""
    tex = (tex.replace(r"5\\times5", r"5\times5").replace(r"\\%", r"\%")
           .replace("(AUC)\\\\\\\\", "(AUC)\\\\").replace("CI]}\\\\\\\\", "CI]}\\\\"))
    tex = tex.replace("\\begin{tabular}", "\\resizebox{\\textwidth}{!}{%\n\\begin{tabular}").replace("\\end{tabular}", "\\end{tabular}}")
    open(os.path.join(OUT, "bench.tex"), "w").write(tex)


FORMULAS = {
    "M1": "a+b", "M2": "ab", "mM1": "a^{-3}+b^{-3}",
    "mM2": "1/(ab)", "F": "a^2+b^2", "R": "(ab)^{-1/2}", "SCI": "(a+b)^{-1/2}",
    "H": "2/(a+b)", "GA": "2\\sqrt{ab}/(a+b)", "AG": "(a+b)/(2\\sqrt{ab})",
    "ABC": "\\sqrt{(a+b-2)/(ab)}", "ABS": "\\sqrt{(a+b-2)/(a+b)}",
    "AZI": "\\bigl(ab/(a+b-2)\\bigr)^3\\;(0\\text{ on }K_2)", "SO": "\\sqrt{a^2+b^2}",
    "SO_red": "\\sqrt{(a-1)^2+(b-1)^2}", "Alb": "|a-b|", "Sigma": "(a-b)^2",
    "redM1": "(a-1)^2/a+(b-1)^2/b",
}
LABEL = {"SO_red": "SO$_{\\mathrm{red}}$", "Sigma": "$\\sigma$", "redM1": "redM$_1$",
         "mM1": "$mM_1$", "mM2": "$mM_2$", "M1": "$M_1$", "M2": "$M_2$"}


def fb_table():
    import importlib.util
    spec = importlib.util.spec_from_file_location("cert", os.path.join(P, "scripts", "70_bid_span_certificate.py"))
    cert = importlib.util.module_from_spec(spec); spec.loader.exec_module(cert)
    F = cert.f_matrix(cert.BID18)
    head = " & ".join(f"$({i},{j})$" for i, j in cert.P4)
    rows = []
    for k, name in enumerate(cert.BID18):
        vals = " & ".join(("$0$" if abs(x) < 1e-12 else f"${x:.4g}$") for x in F[k])
        rows.append(f"{LABEL.get(name, name)} & ${FORMULAS[name]}$ & {vals}\\\\")
    sv = np.linalg.svd(F, compute_uv=False)
    tex = r"""\begin{table}[h]
\centering
\caption{The value matrix $\\Phi_B$ of the $18$ degree-based baseline indices on the ten degree pairs $(i,j)\\in P_4$: entry $f_b(i,j)$, obtained by evaluating each index on $K_{i,j}$ and dividing by its $ij$ edges. Vertex-degree sums $\\sum_v g(d_v)$ are written as $f(a,b)=g(a)/a+g(b)/b$. The matrix has rank $10$""" + f" (smallest/largest singular value ${sv[-1]/sv[0]:.1e}$)" + r""".}
\label{tab:phi}
\resizebox{\\textwidth}{!}{%
\begin{tabular}{ll""" + "r" * 10 + r"""}
\toprule
Index & $f(a,b)$ & """ + head + r"""\\\\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}}
\end{table}
"""
    tex = tex.replace("\\\\\\\\", "\\\\")
    for a in ["Phi", "in P", "textwidth", "sum_v"]:
        tex = tex.replace("\\\\" + a, "\\" + a)
    open(os.path.join(OUT, "phi.tex"), "w").write(tex)


def chem_table():
    """Chemistry-aware screen (scripts/82: curated data, graph-grouped folds, in-fold
    column cleaning); verdict counts against RDKit alone and RDKit + topology."""
    sp = os.path.join(R, "screen_vs_chemistry_grouped_summary.csv")
    if not os.path.exists(sp):
        return
    S = pd.read_csv(sp).set_index("dataset")
    CV = pd.read_csv(os.path.join(R, "screen_vs_chemistry_grouped_cv.csv"))
    cv = {(r.dataset, r.feature_set): r.value for r in CV.itertuples()}
    rows = []
    for d in DS:
        r = S.loc[d]
        rows.append(f"{NICE[d]} & ${int(r.n)}$ & ${int(r.n_groups)}$ & ${int(r.n_rdkit)}$ & ${int(r.n_both)}$ & "
                    f"${cv[(d,'topo')]:.3f}$ & ${cv[(d,'rdkit')]:.3f}$ & ${cv[(d,'rdkit+topo')]:.3f}$ & ${cv[(d,'rdkit+topo+cands')]:.3f}$ & "
                    f"${int(r.rdkit_PASS)}$ / ${int(r.rdkit_INCONCLUSIVE)}$ / ${int(r.rdkit_NEGLIGIBLE)}$ & "
                    f"${int(r.both_PASS)}$ / ${int(r.both_INCONCLUSIVE)}$ / ${int(r.both_NEGLIGIBLE)}$\\\\")
    tex = r"""\begin{table}[h]
\centering
\footnotesize
\caption{The target-aware screen against a chemistry-aware baseline, on the curated data. $n$: molecules; groups: isomorphism classes; $p_{\\mathrm{RD}}$, $p_{\\mathrm{both}}$: RDKit 2D columns, and columns of the union with the $30$ topological indices, after removing non-finite, constant and collinear columns on the full sample (the $18$ degree-based indices contribute ten independent columns). CV: $5$-fold cross-validation grouped by isomorphism class, column cleaning, scaling and penalty selection inside the training folds; RMSE for the regression sets (ridge), ROC-AUC for BBBP ($\\ell_2$-logistic). Candidates: interval verdicts (pass / inconclusive / negligible; the remaining candidate is in span) given RDKit alone and given the union.}
\label{tab:chem}
\setlength{\\tabcolsep}{3pt}
\resizebox{\\textwidth}{!}{%
\begin{tabular}{lrrrrcccccc}
\toprule
 & & & & & \\multicolumn{4}{c}{grouped CV RMSE / AUC} & \\multicolumn{2}{c}{candidates pass/inc./negl.}\\\\
Dataset & $n$ & groups & $p_{\\mathrm{RD}}$ & $p_{\\mathrm{both}}$ & topo & RD & RD+topo & RD+topo+cand & given RD & given RD+topo\\\\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}}
\end{table}
"""
    tex = tex.replace("\\\\\\\\", "\\\\")
    for a in ["mathrm", "ell", "multicolumn", "tabcolsep", "textwidth"]:
        tex = tex.replace("\\\\" + a, "\\" + a)
    open(os.path.join(OUT, "chem.tex"), "w").write(tex)


def generality_table():
    g = pd.read_csv(os.path.join(R, "certificate_generality.csv"))
    order = ["esol", "freesolv", "lipophilicity", "bbbp", "bace", "clintox", "sider", "tox21", "toxcast_data", "qm8", "qm9", "muv", "HIV"]
    nice = dict(NICE, bace="BACE", clintox="ClinTox", sider="SIDER", tox21="Tox21", toxcast_data="ToxCast", qm8="QM8", qm9="QM9", muv="MUV", HIV="HIV")
    gi = g.set_index("dataset")
    g = gi.loc[[o for o in order if o in gi.index]]
    rows = []
    for d, r in g.iterrows():
        cert = "yes" if r.certificate_holds else "\\textbf{no}"
        extra = ""
        if "restricted_n" in g.columns and pd.notna(r.get("restricted_n", np.nan)):
            extra = f" (on the $\\Delta\\le4$ subset, $n={int(r.restricted_n)}$: ranks ${int(r.restricted_rank_mij)}$/${int(r.restricted_rank_BID18)}$, " + ("yes" if r.restricted_certificate_holds else "no") + ")"
        rows.append(f"{nice.get(d, d)} & ${int(r.n_graphs)}$ & ${int(r.n_maxdeg_gt4)}$ & ${int(r.max_degree)}$ & ${int(r.degree_pairs_realised)}$ & "
                    f"${int(r.rank_mij)}$ & ${int(r.rank_BID18)}$ & {sci(r.max_rel_resid_mij_on_BID18)} & {cert}{extra}\\\\")
    tex = r"""\begin{table}[h]
\centering
\footnotesize
\caption{The span certificate on thirteen MoleculeNet SMILES datasets (no targets used). $n$: parsed molecules; $n_{>4}$: molecules containing an atom of degree above four; $\\Delta$: largest degree; $|P_D|$: realised degree pairs; ranks of the count matrix $M$ and of the $18$-index baseline block $B$; largest relative residual of an $m_{ij}$ column on $B$; whether $\\operatorname{rank}[\\mathbf 1,B]=\\operatorname{rank}[\\mathbf 1,M]$ with residuals below $10^{-8}$.}
\label{tab:general}
\setlength{\\tabcolsep}{3pt}
\resizebox{\\textwidth}{!}{%
\begin{tabular}{lrrrrrrcl}
\toprule
Dataset & $n$ & $n_{>4}$ & $\\Delta$ & $|P_D|$ & rk $M$ & rk $B$ & max resid. & certified\\\\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}}
\end{table}
"""
    tex = tex.replace("\\\\\\\\", "\\\\")
    for a in ["textbf", "Delta", "le4", "operatorname", "mathbf", "tabcolsep", "textwidth"]:
        tex = tex.replace("\\\\" + a, "\\" + a)
    open(os.path.join(OUT, "general.tex"), "w").write(tex)


def numbers_macros():
    """Key numbers of the new analyses as LaTeX macros (docs/tables_v2/numbers.tex)."""
    L = []
    def m(name, val):
        L.append(f"\\newcommand{{\\{name}}}{{{val}}}")
    o = pd.read_csv(os.path.join(R, "graph_oracle_floor.csv")).set_index("dataset")
    for d, tag in [("esol", "Esol"), ("freesolv", "Fsv"), ("lipophilicity", "Lipo")]:
        m(f"oracleRmse{tag}", f"{o.loc[d, 'oracle_in_sample_rmse']:.2f}")
        m(f"oracleRsq{tag}", f"{o.loc[d, 'oracle_in_sample_r2']:.2f}")
        m(f"looRmse{tag}", f"{o.loc[d, 'oracle_loo_rmse']:.2f}")
    m("oracleAucBbbp", f"{o.loc['bbbp', 'oracle_in_sample']:.3f}")
    sm = pd.read_csv(os.path.join(R, "benchmark_v2_summary.csv")); cur = pd.read_csv(os.path.join(R, "benchmark_v2_curation.csv")).set_index("dataset")
    for d, tag in [("esol", "Esol"), ("freesolv", "Fsv"), ("lipophilicity", "Lipo")]:
        rd = sm[(sm.dataset == d) & (sm.model == "f_rdkit2d_rf")]["mean"].iloc[0]
        m(f"rdkitRsq{tag}", f"{1 - (rd / cur.loc[d, 'overall_target_sd']) ** 2:.2f}")
    v = pd.read_csv(os.path.join(R, "bid_variant_screen.csv")); pb = v[v.family == "published BID"]
    inst = pb["index"].unique(); base = {k.split(" (")[0].split(", p=")[0].split("_lambda")[0] for k in inst}
    m("nCensusInst", str(len(inst))); m("nCensusIdx", str(len(base)))
    m("nSomborInst", str(sum(1 for k in inst if "Sombor" in k or "SO" in k.split(" (")[0])))
    gs = pd.read_csv(os.path.join(R, "bid_generic_subsets.csv")).iloc[0]
    m("nGenericSubsets", f"{int(gs.full_rank_subsets):,}".replace(",", "\\,")); m("nAllSubsets", f"{int(gs.total_subsets):,}".replace(",", "\\,"))
    cp = os.path.join(R, "nonlinear_nullC_summary.csv")
    if os.path.exists(cp):
        c3 = pd.read_csv(cp)
        for nl, tag in [("C_nonlinear_inspan", "C"), ("D_entropy_type", "D")]:
            for d, dt in [("esol", "Esol"), ("freesolv", "Fsv"), ("lipophilicity", "Lipo"), ("bbbp", "Bbbp")]:
                rr = c3[(c3.null == nl) & (c3.dataset == d)]
                if len(rr): m(f"null{tag}{dt}", f"{rr.p99_abs.iloc[0]:.2f}")
            m(f"nullK{tag}", str(int(c3[c3.null == nl].K.iloc[0])))
        cc = pd.read_csv(os.path.join(R, "nonlinear_nullC_candidates.csv"))
        for d, k, tag in [("esol", "InfoH_deg", "InfoEsol"), ("bbbp", "FourCyc", "FourBbbp")]:
            rr = cc[(cc.dataset == d) & (cc.candidate == k)]
            if len(rr): m("pC" + tag, f"{rr.p_emp_nullC.iloc[0]:.3f}"); m("pD" + tag, f"{rr.p_emp_nullD.iloc[0]:.3f}")
    g = pd.read_csv(os.path.join(R, "certificate_generality.csv"))
    m("nGenDatasets", str(len(g))); m("nGenCertified", str(int(g.certificate_holds.sum())))
    m("nGenMolecules", f"{int(g.n_graphs.sum()):,}".replace(",", "\\,"))
    hiv = g[g.dataset == "HIV"].iloc[0] if (g.dataset == "HIV").any() else None
    if hiv is not None:
        m("hivHyper", str(int(hiv.n_maxdeg_gt4))); m("hivPairs", str(int(hiv.degree_pairs_realised))); m("hivRankM", str(int(hiv.rank_mij))); m("hivRankB", str(int(hiv.rank_BID18)))
        if "restricted_n" in g.columns and pd.notna(hiv.get("restricted_n", np.nan)):
            m("hivRestrictedN", f"{int(hiv.restricted_n):,}".replace(",", "\\,")); m("hivRestrictedResid", sci(hiv.restricted_max_rel_resid).strip("$"))
    m("genMaxResid", sci(g[g.certificate_holds].max_rel_resid_mij_on_BID18.max()).strip("$"))
    fp = os.path.join(R, "nonlinear_noise_floor_summary.csv")
    if os.path.exists(fp):
        s = pd.read_csv(fp); b = s[s.null == "B_inspan"]; a = s[s.null == "A_noise"]
        m("nullKB", str(int(b.K.iloc[0]))); m("nullKA", str(int(a.K.iloc[0])))
        m("nullPmin", f"{1.0 / (int(b.K.iloc[0]) + 1):.3f}")
        for d, tag in [("esol", "Esol"), ("freesolv", "Fsv"), ("lipophilicity", "Lipo"), ("bbbp", "Bbbp")]:
            m(f"nullB{tag}", f"{b[b.dataset == d].p99_abs.iloc[0]:.2f}"); m(f"nullA{tag}", f"{a[a.dataset == d].p99_abs.iloc[0]:.2f}")
        c = pd.read_csv(os.path.join(R, "nonlinear_noise_floor_candidates.csv"))
        m("nullSurvivors", str(int((c.q_bh < 0.05).sum())))
        for d, k, tag in [("esol", "InfoH_deg", "InfoEsol"), ("bbbp", "FourCyc", "FourBbbp"), ("freesolv", "InfoH_deg", "InfoFsv")]:
            rr = c[(c.dataset == d) & (c.candidate == k)]
            if len(rr):
                m("q" + tag, f"{rr.q_bh.iloc[0]:.2f}"); m("p" + tag, f"{rr.p_emp_nullB.iloc[0]:.3f}")

    S = pd.read_csv(os.path.join(R, "screen_vs_chemistry_summary.csv")).set_index("dataset")
    m("chemTopoBeyond", str(int(S.topo_pass_given_rdkit.sum())))
    gp = os.path.join(R, "screen_vs_chemistry_grouped_summary.csv")
    if os.path.exists(gp):
        G = pd.read_csv(gp).set_index("dataset"); CV = pd.read_csv(os.path.join(R, "screen_vs_chemistry_grouped_cv.csv"))
        cv = {(r.dataset, r.feature_set): r.value for r in CV.itertuples()}
        m("chemPassTotal", str(int(G.both_PASS.sum()))); m("chemPassRdTotal", str(int(G.rdkit_PASS.sum())))
        for d, tag in [("esol", "Esol"), ("freesolv", "Fsv"), ("lipophilicity", "Lipo"), ("bbbp", "Bbbp")]:
            m(f"chemPassRd{tag}", str(int(G.loc[d, "rdkit_PASS"]))); m(f"chemInc{tag}", str(int(G.loc[d, "both_INCONCLUSIVE"])))
        deltas = [cv[(d, "rdkit")] - min(cv[(d, "rdkit+topo")], cv[(d, "rdkit+topo+cands")]) for d in ("esol", "freesolv", "lipophilicity")]
        m("chemCvGainMax", f"{max(deltas):.3f}")
        m("chemCvAucGain", f"{max(cv[('bbbp','rdkit+topo')], cv[('bbbp','rdkit+topo+cands')]) - cv[('bbbp','rdkit')]:.3f}")
        m("chemCvLossLipo", f"{cv[('lipophilicity','rdkit+topo')] - cv[('lipophilicity','rdkit')]:.3f}")
    open(os.path.join(OUT, "numbers.tex"), "w").write("\n".join(L) + "\n")


if __name__ == "__main__":
    chem_table(); generality_table(); numbers_macros()
    fb_table()
    data_table()
    bench_table()
    cert_table(); census_table(); nonbid_table()
    print("tables written to", OUT)
