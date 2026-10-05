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
\caption{Span certificate on the four datasets. $|P_D|$: realised degree pairs; ranks are numerical ranks of the column-centred, column-normalised matrices (relative tolerance $10^{-9}$). ``Resid.'' columns give the largest relative residual $\|r\|/\|z-\bar z\|$ on the baseline of the ten counts $m_{ij}$, of the $14$ published BID indices of Table~\ref{tab:census}, and of the $100$ Loyola grid points.}
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
\caption{Census of published BID indices ($f(a,b)$ with $a,b$ the end degrees of an edge) and two non-BID Sombor-type controls. Columns 3--6: largest $|r|$ with any of the $30$ baseline indices on each dataset (the conventional pairwise screen). Last column: largest relative residual on the baseline over the four datasets . For KG-Sombor $e=a+b-2$ is the edge degree. Every BID index is in the span to machine precision; the controls are not.}
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
    a, b = rng("REF_invNirmala"); c1, c2 = rng("REF_noise"); l1, l2 = rng("REF_leak")
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
        rows.append(f"{NICE[d]} & {target[d]} & ${int(r.n_raw)}$ & ${int(r.n_invalid_smiles)}$ & "
                    f"${int(r.n_multifragment_stripped)}$ & ${int(r.n_rows_in_dup_groups)}$ & "
                    f"${int(r.n_conflicting_dup_groups_dropped)}$ & ${int(r.n_final)}$ & "
                    f"${int(r.n_distinct_graphs)}$ & ${r.pct_mols_graph_shared:.1f}$\\\\")
    tex = r"""\begin{table}[h]
\centering
\small
\caption{Datasets and curation. Raw rows of the MoleculeNet files; invalid SMILES; multi-fragment SMILES reduced to the largest fragment; rows in exact canonical-SMILES duplicate groups; duplicate groups dropped for conflicting labels; molecules after curation (also dropping graphs with fewer than two heavy atoms); distinct hydrogen-suppressed unlabelled graphs (Weisfeiler--Lehman hash); and the percentage of curated molecules whose graph is shared with another molecule.}
\label{tab:data}
\setlength{\tabcolsep}{3pt}
\begin{tabular}{llrrrrrrrr}
\toprule
Dataset & Target & raw & invalid & multi-frag. & dup.\ rows & conflicts & curated & graphs & shared (\%)\\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}
\end{table}
"""
    open(os.path.join(OUT, "data.tex"), "w").write(tex)


if __name__ == "__main__":
    data_table()
    cert_table(); census_table(); nonbid_table()
    print("tables written to", OUT)
