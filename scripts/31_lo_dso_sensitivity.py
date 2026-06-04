"""DSO-style sensitivity analysis for the Loyola (LO) index paper.

Honest, reproducible, publication-safe. NO post-hoc decimal-tuned triples.
Only the classical indices LO reduces to, plus the two PRE-SPECIFIED pure-gamma
points LO(0,0,1), LO(0,0,2). Modeled on the four subsections of the DSO paper
(Movahedi-Gutman-Redzepovic-Furtula, MATCH 95 (2026) 141-162).

Graph engine verified: parser SO matches octane_descriptors.csv SO to ~1e-15.
Property data: NIST WebBook columns in octane_descriptors.csv (documented; DSO
used a different source [moleculardescriptors.eu], so DH_vap / entropy
correlations may differ by source -- NOT an index-computation error).

Indices (consistent order everywhere). Edge count m = LO(0,0,0) is EXCLUDED
(constant on octanes -> undefined correlation).
"""
from __future__ import annotations
import csv, math, os, re
import numpy as np
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJ=os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),".."))
RES=os.path.join(PROJ,"results"); FIG=os.path.join(PROJ,"figures")
PROPS=["T_B","dHf","dHvap","S","omega"]
PROP_LABEL={"T_B":r"$T_B$","dHf":r"$\Delta H_f$","dHvap":r"$\Delta H_{\rm vap}$","S":r"$S$","omega":r"$\omega$"}

# index name -> (display label, edge-pair function f(a,b)), in fixed order
def _ga(a,b): return 2*np.sqrt(a*b)/(a+b)
def _ag(a,b): return (a+b)/(2*np.sqrt(a*b))
INDICES=[
 ("M1",   r"$M_1$",        lambda a,b: a+b),
 ("M2",   r"$M_2$",        lambda a,b: a*b),
 ("HM",   r"$HM$",         lambda a,b: (a+b)**2),
 ("mM2",  r"${}^mM_2$",    lambda a,b: 1.0/(a*b)),
 ("R",    r"$R$",          lambda a,b: 1.0/np.sqrt(a*b)),
 ("chi",  r"$\chi$",       lambda a,b: 1.0/np.sqrt(a+b)),
 ("H/2",  r"$H/2$",        lambda a,b: 1.0/(a+b)),
 ("ISI",  r"$ISI$",        lambda a,b: (a*b)/(a+b)),
 ("GA",   r"$GA$",         _ga),
 ("AG",   r"$AG$",         _ag),
 ("LO(0,0,1)", r"$LO(0,0,1)$", lambda a,b: np.exp(np.abs(a-b)/(a+b))),
 ("LO(0,0,2)", r"$LO(0,0,2)$", lambda a,b: np.exp(2*np.abs(a-b)/(a+b))),
]
NAMES=[n for n,_,_ in INDICES]; LABELS=[l for _,l,_ in INDICES]

TWO={'Cl','Br'}
def heavy_adj(smi):
    adj={}; idx=[-1]
    def atom(): idx[0]+=1; adj[idx[0]]=set(); return idx[0]
    prev=[None];stack=[];rings={};nf=[False];i=0;n=len(smi)
    def bond(a,b):
        if a is not None and b is not None and a!=b: adj[a].add(b);adj[b].add(a)
    while i<n:
        ch=smi[i]
        if ch=='[':
            j=smi.index(']',i);tok=smi[i+1:j];i=j+1
            m=re.match(r'\d*([A-Z][a-z]?|[bcnops])',tok);el=m.group(1) if m else 'X'
            if el=='H':continue
            v=atom()
            if not nf[0]:bond(prev[0],v)
            prev[0]=v;nf[0]=False;continue
        if ch in '-=#:/\\':i+=1;continue
        if ch=='(':stack.append(prev[0]);i+=1;continue
        if ch==')':prev[0]=stack.pop();i+=1;continue
        if ch=='.':nf[0]=True;i+=1;continue
        if ch=='%':
            num=smi[i+1:i+3];i+=3
            bond(prev[0],rings.pop(num)) if num in rings else rings.__setitem__(num,prev[0]);continue
        if ch.isdigit():
            bond(prev[0],rings.pop(ch)) if ch in rings else rings.__setitem__(ch,prev[0]);i+=1;continue
        two=smi[i:i+2]
        if two in TWO:i+=2
        elif ch in 'BCNOPSFI' or ch in 'bcnops':i+=1
        else:i+=1;continue
        v=atom()
        if not nf[0]:bond(prev[0],v)
        prev[0]=v;nf[0]=False
    return adj
def pairs_from_smiles(smi):
    adj=heavy_adj(smi);seen=set();out=[]
    for u in adj:
        for v in adj[u]:
            e=(min(u,v),max(u,v))
            if e in seen:continue
            seen.add(e);out.append((len(adj[u]),len(adj[v])))
    return out
def pairs_from_graph(g):
    return [(g.degree(u),g.degree(v)) for u,v in g.edges()]
def index_value(pairs,f):
    return float(sum(f(a,b) for a,b in pairs))
def corr(x,y):
    x=x-x.mean();y=y-y.mean();s=np.linalg.norm(x)*np.linalg.norm(y)
    return float('nan') if s==0 else float(x@y/s)

def main():
    os.makedirs(FIG,exist_ok=True)
    # ---------- load octanes ----------
    oct=list(csv.DictReader(open(os.path.join(RES,"octane_descriptors.csv"))))
    opairs=[pairs_from_smiles(r["smiles"]) for r in oct]
    Y={p:np.array([float(r[p]) for r in oct]) for p in PROPS}
    OV={n:np.array([index_value(p,f) for p in opairs]) for n,_,f in INDICES}

    # ===== save octane descriptors (12 indices + 5 props) =====
    with open(os.path.join(RES,"lo_sensitivity_octane_descriptors.csv"),"w",newline="") as fh:
        w=csv.writer(fh); w.writerow(["name","smiles"]+NAMES+PROPS)
        for k,r in enumerate(oct):
            w.writerow([r["name"],r["smiles"]]+[f"{OV[n][k]:.6f}" for n in NAMES]+[r[p] for p in PROPS])

    # ===== 5.1 prediction correlations =====
    C=np.array([[corr(OV[n],Y[p]) for p in PROPS] for n in NAMES])   # (12,5)
    best=[int(np.nanargmax(np.abs(C[:,j]))) for j in range(len(PROPS))]
    with open(os.path.join(RES,"lo_sensitivity_prediction_correlations.csv"),"w",newline="") as fh:
        w=csv.writer(fh); w.writerow(["index"]+PROPS+["best_in"])
        for i,n in enumerate(NAMES):
            w.writerow([n]+[f"{C[i,j]:+.4f}"+("*" if best[j]==i else "") for j in range(len(PROPS))]
                       +["; ".join(PROPS[j] for j in range(len(PROPS)) if best[j]==i)])
    # LaTeX table (bold best per column)
    with open(os.path.join(RES,"lo_sensitivity_prediction_table_latex.tex"),"w") as fh:
        fh.write("\\begin{table}[H]\n\\centering\n")
        fh.write("\\caption{Pearson correlation $r$ between selected Loyola reductions and pre-specified pure-$\\gamma$ Loyola points and the five physicochemical properties of the 18 octane isomers (NIST WebBook). Bold = best $|r|$ per property. Edge count $m=LO(0,0,0)$ is omitted (constant on octanes).}\n")
        fh.write("\\label{tab:lo_octane}\n\\small\n\\begin{tabular}{lrrrrr}\n\\toprule\n")
        fh.write("Index & "+" & ".join(PROP_LABEL[p] for p in PROPS)+"\\\\\n\\midrule\n")
        for i,(n,lab,_) in enumerate(INDICES):
            cells=[]
            for j in range(len(PROPS)):
                s=f"{C[i,j]:+.3f}"
                cells.append("\\textbf{"+s+"}" if best[j]==i else s)
            fh.write(lab+" & "+" & ".join(cells)+"\\\\\n")
        fh.write("\\bottomrule\n\\end{tabular}\n\\end{table}\n")

    # ===== 5.2 intercorrelation =====
    IC=np.array([[corr(OV[a],OV[b]) for b in NAMES] for a in NAMES])
    with open(os.path.join(RES,"lo_sensitivity_intercorrelation_matrix.csv"),"w",newline="") as fh:
        w=csv.writer(fh); w.writerow([""]+NAMES)
        for i,n in enumerate(NAMES): w.writerow([n]+[f"{IC[i,j]:+.4f}" for j in range(len(NAMES))])
    # how collinear are the pure-gamma points with the classical reductions?
    classic_idx=[i for i,n in enumerate(NAMES) if not n.startswith("LO(")]
    gamma_collinear={}
    for n in ["LO(0,0,1)","LO(0,0,2)"]:
        gi=NAMES.index(n); mx=max(abs(IC[gi,c]) for c in classic_idx)
        gamma_collinear[n]=mx

    # ===== 5.3 degeneracy on order-10 trees =====
    trees=list(nx.nonisomorphic_trees(10)); Nt=len(trees)
    tpairs=[pairs_from_graph(g) for g in trees]
    deg={}
    for n,_,f in INDICES:
        vals=np.array([index_value(p,f) for p in tpairs])
        distinct=len(set(np.round(vals,6)))
        deg[n]=100.0*(1-distinct/Nt)
    with open(os.path.join(RES,"lo_sensitivity_degeneracy_order10_trees.csv"),"w",newline="") as fh:
        w=csv.writer(fh); w.writerow(["index","N_trees","distinct_values","degeneracy_pct"])
        for n,_,f in INDICES:
            vals=np.array([index_value(p,f) for p in tpairs]); d=len(set(np.round(vals,6)))
            w.writerow([n,Nt,d,f"{deg[n]:.2f}"])

    # ===== 5.4 structure sensitivity on decanes =====
    dec=list(csv.DictReader(open(os.path.join(RES,"decane_smiles.csv"))))
    dpairs=[pairs_from_smiles(r["smiles"]) for r in dec]; Nd=len(dec)
    SSrows={}
    for n,_,f in INDICES:
        vals=np.array([index_value(p,f) for p in dpairs]); mu=vals.mean()
        SS=vals.std(ddof=0)/mu; Abr=(vals.max()-vals.min())/mu; SA=SS/Abr if Abr else float('nan')
        SSrows[n]=(SS,Abr,SA)
    with open(os.path.join(RES,"lo_sensitivity_structure_sensitivity_decanes.csv"),"w",newline="") as fh:
        w=csv.writer(fh); w.writerow(["index","SS","Abr","SA"])
        for n,_,f in INDICES:
            SS,Abr,SA=SSrows[n]; w.writerow([n,f"{SS:.5f}",f"{Abr:.5f}",f"{SA:.5f}"])

    # ================= FIGURES =================
    # Fig 1: prediction correlation heatmap
    fig,ax=plt.subplots(figsize=(6.2,5.4))
    im=ax.imshow(C,cmap="RdBu_r",vmin=-1,vmax=1,aspect="auto")
    ax.set_xticks(range(len(PROPS))); ax.set_xticklabels([PROP_LABEL[p] for p in PROPS])
    ax.set_yticks(range(len(NAMES))); ax.set_yticklabels(LABELS,fontsize=9)
    for i in range(len(NAMES)):
        for j in range(len(PROPS)):
            txt=f"{C[i,j]:+.2f}"
            ax.text(j,i,txt,ha="center",va="center",fontsize=7,
                    color="white" if abs(C[i,j])>0.6 else "black",
                    fontweight="bold" if best[j]==i else "normal")
            if best[j]==i: ax.text(j,i-0.34,"$\\ast$",ha="center",va="center",fontsize=10,color="black")
    fig.colorbar(im,ax=ax,fraction=0.046,pad=0.04,label="Pearson $r$")
    ax.set_title("Correlation with octane properties\n(best $|r|$ per property marked $\\ast$)",fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(FIG,"fig_lo_prediction_correlation_heatmap.png"),dpi=200)
    fig.savefig(os.path.join(FIG,"fig_lo_prediction_correlation_heatmap.pdf")); plt.close(fig)

    # Fig 2: intercorrelation heatmap
    fig,ax=plt.subplots(figsize=(6.8,5.8))
    im=ax.imshow(IC,cmap="RdBu_r",vmin=-1,vmax=1,aspect="auto")
    ax.set_xticks(range(len(NAMES))); ax.set_xticklabels(LABELS,rotation=90,fontsize=8)
    ax.set_yticks(range(len(NAMES))); ax.set_yticklabels(LABELS,fontsize=8)
    fig.colorbar(im,ax=ax,fraction=0.046,pad=0.04,label="Pearson $r$")
    ax.set_title("Intercorrelation on the 18 octane isomers",fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(FIG,"fig_lo_intercorrelation_heatmap.png"),dpi=200)
    fig.savefig(os.path.join(FIG,"fig_lo_intercorrelation_heatmap.pdf")); plt.close(fig)

    # Fig 3: degeneracy bars
    fig,ax=plt.subplots(figsize=(7.2,4.2))
    vals=[deg[n] for n in NAMES]
    ax.bar(range(len(NAMES)),vals,color="#4477AA",edgecolor="black",linewidth=0.4)
    ax.set_xticks(range(len(NAMES))); ax.set_xticklabels(LABELS,rotation=45,ha="right",fontsize=8)
    ax.set_ylabel("Degeneracy (\\%)"); ax.set_title("Degeneracy on the 106 order-10 trees (lower = better discrimination)",fontsize=9)
    for i,v in enumerate(vals): ax.text(i,v+0.4,f"{v:.1f}",ha="center",fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(FIG,"fig_lo_degeneracy_order10_trees.png"),dpi=200)
    fig.savefig(os.path.join(FIG,"fig_lo_degeneracy_order10_trees.pdf")); plt.close(fig)

    # Fig 4: structure sensitivity 3 panels
    fig,axes=plt.subplots(1,3,figsize=(12,4.2))
    for ax,k,title in zip(axes,[0,1,2],["SS","Abr","SA"]):
        vv=[SSrows[n][k] for n in NAMES]
        ax.bar(range(len(NAMES)),vv,color="#4477AA",edgecolor="black",linewidth=0.4)
        ax.set_xticks(range(len(NAMES))); ax.set_xticklabels(LABELS,rotation=90,fontsize=7)
        ax.set_title(title,fontsize=11)
    fig.suptitle("Structure-sensitivity metrics on the 75 decane isomers (standard formulas; see methods note)",fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(FIG,"fig_lo_structure_sensitivity_decanes.png"),dpi=200)
    fig.savefig(os.path.join(FIG,"fig_lo_structure_sensitivity_decanes.pdf")); plt.close(fig)

    # ================= console summary =================
    print("OCTANE prediction (best |r| per property):")
    for j,p in enumerate(PROPS):
        i=best[j]; print(f"  {p:7}: {NAMES[i]:10} r={C[i,j]:+.3f}")
    print("\nGamma pure-point collinearity with classical reductions (max |r|):")
    for n,mx in gamma_collinear.items(): print(f"  {n}: max|r| with a classic = {mx:.4f}")
    print(f"\nDegeneracy (order-10 trees, N={Nt}):")
    for n in NAMES: print(f"  {n:10}: {deg[n]:.2f}%")
    print(f"\nStructure sensitivity (decanes, N={Nd}) [SS, Abr, SA]:")
    for n in NAMES: s=SSrows[n]; print(f"  {n:10}: {s[0]:.4f} {s[1]:.4f} {s[2]:.4f}")
    print("\nwrote 6 CSV/TeX files + 4 figures (png+pdf).")

if __name__=="__main__": main()
