# Figure captions

Manuscript-ready captions for the figures in `figures/`. Captions are
phrased to be self-contained: they state what is plotted, the
dataset, the sample size, the metric, and a one-sentence
interpretation. Each caption ends with a brief warning against
overclaiming where appropriate.

---

**Figure 2.** Cumulative explained variance of principal-component
analysis on the 30-index baseline, for three QSPR benchmarks (ESOL,
$n=1127$; FreeSolv, $n=639$; Lipophilicity, $n=4200$). Indices were
standardized to unit variance before PCA. Three components capture
95% of the variance on every dataset, and five components capture
99%. The effective rank of the 30-index baseline is therefore much
smaller than 30; *this does not mean classical indices are
uninformative*, only that on these chemistry benchmarks they live
on a low-dimensional manifold.

**Figure 3.** Number of redundant index pairs out of the
$\binom{30}{2} = 435$ baseline pairs, at two correlation thresholds
($\lvert r\rvert \geq 0.90$ and $\lvert r\rvert \geq 0.95$), for
ESOL ($n=1127$), FreeSolv ($n=639$), and Lipophilicity ($n=4200$).
Roughly 37% of all pairs of classical indices are above the 0.95
threshold on every dataset. The fraction of redundant pairs is
remarkably stable across the three independent benchmarks.

**Figure 4.** Signed Pearson correlation $r$ between each of 34
indices (the Wuzi family at eight canonical parameter points plus
the 30-index classical baseline, less the Padmakar--Ivan index $PI$
which is constant on the 18 octanes) and each of the five octane
physicochemical properties ($T_B$, $\Delta H_f$, $\Delta H_{\mathrm{vap}}$,
$S$, $\omega$) on the 18 isomers. Rows are sorted by maximum
$\lvert r\rvert$. Color encodes $\lvert r\rvert$ while the printed
value carries the sign. The $PI$ degeneracy on octane trees is
$\sigma = 0$, mean $42$. Strong performance on this small dataset
is *necessary but not sufficient* for general QSPR utility (cf.
Figure 8 and the discussion in Section 6).

**Figure 5b** (the only fig5 file shipped; the earlier ESOL-only
Figure 5 was archived during round 4 and is replaced by 5b below).
Wuzi-family redundancy screen extended to all three QSPR benchmarks.
Each row corresponds to one dataset (ESOL, FreeSolv,
Lipophilicity); each column to one $\gamma$ slice. Only one of the
300 grid points falls below the $\lvert r\rvert = 0.95$ threshold
— the FreeSolv cell at $(\alpha, \beta, \gamma) = (-1, -1, 1)$,
marked with a star, where $\max\lvert r\rvert = 0.943$. Its partial
correlation with $\Delta G_{\mathrm{hyd}}$ after controlling for
the 30 baseline indices is $\approx 0.036$, so passing the
correlation threshold here does not translate into meaningful
independent QSPR signal.

**Figure 6.** Percent degeneracy on the 106 non-isomorphic trees of
order 10, for each of the 30 baseline indices and a selection of
Wuzi parameter points. Lower bars correspond to indices that better
discriminate among the 106 distinct tree structures. Spectral and
distance-based indices (Balaban $J$, Estrada, Graph Energy,
Spectral Radius, Hyper-Wiener) are the most discriminating;
several Wuzi parameter points sit in the same regime as classical
degree-based indices like Randić and geometric-arithmetic. This
plot reflects discrimination on a single small graph class and
should not be over-interpreted as a general ranking.

**Figure 7.** Structure-sensitivity metrics on the 75 decane
isomers, top 20 indices by $SS = \sigma / \mu$. $Abr$ is the
abruptness $(\max - \min) / \mu$, and $SA = SS / Abr$. Higher $SS$
indicates more variation across isomers relative to the mean,
which is one (necessary but not sufficient) precondition for
isomer-discriminating descriptors. Following the convention of
Movahedi et al. (2026) §5.4.

**Figure 8.** Pairwise Pearson correlation matrix of the 30
baseline indices on ESOL ($n=1127$). Of the $\binom{30}{2} = 435$
distinct pairs, 159 have $\lvert r\rvert \geq 0.95$, visible as
the dense red core of the matrix. Degree-based indices form a
particularly tightly correlated block; distance-based indices are
mutually correlated as well but weakly correlated with the
degree-based block. Balaban $J$ sits apart from both groups. The
visualization motivates the use of partial-correlation and
PCA-based screening (Figures 2–3) rather than relying on simple
counts of redundant pairs.

**Figure 9.** Downstream ML benchmark: RandomForest performance
(top row) and feature count (bottom row) under four feature
configurations on four datasets (ESOL, FreeSolv, Lipophilicity
regression; BBBP classification, $n=2050$). The primary metric is
RMSE on the three regression panels and ROC-AUC on the BBBP panel.
The four configurations are: `full` (all 30 baseline indices),
`PCA-95` (top-$k$ principal components capturing $\ge 95\%$ of
the baseline variance), `pair` (greedy redundancy filter at
$|r| < 0.95$), and `combined` (pairwise filter plus
$|\mathrm{pcor}(z, y \mid X)| \ge 0.10$). The pairwise screen
matches or marginally improves the full-baseline RandomForest on
every dataset with $4$--$5\times$ fewer features. On BBBP the
pairwise-pruned model attains ROC-AUC $0.860$, marginally exceeding
the full-baseline $0.846$. The strict combined screen removes
every feature on Lipophilicity (no pairwise-pruned index has
$|\mathrm{pcor}|\ge 0.10$ with $\log D$); this is reported as
"killed" in the figure and discussed in the manuscript as an
honest failure mode of the strict combined criterion on
weakly-target-correlated regression tasks.

**Figure 10.** Maximizers of the Wuzi index
$W(T; -1, -1, 1) = \sum_{uv \in E(T)} (d_u d_v)^{-1} (d_u + d_v)^{-1}
\exp(|d_u - d_v|/(d_u + d_v))$ over all non-isomorphic trees of order
$n$, for $n = 10, 11, 12$, drawn as rooted trees with the
highest-degree vertex at the top. Internal vertices are coloured
green and leaves cream. None of the three trees is a caterpillar
(the non-leaf vertices do not form a path). In the language of
Conjecture conj:caterpillar_max (Section 7): $n = 11$ is the spider
$S^{(2)}_5$ with hub degree $(11-1)/2 = 5$ and length-$2$ spokes;
$n = 10$ is the balanced double-spider $DS(3, 3)$ (two hubs of
degree $3$, joined by an edge); $n = 12$ is $DS(4, 3)$. The
\texttt{graph6} strings above each panel uniquely identify each
tree up to isomorphism (degree sequences alone do not, for
$n \ge 6$). Computed from
\texttt{results/wuzi\_caterpillar\_max\_n5\_to\_20.csv}; rendering
script lives in commit \texttt{cbed908}, dataset extended at
\texttt{7f6c44a}. The same pattern (spider for odd $n$, balanced
double-spider for even $n$) holds at every enumerated $n \ge 9$
through $n = 20$.
