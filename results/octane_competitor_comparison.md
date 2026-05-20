# Octane competitor comparison

Pearson correlation $r$ between each index value and each physicochemical property of the $n = 18$ octane isomers.

Indices in the **classical** family are computed by the standard formulas implemented in `src/standard_indices.py`. Indices in the **Wuzi** family are computed from the closed form in `src/wuzi_index.py`. All correlations are computed on the same 18-row octane table (`results/octane_descriptors.csv`) so the comparison is head-to-head and not affected by differences in dataset between sources.


## Combined comparison table

Sign convention: a positive $r$ means the index is positively correlated with the property; a negative $r$ means negatively correlated. What matters chemically is $|r|$ (correlation strength); the sign depends on the index definition.

| Index | Family | $r$ with T_B | $r$ with Delta H_f | $r$ with Delta H_vap | $r$ with S | $r$ with omega |
|---|---|---:|---:|---:|---:|---:|
| M1 | classical | -0.718 | -0.762 | -0.912 | -0.973 | -0.970 |
| M2 | classical | -0.503 | -0.553 | -0.775 | -0.942 | -0.987 |
| R | classical | +0.821 | +0.848 | +0.956 | +0.935 | +0.901 |
| SO | classical | -0.748 | -0.793 | -0.926 | -0.969 | -0.957 |
| GA | classical | +0.823 | +0.856 | +0.959 | +0.942 | +0.909 |
| H | classical | +0.822 | +0.846 | +0.959 | +0.935 | +0.904 |
| ABC | classical | -0.869 | -0.892 | -0.944 | -0.858 | -0.790 |
| AZI | classical | +0.936 | +0.922 | +0.929 | +0.734 | +0.656 |
| W(1, 0, 0) | Wuzi | -0.503 | -0.553 | -0.775 | -0.942 | -0.987 |
| W(-0.5, 0, 0) | Wuzi | +0.821 | +0.848 | +0.956 | +0.935 | +0.901 |
| W(0, -0.5, 0) | Wuzi | +0.801 | +0.829 | +0.953 | +0.950 | +0.927 |
| W(0, -1, 0) | Wuzi | +0.822 | +0.846 | +0.959 | +0.935 | +0.904 |
| W(0, 0, 1) | Wuzi | -0.828 | -0.831 | -0.966 | -0.931 | -0.930 |
| W(0, 0, 2) | Wuzi | -0.830 | -0.844 | -0.967 | -0.939 | -0.926 |
| W(1, 1, 1) | Wuzi | -0.615 | -0.677 | -0.840 | -0.959 | -0.981 |
| W(-1, -1, 1) | Wuzi | +0.763 | +0.830 | +0.747 | +0.632 | +0.494 |

## Best-correlated index per property (by $|r|$)

| Property | Best |$r$|  classical | Best |$r$|  Wuzi |
|---|---|---|
| T_B | AZI (+0.936) | W(0, 0, 2) (-0.830) |
| Delta H_f | AZI (+0.922) | W(-0.5, 0, 0) (+0.848) |
| Delta H_vap | H (+0.959) | W(0, 0, 2) (-0.967) |
| S | M1 (-0.973) | W(1, 1, 1) (-0.959) |
| omega | M2 (-0.987) | W(1, 0, 0) (-0.987) |

## Interpretation

The table is sensitivity-only: it shows how strongly each index covaries with each physicochemical property on the octane benchmark. It is NOT a QSPR utility certification (see the companion methodology paper for that): a single scalar index that correlates strongly with one property on 18 isomers can still be redundant with the classical baseline in a multivariate sense.

Three honest readings of the table:

1. The Wuzi family contains, at specific parameter settings, the classical indices $M_2$ (at $(1, 0, 0)$), $R$ (at $(-1/2, 0, 0)$), SCI (at $(0, -1/2, 0)$), and $H/2$ (at $(0, -1, 0)$). At these points the Wuzi correlations agree with the corresponding classical values up to floating point.
2. The strongest absolute correlations on most properties are achieved by classical indices well established in the math-chem literature ($mM_2$, $AZI$, $ABC$ on octanes), not by the Wuzi family at any of the sampled parameter triples. We do not claim Wuzi to be a superior octane-property predictor.
3. The Wuzi family at the mixed-sign triple $(-1, -1, 1)$ produces correlation magnitudes in the same general regime as established BID indices on this benchmark, supporting the use of the family as a worked case study of a parametric BID generalization without claiming predictive supremacy.

This table is the entirety of the predictive evidence reported in the Wuzi paper. The multi-dataset orthogonality screening of the family and the downstream RandomForest benchmark live in the companion methodology paper (Paper 2).
