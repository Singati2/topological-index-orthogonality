# Why decimal-tuned parameter optimization was excluded

The Loyola index $LO(G;\alpha,\beta,\gamma)$ is a three-parameter family. One could,
after computing the correlations of the classical indices with the octane properties,
search a fine grid of $(\alpha,\beta,\gamma)$ and report whichever triple maximises the
correlation for each property (e.g. `LO(1.05,0,0)` instead of `LO(1,0,0)=M2`).

**This was deliberately *not* done, because it is in-sample overfitting.** The 18 octane
isomers are simultaneously the tuning set and the evaluation set. Optimising three
continuous parameters against the same 18 points is guaranteed to match or beat any
fixed index — the same "win" is obtained for *any* target, including random noise — so
it measures nothing about the index. We verified this directly: an in-sample decimal
tweak "beats" the best classic on all five properties, but under leakage-free
leave-one-out validation the apparent advantage disappears on 4 of 5 properties.

Accordingly, the sensitivity analysis uses **only**:

1. the classical indices that $LO$ reduces to **exactly or up to a constant factor**
   ($M_1, M_2, HM, {}^mM_2, R, \chi, H/2, ISI, GA, AG$), and
2. the **pre-specified** pure-$\gamma$ Loyola points $LO(0,0,1)$ and $LO(0,0,2)$.

No parameter value was chosen after inspecting the property correlations. The edge
count $m=LO(0,0,0)$ is excluded because it is constant on the octanes (undefined
correlation).

A separate, leakage-free large-$n$ benchmark (tuned vs. tuned, on MoleculeNet datasets)
is reported elsewhere and is the appropriate place for any tuning-based claim; it is
kept out of this section.
