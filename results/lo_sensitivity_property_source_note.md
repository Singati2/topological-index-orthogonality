# Octane property-source note

The five octane physicochemical properties used in Section 5.1 — boiling point $T_B$,
enthalpy of formation $\Delta H_f$, heat of vaporization $\Delta H_{\rm vap}$, entropy
$S$, acentric factor $\omega$ — are taken from the **NIST WebBook** columns already
present in `results/octane_descriptors.csv`.

**Important:** the DSO paper (Movahedi–Gutman–Redžepović–Furtula, *MATCH* 95 (2026)
141–162) took its octane properties from a *different* source — reference [14] there is
the moleculardescriptors.eu / Milano Chemometrics octane dataset
(`web.archive.org/.../moleculardescriptors.eu/dataset/dataset.htm`).

A cross-check on the plain Sombor index `SO` (identical formula in both works) shows our
correlations **match the DSO paper exactly on $T_B$, $\Delta H_f$, and $\omega$**
(differences < 0.003) but **differ on $\Delta H_{\rm vap}$ and $S$** (ours
$-0.926/-0.969$, DSO $-0.903/-0.947$). Because the `SO` *index values* match on three of
five properties, our index computation is correct (a code error would shift all five);
the difference on $\Delta H_{\rm vap}$ and $S$ therefore reflects the **property-data
source**, not an index-computation error.

We did **not** reconstruct, infer, or approximate the moleculardescriptors.eu property
values. We therefore make **no** claim that our table is "aligned to DSO" or that "the
SO row matches DSO." If the actual Milano property values are obtained and verified, the
table can be recomputed on that source for direct comparability; until then it is on
NIST, documented here.
