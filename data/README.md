# Datasets

Four MoleculeNet benchmarks (Wu et al., *Chem. Sci.* 9 (2018) 513–530),
redistributed here unchanged from the public DeepChem/MoleculeNet mirror so
that every result in the paper can be reproduced offline. The terms of the
original sources apply; please cite them (see below) if you use the data.

| File | Source URL (DeepChem mirror) | Rows | SMILES column | Target column | Original source |
|---|---|---|---|---|---|
| `esol.csv` | `.../datasets/delaney-processed.csv` | 1128 | `smiles` | `measured log solubility in mols per litre` | Delaney, *J. Chem. Inf. Comput. Sci.* 44 (2004) 1000–1005 |
| `freesolv.csv` | `.../datasets/SAMPL.csv` | 642 | `smiles` | `expt` (experimental hydration free energy, kcal/mol) | Mobley & Guthrie, *J. Comput.-Aided Mol. Des.* 28 (2014) 711–720 |
| `lipophilicity.csv` | `.../datasets/Lipophilicity.csv` | 4200 | `smiles` | `exp` (log D at pH 7.4) | ChEMBL document CHEMBL3301361 (AstraZeneca) |
| `bbbp.csv` | `.../datasets/BBBP.csv` | 2050 | `smiles` | `p_np` (0/1) | Martins et al., *J. Chem. Inf. Model.* 52 (2012) 1686–1697 |

Mirror base URL: `https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/`.

## Integrity

`SHA256SUMS` lists the SHA-256 of each raw file. On 2026-10-05 each file was
re-downloaded from the mirror and found byte-identical
(`python scripts/76_verify_datasets.py --check-remote`; report in
`results/dataset_verification.md`). Check locally with
`cd data && shasum -a 256 -c SHA256SUMS`.

Verified properties: no missing SMILES or targets; BBBP labels are 0/1;
11 BBBP SMILES do not parse in RDKit; 105 BBBP SMILES and 1 Lipophilicity
SMILES contain several fragments (the largest is kept).

## Curated data (`curated/`)

`<name>_curated.csv` has columns `canonical_smiles`, `target`,
`graph_wl_hash` (Weisfeiler–Lehman hash, 4 iterations, of the
hydrogen-suppressed, atom- and bond-unlabelled graph). Curation rule
(`scripts/74_dedup_sensitivity.py`, identical to `scripts/72_benchmark_v2.py`):
keep the largest fragment; RDKit canonical SMILES; average regression
duplicates whose labels agree within 0.5 target units and drop duplicate
groups that disagree by more (or carry conflicting BBBP labels); drop graphs
with fewer than two heavy atoms.

| Dataset | Curated molecules | Distinct graphs | Molecules sharing a graph |
|---|---|---|---|
| ESOL | 1115 | 628 | 55.5% |
| FreeSolv | 639 | 236 | 74.3% |
| Lipophilicity | 4200 | 3472 | 28.2% |
| BBBP | 1955 | 1693 | 22.6% |
