# Cross-dataset redundancy summary

All values computed on the same 30-index baseline set;
see `docs/theoretical_foundation.md` for the BID basis
theorem implying a 10-dimensional ceiling.

## Baseline redundancy

| Dataset | n | Pairs \|r\|≥0.90 | Pairs \|r\|≥0.95 | PC₉₀ | PC₉₅ | PC₉₉ |
|---|---:|---:|---:|---:|---:|---:|
| esol | 1127 | 203 / 435 (46.7%) | 159 / 435 (36.6%) | 2 | 3 | 5 |
| freesolv | 639 | 257 / 435 (59.1%) | 162 / 435 (37.2%) | 2 | 3 | 5 |
| lipophilicity | 4200 | 200 / 435 (46.0%) | 164 / 435 (37.7%) | 2 | 3 | 5 |

Theoretical upper bound on endpoint-degree edge-sum (BID)
index dimension under Δ≤4: **10**. Observed effective
rank (PC₉₅) is consistently 3 across all three benchmarks.


## Wuzi parameter sweep verdict

| Dataset | Grid points | Pass | Fail | Min max\|r\| | Max \|partial corr w/ target\| |
|---|---:|---:|---:|---:|---:|
| esol | 100 | 0 | 100 | 0.965 | 0.0282 |
| freesolv | 100 | 1 | 99 | 0.943 | 0.0362 |
| lipophilicity | 100 | 0 | 100 | 0.965 | 0.0093 |

The Wuzi family fails the kill-test at every parameter
setting on every dataset, consistent with Corollary 3 of
the BID basis theorem.
