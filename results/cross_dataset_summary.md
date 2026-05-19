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


## Wuzi parameter sweep — redundancy-screen verdict

Verdict is **PASS** if max |r| with every classical baseline
is strictly below 0.95, **FAIL** otherwise. The threshold
|r| = 0.95 is conventional; the substantive question is
whether the partial correlation with the target property
is meaningfully different from zero after controlling for
the baseline (rightmost column).

| Dataset | Grid points | Pass | Fail | Min max\|r\| | Max \|partial corr w/ target\| |
|---|---:|---:|---:|---:|---:|
| esol | 100 | 0 | 100 | 0.965 | 0.0282 |
| freesolv | 100 | 1 | 99 | 0.943 | 0.0362 |
| lipophilicity | 100 | 0 | 100 | 0.965 | 0.0093 |

On ESOL and Lipophilicity, every one of the 100 grid
points is highly correlated with at least one classical
baseline (|r| ≥ 0.95). On FreeSolv, a single grid point
falls just below the threshold, but its partial correlation
with the target after controlling for the 30 baselines is
≈ 0.04, so it does not provide useful independent QSPR signal.
This behavior is consistent with the BID-basis observation
(`docs/theoretical_foundation.md`): every endpoint-degree
edge-sum index lies in a vector space of dimension at most
10 on Δ ≤ 4 graphs, so high empirical redundancy with the
30-index baseline is the expected outcome.
