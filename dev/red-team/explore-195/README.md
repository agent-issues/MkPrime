# #195 exploration: single-topology posteriors and the assessment phase

Evidence for agent-issues/MkPrime#195; not package code. Scripts load the
renamed agent build (`MKP_LIB`, `MKP_PKG` in `common.R`). Each run is short
and guarded by `maxTime`, `setTimeLimit()` and a shell `timeout`.

| Script | Question | Output |
|---|---|---|
| `q1-topology-counts.R` | Do real datasets give one-topology posteriors? | `q1-results.csv` |
| `q2-assessment-phase.R` | Martin's assessment phase on congruent / stuck / mixing chains | `q2-results-Longrich2010.csv` |
| `q2b-sim8.R` | Same, on a concentrated 8-tip posterior, two budgets | `q2b-results.csv` |
| `q3-reachability.R` | Today's stop reason for the RT9-D configurations | `q3-results.csv` |
| `q4-bandit.R` | Does `.MinEssRate()` still reward frozen topology? | `q4-rates.csv`, `q4-h2h.csv` |

The parallel row of `q3-results.csv` is void: callr workers cannot load the
renamed build (`could not find function ".RunMkPrimeSingleRun"`).
