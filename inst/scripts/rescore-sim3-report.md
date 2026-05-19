# Sim 3 v4 / v4-cross — root-invariant rescoring report

## 1. The fix

`hasBipart()` in every `run_v4*.R` and `run_v4cross*.R` script used
`ape::prop.part`, which returns the *rooted* partitions of a tree.
When MCMC trees are rooted in different places across the sample, the
same unrooted bipartition appears in `prop.part`'s output as the
target tip-set on some trees and its complement on others — so
`setequal(p, target)` returns FALSE for every tree whose root sits
inside the target.

`inst/simulations/ecology/sim3-scoring.R` replaces `hasBipart` with
`HasBipartSplits`, which enumerates unrooted bipartitions via
`TreeTools::as.Splits` and accepts a split as present if some
Splits row matches either the target *or* its complement.

Unit test (`((A1,A2),(A3,(B1,B2)))` re-rooted on B1, asking about
`{B1,B2}`): legacy returns FALSE, corrected returns TRUE.

## 2. Bug bite varies enormously by run

Across the 22 rescored chains, the saved tree lists contain 3–20
distinct root configurations. Wherever a target bipartition's tip
set ever contained the root, the legacy scorer under-counted that
split. CID is root-invariant; corrected and legacy CIDs match to
machine epsilon for all 22 runs (max |Δ| ≈ 0).

## 3. v4cross-b-pt5 rep 01 — the original symptom

|                   | Blind, legacy | Blind, corr. | Aware, legacy | Aware, corr. |
|-------------------|--------------:|-------------:|--------------:|-------------:|
| trueClade_A       |         0.973 |        0.973 |         0.552 |        0.984 |
| trueClade_B       |         0.448 |        0.995 |         0.967 |        1.000 |
| trueClade_C       |         1.000 |        1.000 |         0.995 |        1.000 |
| trueClade_D       |         0.885 |        1.000 |         0.770 |        1.000 |
| trueSister_AC     |         0.978 |        1.000 |         0.437 |        1.000 |
| trueSister_BD     |         0.148 |        1.000 |         0.607 |        1.000 |
| falseInner        |             0 |            0 |             0 |            0 |
| falseClade_AB     |             0 |            0 |             0 |            0 |

The legacy P(AC)=0.978 vs P(BD)=0.148 asymmetry collapses: both go
to 1.000 once root-invariance is restored. The "aware clade-C
collapse" (legacy P(clC)=0.219 in v4cross-b-pt-aware) likewise
vanishes — corrected = 1.000.

## 4. v4cross-b-pt5 PT five-rep means (blind vs aware, post-burnin)

|                | Blind, corr. mean (sd) | Aware, corr. mean (sd) | Blind, legacy mean | Aware, legacy mean |
|----------------|----------------------:|----------------------:|-------------------:|-------------------:|
| trueClade_A    | 0.995 (0.012)         | 0.997 (0.007)         | 0.981              | 0.776              |
| trueClade_B    | 0.999 (0.002)         | 1.000 (0.000)         | 0.787              | 0.892              |
| trueClade_C    | 1.000 (0.000)         | 1.000 (0.000)         | 0.885              | 0.909              |
| trueClade_D    | 1.000 (0.000)         | 1.000 (0.000)         | 0.785              | 0.814              |
| trueSister_AC  | 0.999 (0.002)         | 0.996 (0.008)         | 0.750              | 0.479              |
| trueSister_BD  | 0.999 (0.002)         | 0.996 (0.008)         | 0.374              | 0.578              |
| falseInner     | 0                     | 0                     | 0                  | 0                  |
| falseClade_AB  | 0                     | 0.001                 | 0                  | 0.001              |

**The "aware loses on v4-cross" headline does NOT survive correction.**
Once the rooting artefact is removed, both chains recover the true
topology at probability ~1 across all five reps. The legacy
ranking that put blind ahead of aware on AC/BD was driven entirely
by chain-specific root-placement preferences interacting with the
buggy scorer, not by any genuine topological difference.

## 5. False-clade support

P(falseInner) = 0 across all 22 chains, corrected and legacy alike.
P(falseClade_AB) = 0 everywhere except v4cross-b-pt5-rep05-aware
(0.006, equal in both scorings). The bug did not hide active false-
clade support; the unrooted scoring confirms that the ecology-driven
spurious clades never received non-trivial posterior support in any
saved run.

## 6. v4 (non-cross) runs

| Run            | Metric         | Corrected | Legacy |
|----------------|----------------|----------:|-------:|
| v4-v4a-blind   | trueClade_A    |     0.906 |  0.906 |
|                | trueClade_C    |     0.813 |  0.813 |
|                | trueSister_AC  |     0.912 |  0.819 |
| v4-v4b-blind   | trueSister_AC  |     1.000 |  0.574 |
| v4-v4c-blind   | trueSister_AC  |     1.000 |  0.911 |
| v4-v4c-aware   | trueClade_C    |     1.000 |  0.005 |
|                | trueSister_AC  |     1.000 |  0.000 |

v4a is barely affected (low rooting variability). v4b–v4c are
materially affected: notably v4-v4c-aware legacy reported
P(trueSister_AC) = 0 (catastrophic), corrected = 1.000. Any v4
narrative quoted from legacy CSVs needs revisiting.

## 7. Files committed

- `inst/simulations/ecology/sim3-scoring.R` — root-invariant helpers.
- `inst/scripts/rescore-sim3.R` — re-scoring driver (Hamilton-ready).
- `inst/scripts/rescore-sim3-results.{rds,csv}` — tidy scores
  (corrected and legacy side by side) for all 22 chains.
- `inst/scripts/rescore-sim3-root-dist.csv` — per-chain root-position
  distributions (3–20 distinct rootings observed per chain).
- `inst/scripts/rescore-sim3-report.md` — this report.

The original `run_v4*.R` and `run_v4cross*.R` scripts are deliberately
unchanged (historical record); future v4-style scripts should
`source("inst/simulations/ecology/sim3-scoring.R")` and call
`ScoreTreesUnrooted` instead of inlining `hasBipart`.
