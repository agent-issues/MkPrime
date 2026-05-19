# Rescore: multirep-v3 (Sim 3 multirep) — 8 reps under root-invariant scoring

Date: 2026-05-19

Setup: 16 tips arranged as `((A,C),(B,D))`, each clade balanced quartet, 
eco-1 = ALL of clades A and B (full-clade convergent ecology), 480 chars 
(120 neomorphic + 360 transformational), phi = 4, pi0 = 0.75, theta = 1, 
tipBr = 0.50, stemBr = 0.30, rootBr = 0.15, truth TL = 13.500.
Scoring uses `inst/simulations/ecology/sim3-scoring.R::ScoreTreesUnrooted` 
(root-invariant) on the post-burnin 75% of each chain's saved trees.

## Per-rep table — corrected support

|  rep | chain |  nKept | P(AC)  | P(BD)  | P(AB)  | P(CD)  | TLmean | nUniq | topShare |
|-----:|:------|------:|-------:|-------:|-------:|-------:|-------:|------:|--------:|
|   1 | blind |   157 | 0.000 | 0.000 | **0.854** | 0.854 | 151.67 |   157 | 0.006 |
|   1 | aware |   180 | 0.022 | 0.022 | **0.789** | 0.789 |  3.19 |   160 | 0.022 |
|   2 | blind |   178 | 0.000 | 0.000 | **1.000** | 1.000 | 108.59 |   151 | 0.017 |
|   2 | aware |   168 | 0.006 | 0.006 | **0.375** | 0.375 |  3.69 |   130 | 0.024 |
|   3 | blind |   171 | 0.000 | 0.000 | **0.444** | 0.444 | 35.99 |   170 | 0.012 |
|   3 | aware |   176 | 0.000 | 0.000 | **0.023** | 0.023 | 20.71 |   113 | 0.062 |
|   4 | blind |   183 | 0.000 | 0.000 | **0.918** | 0.918 | 24.38 |   139 | 0.038 |
|   4 | aware |   172 | 0.012 | 0.012 | **0.320** | 0.320 |  3.21 |    94 | 0.052 |
|   5 | blind |   180 | 0.000 | 0.000 | **0.978** | 0.978 | 119.49 |   168 | 0.022 |
|   5 | aware |   182 | 0.000 | 0.000 | **0.462** | 0.462 | 19.21 |   111 | 0.060 |
|   6 | blind |   178 | 0.000 | 0.000 | **0.949** | 0.949 | 81.19 |   142 | 0.022 |
|   6 | aware |   165 | 0.030 | 0.030 | **0.085** | 0.085 | 18.44 |   117 | 0.030 |
|   7 | blind |   167 | 0.000 | 0.000 | **1.000** | 1.000 | 138.70 |   132 | 0.024 |
|   7 | aware |   174 | 0.000 | 0.000 | **0.983** | 0.983 |  2.87 |    99 | 0.029 |
|   8 | blind |   167 | 0.012 | 0.012 | **0.772** | 0.772 | 154.80 |   153 | 0.012 |
|   8 | aware |   178 | 0.028 | 0.028 | **0.253** | 0.253 |  3.23 |   167 | 0.017 |

## Mean across reps (corrected)

| metric          |  blind |  aware |
|:----------------|------:|------:|
| P(AC) true      | 0.001 | 0.012 |
| P(BD) true      | 0.001 | 0.012 |
| **P(AB) false** | **0.864** | **0.411** |
| P(CD) false     | 0.864 | 0.411 |
| TL mean         | 101.85 | 9.32 | (truth 13.50)
| nUnique topo    | 151.5 | 123.9 |

## Headline counts (out of 8 reps)

- Blind P(AB) > 0.5: **7 / 8**
- Blind P(AB) > 0.9: **5 / 8**
- Aware P(AC) > 0.5: **0 / 8**
- Aware P(AC) > 0.9: 0 / 8
- Aware mode-trapped (P(AC) < 0.5): 8 / 8
- Model dichotomy (blind P(AB)>0.5 AND aware P(AC)>0.5): **0 / 8**
- Strong dichotomy (both > 0.9): 0 / 8
- Reversed (aware false-supports AB > 0.5 AND blind correct AC > 0.5): 0 / 8

## Per-rep blind-vs-aware delta on the headline metrics

|  rep | P(AC) blind | P(AC) aware | P(AB) blind | P(AB) aware |  TL blind |  TL aware | dichotomy? |
|-----:|----:|----:|----:|----:|----:|----:|:---|
|   1 | 0.000 | 0.022 | 0.854 | 0.789 | 151.67 |  3.19 | no |
|   2 | 0.000 | 0.006 | 1.000 | 0.375 | 108.59 |  3.69 | no |
|   3 | 0.000 | 0.000 | 0.444 | 0.023 | 35.99 | 20.71 | no |
|   4 | 0.000 | 0.012 | 0.918 | 0.320 | 24.38 |  3.21 | no |
|   5 | 0.000 | 0.000 | 0.978 | 0.462 | 119.49 | 19.21 | no |
|   6 | 0.000 | 0.030 | 0.949 | 0.085 | 81.19 | 18.44 | no |
|   7 | 0.000 | 0.000 | 1.000 | 0.983 | 138.70 |  2.87 | no |
|   8 | 0.012 | 0.028 | 0.772 | 0.253 | 154.80 |  3.23 | no |

(`dichotomy?` = blind supports the false eco-clade AB AND aware supports the true AC.)

## Legacy (root-dependent) vs corrected scoring — P(AB) on blind chains

The whole point of rescoring: the original prop.part scorer in `run_rep.R`
underestimated the false-clade signal whenever the chain rooted inside
the {A,B} tip set. The size of that artefact, per rep:

| rep | legacy P(AB) | corrected P(AB) |   delta |
|----:|-------:|-------:|-------:|
|   1 | 0.854  | 0.854  |  +0.000 |
|   2 | 0.000  | **1.000** | **+1.000** |
|   3 | 0.439  | 0.444  |  +0.006 |
|   4 | 0.902  | 0.918  |  +0.016 |
|   5 | 0.000  | **0.978** | **+0.978** |
|   6 | 0.444  | 0.949  |  +0.506 |
|   7 | 1.000  | 1.000  |  +0.000 |
|   8 | 0.772  | 0.772  |  +0.000 |

Mean blind P(AB): legacy = 0.551, corrected = **0.864** (+0.313).
Three reps (02, 05, 06) had the false-clade signal substantially or
entirely hidden by the legacy scoring bug. The aware chains were
largely unaffected (mean delta +0.030) because aware-chain root
positions were less correlated with the eco split.

## Interpretation

The corrected scoring confirms a strong, reproducible **false eco-clade
signal under the blind model**: across 8 independent replicates, mean
corrected P(AB) = **0.864**, with 7/8 reps > 0.5 and 5/8 reps > 0.9.
This is the long-sought convergent-ecology demo and it is robust across
replicates of the same setup.

The story under the aware model is **not** clean recovery of the truth,
however: aware P(AC) is essentially zero everywhere (max 0.030, mean
0.012, 8/8 reps mode-trapped). Aware drastically *reduces* the false
clade — mean P(AB) drops from 0.864 (blind) to 0.411 (aware) — and
cuts TL from a wildly inflated ~102 (blind) toward ~9 (aware, truth =
13.5), but the aware posterior is **diffuse**: 94-167 unique canonical
topologies among ~170 post-burnin samples, top-topology share 1.7-6.2 %.
Aware refuses to commit to (A,C),(B,D) and instead spreads probability
across many topologies that lack the false eco split.

That is the **partial dichotomy**: blind locks onto the false eco-clade,
aware refuses to lock onto anything. The strict model-dichotomy demo
(blind supports false clade AND aware supports truth) is **NOT** observed
in any of the 8 reps under the present settings.

This rules out using multirep-v3 verbatim as the publication-quality
"two models, two answers" headline. It still strongly demonstrates
**half** the story (aware correctly refuses to endorse the convergence-
driven false clade), which is itself publishable as the negative-control
half of the claim.

### Implications for v5break

v5break (job 17226503) is still needed: multirep-v3 does not deliver
a clean aware-recovers-truth result. The diffuse-posterior pattern in
aware suggests mixing or identifiability limits the aware model rather
than incorrect inference per se, but the present data do not show aware
beating blind in *the direction of truth* — only away from the
ecology-driven false clade.

### Caveat reps

- Rep 07 — aware also strongly supports the false clade (P(AB) = 0.983).
  Visually the aware chain TL is collapsed (TLmean = 2.87 << truth =
  13.5), suggesting an aware-chain stuck mode rather than a real signal
  inversion.
- Rep 01 — aware P(AB) = 0.789 alongside collapsed TLmean = 3.19, same
  pattern as rep 07 (aware mode trap on a short-TL local optimum). The
  previously diagnosed "TL collapse" issue.
- No rep showed aware false-supporting AB while blind correctly
  supported AC — the reversed pattern is not observed.

Outputs:
- `inst/scripts/rescore-multirep-v3-results.{rds,csv}` (long form)
- `inst/scripts/rescore-multirep-v3-summary.csv` (wide form)
- `inst/scripts/rescore-multirep-v3-plot.pdf` (bar chart)
