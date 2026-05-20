# Aware multirep-v3 posterior shape characterisation

Date: 2026-05-19

Goal: is the multirep-v3 aware chain doing the methodologically RIGHT thing 
(broad honest uncertainty, including truth in the credibility set) or the WRONG 
thing (failing to find truth at all)?

Setup: 16 tips, full-clade eco-1 (clades A and B as the 8-tip false eco clade), 
truth = ((A,C),(B,D)); tipBr=0.50, stemBr=0.30, rootBr=0.15, truth TL=13.500. 
Scoring is root-invariant via `inst/ecology/simulations/sim3-scoring.R`.

## Per-rep x chain posterior summary

|  rep | chain |    n | P(AC) | P(AB) | CID mean | CID 05 | CID 95 | nUniq | |CS95| | truth in CS? |
|-----:|:------|----:|----:|----:|------:|-----:|-----:|----:|----:|:----|
|   1 | aware |  180 | 0.022 | 0.789 | 0.446 | 0.332 | 0.545 |  160 |  151 | no |
|   1 | blind |  157 | 0.000 | 0.854 | 0.519 | 0.401 | 0.600 |  157 |  150 | no |
|   2 | aware |  168 | 0.006 | 0.375 | 0.460 | 0.317 | 0.583 |  130 |  122 | no |
|   2 | blind |  178 | 0.000 | 1.000 | 0.464 | 0.368 | 0.536 |  151 |  143 | no |
|   3 | aware |  176 | 0.000 | 0.023 | 0.502 | 0.373 | 0.657 |  113 |  105 | no |
|   3 | blind |  171 | 0.000 | 0.444 | 0.503 | 0.393 | 0.618 |  170 |  162 | no |
|   4 | aware |  172 | 0.012 | 0.320 | 0.404 | 0.273 | 0.506 |   94 |   86 | no |
|   4 | blind |  183 | 0.000 | 0.918 | 0.496 | 0.398 | 0.582 |  139 |  130 | no |
|   5 | aware |  182 | 0.000 | 0.462 | 0.473 | 0.393 | 0.605 |  111 |  102 | no |
|   5 | blind |  180 | 0.000 | 0.978 | 0.475 | 0.378 | 0.550 |  168 |  159 | no |
|   6 | aware |  165 | 0.030 | 0.085 | 0.350 | 0.235 | 0.471 |  117 |  109 | no |
|   6 | blind |  178 | 0.000 | 0.949 | 0.418 | 0.322 | 0.520 |  142 |  134 | no |
|   7 | aware |  174 | 0.000 | 0.983 | 0.399 | 0.280 | 0.553 |   99 |   91 | no |
|   7 | blind |  167 | 0.000 | 1.000 | 0.521 | 0.397 | 0.606 |  132 |  124 | no |
|   8 | aware |  178 | 0.028 | 0.253 | 0.417 | 0.322 | 0.511 |  167 |  159 | no |
|   8 | blind |  167 | 0.012 | 0.772 | 0.472 | 0.372 | 0.543 |  153 |  145 | no |

Notes:
- `P(AC)` = posterior probability that the TRUE sister-of-clades bipartition `((A,C) vs rest)` is present.
- `P(AB)` = posterior probability of the FALSE (eco-driven) bipartition `((A,B) vs rest)`.
- `CID` = normalised Clustering Information Distance to truth, lower is better; 0 = identical.
- `|CS95|` = number of unique canonical topologies in the 95% credibility set.

## Aggregate aware-vs-blind comparison

Mean across the 8 reps (corrected, root-invariant):

| metric                       |   blind |   aware |  delta (aware - blind) |
|:-----------------------------|--------:|--------:|----------------------:|
| mean P(AC) true              | 0.001 | 0.012 |  +0.011 |
| mean P(AB) false             | 0.864 | 0.411 |  -0.453 |
| mean CID-to-truth            | 0.484 | 0.432 |  -0.052 |
| mean nUnique topologies      | 151.5 | 123.9 |  -27.6 |
| mean |CS95|                  | 143.4 | 115.6 |  -27.8 |
| reps with truth in CS95      | 0 / 8 | 0 / 8 | |
| reps with truth in posterior | 0 / 8 | 0 / 8 | |

Paired comparison of mean CID-to-truth across reps:
- Reps where aware's mean CID is LOWER than blind's: **8 / 8**
- Reps where aware's mean CID is HIGHER than blind's: 0 / 8
- Paired Wilcoxon (aware - blind, two-sided): V = 0.0, p = 0.0143

Per-rep CID-to-truth comparison:

|  rep | blind mean | aware mean | delta | blind median | aware median | blind 95% | aware 95% |
|----:|----:|----:|----:|----:|----:|----:|----:|
|   1 | 0.519 | 0.446 | -0.073 | 0.524 | 0.452 | 0.600 | 0.545 |
|   2 | 0.464 | 0.460 | -0.004 | 0.473 | 0.467 | 0.536 | 0.583 |
|   3 | 0.503 | 0.502 | -0.001 | 0.505 | 0.507 | 0.618 | 0.657 |
|   4 | 0.496 | 0.404 | -0.092 | 0.500 | 0.436 | 0.582 | 0.506 |
|   5 | 0.475 | 0.473 | -0.002 | 0.480 | 0.468 | 0.550 | 0.605 |
|   6 | 0.418 | 0.350 | -0.069 | 0.418 | 0.350 | 0.520 | 0.471 |
|   7 | 0.521 | 0.399 | -0.122 | 0.530 | 0.399 | 0.606 | 0.553 |
|   8 | 0.472 | 0.417 | -0.054 | 0.476 | 0.425 | 0.543 | 0.511 |

## Top-5 canonical topologies per rep x chain

`hasAC`, `hasAB` are unrooted-bipartition presence flags for the top topology of that mass.

### rep 1 / blind — top 5 (of 157 unique; |CS95|=150, CS mass=0.955)

| rank | prob  | hasAC | hasAB |   CID |
|----:|-----:|:---|:---|----:|
|   1 | 0.006 | . | . | 0.467 |
|   2 | 0.006 | . | . | 0.548 |
|   3 | 0.006 | . | . | 0.562 |
|   4 | 0.006 | . | AB | 0.564 |
|   5 | 0.006 | . | AB | 0.533 |

### rep 1 / aware — top 5 (of 160 unique; |CS95|=151, CS mass=0.950)

| rank | prob  | hasAC | hasAB |   CID |
|----:|-----:|:---|:---|----:|
|   1 | 0.022 | . | AB | 0.507 |
|   2 | 0.017 | . | AB | 0.500 |
|   3 | 0.017 | . | AB | 0.547 |
|   4 | 0.011 | . | AB | 0.511 |
|   5 | 0.011 | . | . | 0.420 |

### rep 2 / blind — top 5 (of 151 unique; |CS95|=143, CS mass=0.955)

| rank | prob  | hasAC | hasAB |   CID |
|----:|-----:|:---|:---|----:|
|   1 | 0.017 | . | AB | 0.530 |
|   2 | 0.017 | . | AB | 0.513 |
|   3 | 0.017 | . | AB | 0.492 |
|   4 | 0.017 | . | AB | 0.402 |
|   5 | 0.017 | . | AB | 0.451 |

### rep 2 / aware — top 5 (of 130 unique; |CS95|=122, CS mass=0.952)

| rank | prob  | hasAC | hasAB |   CID |
|----:|-----:|:---|:---|----:|
|   1 | 0.024 | . | . | 0.507 |
|   2 | 0.024 | . | AB | 0.496 |
|   3 | 0.018 | . | . | 0.405 |
|   4 | 0.018 | . | . | 0.430 |
|   5 | 0.018 | . | AB | 0.343 |

### rep 3 / blind — top 5 (of 170 unique; |CS95|=162, CS mass=0.953)

| rank | prob  | hasAC | hasAB |   CID |
|----:|-----:|:---|:---|----:|
|   1 | 0.012 | . | . | 0.474 |
|   2 | 0.006 | . | . | 0.637 |
|   3 | 0.006 | . | . | 0.513 |
|   4 | 0.006 | . | . | 0.498 |
|   5 | 0.006 | . | . | 0.506 |

### rep 3 / aware — top 5 (of 113 unique; |CS95|=105, CS mass=0.955)

| rank | prob  | hasAC | hasAB |   CID |
|----:|-----:|:---|:---|----:|
|   1 | 0.062 | . | . | 0.439 |
|   2 | 0.040 | . | . | 0.444 |
|   3 | 0.034 | . | . | 0.465 |
|   4 | 0.034 | . | . | 0.533 |
|   5 | 0.028 | . | . | 0.392 |

### rep 4 / blind — top 5 (of 139 unique; |CS95|=130, CS mass=0.951)

| rank | prob  | hasAC | hasAB |   CID |
|----:|-----:|:---|:---|----:|
|   1 | 0.038 | . | AB | 0.493 |
|   2 | 0.022 | . | AB | 0.481 |
|   3 | 0.022 | . | AB | 0.500 |
|   4 | 0.022 | . | AB | 0.526 |
|   5 | 0.022 | . | AB | 0.528 |

### rep 4 / aware — top 5 (of 94 unique; |CS95|=86, CS mass=0.953)

| rank | prob  | hasAC | hasAB |   CID |
|----:|-----:|:---|:---|----:|
|   1 | 0.052 | . | . | 0.436 |
|   2 | 0.041 | . | . | 0.287 |
|   3 | 0.029 | . | . | 0.486 |
|   4 | 0.029 | . | . | 0.487 |
|   5 | 0.029 | . | . | 0.429 |

### rep 5 / blind — top 5 (of 168 unique; |CS95|=159, CS mass=0.950)

| rank | prob  | hasAC | hasAB |   CID |
|----:|-----:|:---|:---|----:|
|   1 | 0.022 | . | AB | 0.372 |
|   2 | 0.011 | . | AB | 0.506 |
|   3 | 0.011 | . | AB | 0.498 |
|   4 | 0.011 | . | AB | 0.392 |
|   5 | 0.011 | . | AB | 0.550 |

### rep 5 / aware — top 5 (of 111 unique; |CS95|=102, CS mass=0.951)

| rank | prob  | hasAC | hasAB |   CID |
|----:|-----:|:---|:---|----:|
|   1 | 0.060 | . | AB | 0.468 |
|   2 | 0.060 | . | . | 0.489 |
|   3 | 0.027 | . | . | 0.477 |
|   4 | 0.022 | . | AB | 0.468 |
|   5 | 0.022 | . | . | 0.428 |

### rep 6 / blind — top 5 (of 142 unique; |CS95|=134, CS mass=0.955)

| rank | prob  | hasAC | hasAB |   CID |
|----:|-----:|:---|:---|----:|
|   1 | 0.022 | . | AB | 0.403 |
|   2 | 0.022 | . | AB | 0.372 |
|   3 | 0.022 | . | AB | 0.461 |
|   4 | 0.022 | . | AB | 0.389 |
|   5 | 0.017 | . | AB | 0.452 |

### rep 6 / aware — top 5 (of 117 unique; |CS95|=109, CS mass=0.952)

| rank | prob  | hasAC | hasAB |   CID |
|----:|-----:|:---|:---|----:|
|   1 | 0.030 | . | . | 0.430 |
|   2 | 0.024 | . | . | 0.332 |
|   3 | 0.024 | AC | . | 0.285 |
|   4 | 0.018 | . | . | 0.471 |
|   5 | 0.018 | . | . | 0.443 |

### rep 7 / blind — top 5 (of 132 unique; |CS95|=124, CS mass=0.952)

| rank | prob  | hasAC | hasAB |   CID |
|----:|-----:|:---|:---|----:|
|   1 | 0.024 | . | AB | 0.408 |
|   2 | 0.018 | . | AB | 0.544 |
|   3 | 0.018 | . | AB | 0.587 |
|   4 | 0.018 | . | AB | 0.487 |
|   5 | 0.018 | . | AB | 0.332 |

### rep 7 / aware — top 5 (of 99 unique; |CS95|=91, CS mass=0.954)

| rank | prob  | hasAC | hasAB |   CID |
|----:|-----:|:---|:---|----:|
|   1 | 0.029 | . | AB | 0.434 |
|   2 | 0.029 | . | AB | 0.326 |
|   3 | 0.023 | . | AB | 0.458 |
|   4 | 0.023 | . | AB | 0.553 |
|   5 | 0.023 | . | AB | 0.553 |

### rep 8 / blind — top 5 (of 153 unique; |CS95|=145, CS mass=0.952)

| rank | prob  | hasAC | hasAB |   CID |
|----:|-----:|:---|:---|----:|
|   1 | 0.012 | . | AB | 0.514 |
|   2 | 0.012 | . | AB | 0.506 |
|   3 | 0.012 | . | AB | 0.403 |
|   4 | 0.012 | . | . | 0.477 |
|   5 | 0.012 | . | . | 0.524 |

### rep 8 / aware — top 5 (of 167 unique; |CS95|=159, CS mass=0.955)

| rank | prob  | hasAC | hasAB |   CID |
|----:|-----:|:---|:---|----:|
|   1 | 0.017 | . | . | 0.439 |
|   2 | 0.011 | . | . | 0.405 |
|   3 | 0.011 | . | . | 0.437 |
|   4 | 0.011 | . | . | 0.366 |
|   5 | 0.011 | . | AB | 0.394 |

## Truth-in-posterior summary

Each chain's posterior is a finite sample; a canonical topology hash either appears in the post-burnin draws or does not.

| rep | blind: truth in post? | blind: truth prob | aware: truth in post? | aware: truth prob |
|----:|:---|----:|:---|----:|
|   1 | no | 0.000 | no | 0.000 |
|   2 | no | 0.000 | no | 0.000 |
|   3 | no | 0.000 | no | 0.000 |
|   4 | no | 0.000 | no | 0.000 |
|   5 | no | 0.000 | no | 0.000 |
|   6 | no | 0.000 | no | 0.000 |
|   7 | no | 0.000 | no | 0.000 |
|   8 | no | 0.000 | no | 0.000 |

## Verdict criteria

Interpretation grid:
- **Broad honest uncertainty**: aware's CID-to-truth distribution is LOWER (closer to truth) than blind's on most reps, AND the credibility set is broader (larger |CS95|), AND truth is occasionally inside CS95.
- **Failing to find truth**: aware's posterior never includes AC and CID-to-truth is no better than blind's.
- **Mode-trapping into wrong neighbourhood**: aware concentrates probability on a NON-AB, NON-AC topology that is itself far from truth, so neither bipartition flag fires but CID is poor.

See report Section 'Aggregate aware-vs-blind comparison' above for the numeric verdict inputs.

Outputs:
- `inst/ecology/scripts/aware-multirep-v3-posterior-shape.{rds,csv}`
- `inst/ecology/scripts/aware-multirep-v3-posterior-shape-credset.csv`
- `inst/ecology/scripts/aware-multirep-v3-posterior-shape-plot.pdf`
