# Profiling focus areas — mkp

Hot paths inside the C++ MCMC engine. Ranked by (per-call cost) × (call frequency on the inner loop).

**Baseline from `data-raw/step6-profile-result.rds` (2026-05-18):** 1.57 min wall, `.Call → run_mcmc_batch_cpp` = **99.27 % self.pct**. The Rcpp boundary contributes essentially zero overhead, so all findings will be C++-internal (`[Optimise]` or `[AT-LIMIT]`); no `[Port]` candidates remain. The 0.5 % "bytecode compiler" tail (`cmpfun`, `genCode`, `tryInline`) is fixed startup cost and not addressable.

Methodology: profvis surfaces only the top-level R wrapper because by-self time inside `.Call` is unobservable from R. The rotation therefore drives a representative workload through each focus area and uses VTune (or `perf` on Hamilton) for inner-loop attribution.

| # | Area | Files | Why hot (1 line) | Baseline cost | Last profiled | Status |
|---|------|-------|------------------|---------------|---------------|--------|
| 1 | Felsenstein pruning (CL accumulation) | `src/mcmc_likelihood.cpp`, `src/node_cl_cache.h` | Called per evaluated tree per partition; dominant cost in a typical MCMC iteration | 99.8 s / 800 iter (Sun2018, EG, prod move weights) | 2026-05-18 | PROFILED |
| 2 | Partial CL invalidation under SPR | `src/node_cl_cache.h` (`find_dirty_spr`, partial-eval walks), `src/tree_moves.cpp` | M-158 (be2f67b) introduced partial-CL evaluation for SPR; the dirty-set walk is the new hot loop for any accepted tree move | — | — | NEW |
| 3 | Per-character k′ Gibbs sweep | `src/mcmc_likelihood.cpp` (Gibbs kPrime), `src/mcmc.cpp` (cases dispatching to it) | M-154+M-155+M-172 (e8819a1, 224bda8, 2b9b4c9) — block + within-partition pattern compression; called every iteration for every transformational character | — | — | NEW |
| 4 | Rate-matrix exponentiation P(t) = exp(Qt) | `src/rate_matrix.cpp` | Analytical JC; called per (k′, branch-length) combo; if cache miss this is per-edge per-iteration | — | — | NEW |
| 5 | Ascertainment correction | `src/ascertainment.cpp` | Constant-site / singleton corrections; called per evaluated likelihood; M-145+M-143 cache invalidation may force re-eval | — | — | NEW |
| 6 | Tree-move proposal cost (SPR/NNI/TBR) | `src/tree_moves.cpp`, `src/proposals.cpp` | Move proposal + Hastings ratio; per-iteration; M-160 reverted DR-NNI but TBR/pSPR still in dispatch | — | — | NEW |
| 7 | ACRV (across-character rate variation) | `src/acrv.cpp` | Lognormal ACRV; sweeping the rate categories per character per evaluation when `qHet=TRUE` | — | — | NEW |
| 8 | Relabelling correction for Mk′ | `src/corrections.cpp` | Called per likelihood eval for transformational characters; tiny file but factorial-style calc | — | — | NEW |
| 9 | MH dispatch + adaptive tuning (R-side hot) | `R/RunMkPrime.R` (.RunMkPrimeSingleRun batch loop, .StateToRow), inside-batch R callbacks | 0.7 % residual after `.Call`; not large but the only R on the hot path — confirm no allocation surprises | — | — | NEW |
| 10 | sample-row write (StateToRow + stream buffer) | `R/streaming.R`, `R/RunMkPrime.R` (.StateToRow), `src/mcmc.cpp:4845-4881` | Each saved sample flushes through R; STREAM-002 + STREAM-003 already touched this; check allocation churn | — | — | NEW |

**Parked / SKIPPED unless explicit:** setup phases (`R/MkPrimeData.R`, `R/MkPrimeModel.R`), CLI parsing, file IO, the prior-density paths (`R/MkPrimeModel.R::LogPrior*`) — they run once per session or once per chain init.

## Notes on the ranking

- **Area 1 (Felsenstein pruning)** is mechanically the dominant cost (cost ∝ nTip × nChar × nState per partition per eval). Even after the M-161 per-partition CL cache (ee24fbb), every accepted tree move dirties some subtree and triggers re-eval. Confirm with a one-shot VTune run before ranking sub-areas inside it.
- **Area 2 (partial CL invalidation)** is fresh code surface (M-158 + M-161 + M-162B + M-164) — likely opportunities for cache-locality and dirty-set-pruning wins.
- **Area 3 (Gibbs kPrime)** has had three rounds of optimisation already (M-154/M-155/M-172) and may be `AT-LIMIT`; first VTune round should confirm.
- **Area 4 (rate matrix)** is analytical (closed-form for JC); probably negligible exclusive time but worth confirming.
- **Areas 9 and 10 (R-side residuals)** are tagged here for completeness — given the 99.27 % `.Call` dominance, expect both to be filed `AT-LIMIT` after first inspection.

## Rotation rule

Pick the next area where `status ∈ {NEW, PROFILED}` AND (`last_profiled` is empty OR `last_profiled < last_code_change` for that area's files). After scaffolding, area 1 (Felsenstein pruning) is up first.
