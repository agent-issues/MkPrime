# Plan: Benchmark and Optimize Tree ESS Calculation

**Agent:** F  
**Worktree:** `mkp-tree-ess` (feature/tree-ess)  
**Date:** 2026-03-28

---

## Context

Tree-topology ESS is computed via `ConvergenceDiagnostics(posterior, trees = TRUE)`.
The pipeline has three stages:

```
.ComputeTreeEss()
  └→ .TreeESS(chain)                     [R/treeESS.R]
       ├→ TreeDist::RobinsonFoulds(trees) → as.matrix()   [Stage 1: O(n² · t) pairwise RF]
       ├→ .FrechetCorrelationESS(dmat)    → C++ inner loop [Stage 2: O(n²/2 + L)]
       └→ .MedianPseudoESS(dmat)          → coda per row   [Stage 3: O(n · ESS)]
```

Where **n** = number of trees (capped at 1000) and **t** = number of tips.

The C++ Fréchet ESS inner loop (Stage 2) has already been optimized from
O(n²+nL) to O(n²/2+L), with matrix symmetry exploitation and precomputed
partial sums. The existing `benchmark/ab_frechet_ess.R` script covers this
component in isolation.

**What's missing:** No profiling of the full end-to-end pipeline. We don't
know how time distributes across the three stages, so we don't know where
to focus optimization effort.

---

## Approach: Measure → Identify → Optimize → Verify

Following the profiling skill's guidance: never optimize without a profile,
and always A/B benchmark changes.

---

## Phase 1: End-to-End Profiling

**Goal:** Determine what fraction of `.TreeESS()` wall time is spent in each
stage. This tells us where optimization effort has the highest leverage.

### Step 1.1: Build a profiling harness

Create `benchmark/profile-tree-ess.R` that:

1. Generates realistic tree samples (autocorrelated NNI chains from random
   starting trees, varying **n** from 100 to 1000 and **t** from 20 to 100).
2. Times each stage independently:
   - `TreeDist::RobinsonFoulds(trees)` (distance computation)
   - `as.matrix()` on the dist object (conversion)
   - `dmat * dmat` (element-wise squaring — creates a copy)
   - `frechet_correlation_ess_cpp()` (C++ inner loop)
   - `apply(dmat, 1, coda::effectiveSize)` (median pseudo-ESS)
3. Reports wall-time breakdown as a table: `n`, `t`, `stage`, `median_ms`.

No renamed package copies needed — this is a pure profiling pass using the
installed release build, not an A/B comparison.

### Step 1.2: Run the profiling harness

Execute from a clean `Rscript` subprocess (not RStudio — avoids load_all
pitfall). Build with `-O2` (release).

**Expected outcome:** Hypothesis is that Stage 1 (RF distances) dominates
overwhelmingly at large n, with Stage 3 (median pseudo-ESS) a distant second.
Stage 2 (C++ inner loop) is likely negligible after the O(n²/2+L) optimization.

---

## Phase 2: Targeted Optimization (contingent on Phase 1 results)

### If Stage 1 dominates (RF distance computation):

RF distance computation is external (`TreeDist::RobinsonFoulds`). Options:

**(a) Check for faster TreeDist API.**
`TreeDist::RobinsonFoulds` with `reportMatching = FALSE` (default) should
already be the fast path. Verify we're not passing unnecessary arguments.

**(b) Integer distance matrix — avoid double conversion.**
RF distances are always even integers. If we can get them as integers and
square in-place, we save the `as.matrix()` copy and the `dmat * dmat` copy.
Check whether `TreeDist::RobinsonFoulds` returns a `dist` of integers.

**(c) Fuse squaring into the C++ call.**
Currently `.FrechetCorrelationESS()` does `dmat_sq <- dmat * dmat` (R-level
copy of the full n×n matrix) then passes to C++. Instead, accept the
unsquared matrix in C++ and square on-the-fly during the upper-triangle
traversal — eliminates one O(n²) allocation+copy.

**(d) Parallelise RF distance computation.**
`TreeDist::RobinsonFoulds` is already vectorised internally; check whether
it supports parallel execution or whether we should compute the upper
triangle in parallel ourselves.

**(e) Reduce n more aggressively (statistical argument).**
The current cap is 1000 trees. For tree ESS, subsampling to 500 or even 250
may give sufficiently accurate estimates with 4× less work. We could
benchmark the accuracy–speed tradeoff by computing ESS at various subsample
sizes and checking convergence of the estimate.

### If Stage 3 is significant (median pseudo-ESS):

**(f) Batch ESS via spectral method.**
`coda::effectiveSize()` uses AR fitting per series. For n rows of length n,
this is n independent fits. A spectral (FFT-based) batch approach could
compute all n ESS values in one pass. Alternatively, move the per-row ESS
to C++ using the same autocorrelation-sum approach as Fréchet ESS.

**(g) Use only a row subsample.**
For medianPseudoESS, we only need the *median* of n ESS values. We could
subsample rows (e.g. every 4th row) and still get a good median estimate.

### If Stage 2 matters (unlikely, but for completeness):

**(h) SIMD vectorisation hints.**
The column-contiguous traversal in `tree_ess.cpp` is already friendly to
auto-vectorisation. If profiling shows it matters, add `#pragma omp simd`
or manual intrinsics for the inner accumulation loops.

---

## Phase 3: A/B Benchmarking of Changes

For each optimization applied in Phase 2:

1. Install the current code as the **reference** build (renamed `mkpRef`).
2. Apply the change.
3. Install as the **dev** build (renamed `mkpDev`).
4. Run `benchmark/ab_frechet_ess.R` (for Stage 2 changes) or a new
   end-to-end A/B script (for Stage 1/3 changes) in a single subprocess.
5. Include a canary benchmark on unrelated code.
6. Require >5% improvement with canary <3% to accept the change.

Use the existing `rename_and_install()` recipe from `PROFILING.md`.

---

## Phase 4: Correctness Verification

After any code change:

1. Run `tests/testthat/test-treeESS.R` — reference values validated against
   treess v1.0.1 must still match within tolerance.
2. For changes to `.FrechetCorrelationESS()` or the C++ code: verify exact
   numerical agreement between ref and dev builds on all test matrix sizes
   (the A/B script already checks this with `|r-v|/max(|r|,1) > 1e-6`).
3. For subsampling changes: document the accuracy tradeoff (ESS estimate
   variance vs. sample size).

---

## Deliverables

| Artifact | Description |
|----------|-------------|
| `benchmark/profile-tree-ess.R` | End-to-end profiling harness |
| Profile results table | Stage-by-stage timing at various (n, t) |
| Code changes (if any) | Targeted optimizations informed by profiling |
| A/B benchmark results | Before/after comparison with canary |
| Updated `benchmark/PROFILING.md` | Document tree-ESS profiling findings |

---

## Non-Goals

- **VTune profiling of tree_ess.cpp specifically.** The C++ inner loop is
  likely not the bottleneck. VTune is more useful for the MCMC hot path
  (already documented in PROFILING.md). If Phase 1 reveals the C++ code
  matters, we'll add a VTune pass.
- **Changing the TreeDist dependency.** We won't replace Robinson-Foulds
  with a different tree distance metric — that's a scientific decision,
  not a performance one.
- **Optimizing `.ComputeTreeEss()` orchestration.** The per-run loop and
  subsampling logic is trivial overhead compared to the distance computation.

---

## Execution Order

1. **Phase 1** (Step 1.1 + 1.2): Write and run the profiling harness → ~30 min
2. **Interpret results** and select 1–2 optimizations from Phase 2 → present to user
3. **Phase 3**: Implement, A/B benchmark, iterate
4. **Phase 4**: Verify correctness, commit
