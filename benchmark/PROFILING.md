# MkPrime — Profiling & Performance Opportunities

## Status

Phase 8 hot-path work (M-061–M-065) is complete:
- C++ MCMC inner loop (`run_mcmc_batch_cpp`)
- C++ proposals (NNI, SPR, BetaSimplex)
- Pre-allocated `ClWorkspace` flat buffer (M-063)
- Per-partition partial likelihood cache (M-064)
- Parent/child vector arguments (no edge-matrix round-trips) (M-065)

This document records the remaining opportunities identified by static
analysis (2026-03-27), ordered by expected impact.

---

## Build for profiling (Windows)

### 1. Add debug symbols to `src/Makevars.win`

```makefile
PKG_CXXFLAGS = $(SHLIB_CXXFLAGS) -g -fno-omit-frame-pointer
```

**Remove this file after profiling** — `-g` inflates the DLL and `-fno-omit-frame-pointer`
reserves a register at the cost of a marginal performance penalty.

### 2. Install with symbols (strips `-s` from DLLFLAGS)

```r
vtune_lib <- ".vtune-lib"
old_mf <- Sys.getenv("MAKEFLAGS")
stale <- Sys.glob(file.path("src", c("*.o", "*.dll", "*.so")))
if (length(stale)) file.remove(stale)
Sys.setenv(MAKEFLAGS = "DLLFLAGS=-static-libgcc")
install.packages(".", lib = vtune_lib, repos = NULL, type = "source",
                 INSTALL_opts = "--no-multiarch")
Sys.setenv(MAKEFLAGS = old_mf)
```

### 3. Collect hotspots

```bash
# Hardware sampling (requires Intel PMU driver)
vtune -collect hotspots -result-dir vtune-out -- Rscript benchmark/vtune-driver.R

# Software sampling (always works)
vtune -collect hotspots -knob sampling-mode=sw -result-dir vtune-out -- Rscript benchmark/vtune-driver.R

# Filter to package
vtune -report hotspots -result-dir vtune-out -filter "module=mkp.dll"
```

---

## A/B renaming recipe

To load both reference and dev builds in one process (avoids DLL-lock and
cross-session thermal drift):

```r
rename_and_install <- function(new_name, source_dir, lib_dir) {
  pkg <- file.path(tempdir(), paste0(new_name, "_src"))
  if (dir.exists(pkg)) unlink(pkg, recursive = TRUE)
  dir.create(pkg, recursive = TRUE)
  for (d in c("R", "src", "inst", "man", "vignettes")) {
    src_d <- file.path(source_dir, d)
    if (dir.exists(src_d)) file.copy(src_d, pkg, recursive = TRUE)
  }
  for (f in c("DESCRIPTION", "NAMESPACE", "LICENSE", "LICENSE.md")) {
    src_f <- file.path(source_dir, f)
    if (file.exists(src_f)) file.copy(src_f, pkg)
  }
  rewrite <- function(path, pattern, replacement) {
    lines <- readLines(path, warn = FALSE)
    writeLines(gsub(pattern, replacement, lines), path)
  }
  rewrite(file.path(pkg, "DESCRIPTION"), "^Package: mkp$",
          paste0("Package: ", new_name))
  rewrite(file.path(pkg, "NAMESPACE"), "\\bmkp\\b", new_name)
  rewrite(file.path(pkg, "R", "RcppExports.R"),
          "_mkp_", paste0("_", new_name, "_"))
  rewrite(file.path(pkg, "src", "RcppExports.cpp"),
          "_mkp_", paste0("_", new_name, "_"))
  rewrite(file.path(pkg, "src", "RcppExports.cpp"),
          "R_init_mkp", paste0("R_init_", new_name))
  stale <- Sys.glob(file.path(pkg, "src", c("*.o", "*.dll", "*.so")))
  if (length(stale)) file.remove(stale)
  dir.create(lib_dir, recursive = TRUE, showWarnings = FALSE)
  install.packages(pkg, lib = lib_dir, repos = NULL, type = "source",
                   INSTALL_opts = "--no-multiarch")
}

# Install reference (current HEAD on main, or tagged commit):
rename_and_install("mkpRef", ".", ".bench-lib/ref")
# Make your change, then:
rename_and_install("mkpDev", ".", ".bench-lib/dev")
# Then run: Rscript benchmark/bench-ab.R
```

---

## Static-analysis opportunities

### Tier 1 — Algorithmic (high expected impact)

#### OPP-1: JC matrix-vector product O(k²) → O(k)

**Files:** `mcmc_likelihood.cpp` (`pruning_jc_flat`, `pruning_jc_acrv_flat`),
           `likelihood.cpp` (`pruning_jc`), `acrv.cpp` (`pruning_jc_acrv`)

**Current code (inner loop per edge):**
```cpp
for (int c = 0; c < nChar; ++c) {
  int offset = c * kStates;
  for (int i = 0; i < kStates; ++i) {
    double sum = 0.0;
    for (int j = 0; j < kStates; ++j)
      sum += ((i == j) ? p_same : p_diff) * clCh[offset + j];
    clPar[offset + i] = sum;
  }
}
```

**Problem:** O(k²) FMAs per character per edge, plus a branch (`i==j`) that
defeats auto-vectorisation. For the JC model, the matrix has only two distinct
values (p_same on diagonal, p_diff everywhere else), so the product simplifies:

```
new_cl[i] = p_same * cl[i] + p_diff * (sum_cl − cl[i])
           = p_diff * sum_cl + (p_same − p_diff) * cl[i]
```

**Replacement:**
```cpp
double diff_coeff = p_same - p_diff;  // precomputed outside c-loop
for (int c = 0; c < nChar; ++c) {
  int offset = c * kStates;
  double sum_cl = 0.0;
  for (int j = 0; j < kStates; ++j) sum_cl += clCh[offset + j];
  for (int i = 0; i < kStates; ++i)
    clPar[offset + i] = p_diff * sum_cl + diff_coeff * clCh[offset + i];
}
```

**Expected gain:** For k=2 neutral on FMA count but removes the branch and
enables auto-vectorisation. For k=4: ~4× fewer FMAs in the inner product.
For k=7: ~7× fewer. This affects every edge traversal for transformational
characters; impact scales with mean kPrime.

The multiply+update path (`*=`) is analogous — precompute `sum` from `clCh`,
then write the product into `clPar[offset+i] *= ...` after the fact.

#### OPP-2: Cache `maxNode` in McmcState

**File:** `mcmc_likelihood.cpp` (`cpp_partition_log_likelihood`)

**Current:** Every partition call scans all `nEdge` pairs to find `maxNode`.
Called once per partition per MCMC iteration (and once more in the flat
workspace path in `pruning_jc_flat`).

**Fix:** Store `maxNode` in `McmcState` (or `McmcData`). Update it once at
`fill_partition_cache` time and again if topology changes are accepted.

**Expected gain:** Eliminates O(nPartitions × nEdge) scan per iteration at
zero cost. For hyoliths (63 tips, ~122 edges, ~3 partitions): ~366
comparisons saved per iteration.

#### OPP-3: Precompute `qnorm` midpoints in McmcData

**File:** `mcmc_likelihood.cpp` (`cpp_acrv_rates`)

**Current:**
```cpp
double z = R::qnorm(mid, 0.0, 1.0, 1, 0);  // called nCat times, every partition
```

The midpoints `mid = (i + 0.5) / nCat` depend only on `nCat`, not on `rateLogSd`.
The `z` values are therefore constant for the lifetime of the MCMC run.

**Fix:** Add `std::vector<double> acrvZ` to `McmcData`; populate in
`prepare_mcmc_data()`. Then `cpp_acrv_rates` just computes
`exp(mu + rateLogSd * acrvZ[i])` — no transcendental calls.

**Expected gain:** Removes `nCat` `qnorm` calls per partition per iteration.
For nCat=4 and 3 partitions: 12 `qnorm` calls eliminated every iteration.

---

### Tier 2 — Memory / allocation

#### OPP-4: `ascertainment.cpp` uses jagged CL heap allocation

**Files:** `ascertainment.cpp` (`constant_site_prob_jc`, `constant_site_prob_mkn`,
           `singleton_site_prob_jc`, `singleton_site_prob_mkn`)

**Current:** All four functions allocate `std::vector<std::vector<double>> CL`
on every call. These are invoked on every `coding != 0` partition evaluation
(always true for variable coding). The flat-buffer optimisation (M-063)
was not applied here.

**Fix:** Thread `ClWorkspace*` through `cpp_partition_log_likelihood` to the
ascertainment functions, or use a small stack-allocated array (the ascertainment
pruning uses only `nChar = kStates` pseudo-characters, so the CL buffer is tiny
and suitable for `alloca`/VLAs or a compile-time bounded stack array).

**Expected gain:** Removes 4 heap alloc/free round-trips per partition per
iteration when coding correction is active.

#### OPP-5: `beta_simplex_proposal` returns `Rcpp::List`

**File:** `proposals.cpp`

**Current:** Returns `List::create(_["value"] = ..., _["logHastings"] = ...)`.
The caller in `do_move_impl` immediately pulls both fields out. Constructing
and destroying an Rcpp `List` (with string key hashing) inside the hot loop
is avoidable overhead.

**Fix:** Change signature to take output references:
```cpp
void beta_simplex_proposal(NumericVector& x, int index, double tuning,
                           double& logHastings);
```
or return a simple `std::pair<NumericVector, double>`.

**Expected gain:** Eliminates one Rcpp `List` alloc/dealloc per BetaSimplex
proposal. BetaSimplex is the most frequent branch-length move.

#### OPP-6: NNI/SPR rollback clones entire vectors

**File:** `mcmc.cpp` (`do_move_impl`, cases 5 and 6)

**Current:**
```cpp
oldParent = clone(state->parent);
oldChild  = clone(state->child);
oldRelBr  = clone(state->relBrLengths);
```
Clones all `nEdge` parent, child, and relBrLength values on every topology
proposal — even when the move will be rejected.

**For NNI**, only 4 edges change (the 4 edges incident to the chosen internal
edge). The rollback could save just those indices+values.

**Expected gain:** For a 63-tip tree (122 edges), each clone copies ~122
integers/doubles. Replacing with targeted saves of 4–8 values cuts the
rollback cost by ~15–30×. NNI/SPR together make up a substantial fraction
of all proposals.

#### OPP-7: Transformational partition sub-matrix extraction

**File:** `mcmc_likelihood.cpp` (`cpp_partition_log_likelihood`, trans branch)

**Current:** On every kPrime-affecting call, allocates `IntegerMatrix sub(nTip, nSub)`
for each distinct kPrime group in the partition. For nTip=63 and nSub=O(nChar),
this is a non-trivial allocation inside the hot path.

**Fix:** Pre-allocate a scratch `IntegerMatrix` of maximum possible size in
`ClWorkspace` or `McmcState` and reuse it. Alternatively, pass column
indices and use pointer arithmetic directly against `part.tipStates`'s raw
data without constructing a new matrix.

---

### Tier 3 — Profiling-dependent

#### OPP-8: `edgeSamples` List growth

**File:** `mcmc.cpp` (`run_mcmc_batch_cpp`)

Each thinned sample pushes an `IntegerMatrix(nEdge, 2)` into an Rcpp `List`.
For large trees and small `thin`, this becomes a significant allocation
sequence. The current approach is fine for typical use (thin=10+), but for
dense sampling it may show up in profiles.

#### OPP-9: Numerical underflow (correctness + future robustness)

No log-rescaling is applied during the pruning traversal. For >100-taxon
trees with many characters and small branch lengths, products of conditional
likelihoods can approach floating-point zero before the root sum. Currently
safe for the intended use case, but worth tracking as tree sizes grow.

---

## Benchmark interpretation guide

| Diff vs reference | Interpretation |
|---|---|
| < 3% | Noise — cannot distinguish from thermal/scheduler variance |
| 3–5% | Possibly real; run again to confirm direction |
| > 5% | Likely real; confirm with a second independent run |
| Canary moves > 3% | Discard entire run — environment was unstable |

## VTune results (2026-03-28)

**Workload:** Sun2018 hyoliths (54 taxa, 225 chars, all transformational),
1 000 000 iterations, default move schedule (Gibbs SPR + Gibbs subtree
swap + TBR + NNI + SPR), no parallel tempering, no warmup.

**Collection:** User-mode sampling (software); hardware PMU driver not loaded.
Intel i7-10700, `-O2 -g -fno-omit-frame-pointer`, no stripping.

### Module breakdown (7.1 s CPU time)

| Module | CPU time | % of total | Interpretation |
|--------|----------|------------|----------------|
| R.dll | 4.58 s | 64.5% | R interpreter overhead (variable lookup, environment ops, pairlist construction) |
| ucrtbase.dll | 1.69 s | 23.8% | `malloc_base` (8.0%), `free_base` (5.4%), `strcmp`/`strcoll` (2.0%), other |
| MkPrime.dll | 0.68 s | 9.5% | C++ MCMC engine |
| Rcpp.dll | 0.09 s | 1.3% | `Rcpp_precious_preserve` etc. |

**R-side allocation overhead (R.dll + ucrtbase.dll):** `Rf_allocVector3`
(7.8%) + `malloc_base` (8.0%) + `free_base` (5.4%) + `Rf_cons` (5.2%) =
**26.4%** of total CPU on memory allocation and R object construction.

### Top functions overall (> 1% CPU)

| Function | Module | CPU time | % total |
|----------|--------|----------|---------|
| `func@0x2bb4ea9b0` (R internal) | R.dll | 0.572 s | 8.1% |
| `malloc_base` | ucrtbase.dll | 0.570 s | 8.0% |
| `Rf_allocVector3` | R.dll | 0.552 s | 7.8% |
| `free_base` | ucrtbase.dll | 0.383 s | 5.4% |
| `Rf_cons` | R.dll | 0.370 s | 5.2% |
| `func@0x2bb531da0` (R internal) | R.dll | 0.307 s | 4.3% |
| `TreeTools::postorder_order` | MkPrime.dll | 0.230 s | 3.2% |
| `Rf_findVarInFrame3` | R.dll | 0.198 s | 2.8% |
| Unresolved (static) | MkPrime.dll | 0.157 s | 2.2% |
| `func@0x2bb536860` (R internal) | R.dll | 0.154 s | 2.2% |
| `SETCAR` | R.dll | 0.152 s | 2.1% |
| `SETCDR` | R.dll | 0.139 s | 2.0% |
| `nni_proposal_impl` | MkPrime.dll | 0.107 s | 1.5% |
| `Rf_NewEnvironment` | R.dll | 0.106 s | 1.5% |
| `spr_proposal_impl` | MkPrime.dll | 0.093 s | 1.3% |

### Within MkPrime.dll (0.68 s, 100% = MkPrime)

| Function | CPU time | % MkPrime | % total | Notes |
|----------|----------|-----------|---------|-------|
| `TreeTools::postorder_order` | 0.230 s | 34% | 3.2% | Called every topology proposal to maintain preorder invariant |
| Unresolved static (addr in `block_gibbs_branch_sweep_impl` range) | 0.157 s | 23% | 2.2% | Static symbol; blockGibbsBranch=FALSE so likely mis-attributed or compiler-merged code |
| `nni_proposal_impl` | 0.107 s | 16% | 1.5% | NNI proposals |
| `spr_proposal_impl` | 0.093 s | 14% | 1.3% | SPR proposals |
| `[MkPrime.dll]` (inlined/unresolved) | 0.029 s | 4% | 0.4% | |
| `run_mcmc_batch_cpp` | 0.016 s | 2% | 0.2% | Batch loop overhead itself |
| `Rcpp::Vector<14>::Vector` (clone) | 0.016 s | 2% | 0.2% | Copy constructor |
| `cpp_partition_log_likelihood` (internal) | 0.015 s | 2% | 0.2% | Partition-level likelihood |
| `get_mcmc_state` | 0.015 s | 2% | 0.2% | State extraction |

### Key findings

1. **R↔C++ boundary dominates, not the C++ hot path.** The hardcoded
   `batchSize = 200` in `RunMkPrime.R` means 5 000 R↔C++ round-trips per
   million iterations. Each round-trip constructs an Rcpp `List`, crosses
   the R/C++ boundary, and the R side processes counts, stores samples,
   and rebuilds `scaleTunings`. This R-side overhead is **90%+ of CPU**.

2. **Felsenstein pruning is NOT the bottleneck.** `pruning_jc_flat`,
   `pruning_jc_acrv_flat`, `pruning_mkn_flat`, and all ascertainment
   functions (`constant_site_prob_jc`, `singleton_site_prob_jc`, etc.)
   are completely absent from the profile. The flat-buffer workspace
   optimization (M-063) was highly effective.

3. **Within the C++ path, tree reordering dominates.**
   `TreeTools::postorder_order` (34% of MkPrime time) is called every
   iteration to maintain the preorder invariant after topology changes.

4. **Clone overhead (OPP-6) is negligible** — only 0.2% of total CPU.
   The NNI/SPR vector clones don't register as material.

5. **OPP-1 (JC O(k²)→O(k)) is not the current bottleneck.** The pruning
   inner loop doesn't even appear in the profile at this dataset size.
   Still worth implementing for larger k and algorithmic correctness.

### Revised priority ranking

Based on VTune measurements (from highest to lowest expected impact):

1. **NEW: Adaptive batch size** (M-106, promoted to P1). Reducing 5000
   R↔C++ round-trips to 500 or fewer would eliminate ~90% of the dominant
   overhead. Expected: 2–5× wall-clock speedup.

2. **NEW: Cache/patch postorder after NNI** (M-108). NNI changes 4
   edges; patching the traversal order instead of recomputing from
   scratch removes the #1 C++ hotspot (3.2% of total).

3. **OPP-1 (JC O(k) product).** Not the current bottleneck but scales
   with mean kPrime. Implement when tackling datasets with k > 4.

4. **OPP-2 + OPP-3 (cache maxNode / qnorm midpoints).** Trivial to
   implement, negligible measured impact, still good housekeeping.

5. **OPP-4–7 (C++ allocation optimizations).** All negligible in the
   profile. Defer unless a future profile (after batch size fix) shows
   them appearing.

## Suggested profiling order (revised)

1. ~~Run VTune hotspot collection~~ — **DONE** (2026-03-28).
2. ~~Implement M-106 (adaptive batch size)~~ — **DONE**. Marginal gain (<1%)
   because C++ dominates at realistic acceptance rates. The VTune "90% R
   overhead" finding was an artifact of 0% acceptance rates in the profiling
   workload (all proposals rejected → no likelihood evaluation → trivial C++).
3. ~~Re-profile with Rprof~~ — **DONE**. Confirmed 90% of wall time is in
   `run_mcmc_batch_cpp`. R orchestration is <2%.
4. ~~Implement OPP-1 (JC O(k) product in ACRV path)~~ — **DONE**. 37% speedup.
5. ~~Tip-init hoist in ACRV pruning~~ — **DONE**. Additional 15% on top of OPP-1.
6. ~~Fast-path for uniform kPrime partitions~~ — **DONE**. Additional 9%.
7. ~~O(1) kPrime rollback~~ — **DONE**. Additional 4%.

**Cumulative speedup: 1.80×** (1055 → 1895 iter/s on Sun2018 hyoliths, 20k iter).

## Optimization round 2 (2026-03-28)

**Workload:** Sun2018 hyoliths (54 taxa, 225 chars, all transformational),
20 000 iterations, no Gibbs/weighted moves, TBR+NNI+SPR+BetaSimplex+kPrime.

| Optimization | iter/s | Δ vs baseline | Cumulative |
|---|---|---|---|
| Baseline (batchSize=200, no OPP-1) | 1055 | — | 1.00× |
| M-106: adaptive batch (500/5000) | 1055 | +0% | 1.00× |
| OPP-1: JC O(k) in ACRV | 1450 | +37% | 1.37× |
| Tip-init hoist | 1675 | +15% | 1.59× |
| Fast-path kPrime partitions | 1824 | +9% | 1.73× |
| O(1) kPrime rollback | 1895 | +4% | 1.80× |

Rprof after all optimizations: 90% `.Call`/`run_mcmc_batch_cpp`, 7% model
setup (one-time), <2% R orchestration. Ascertainment correction is ~12% of
per-evaluation cost. Remaining C++ time is dominated by the Felsenstein
pruning traversal (O(nEdge × nChar × k × nCat) per likelihood evaluation).

### Next optimization targets (diminishing returns)

- ~~M-108 (cache/patch postorder after NNI)~~: superseded by OPP-6b in-place
  NNI which eliminates `preorder_weighted_impl` entirely for NNI moves.
- Convert inner pruning functions from `NumericVector` to `const double*`
  to eliminate R allocation from the hot path (~0.1-0.2% estimated).
- Partial likelihood caching across proposals (major architectural change).

## Optimization round 3 (2026-03-28)

**Workload:** Sun2018 hyoliths (54 taxa, 225 chars, all transformational),
20 000 iterations, NNI+SPR+BetaSimplex+kPrime+scalar moves.

| Optimization | iter/s | Δ vs round 2 | Cumulative |
|---|---|---|---|
| Round 2 baseline | 1895 | — | 1.80× |
| OPP-6b: in-place NNI + O(1) BetaSimplex rollback + ascertainment tip-init hoist | 2611 | +38% | 2.47× |

**OPP-6b: In-place NNI** (dominant contributor). NNI only changes 2 parent
assignments. The reversed edge ordering remains a valid postorder for
Felsenstein pruning because the moved subtrees stay at positions after their
new parent's edge in the original preorder. This eliminates per-NNI:
- `preorder_weighted_impl` call (VTune's #1 C++ hotspot, 3.2% of total)
- 2× vector clone (parent + child, ~106 ints each)
- 1× Rcpp::List construction + destruction
- absLen vector + relBr recomputation

Rollback: save/restore 2 parent values (O(1)).

**O(1) BetaSimplex rollback.** `beta_simplex_impl` now outputs old values
via output params; cases 4 and 12 save 2 elements instead of cloning the
full relBrLengths vector.

**Ascertainment tip-init hoist.** All 4 ascertainment functions
(`constant_site_prob_jc/mkn`, `singleton_site_prob_jc/mkn`) now initialize
tips once before the ACRV category loop, matching the main pruning functions.

### Next optimization targets

- SPR/TBR still call `preorder_weighted_impl` (harder to avoid — topology
  changes are more complex than NNI's 2-parent swap).
- Partial likelihood caching across proposals (major architectural change).
- Convert inner pruning `NumericVector` args to `const double*` (~0.1-0.2%).
