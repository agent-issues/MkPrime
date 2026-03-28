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

## Suggested profiling order

1. Run VTune hotspot collection (see above) on current build to rank
   actual CPU time by function.
2. Implement OPP-1 (JC O(k) product) — highest algorithmic leverage.
3. Benchmark with `bench-ab.R` to quantify the gain before moving on.
4. Implement OPP-2 and OPP-3 together (both are trivial and additive).
5. Profile again to confirm remaining hotspots and decide on Tier 2 items.

---

## VTune hotspot results — Gibbs/weighted moves (2026-03-28)

**Build:** `-O2 -g -fno-omit-frame-pointer`, symbols preserved (no `-s`)
**CPU:** Intel i7-10700 (10th gen), 2.904 GHz, 16 logical cores
**Dataset:** Sun2018 (54 taxa, 225 characters, 3 partition types)
**Collection:** User-mode sampling (software), no hardware PMU driver

### Top hotspots by config (MkPrime.dll only, % of total CPU)

| Function | Config b (Gibbs) | Config c (WeightedSPR) | Config d (BlockGibbs) |
|----------|:---:|:---:|:---:|
| `pruning_jc_acrv` | 42.0% (11.7s) | 41.5% (8.9s) | 40.2% (9.4s) |
| `_expl_internal` (exp) | 15.2% (4.2s) | 13.0% (2.8s) | 14.4% (3.4s) |
| `std::fill` (CL zeroing) | 6.3% (1.8s) | 6.5% (1.4s) | 5.9% (1.4s) |
| `constant_site_prob_jc` | 4.2% (1.2s) | 4.3% (0.9s) | 4.4% (1.0s) |
| `malloc_base` (ucrtbase) | 4.3% (1.2s) | 4.5% (1.0s) | 4.3% (1.0s) |
| `vector<vector<double>> ctor` | 1.6% | 1.0% | — |
| `vector<vector<double>> dtor` | 1.0% | 1.9% | — |
| `Rcpp::Matrix::offset` | 1.7% | 2.1% | — |
| `preorder_weighted_impl` | <0.1% | — | — |
| `clone`-related | <0.1% | — | — |
| R math (qbeta/rbeta/dbeta) | not visible | not visible | not visible |

### Key findings

1. **Felsenstein pruning dominates uniformly.** `pruning_jc_acrv` + `exp()` =
   ~55–57% of total CPU time across all configs. This is the O(k²) inner
   loop (OPP-1) plus `exp(Qt)` per edge. The profile is remarkably stable
   regardless of which Gibbs/weighted moves are enabled.

2. **CL workspace zeroing (`std::fill`) costs 6–7%.** The flat CL buffer is
   zeroed before each pruning pass. For 54 taxa with 6 ACRV categories,
   this is ~650 × nChar × kStates doubles zeroed per full LL evaluation.
   Could be reduced by lazy zeroing (only clear what's used) or by tracking
   which entries are stale.

3. **Ascertainment correction (`constant_site_prob_jc`) at 4–5%.** This uses
   the old jagged `vector<vector<double>>` allocation (OPP-4). The
   ctor+dtor overhead is ~2.5% combined. Threading `ClWorkspace*` through
   would eliminate this.

4. **`clone()` and `preorder_weighted_impl` are negligible.** M-104 (reduce
   clone overhead) has much lower priority than expected. The clone
   allocations are tiny (~nEdge integers/doubles) compared to the pruning
   cost per candidate. `preorder_weighted_impl` is under 0.1% in all
   configs.

5. **`malloc_base` at 4.3%** is distributed across all allocations (Rcpp
   vector construction, jagged CL in ascertainment, etc.). No single
   allocation source dominates.

6. **R distribution functions are not visible.** `qbeta`, `rbeta`, `dbeta`
   do not appear in the weighted-move profiles. The BranchBins precomputation
   means `qbeta` is only called once at init, and the per-iteration Beta
   distribution draws are negligible.

### Revised optimization priority

Based on measured hotspots, the priority order for Gibbs/weighted move
performance is:

1. **OPP-1: JC O(k) product** — 42% of CPU. Highest leverage, helps all
   move types equally.
2. **M-105: Partial likelihood reuse** — reduces the *number* of
   `pruning_jc_acrv` calls per Gibbs/weighted move. Currently each
   candidate in Gibbs SPR does a full-tree pruning; reusing CLs for
   unchanged subtrees would cut this to O(depth) per candidate.
3. **OPP-4: Ascertainment CL allocation** — 4–5% direct + 2.5% ctor/dtor.
   Thread ClWorkspace through ascertainment functions.
4. **CL zeroing** — 6% from `std::fill`. Lazy zeroing or dirty-flag
   tracking.
5. ~~M-104: clone() reduction~~ — deprioritised; negligible in profile.
6. **OPP-2/OPP-3** — trivial gains; implement opportunistically.
