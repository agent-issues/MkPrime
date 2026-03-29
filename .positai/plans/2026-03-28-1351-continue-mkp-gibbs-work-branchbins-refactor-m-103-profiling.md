# BranchBins refactor commit + Gibbs/Weighted hot-path profiling

## Overview

Two pieces of work:

1. **Complete the BranchBins refactor** — compile, test, commit the
   uncommitted changes that move `BranchBins` from a static global with
   lazy init into `McmcData` (precomputed at MCMC setup).

2. **Profile the Gibbs/weighted move hot path** — complementary to
   Agent C's M-103 (which profiles the standard move schedule on `main`).
   This profiles the *Gibbs/weighted-specific* code paths to determine
   where the 7–60× per-iteration cost premium (M-102) comes from.

## Part 1: BranchBins refactor

### What changed

- `mcmc_state.h`: New `BranchBins` struct with `init()` method; added
  `branchBins` member to `McmcData`.
- `mcmc.cpp`: Removed static `s_branchBins` global and `get_branch_bins()`.
  All five weighted/block-Gibbs functions now read `data->branchBins` and
  `bins.concentration` instead of calling a lazy-init function and
  recomputing `2.0 * nBins` inline.
- `mcmc_likelihood.cpp`: `set_branch_bins()` now also calls
  `d->branchBins.init(nBins)`.
- `R/RcppExports.R`: Line-ending change only (no signature change).

### Steps

1. Clean stale objects: `rm -f src/*.o src/*.dll`
2. Build: `Rscript -e "pkgbuild::compile_dll(debug = FALSE); devtools::load_all()"`
3. Run targeted tests (Gibbs/weighted move tests):
   `Rscript -e "devtools::load_all(); testthat::test_file('tests/testthat/test-gibbs-moves.R')"`
   and/or the full test suite
4. Commit on `feature/gibbs-weighted-moves`

## Part 2: Gibbs/weighted hot-path profiling (VTune)

### Motivation

M-102 showed Gibbs/weighted moves cost 7–60× more per iteration than
standard moves (-O2 build). The question is *where* that time goes:

| Candidate | Why it might dominate |
|-----------|---------------------|
| `clone()` allocations | gibbs_spr: 2+nCand clones; weighted_spr: 5+2×nCand+nBins×nCand |
| Felsenstein pruning via `compute_full_loglik_at` | Called once per candidate (Gibbs) or once per candidate×bin (weighted) |
| `preorder_weighted_impl` | Called per candidate (Gibbs SPR) or per candidate×bin (weighted SPR) |
| `R::qbeta`/`R::rbeta`/`R::dbeta` | Per-bin distribution function calls |
| `swap_subtrees_impl` | Called per partner in `gibbs_subtree_swap_impl` |

VTune will rank these by actual CPU time.

### Approach

**Driver script:** Create `benchmark/vtune-driver-gibbs.R` — a variant of
the existing `vtune-driver.R` that exercises Gibbs/weighted move
configurations. Run four workloads (matching M-091/M-102 configs):

| Config | MCMC flags | Rationale |
|--------|-----------|-----------|
| (b) +Gibbs | `gibbsSpr=T, gibbsSubtreeSwap=T` | Profiles clone + pruning cost |
| (c) +Weighted | `weightedSpr=T` | Profiles bin-marginalisation overhead |
| (d) +BlockGibbs | (b) + `blockGibbsBranch=T` | Profiles sweep cost |
| (e) +WeightedBranch | `weightedBranchScale=T` | Profiles bin-based branch scaling |

Each config runs enough iterations to give VTune ~30s of CPU time.
Configs (b)–(e) can share setup and just vary the MCMC flags.

### Build steps

1. Create temporary `src/Makevars.win`:
   ```makefile
   PKG_CXXFLAGS = $(SHLIB_CXXFLAGS) -g -fno-omit-frame-pointer
   ```
2. Clean stale objects
3. Install into `.vtune-lib` with `DLLFLAGS` override (no `-s` strip)
4. Remove `src/Makevars.win`

### VTune collection

```bash
VTUNE="C:/Program Files (x86)/Intel/oneAPI/vtune/latest/bin64/vtune.exe"

# One collection per config (b, c, d, e):
"$VTUNE" -collect hotspots -result-dir vtune-gibbs-b \
    -- Rscript benchmark/vtune-driver-gibbs.R b

"$VTUNE" -collect hotspots -result-dir vtune-gibbs-c \
    -- Rscript benchmark/vtune-driver-gibbs.R c

# etc.
```

If hardware PMU fails, fall back to `-knob sampling-mode=sw`.

### Analysis

```bash
"$VTUNE" -report hotspots -result-dir vtune-gibbs-b -filter "module=mkp.dll"
```

For each config, identify top-5 hotspots by self-time %. Cross-reference
with the candidate table above. Key questions:

1. What fraction of time is `compute_full_loglik_at` (i.e., pruning)?
   If >60%, then M-105 (partial likelihood reuse) is the top priority.
2. What fraction is `clone()` / heap allocation? If >10%, M-104 is
   worth pursuing.
3. What fraction is `preorder_weighted_impl`? If significant, caching
   the preorder traversal between candidates matters.
4. Are `R::qbeta`/`R::rbeta` material? If so, consider C++ replacements.

### Deliverables

- `benchmark/vtune-driver-gibbs.R` — Gibbs/weighted VTune driver script
- Hotspot results appended to `benchmark/PROFILING.md` (new section)
- Updated priority ranking for M-104, M-105 based on measured data
- New tasks filed in `to-do.md` if unexpected hotspots appear

### Cleanup

- Remove `vtune-gibbs-*` result directories (or add to `.gitignore`)
- Remove `.vtune-lib/` (or keep for future; add to `.gitignore`)
- Ensure `src/Makevars.win` is removed

## Execution order

1. Part 1 first (BranchBins refactor) — ensures profiling runs on clean,
   committed code
2. Part 2 (VTune profiling) — uses the freshly committed code
