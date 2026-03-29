# M-103: Profile the C++ MCMC Hot Path (VTune)

## Goal

Identify the top CPU hotspots in the compiled MCMC engine using Intel VTune
hardware sampling. This determines which optimization (M-104, M-105, M-106,
OPP-4 through OPP-9) yields the most gain.

## Context

- **VTune IS available** on this system (Intel i7-10700, 10th gen — hardware
  sampling works). The `r-package-profiling` skill documents the full workflow.
- `benchmark/vtune-driver.R` already exists — it loads Sun2018 (a hyolith dataset) (54 taxa,
  225 chars, 3 partition types) and runs 200 000 iterations under the
  standard move schedule.
- `benchmark/PROFILING.md` documents the build recipe and A/B renaming pattern.
- OPP-1 through OPP-6 were identified by static analysis. VTune will confirm
  which of these actually matter at runtime and whether anything was missed.

## Steps

### Step 1: Build with symbols

Create a temporary `src/Makevars.win` with profiling flags:
```makefile
PKG_CXXFLAGS = $(SHLIB_CXXFLAGS) -g -fno-omit-frame-pointer
```

Then install into `.vtune-lib` with the `DLLFLAGS` override (per the
`r-package-profiling` skill):

```r
vtune_lib <- file.path(getwd(), ".vtune-lib")
stale <- Sys.glob(file.path("src", c("*.o", "*.dll", "*.so")))
if (length(stale)) file.remove(stale)
old_mf <- Sys.getenv("MAKEFLAGS")
Sys.setenv(MAKEFLAGS = "DLLFLAGS=-static-libgcc")
install.packages(".", lib = vtune_lib, repos = NULL, type = "source",
                 INSTALL_opts = "--no-multiarch")
Sys.setenv(MAKEFLAGS = old_mf)
```

Clean up: remove `src/Makevars.win` after the install succeeds.

### Step 2: Collect VTune hotspots

Run from the `mkp/` root:
```bash
vtune -collect hotspots -result-dir vtune-out \
      -- Rscript benchmark/vtune-driver.R
```

If hardware PMU driver is not loaded, fall back to software sampling:
```bash
vtune -collect hotspots -knob sampling-mode=sw -result-dir vtune-out \
      -- Rscript benchmark/vtune-driver.R
```

### Step 3: Read the report

Filter to the package DLL:
```bash
vtune -report hotspots -result-dir vtune-out -filter "module=mkp.dll"
```

This gives a ranked list of functions by CPU time percentage. Expected
candidates (from static analysis):

| Function | Static-analysis prediction |
|----------|---------------------------|
| `pruning_jc_flat` / `pruning_jc_acrv_flat` | Dominant — O(k²) inner loop (OPP-1) |
| `pruning_mkn_flat` / `pruning_mkn_acrv_flat` | Moderate — full 2×2 matrix multiply |
| `constant_site_prob_jc` / `singleton_site_prob_jc` | Per-partition heap alloc (OPP-4) |
| `do_move_impl` | Dispatch + clone overhead (OPP-6) |
| `nni_proposal_impl` / `spr_proposal_impl` | Proposal generation |
| `preorder_weighted_impl` | Post-proposal tree reordering |
| `beta_simplex_impl` | Branch-length proposals |
| R runtime (`R::qbeta`, `R::rbeta`, etc.) | R distribution function cost |

### Step 4: Interpret and cross-reference

For each hotspot:

1. **What % of total CPU time?** (VTune's self-time column)
2. **Does it match a known OPP?** If so, confirm or adjust expected gain.
3. **Is it a surprise?** If a function not in OPP-1–9 appears, investigate.
4. **Is R runtime significant?** If `libR.dll` or `Rmath.dll` functions
   appear in the top 10, it means R↔C++ boundary or distribution functions
   are material — this would motivate caching or C++ reimplementations.

### Step 5: Update PROFILING.md and file tasks

- Append a "VTune results" section to `benchmark/PROFILING.md` with the
  hotspot table (function, % CPU, interpretation).
- Re-rank OPP-4 through OPP-9 based on actual measurements.
- File new `M-nnn` tasks in `to-do.md` for any newly-identified bottlenecks.
- Update M-104/M-105 descriptions if the profile changes their priority.

### Step 6: Cleanup

- Remove `vtune-out/` directory (or `.gitignore` it) — it's large and
  machine-specific.
- Confirm `src/Makevars.win` was removed (Step 1 cleanup).
- Remove `.vtune-lib/` or leave for future runs (already gitignored or
  should be added to `.gitignore`).

## Out of scope

- Implementing any optimizations (those are M-104, M-105, M-106).
- Gibbs/weighted move profiling (M-102 already benchmarked those; they're
  default OFF and not in the standard move schedule).
- A/B benchmarking (that's for measuring the impact of changes, not for
  identifying hotspots).
