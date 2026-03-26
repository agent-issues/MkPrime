# MkPrime — Strategic Coordination

Last updated: 2026-03-26

## Project State

**Phase:** 5 (Parallel tempering + convergence) — COMPLETE. Ready for Phase 6 (Progress display + GUI hooks).

MkPrime is a new R package for Bayesian phylogenetic inference under the
Mk' model. The architecture follows StratoBayes (C++ hot loop via Rcpp,
R orchestration). The full design plan is in `../.positai/plans/mkprime-r.md`.

## Phase Plan

### Phase 1: Foundation (package skeleton + data layer)
**Status:** COMPLETE (2026-03-26)
**Goal:** Working package that accepts `phyDat` input, classifies characters
into neomorphic/transformational/known types, partitions by kObs, and produces
the internal data structures the C++ likelihood engine will consume.

**Exit criteria:**
- `R CMD build` succeeds
- `MkPrimeData()` correctly classifies and partitions test datasets
- Edge cases handled (invariant chars, missing data, single-state chars)
- Unit tests pass

### Phase 2: Likelihood engine (C++)
**Status:** COMPLETE (2026-03-26)
**Goal:** Correct likelihood computation for all three character types,
including ACRV and ascertainment bias correction.

**Exit criteria:**
- JC(k') and MkN P(t) matrices correct (validated against analytical values)
- Felsenstein pruning produces correct likelihoods on hand-crafted examples
- Ascertainment correction and relabelling correction implemented
- ACRV integration working
- **Likelihoods match RevBayes to within ~1e-6** on reference datasets

### Phase 3: Basic MCMC (fixed topology)
**Status:** COMPLETE (2026-03-26)
**Goal:** Single-chain MCMC on continuous parameters + k'_i with fixed tree.

**Exit criteria:**
- MH proposals for continuous parameters (Scale, BetaSimplex) ✓
- BoundedIntegerWalk for k'_i ✓
- Produces valid posterior samples ✓
- Console progress bar ✓
- MkPosterior result object with basic summaries ✓
- Integration validation: posterior recovers simulated parameters ✓

### Phase 4: Tree search
**Status:** COMPLETE (2026-03-26)
**Goal:** NNI and SPR topology proposals, tree logging, parsimony-informed
priors, integration validation. Tasks M-022 through M-028.

### Phase 5: Parallel tempering + convergence
**Status:** COMPLETE (2026-03-26)
**Goal:** Multiple chains with adaptive temperature spacing, convergence
monitoring, stopping rules, checkpointing.

**Task breakdown:**
- M-029: Parallel tempering core (multi-chain + temperature ladder)
- M-030: Chain swap proposals between adjacent temperatures
- M-031: Adaptive temperature tuning
- M-032: Independent runs (nRuns outer loop)
- M-033: Convergence monitoring (ESS + PSRF)
- M-034: Stopping rules (ESS/PSRF/time/iter thresholds)
- M-035: Checkpointing (save/restore)
- M-036: MkPosterior multi-run updates

**Design decisions:**
- Geometric temperature ladder: β_i = heat^((i-1)/(nChains-1)). Cold
  chain β=1, hottest β=heat (default 0.2).
- Heated acceptance: only likelihood is tempered (β×logLik + logPrior).
  Prior is unheated to maintain proper support.
- State stores unheated logLik/logPrior. Heated posterior computed on the
  fly for MH acceptance.
- Per-chain independent tuning: heated chains need different proposal widths.
- Independent runs for PSRF: nRuns≥2 enables Gelman–Rubin convergence
  diagnostics across cold chains from different runs.
- Sequential execution of runs. Parallel (future/furrr) deferred.

**Exit criteria:**
- Parallel tempering improves mixing vs single chain (measured by ESS/iter)
- Chain swap acceptance 23–30% between adjacent pairs
- PSRF < 1.05 on simulated data with sufficient iterations
- Checkpoint save/restore produces identical continuation
- MkPosterior reports combined diagnostics across runs

### Phase 6: Progress display + GUI hooks
**Status:** Not started
**Goal:** Live MCMC traces, PNG-based Shiny integration, TreeSearch GUI module.

### Phase 7: Extensions (future)
**Status:** Not started
**Goal:** `"informative"` coding, partition rates, siteMatrices, model
comparison, TBR, HMC.

---

## Agent Allocation

Currently single-agent. Multi-agent (2 agents) planned from Phase 4 onward,
when tree moves, tempering, and convergence monitoring become independent
workstreams.

## Related Projects

| Project | Relationship |
|---------|-------------|
| `mkprime/` | RevBayes-scripted Mk' inference — reference implementation and results |
| `revbayes-ms/` | Native C++ `dnMkPrime` for RevBayes — mathematical reference |
| `neotrans/` | Best-practice RevBayes phylogenetics — model specification reference |
| `StratoBayes/` | Architectural template — MCMC engine, tempering, progress display |
| `TreeSearch` | Future GUI integration via `Suggests: MkPrime` |
