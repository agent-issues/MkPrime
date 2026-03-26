# MkPrime — Strategic Coordination

Last updated: 2026-03-26

## Project State

**Phase:** 3 (Basic MCMC) — COMPLETE. Ready for Phase 4 (Tree search).

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
**Status:** Not started
**Goal:** SPR and NNI topology proposals, compound Dirichlet branch length
prior, tree logging.

**Exit criteria:**
- Topology moves preserve tree invariants
- Branch length prior correctly implemented
- Tree logging to file (Newick)
- Recovers known tree from simulated data

### Phase 5: Parallel tempering + convergence
**Status:** Not started
**Goal:** Multiple chains with adaptive temperature spacing, convergence
monitoring, stopping rules, checkpointing.

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
