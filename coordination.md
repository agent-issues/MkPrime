# MkPrime — Strategic Coordination

Last updated: 2026-03-30 10:20

## Project State

**Phase:** All core phases (1-10) complete. Phase 6b active: M-080 (TreeSearch hook) assigned to C.

**Recent milestones (since 2026-03-29):**
- Post-phase optimization complete: Bactrian proposals (M-118), pSPR (M-119), 2D joint Bactrian (M-120), node-level CL dirty flags (M-121), block Dirichlet branch proposal (M-125).
- Streaming/UI polish: buffer flush guards (M-132), warmup trace fix (M-133), treeThin (M-134), ESS label cleanup (M-136), auto-derive treeFile (M-137), test warning suppression (M-138), ticker simplification (M-139).
- **Dropped `coda` dependency**: native rank-normalized R-hat (Vehtari et al. 2021) and FFT-based ESS replace all `coda` calls. Convergence thresholds tightened (R-hat ≤ 1.01). 40 new unit tests.

**Open tasks:** M-131 (1 specific) + M-080 assigned (C). M-124 completed 2026-03-30 (PR draft).

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
- Independent runs for R-hat: nRuns≥2 enables rank-normalized split-R-hat
  (Vehtari et al. 2021) across cold chains from different runs. Replaces
  classical PSRF (Gelman–Rubin) as of 2026-03-30.
- Sequential execution of runs. Parallel execution implemented in Phase 10 (M-093).

**Exit criteria:**
- Parallel tempering improves mixing vs single chain (measured by ESS/iter)
- Chain swap acceptance 23–30% between adjacent pairs
- R-hat < 1.01 on simulated data with sufficient iterations
- Checkpoint save/restore produces identical continuation
- MkPosterior reports combined diagnostics across runs

### Phase 6: Progress display + GUI hooks
**Status:** Core COMPLETE (2026-03-26). Shiny module deferred to 6b.
**Goal:** Live MCMC traces, PNG-based Shiny integration, TreeSearch GUI module.

**Completed (M-037–M-039):**
- Progress callback infrastructure (plot_every, progress_fn in MkPrimeMCMC)
- mkp_trace_plot(): multi-panel base-R trace plots (per-run, warmup-aware)
- mkp_png_progress(): factory fn writing atomic PNG + JSON for Shiny polling

- M-040: EasyMkPrime() standalone Shiny app — bslib UI, callr background
  MCMC, PNG progress polling, results tabs (traces, summary, consensus).

**Phase 6b: TreeSearch GUI integration (tasks defined, not started)**

Architecture: MCMC runs as a detached OS process (survives Shiny timeout/browser
close); app is a log-file viewer that can connect/reconnect at will. Relies on
M-075 disk logging and M-076 ReadMkLog. Tasks are sequential:
- M-081: Cancel-file support in RunMkPrime (prerequisite for clean stop)
- M-078: `MkBayesianUi()` / `MkBayesianServer()` — reusable Shiny module
- M-079: Refactor `app.R` to use `MkBayesianServer`
- M-080: TreeSearch-side `mod_bayesian.R` + "Bayesian (Mk')" tab *(TreeSearch-a)*

### Phase 7: Extensions
**Status:** COMPLETE (2026-03-28). All items done including M-052 (Q-het) and M-053 (TBR).
**Goal:** `"informative"` coding, partition rates, siteMatrices, model
comparison, TBR, HMC.

**Completed (M-044–M-051):**
- `coding = "informative"` ascertainment correction (singleton site probs)
- Partition rate scalar `rate_neo` (neomorphic-specific rate multiplier)
- Stepping-stone marginal likelihood estimation (`mkp_stepping_stone()`)

**Completed (later phases):**
- ~~Beta-distributed Q-matrix heterogeneity (siteMatrices)~~ -> M-052 (done)
- ~~TBR moves~~ -> M-053 (done)
- ~~HMC for branch lengths~~ → replaced by M-054 block Gibbs branch sweep (done)

---


### Phase 8: MCMC performance — C++ inner loop
**Status:** COMPLETE (2026-03-27). All sub-phases done (8a-8e). ~85% cumulative speedup.
**Goal:** Port MCMC hot path from R to C++, eliminating per-iteration R
overhead. Use TreeTools C++ headers for tree manipulation.

**Sub-phases:**
- **8a (P1):** Quick R-side wins ? eliminate redundant `ape::reorder.phylo`
  from likelihood, pre-compute move weights, bypass validation in hot path.
- **8b (P2):** Port NNI/SPR/BetaSimplex proposals to C++, using TreeTools
  C++ headers for tree reordering and descendant-finding.
- **8c (P2):** C++ MCMC inner loop ? `do_move()` propose/evaluate/accept
  cycle, then the outer iteration loop with R callbacks for progress/sampling.
- **8d (P3):** Partial likelihood recalculation ? cache partition likelihoods,
  only recompute affected partitions per move type.

**Design notes:**
- TreeTools provides C++ headers (`preorder_edges_and_nodes`,
  `descendant_edges`, `postorder_order`) via `inst/include/TreeTools/`.
  Link via `LinkingTo: TreeTools` in DESCRIPTION.
- `TreeTools::Preorder()` is ~1.6x faster than `ape::reorder.phylo()`.
- The conditional-likelihood workspace (`std::vector<std::vector<double>>`)
  is currently reallocated per likelihood call ? should be pre-allocated
  once in a C++ state struct and reused.
- R remains responsible for: initialization, progress display, sample
  storage, checkpointing, adaptation logic, and result assembly.

### Phase 9: Advanced tree moves — Gibbs & Weighted proposals
**Status:** COMPLETE (2026-03-29). `feature/gibbs-weighted-moves` merged to main (M-105/M-109/M-111 optimization stack: 3.9× SPR, 1.9× swap, 2.1× combined Gibbs speedup).
**Goal:** Add five new topology/branch-length moves from RevBayes's
`mcmc_tree_moves` branch, adapted for MkPrime's unrooted representation
and flat-buffer likelihood engine. Add an adaptive move scheduler that
learns optimal move frequencies during warmup.

**New moves:**
| Move | Type | Cost | Key idea |
|------|------|------|----------|
| GibbsSPR | Topology | O(N) evals | Sample reattachment ∝ posterior weight |
| GibbsSubtreeSwap | Topology | O(N) evals | Sample swap partner ∝ posterior weight |
| WeightedBranchLengthScale | Branch | O(B) evals | Integrate over discretised branch fractions |
| WeightedSPR | Topology+Branch | O(N×B) evals | Marginalise topology × branch fraction jointly |
| WeightedSubtreeSwap | Topology+Branch | O(N×B) evals | Same as WeightedSPR with swap operation |

**Adaptive scheduler (M-092):** Track `accept_rate / mean_cost` per move
during warmup; reweight via softmax with annealing temperature; freeze
after warmup. User can pin specific weights via `moveWeights` in
`MkPrimeMCMC()`.

**Sub-phases:**
- **9a (P1 blockers):** M-083 (`compute_full_loglik()` helper), M-084 (subtree-swap op)
- **9b (Gibbs):** M-085 (GibbsSPR), M-086 (GibbsSubtreeSwap)
- **9c (Weighted):** M-087 (WeightedBranchLengthScale), M-088 (WeightedSPR), M-089 (WeightedSubtreeSwap)
- **9d (Integration):** M-090 (dispatch wiring), M-092 (adaptive scheduler), M-091 (mixing validation)

**Plan file:** `.positai/plans/2026-03-27-1100-gibbsweighted-tree-moves-for-unrooted-trees.md`

**Design decisions:**
- Gibbs/Weighted moves use `exp(β × logLik + logPrior)` weights — compatible
  with parallel tempering (β from chain's heat; prior unheated).
- Weighted moves default OFF (`weightedSpr = FALSE`) — too expensive for small
  datasets; user opts in for large trees where O(N×B) is worth the gain.
- Adaptive scheduler freezes weights at warmup end to preserve detailed balance.

**M-091 mixing validation** (2026-03-28, Agent B, Vinther 2008 hyoliths,
23 taxa / 54 chars, 20k iter, debug build):

| Config | Wall (s) | Slowdown | ESS(logP) | ESS(TL) | ESS/s(logP) | ESS/s(TL) | SPR acc | Gibbs SPR acc |
|--------|----------|----------|-----------|---------|-------------|-----------|---------|---------------|
| (a) Baseline | 3.1 | 1× | 3.5 | 26.4 | 1.13 | 8.48 | 2.2% | — |
| (b) +Gibbs | 23.4 | 7.5× | 8.7 | 37.7 | 0.37 | 1.61 | 2.1% | 31.6% |
| (c) +Weighted | 74.4 | 23.9× | 5.6 | 15.6 | 0.08 | 0.21 | 2.7% | 35.1% |

Findings: (1) Gibbs SPR acceptance is ~15× higher than standard SPR
(32% vs 2%), but each Gibbs proposal costs O(N) likelihood evaluations.
On this small dataset the per-iteration cost increase (7.5×) outweighs
the mixing gain, so **ESS/second favours the baseline for small trees**.
(2) WeightedSPR adds branch-fraction integration (O(N×B)) on top, making
it 24× slower per iteration with no net mixing benefit here. (3) All
configs show low absolute ESS — none mixed well in 20k iterations from a
random starting tree. (4) Caveats: debug build (-O0); the C++ inner
loops should be substantially faster at -O2. Gibbs/Weighted moves may
show positive ESS/s on larger trees where standard SPR acceptance is
near zero. A production-build comparison on a 50+ taxon dataset is
needed before concluding.  (5) Bug found during validation: the adaptive
scheduler floor (`wMin`) was divided by nFree, giving 0.7% per-move
floor instead of the intended 5%. Fixed in `148d31c`.

---

### Phase 10: Run-level parallelism via `future`
**Status:** COMPLETE (2026-03-27). All three sub-tasks done by Agent E on `feature/parallel-runs`. Branch ready for merge review.
**Goal:** Parallelize independent MCMC runs across CPU cores and HPC
nodes using the `future` package as a backend-agnostic parallelism
layer. The user sets `future::plan()` before calling `RunMkPrime()`;
the package never sets a plan itself (CRAN policy).

**Design decisions:**
- `future` in Suggests (opt-in). Falls back to sequential with a message
  if not installed or `parallel = FALSE` (the default).
- Persistent workers: each run is one long-lived `future::future()` call,
  not a stream of segments. Avoids per-segment XPtr reconstruction and
  HPC job-submission overhead.
- XPtr reconstruction: workers call `.InitMcmcData()` + `init_mcmc_state()`
  internally from serialized R inputs — same pattern as `ResumeMkPrime()`.
- Convergence signaling: orchestrator polls log files via `ReadMkLog()`;
  writes cancel files (M-081) when ESS/PSRF criteria are met.
- Streaming required for mid-run convergence stopping; auto-assigned to
  `tempfile()` if `logFile` is NULL in parallel mode.
- Within-run chain parallelism (OpenMP etc.) is explicitly out of scope.

**Sub-phases:**
- **M-094 (P1):** ✓ Extract `.RunMkPrimeSingleRun()` — refactor + dedup. Fixed multi-run checkpoint regression.
- **M-095 (P1):** ✓ Parallel orchestration (`.RunParallelRuns()`, `.CheckConvergenceFromLogs()`, `future` in Suggests).
- **M-096 (P1):** ✓ Tests (`test-parallel.R`, sequential + multisession) + HPC usage section in `RunMkPrime()` roxygen.
- **Note:** `man/` not regenerated (stale installed package); run `roxygen2::roxygenise()` after installing `feature/parallel-runs`.

**Plan file:** `.positai/plans/2026-03-27-1307-plan.md`

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
