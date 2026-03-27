# MkPrime Task Queue

## How this works

- Tasks are sorted by priority (highest first within each status group).
- An agent claims a task by changing its status to `ASSIGNED (X)`.
- On completion, **delete** the row from this file and append a summary row
  to `completed-tasks.md`.

Task IDs use `M-nnn` prefix (MkPrime) to avoid collision with TreeSearch
`T-nnn` IDs.

Standing tasks (S-RED, S-PROF, S-COORD) are always present. When one is
completed, reset it to OPEN. Their effective priority is dynamic:
- ≥6 OPEN specific tasks → standing tasks are P3
- 3–5 OPEN specific tasks → standing tasks are P2
- <3 OPEN specific tasks → standing tasks are P1

---

## Phase 9: Advanced tree moves — Gibbs & Weighted proposals

All work on this phase happens in the `mkp-gibbs` worktree
(`feature/gibbs-weighted-moves` branch). Plan file:
`.positai/plans/2026-03-27-1100-gibbsweighted-tree-moves-for-unrooted-trees.md`.

### 9a: Infrastructure (blockers for all subsequent moves)

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-083 | P1 | DONE (B) | **Extract `compute_full_loglik()` helper.** All five new Gibbs/Weighted moves need to evaluate the full posterior likelihood from a C++ context — multiple times per proposal. Currently the only entry point is the R-callable `MkpLogLikelihood()` and the hot-loop-internal call in `mcmc_engine.cpp`. Expose a `compute_full_loglik(state, data)` C++ function (and optionally `compute_full_loglik_with_prior()`) callable from within proposal code. This is the P1 blocker for M-085–M-089. Signature: takes the flat-buffer state struct + partition data (pointers already held in `McmcState`) and returns `double` log-likelihood (unheated). No R boundary crossing; pure C++. |
| M-084 | P1 | DONE (B) | **Subtree swap topology operation + `get_valid_swap_partners()`.** NNI and SPR are already in `tree_moves.cpp`. Add: (a) `swap_subtrees(parent, child, relBrLengths, nodeA, nodeB)` — detaches two non-nested, non-sibling subtrees and exchanges their attachment points, then calls `postorder_reorder()`; returns updated parent/child/relBrLengths. (b) `get_valid_swap_partners(parent, child, pruneNode)` — returns all nodes that are non-nested and non-sibling relative to `pruneNode` (valid swap candidates). (c) Unit tests: round-trip swap returns original tree; swap on 5-taxon tree matches hand-computed result; `get_valid_swap_partners` excludes nested/sibling nodes correctly. This operation is used by both M-086 (GibbsSubtreeSwap) and M-089 (WeightedSubtreeSwap). |

### 9b: Gibbs moves

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-085 | P2 | DONE (B) | **GibbsSPR.** Enumerate all valid SPR reattachment candidates (non-nested, non-sibling edges), compute the tempered posterior weight `exp(β × logLik + logPrior)` for each, normalise to a probability vector, sample one candidate proportionally, and apply. Return `logHastings = log(1/w_chosen) - log(1/w_original_position)` so the MH ratio is correct. Cost: O(N) likelihood evals where N = number of candidates. Tempered: `exp(β × logLik)` (β from chain's heat; logPrior unheated). Depends on M-083. |
| M-086 | P2 | DONE (B) | **GibbsSubtreeSwap.** Same pattern as M-085 but uses the subtree-swap operation from M-084. Enumerate valid swap partners, weight by posterior, sample ∝ weight, apply swap, return correct Hastings ratio. Depends on M-083, M-084. |

### 9c: Weighted moves

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-087 | P2 | DONE (B) | **WeightedBranchLengthScale.** Discretize the branch-fraction proposal space into B bins (default B = 10) using Beta(0.25, 0.25) quantiles as break points. For each bin midpoint, evaluate the posterior, weight by `exp(β × logLik) × proposal_density`. Sample a bin from the normalised weights (CDF inverse), then draw the actual branch fraction from a Beta centred on that bin's midpoint. Compute the piecewise Hastings ratio. Cost: O(B) likelihood evals. Depends on M-083. |
| M-088 | P2 | DONE (B) | **WeightedSPR.** Gibbs SPR (M-085) extended to also integrate over branch fractions at each candidate position: for each candidate topology, marginalise over B branch-fraction bins, producing a marginal weight. Sample topology from marginal weights, then sample branch fraction from conditional distribution given chosen topology. Compute combined Hastings ratio. Cost: O(N × B) likelihood evals. Depends on M-083, M-085. |
| M-089 | P3 | OPEN | **WeightedSubtreeSwap.** Same as M-088 but uses subtree-swap topology operation from M-084. Cost: O(N × B) evals. Depends on M-083, M-084, M-086, M-088 (reuse integration logic). |

### 9d: Integration, scheduling & validation

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-090 | P2 | OPEN | **Wire new moves into dispatch + `MkPrimeMCMC()` parameters.** Add the five new moves to the `do_move_impl()` dispatch table in `mcmc_engine.cpp`. Add user-facing parameters to `MkPrimeMCMC()`: `gibbsSpr` (logical, default TRUE), `gibbsSubtreeSwap` (logical, default TRUE), `weightedBranchScale` (logical, default FALSE — expensive), `weightedSpr` (logical, default FALSE), `weightedSubtreeSwap` (logical, default FALSE), `nBranchBins` (integer, default 10L). Update `MkPrimeModel()` docs and `MkPrimeMCMC()` docs. Add move-presence guards so disabled moves are excluded from the move pool entirely. Depends on M-085, M-086, M-087, M-088, M-089. |
| M-092 | P2 | OPEN | **Adaptive move scheduler.** During warmup, track per-move `(n_accepted, n_proposed, wall_time_ns)` in the C++ state struct. Every adaptation window (200 iters), recompute move weights: `score[m] = accept_rate[m] / mean_cost_s[m]`; update weights via softmax with temperature that anneals from 2.0 → 0.5 over warmup, floored at `w_min = 0.05`. Freeze weights at warmup end (detailed balance preserved post-warmup). User override: `moveWeights = list(nni = 0.3, spr = 0.3, ...)` in `MkPrimeMCMC()` pins specific move weights and excludes them from adaptation. Log final scheduled weights to the run header in the TSV log file. Depends on M-090. |
| M-091 | P3 | OPEN | **Mixing validation: Gibbs/Weighted vs standard moves.** On the Vinther hyoliths dataset, compare ESS/wall-time between: (a) baseline (NNI + SPR only), (b) + GibbsSPR + GibbsSubtreeSwap, (c) + WeightedSPR. Use `MkpWatchLog()` output and `ConvergenceDiagnostics()` to compare. Report topology ESS/hour and scalar-parameter ESS/hour. Document findings in a comment in `coordination.md`. Depends on M-090, M-092. |

---

## Phase 8: MCMC performance — C++ inner loop

The MCMC inner loop currently runs entirely in R. Every iteration incurs
R function-call overhead, list allocation/copy, GC pressure, and redundant
computation. This phase ports the hot path to C++, using TreeTools C++
headers where available.

### 8b: Port proposals to C++

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-058 | P2 | DONE (C) | **C++ NNI proposal.** Port `ProposeNni()` (proposals.R:124-177) to C++. Current R implementation does: find internal edges via `which()`, pick random internal edge, find children/siblings via `which()`, swap two edge endpoints, call `ape::reorder.phylo`. In C++: operate directly on integer parent/child vectors, use TreeTools `preorder_edges_and_nodes` header for reordering. Return new edge matrix + relative branch lengths + logHastings. Avoids R list allocation and `ape::reorder.phylo` overhead per NNI proposal. |
| M-059 | P2 | DONE (C) | **C++ SPR proposal.** Port `ProposeSpr()` (proposals.R:202-289) to C++. The R version's `.Descendants()` BFS (lines 294-307) uses `c()` concatenation in a while loop — classic R antipattern. Use TreeTools `descendant_edges` C++ header or write a direct BFS on integer arrays. The rest of the SPR logic (prune node, suppress, regraft, reorder) is straightforward edge-index manipulation. Return new edge matrix + relative branch lengths + logHastings. |
| M-060 | P2 | DONE (C) | **C++ BetaSimplex proposal.** Port `ProposeBetaSimplex()` (proposals.R:24-80) to C++. Currently allocates a new vector per call and computes two `dbeta()` evaluations. For branch-length simplex proposals on large trees (nEdge can be hundreds), this vector allocation is significant. C++ version modifies in-place and returns logHastings. |

### 8c: C++ MCMC inner loop

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-061 | P2 | DONE (C) | **C++ `.DoMove()` — propose/evaluate/accept cycle.** This is the single highest-impact change. Currently `.DoMove()` (RunMkPrime.R:754-829) is called millions of times and each call: (1) copies the full state list (`proposed <- state`), (2) dispatches through R `switch()`, (3) calls the proposal function in R, (4) calls `LogPrior()` in R, (5) constructs a temporary tree and calls `MkpLogLikelihood()`, (6) returns a new R list. On rejection (~70-80% of calls), all allocations are garbage. A C++ `do_move()` would hold state in a C++ struct, call C++ proposals directly, call the existing C++ pruning functions, and return only accept/reject + updated state — eliminating all per-iteration R overhead. Depends on M-058, M-059, M-060. Also requires porting `LogPrior()` to C++. |
| M-062 | P2 | DONE (B) | **C++ outer iteration loop.** Once `.DoMove()` is in C++, port the per-chain/per-run/per-iteration triple loop (RunMkPrime.R:127-198) to C++. The C++ loop handles: move selection (weighted sampling), per-chain state updates, chain swap proposals, acceptance tracking. Call back to R only for: progress display (every N iterations), sample storage (every thin-th iteration post-warmup), adaptation (every 200 iterations during warmup), and checkpointing. This eliminates R interpreter overhead from the hot loop entirely. Depends on M-061. |
| M-063 | P3 | DONE (B) | **C++ state struct with pre-allocated workspace.** |

### 8d: Likelihood computation efficiency

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-064 | P3 | DONE (C) | **Partial likelihood recalculation for parameter-only moves.** When only a scalar parameter changes (tree_length, rate_loss, rate_log_sd, rate_neo, p), the tree topology and tip states are unchanged. For `tree_length` and `rate_neo` moves, only branch lengths change — the transition probability matrices change but the tree structure doesn't. For `rate_log_sd` moves, only the ACRV rate categories change. For `kPrime` moves, only one partition sub-group changes. Currently `MkpLogLikelihood()` recomputes everything from scratch. Implement move-type-specific partial recomputation: cache partition likelihoods and only recompute the affected partition(s). |
| M-065 | P3 | DONE (C) | **Avoid re-extracting parent/child/edgeLength from tree every likelihood call.** `MkpLogLikelihood()` extracts `tree$edge[,1]`, `tree$edge[,2]`, `tree$edge.length` as R vectors, then passes them to C++. In the C++ inner loop (M-061+), these would be held directly in the C++ state struct, eliminating this per-call overhead. This task is partly subsumed by M-063 but is listed separately for tracking. |

---

## Progress display & UX

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-076 | P2 | DONE (C) | **`MkpWatchLog()` — live trace plot from disk log files.** Polling watcher that reads one or more log files (from `logFile` in `MkPrimeMCMC`) and redraws trace plots each cycle. Ctrl+C stops cleanly and returns last data invisibly. `MkLogPaths()` helper expands base log name to per-run paths. Auto-selects key scalar params (excludes br_, log_likelihood). |
| M-075 | P2 | DONE (C) | **Streaming log output (`logFile` / `bufferSize`).** Implement Tracer-compatible TSV output during MCMC. `logFile` writes scalar parameter samples to disk in batches (`bufferSize`), flushed at each buffer boundary. `ReadMkLog()` reads log files back into R. `nRuns > 1` creates per-run log files (`run_1.log`, `run_2.log`, …). Checkpoint includes streaming state for seamless resume. |
| M-074 | P2 | DONE (C) | **Per-parameter progress table during MCMC.** At each `checkEvery` interval (post-warmup), print a convergence table showing ESS and PSRF (when `nRuns >= 2`) for each key parameter — mirroring the format of `print.MkpDiagnostics`. Extend `.CheckConvergence()` to return full per-parameter `ess` and `psrf` vectors (and work for `nRuns = 1`, ESS-only). Add `.PrintProgressTable()` helper in `Convergence.R`. Update `doCheck` blocks in both `RunMkPrime()` and `ResumeMkPrime()`. |
| M-071 | P1 | DONE (C) | **Progress bar redesign.** Replace `{cli::pb_bar} current/total \| accept% ` with `iter \| ESS: N \| logP: N \| PSRF: N`. A bar is misleading because the run works toward a *convergence* condition, not a fixed iteration count. Make `nIter` default to `Inf` in `MkPrimeMCMC()` (rely on `minEss`/`maxPsrf`/`maxTime` stopping criteria). Update `ResumeMkPrime` progress display to match. Update vignette MCMC configs and docs. |
| M-072 | P2 | DONE (C) | **hyoliths.qmd tree-summary-demo chunk fixes.** (a) Add `par(mar = rep(0, 4))` before `plot(consensus, ...)` — do this via a chunk `fig.par` option or an explicit `par()` call. (b) The plot title says "rogues excluded" even when no rogue detection has run; make the title conditional on whether rogues were actually identified and removed. |

## Convergence diagnostics

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-073 | P1 | DONE (C) | **Tree ESS in ConvergenceDiagnostics.** Add `trees = TRUE` argument. If `treess` and `TreeDist` installed and `!isFALSE(trees)`, compute topology ESS using `treess::treess(perRunTrees, TreeDist::RobinsonFoulds, methods = treess::getESSMethods(TRUE))`. Subsample each run to ≤1000 trees; if `interactive()` and total trees > threshold, emit a progress message. Extract `frechetCorrelationESS` and `medianPseudoESS`, sum across runs. `print.MkpDiagnostics` always shows a topology row (NA when not computed). Add `treess` and `TreeDist` to Suggests. |

---

## Posterior & analysis improvements

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-066 | P1 | DONE (C) | **Automatic burnin selection for MkPosterior.** The posterior object should support tunable burnin that maximizes ESS while minimizing PSRF. Add a `burnin` parameter (number of post-warmup samples to discard) with an automatic selection method that optimizes the ESS/PSRF trade-off. The `MkPosterior` accessors (samples, trees, summary, etc.) should default to returning only post-burnin samples from all runs combined. Provide a method to adjust burnin after the fact (e.g. `SetBurnin(posterior, n)` or a `burnin` argument to accessor functions). Currently, only warmup is discarded; additional burnin is common practice when chains take time to find the typical set even after warmup ends. |
| M-067 | P2 | DONE (C) | **Vignette and package documentation references via inst/REFERENCES.bib + Rdpack.** Create `inst/REFERENCES.bib` with all cited references (Sun et al. 2018, Lewis 2001, Xie et al. 2011, etc.). Set up Rdpack for `\insertRef{}` in roxygen docs. Update `hyoliths.qmd` to use `bibliography: ../inst/REFERENCES.bib` (or copy to vignettes/). Check `../TreeTools/` for a working template of this setup. Add `Rdpack` to Imports in DESCRIPTION. |
| M-068 | P2 | DONE (C) | **Rogue taxon suppression in hyoliths.qmd consensus tree.** In the tree summary section of `vignettes/hyoliths.qmd`, use `Rogue::QuickRogue()` to identify rogue taxa, then exclude them from the consensus tree. Display which taxa were identified as rogues and show the cleaned consensus. Add `Rogue` to Suggests in DESCRIPTION. |

---

## Phase 7d: Deferred extensions

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-052 | P3 | OPEN | Beta-distributed Q-matrix heterogeneity (siteMatrices). |
| M-053 | P3 | OPEN | TBR moves (if mixing diagnostics show SPR is insufficient). |
| M-054 | P3 | OPEN | HMC for branch lengths (if MH mixing insufficient). |

---

## EasyMkPrime Shiny GUI

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-077 | P2 | DONE (C) | **EasyMkPrime: auto-detect neomorphic characters from data.** `AutoDetectNeomorphic()` standalone helper (MkPrimeData.R); app.R wired to call it on data load, populate neomorphic text input, and show notification. Per-character override via text input; type breakdown shown in dataInfo. |

---

## Phase 6b: TreeSearch GUI integration

The goal is a **"Bayesian (Mk')"** mode inside EasyTrees (TreeSearch's Shiny GUI).

### Architecture: MCMC as a detached OS process; Shiny as a log viewer

Long-running MCMC (hours to days) must not be tied to the lifetime of the Shiny
session. The design is:

1. **Detached process**: "Run" writes a script to a tempfile and launches it as a
   detached `Rscript` (via `processx::process$new(supervise = FALSE)`), so the
   MCMC survives browser close, session timeout, and OS-sleep resumption.
2. **Job file**: a small RDS written beside the logs records `{ pid, logFiles,
   cancelFile, startTime, status }`, enabling the app to reconnect to a running
   job after restart.
3. **Disk-based progress**: the MCMC already writes Tracer-compatible TSV log
   files every `bufferSize` samples (M-075). The app polls these with
   `ReadMkLog()` (M-076) and renders trace plots / ESS in-place.
4. **Clean stop**: "Stop" calls `file.create(cancelFile)`. The MCMC R callback
   checks for this file at each progress interval and exits cleanly (saves
   checkpoint). This requires adding cancel-file checking to RunMkPrime (M-081).

This replaces the `callr::r_bg(supervise = TRUE)` approach used in the current
`app.R` (which dies with the session).

Tasks are sequential: M-081 (cancel support) → M-078 (module) → M-079 (refactor
app.R) → M-080 (TreeSearch hook).

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-081 | P2 | DONE (C) | **Cancel-file support in RunMkPrime.** `cancelFile` param in `MkPrimeMCMC()`; checked every 200-iter batch in both RunMkPrime and ResumeMkPrime; flushes buffers + saves checkpoint + sets `stop_reason = "cancelled"`. `.FlushAndSaveCheckpoint()` helper eliminates doCheck duplication. `MkCancelPath(jobDir)` exported. 5 tests. |
| M-082 | P2 | OPEN | **Checkpoint integration in EasyMkPrime (detached-process design).** In the new detached-process architecture (M-078/M-079), the MCMC script launched by `processx::process$new()` should write its checkpoint file into the same `logDir` as the TSV log files — no user file-system interaction needed. The job.rds (written by M-078) records the checkpoint path alongside `logFiles` and `cancelFile`, so the "Reconnect" button can resume a crashed/stopped run transparently. Specific tasks: (a) pass `checkpointFile = file.path(logDir, "checkpoint.rds")` when constructing the detached script; (b) "Reconnect" loads job.rds and checks for checkpoint existence before calling `RunMkPrime(overwrite=FALSE)`; (c) "Stop" flush + checkpoint are already guaranteed by M-081's cancel-file mechanism. If the checkpoint approach is not user-transparent, present the user with options (keep checkpoint vs. discard). *(From u.397.)* |
| M-078 | P2 | DONE (E) | **`MkBayesianUi()` / `MkBayesianServer()` — Shiny module (detached-process design).** New `R/BayesianModule.R`. `MkBayesianUi(id)`: accordion with MCMC config inputs (nRuns, nChains, heat, warmup, minEss, maxTime), `logDir` path field, neomorphic-character text field (auto-populated from `AutoDetectNeomorphic()` when dataset changes), Run / Stop / Reconnect buttons, progress area (trace plot + ESS table, polled from log files). `MkBayesianServer(id, dataset, startTree = NULL)`: `dataset` is a reactive phyDat from the host. "Run" writes a script + launches detached `processx::process$new(supervise = FALSE)`; writes job.rds. "Reconnect" loads an existing job.rds to resume monitoring a running job. Progress polling: `invalidateLater(5000)` reads log files via `ReadMkLog()` and redraws traces. "Stop" calls `file.create(cancelFile)`. Returns `list(jobFile = reactive(path\|NULL), trees = reactive(multiPhylo\|NULL), status = reactive(chr))`. Requires `processx` in Suggests. `shiny::testServer` smoke test. |
| M-079 | P2 | OPEN | **Refactor `inst/MkPrime/app.R` to use `MkBayesianServer`.** Replace the inline MCMC-launch / progress-poll / result-read logic in `app.R` with a call to `MkBayesianServer`. The standalone app retains its own data-loading and starting-tree sections. The `trees` reactive from the module feeds the existing Traces / Summary / Consensus tabs. Verifies the standalone app still works end-to-end with the new detached-process design. |
| M-080 | P2 | OPEN | **TreeSearch `mod_bayesian.R` — "Bayesian (Mk')" tab in EasyTrees.** *(TreeSearch-a repo; tracked here for coordination.)* `Suggests: MkPrime` added to TreeSearch DESCRIPTION. New `inst/Parsimony/server/mod_bayesian.R`: `bayesian_ui(id)` wraps `MkPrime::MkBayesianUi(id)` with availability check; `bayesian_server(id, r, HaveData, UpdateAllTrees, ...)` calls `MkPrime::MkBayesianServer(id, dataset = reactive(r$dataset))`. On job completion (status == "done"), reads post-burnin trees from the log directory and inserts into `r$allTrees` (displayed as "N posterior trees (unscored)"). New "Bayesian (Mk')" nav tab in EasyTrees alongside the existing parsimony tabs. |

---

## Standing Tasks

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| S-RED | dyn | OPEN | **Standing: Red-team review.** Review recent code changes for correctness bugs, edge cases, and safety issues. Focus areas: C++ MCMC hot path (likelihood, proposals, flat-buffer pruning, OPP-1–6 changes), MCMC acceptance logic (rollback correctness for each move type), R orchestration (convergence checking, streaming buffers, checkpoint/resume). File any bugs found as `u.nnn` issue files. When completed, record the focus area and outcome in the Notes column and reset to OPEN. Priority: ≥6 open tasks → P3, 3–5 → P2, <3 → P1. | Last run: — |
| S-PROF | dyn | OPEN | **Standing: Performance profiling.** Profile the compiled MCMC hot path using VTune (see `r-package-profiling` skill) or `bench::mark()` microbenchmarks. Identify the current top hotspot after OPP-1–6. Check whether `pruning_jc_flat` / `pruning_jc_acrv_flat` show further vectorisation opportunities, whether chain-swap overhead is visible at scale, or whether R↔C++ boundary crossings dominate for small datasets. File any actionable findings as new `M-nnn` tasks. When completed, record the focus and key finding in Notes and reset to OPEN. Priority: same dynamic rule as S-RED. | Last run: — |
| S-COORD | dyn | OPEN | **Standing: Coordination review.** Review `to-do.md` and `completed-tasks.md` for consistency: stale ASSIGNED statuses, tasks that are done but not archived, priorities that need adjusting. Check agent log files (`agent-*.md`) for blocked work or stale context. Scan `u.nnn` issue files and triage any that have accumulated. Verify the phase-8 OPP commit is clean and push if ready. When completed, record round number and actions taken in Notes and reset to OPEN. Priority: same dynamic rule as S-RED. | Last run: 2026-03-27 round 1 (this session). Actions: added standing-task rules + S-RED/S-PROF/S-COORD rows; triaged u.397 → M-082; fixed `completed-tasks.md` (removed `EOF 2>&1` artefact, merged duplicate M-071 rows); added `test_run*.log` / `test_output.txt` to `.gitignore`. 7 OPEN specific tasks (M-052/053/054/078/079/080/082) → standing tasks P3. |
