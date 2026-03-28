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
| M-083 | P1 | OPEN | **Extract `compute_full_loglik()` helper.** All five new Gibbs/Weighted moves need to evaluate the full posterior likelihood from a C++ context — multiple times per proposal. Currently the only entry point is the R-callable `MkpLogLikelihood()` and the hot-loop-internal call in `mcmc_engine.cpp`. Expose a `compute_full_loglik(state, data)` C++ function (and optionally `compute_full_loglik_with_prior()`) callable from within proposal code. This is the P1 blocker for M-085–M-089. Signature: takes the flat-buffer state struct + partition data (pointers already held in `McmcState`) and returns `double` log-likelihood (unheated). No R boundary crossing; pure C++. |
| M-084 | P1 | OPEN | **Subtree swap topology operation + `get_valid_swap_partners()`.** NNI and SPR are already in `tree_moves.cpp`. Add: (a) `swap_subtrees(parent, child, relBrLengths, nodeA, nodeB)` — detaches two non-nested, non-sibling subtrees and exchanges their attachment points, then calls `postorder_reorder()`; returns updated parent/child/relBrLengths. (b) `get_valid_swap_partners(parent, child, pruneNode)` — returns all nodes that are non-nested and non-sibling relative to `pruneNode` (valid swap candidates). (c) Unit tests: round-trip swap returns original tree; swap on 5-taxon tree matches hand-computed result; `get_valid_swap_partners` excludes nested/sibling nodes correctly. This operation is used by both M-086 (GibbsSubtreeSwap) and M-089 (WeightedSubtreeSwap). |

### 9b: Gibbs moves

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-085 | P2 | OPEN | **GibbsSPR.** Enumerate all valid SPR reattachment candidates (non-nested, non-sibling edges), compute the tempered posterior weight `exp(β × logLik + logPrior)` for each, normalise to a probability vector, sample one candidate proportionally, and apply. Return `logHastings = log(1/w_chosen) - log(1/w_original_position)` so the MH ratio is correct. Cost: O(N) likelihood evals where N = number of candidates. Tempered: `exp(β × logLik)` (β from chain's heat; logPrior unheated). Depends on M-083. |
| M-086 | P2 | OPEN | **GibbsSubtreeSwap.** Same pattern as M-085 but uses the subtree-swap operation from M-084. Enumerate valid swap partners, weight by posterior, sample ∝ weight, apply swap, return correct Hastings ratio. Depends on M-083, M-084. |

### 9c: Weighted moves

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-087 | P2 | OPEN | **WeightedBranchLengthScale.** Discretize the branch-fraction proposal space into B bins (default B = 10) using Beta(0.25, 0.25) quantiles as break points. For each bin midpoint, evaluate the posterior, weight by `exp(β × logLik) × proposal_density`. Sample a bin from the normalised weights (CDF inverse), then draw the actual branch fraction from a Beta centred on that bin's midpoint. Compute the piecewise Hastings ratio. Cost: O(B) likelihood evals. Depends on M-083. |
| M-088 | P2 | OPEN | **WeightedSPR.** Gibbs SPR (M-085) extended to also integrate over branch fractions at each candidate position: for each candidate topology, marginalise over B branch-fraction bins, producing a marginal weight. Sample topology from marginal weights, then sample branch fraction from conditional distribution given chosen topology. Compute combined Hastings ratio. Cost: O(N × B) likelihood evals. Depends on M-083, M-085. |
| M-089 | P3 | OPEN | **WeightedSubtreeSwap.** Same as M-088 but uses subtree-swap topology operation from M-084. Cost: O(N × B) evals. Depends on M-083, M-084, M-086, M-088 (reuse integration logic). |

### 9d: Integration, scheduling & validation

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-090 | P2 | DONE (E) | **Wire new moves into dispatch + `MkPrimeMCMC()` parameters.** Add the five new moves to the `do_move_impl()` dispatch table in `mcmc_engine.cpp`. Add user-facing parameters to `MkPrimeMCMC()`: `gibbsSpr` (logical, default TRUE), `gibbsSubtreeSwap` (logical, default TRUE), `weightedBranchScale` (logical, default FALSE — expensive), `weightedSpr` (logical, default FALSE), `weightedSubtreeSwap` (logical, default FALSE), `nBranchBins` (integer, default 10L). Update `MkPrimeModel()` docs and `MkPrimeMCMC()` docs. Add move-presence guards so disabled moves are excluded from the move pool entirely. Depends on M-085, M-086, M-087, M-088, M-089. |
| M-092 | P2 | OPEN | **Adaptive move scheduler.** During warmup, track per-move `(n_accepted, n_proposed, wall_time_ns)` in the C++ state struct. Every adaptation window (200 iters), recompute move weights: `score[m] = accept_rate[m] / mean_cost_s[m]`; update weights via softmax with temperature that anneals from 2.0 → 0.5 over warmup, floored at `w_min = 0.05`. Freeze weights at warmup end (detailed balance preserved post-warmup). User override: `moveWeights = list(nni = 0.3, spr = 0.3, ...)` in `MkPrimeMCMC()` pins specific move weights and excludes them from adaptation. Log final scheduled weights to the run header in the TSV log file. Depends on M-090. |
| M-091 | P3 | OPEN | **Mixing validation: Gibbs/Weighted vs standard moves.** On the Vinther hyoliths dataset, compare ESS/wall-time between: (a) baseline (NNI + SPR only), (b) + GibbsSPR + GibbsSubtreeSwap, (c) + WeightedSPR. Use `MkpWatchLog()` output and `ConvergenceDiagnostics()` to compare. Report topology ESS/hour and scalar-parameter ESS/hour. Document findings in a comment in `coordination.md`. Depends on M-090, M-092. |

---

## Phase 7d: Deferred extensions

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-052 | P3 | OPEN | Beta-distributed Q-matrix heterogeneity (siteMatrices). |
| M-053 | P3 | OPEN | TBR moves (if mixing diagnostics show SPR is insufficient). |
| M-054 | P3 | OPEN | HMC for branch lengths (if MH mixing insufficient). |

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

M-081 → M-078 → M-079 → M-082/M-093 all complete on main. Remaining: M-080
(TreeSearch-a repo).

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-080 | P2 | OPEN | **TreeSearch `mod_bayesian.R` — "Bayesian (Mk')" tab in EasyTrees.** *(TreeSearch-a repo; tracked here for coordination.)* `Suggests: MkPrime` added to TreeSearch DESCRIPTION. New `inst/Parsimony/server/mod_bayesian.R`: `bayesian_ui(id)` wraps `MkPrime::MkBayesianUi(id)` with availability check; `bayesian_server(id, r, HaveData, UpdateAllTrees, ...)` calls `MkPrime::MkBayesianServer(id, dataset = reactive(r$dataset))`. On job completion (status == "done"), reads post-burnin trees from the log directory and inserts into `r$allTrees` (displayed as "N posterior trees (unscored)"). New "Bayesian (Mk')" nav tab in EasyTrees alongside the existing parsimony tabs. |

---

## Standing Tasks

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| S-RED | dyn | OPEN | **Standing: Red-team review.** Review recent code changes for correctness bugs, edge cases, and safety issues. Focus areas: C++ MCMC hot path (likelihood, proposals, flat-buffer pruning, OPP-1–6 changes), MCMC acceptance logic (rollback correctness for each move type), R orchestration (convergence checking, streaming buffers, checkpoint/resume). File any bugs found as `u.nnn` issue files. When completed, record the focus area and outcome in the Notes column and reset to OPEN. Priority: ≥6 open tasks → P3, 3–5 → P2, <3 → P1. \| Last run: 2026-03-27 round 3 (E). Focus: Phase 10 parallel mode (M-094–M-096). Found 3 bugs, all fixed inline: (1) .BuildResult() used mcmc$logFile not logFilePaths — crashes when parallel auto-assigns log; (2) treeFile shared across workers — concurrent write corruption; (3) checkpoint save in sequential else-branch only — parallel+checkpointFile silently wrote nothing, blocking resume. All 3 have regression tests. |
| S-PROF | dyn | OPEN | **Standing: Performance profiling.** Profile the compiled MCMC hot path using VTune (see `r-package-profiling` skill) or `bench::mark()` microbenchmarks. Identify the current top hotspot after OPP-1–6. Check whether `pruning_jc_flat` / `pruning_jc_acrv_flat` show further vectorisation opportunities, whether chain-swap overhead is visible at scale, or whether R↔C++ boundary crossings dominate for small datasets. File any actionable findings as new `M-nnn` tasks. When completed, record the focus and key finding in Notes and reset to OPEN. Priority: same dynamic rule as S-RED. \| Last run: — |
| S-COORD | dyn | OPEN | **Standing: Coordination review.** Review `to-do.md` and `completed-tasks.md` for consistency: stale ASSIGNED statuses, tasks that are done but not archived, priorities that need adjusting. Check agent log files (`agent-*.md`) for blocked work or stale context. Scan `u.nnn` issue files and triage any that have accumulated. When completed, record round number and actions taken in Notes and reset to OPEN. Priority: same dynamic rule as S-RED. \| Last run: 2026-03-28 round 3. Actions: (1) pruned completed Phase 10 section from to-do.md (M-094/095/096 all DONE and archived); (2) updated coordination.md project state: Phase 10 "planned" → "complete", added logseries branch note, bumped date; (3) agent logs: all four stale (agent-a: Phase 2 era, agent-b: M-062, agent-c: IDLE with logseries branch ready to merge, agent-e: ACTIVE but all tasks done); (4) no u.nnn files; (5) logseries k' prior on separate branch has no task ID — noted in coordination.md for awareness. 14 OPEN specific tasks → standing tasks remain P3. |
