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

---

## Bugs

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| | | | *(No open bugs)* |

## Optimization roadmap (Gibbs/weighted moves)

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-109 | P1 | DONE | **Eliminate per-candidate clone() in Gibbs/weighted moves.** Added `preorder_into()` (lightweight DFS preorder into pre-allocated buffers). Refactored: `gibbs_subtree_swap_impl` (5→0 clones/cand), `gibbs_spr_impl_full` (3→0), `gibbs_spr_impl` step 10 (3→0), `weighted_spr_impl` (2+1/bin→0), `weighted_subtree_swap_impl` (1+1/bin→0). All use in-place save/restore + shared output buffers. |
| M-111 | P2 | OPEN | **Partial CL reuse for Gibbs subtree swap.** Port M-105's `caching_downpass`/`evaluate_candidate` pattern from `gibbs_spr_impl` to `gibbs_subtree_swap_impl`. Swap still does full per-candidate evaluation via `compute_full_loglik_at` (M-109 eliminated the clone overhead but not the O(N×C) per-candidate cost). Benchmark (2026-03-28, Sun2018 3k iter): SPR-only with M-105 = 670 iter/s; swap-only = 367 iter/s; swap is now the dominant Gibbs bottleneck. |
| M-105 | P2 | DONE | **Partial likelihood reuse for Gibbs SPR.** Merged from `feature/gibbs-weighted-moves`. Benchmark confirms 3.9× speedup for Gibbs SPR (670 vs 180 iter/s). |
| M-108 | P4 | OPEN | **~~Cache/patch postorder after NNI.~~** Superseded for NNI by OPP-6b (in-place NNI eliminates `preorder_weighted_impl` entirely). Remaining value only for SPR/TBR which still call `preorder_weighted_impl`. Priority demoted; SPR/TBR postorder cost is small relative to their per-candidate pruning. |

## Misc / UI improvements

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-110 | P3 | OPEN | **Evaluate starting tree strategy.** NJ trees can produce negative branch lengths (clamped to 1e-8 with a warning). Consider (Goloboff-)Wagner trees or uniform branch lengths as alternatives. Benchmark convergence speed from different starting points. (from u.465) |

## Phase 7d: Deferred extensions

| ID | Priority | Status | Description |
|----|----------|--------|-------------|


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
| M-080 | P2 | ASSIGNED (C) | **TreeSearch `mod_bayesian.R` — "Bayesian (Mk')" tab in EasyTrees.** *(TreeSearch-a repo; tracked here for coordination.)* `Suggests: MkPrime` added to TreeSearch DESCRIPTION. New `inst/Parsimony/server/mod_bayesian.R`: `bayesian_ui(id)` wraps `MkPrime::MkBayesianUi(id)` with availability check; `bayesian_server(id, r, HaveData, UpdateAllTrees, ...)` calls `MkPrime::MkBayesianServer(id, dataset = reactive(r$dataset))`. On job completion (status == "done"), reads post-burnin trees from the log directory and inserts into `r$allTrees` (displayed as "N posterior trees (unscored)"). New "Bayesian (Mk')" nav tab in EasyTrees alongside the existing parsimony tabs. |

---

## Standing Tasks

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| S-RED | dyn | OPEN | **Standing: Red-team review.** Review recent code changes for correctness bugs, edge cases, and safety issues. Focus areas: C++ MCMC hot path (likelihood, proposals, flat-buffer pruning, OPP-1–6 changes), MCMC acceptance logic (rollback correctness for each move type), R orchestration (convergence checking, streaming buffers, checkpoint/resume). File any bugs found as `u.nnn` issue files. When completed, record the focus area and outcome in the Notes column and reset to OPEN. Priority: ≥6 open tasks → P3, 3–5 → P2, <3 → P1. \| Last run: 2026-03-27 round 3 (E). Focus: Phase 10 parallel mode (M-094–M-096). Found 3 bugs, all fixed inline: (1) .BuildResult() used mcmc$logFile not logFilePaths — crashes when parallel auto-assigns log; (2) treeFile shared across workers — concurrent write corruption; (3) checkpoint save in sequential else-branch only — parallel+checkpointFile silently wrote nothing, blocking resume. All 3 have regression tests. |
| S-PROF | dyn | OPEN | **Standing: Performance profiling.** Profile the compiled MCMC hot path using VTune (see `r-package-profiling` skill) or `bench::mark()` microbenchmarks. Identify the current top hotspot after OPP-1–6. Check whether `pruning_jc_flat` / `pruning_jc_acrv_flat` show further vectorisation opportunities, whether chain-swap overhead is visible at scale, or whether R↔C++ boundary crossings dominate for small datasets. File any actionable findings as new `M-nnn` tasks. When completed, record the focus and key finding in Notes and reset to OPEN. Priority: same dynamic rule as S-RED. \| Last run: 2026-03-28 round 2 (C). Focus: M-106 optimization round — Rprof + bench::mark profiling of real workload (Sun2018, 20k iter). Key finding: VTune "90% R overhead" was artifact of 0% acceptance rates; real bottleneck is C++ Felsenstein pruning (90% of wall time). Implemented 5 optimizations for 1.80× cumulative speedup. Remaining C++ time dominated by pruning traversal (diminishing returns). |
| S-COORD | dyn | OPEN | **Standing: Coordination review.** Review `to-do.md` and `completed-tasks.md` for consistency: stale ASSIGNED statuses, tasks that are done but not archived, priorities that need adjusting. Check agent log files (`agent-*.md`) for blocked work or stale context. Scan `u.nnn` issue files and triage any that have accumulated. When completed, record round number and actions taken in Notes and reset to OPEN. Priority: same dynamic rule as S-RED. \| Last run: 2026-03-28 round 3. Actions: (1) pruned completed Phase 10 section from to-do.md (M-094/095/096 all DONE and archived); (2) updated coordination.md project state: Phase 10 "planned" → "complete", added logseries branch note, bumped date; (3) agent logs: all four stale (agent-a: Phase 2 era, agent-b: M-062, agent-c: IDLE with logseries branch ready to merge, agent-e: ACTIVE but all tasks done); (4) no u.nnn files; (5) logseries k' prior on separate branch has no task ID — noted in coordination.md for awareness. 14 OPEN specific tasks → standing tasks remain P3. |
