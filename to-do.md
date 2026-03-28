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
| M-100 | P2 | OPEN | **Reject `qHeterogeneity = TRUE` + `coding = "informative"` at model validation time.** `het_singleton_site_prob()` is a stub returning 0, so combining Het with informative coding silently produces wrong ascertainment corrections (and therefore wrong likelihoods). Add a validation check in `MkPrimeModel()` (or at the start of `RunMkPrime()`) that errors with a clear message. Remove the guard when Phase 7's F81 singleton correction is implemented. |

## Optimization roadmap (Gibbs/weighted moves)

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-102 | P1 | ASSIGNED (E) | **Production build comparison (-O2).** The 78× slowdown for block Gibbs includes debug overhead (-O0). Re-run the same Sun 2018 comparison (configs a/b/c/d) at -O2 to get realistic cost ratios before investing in algorithmic work. Low effort, high information value. |
| M-103 | P1 | OPEN | **Profile the C++ hot path (S-PROF).** Use `bench::mark()` or VTune to identify whether time is dominated by Felsenstein pruning, `preorder_weighted_impl`, `clone()` allocations, or R distribution functions (`R::qbeta`, `R::rbeta`, `R::dbeta`). Determines which optimization yields the most. |
| M-104 | P2 | OPEN | **Reduce `clone()` overhead in Gibbs SPR.** Each of the ~100 candidates in `gibbs_spr_impl` allocates 3 new Rcpp vectors (`clone(parent)`, `clone(child)`, `clone(absLen)`) = 300 heap allocations per call. A single reusable working buffer would eliminate this. Blocked on M-103 (profile first to confirm this is material). |
| M-105 | P2 | OPEN | **Partial likelihood reuse for Gibbs/weighted SPR.** When pruning subtree v and evaluating regraft positions, the CLs for all nodes not on the regraft path are unchanged. Precomputing the "pruned tree" CL once and only updating the affected path per candidate could reduce per-candidate cost from O(N×C) to O(D×C) where D is the path depth (typically O(log N)). High effort, potentially large payoff. Blocked on M-103. |
| M-106 | P3 | OPEN | **Adaptive batch sizing.** Currently fixed at 200 iterations per batch. For expensive move configurations (block Gibbs, weighted SPR), a smaller batch (e.g. 50) would improve UI responsiveness. Could auto-tune based on observed per-iteration wall time from the first batch. |

## Misc / UI improvements

| ID | Priority | Status | Description |
|----|----------|--------|-------------|

## Phase 7d: Deferred extensions

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-053 | P3 | DONE | TBR moves — merged. 5.9× ESS/s gain at 180 taxa; default off (`tbr = FALSE`). |

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
| S-PROF | dyn | OPEN | **Standing: Performance profiling.** Profile the compiled MCMC hot path using VTune (see `r-package-profiling` skill) or `bench::mark()` microbenchmarks. Identify the current top hotspot after OPP-1–6. Check whether `pruning_jc_flat` / `pruning_jc_acrv_flat` show further vectorisation opportunities, whether chain-swap overhead is visible at scale, or whether R↔C++ boundary crossings dominate for small datasets. File any actionable findings as new `M-nnn` tasks. When completed, record the focus and key finding in Notes and reset to OPEN. Priority: same dynamic rule as S-RED. \| Last run: — |
| S-COORD | dyn | OPEN | **Standing: Coordination review.** Review `to-do.md` and `completed-tasks.md` for consistency: stale ASSIGNED statuses, tasks that are done but not archived, priorities that need adjusting. Check agent log files (`agent-*.md`) for blocked work or stale context. Scan `u.nnn` issue files and triage any that have accumulated. When completed, record round number and actions taken in Notes and reset to OPEN. Priority: same dynamic rule as S-RED. \| Last run: 2026-03-28 round 3. Actions: (1) pruned completed Phase 10 section from to-do.md (M-094/095/096 all DONE and archived); (2) updated coordination.md project state: Phase 10 "planned" → "complete", added logseries branch note, bumped date; (3) agent logs: all four stale (agent-a: Phase 2 era, agent-b: M-062, agent-c: IDLE with logseries branch ready to merge, agent-e: ACTIVE but all tasks done); (4) no u.nnn files; (5) logseries k' prior on separate branch has no task ID — noted in coordination.md for awareness. 14 OPEN specific tasks → standing tasks remain P3. |
