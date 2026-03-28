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
| M-106 | P1 | OPEN | **Adaptive batch sizing (both directions).** Currently fixed at 200 iterations per batch. VTune (M-103) shows R↔C++ round-trip overhead is **90%+ of CPU** for the default move schedule — the 200-iteration batch causes 5 000 round-trips per million iterations. Increase default batch size to ~2 000–5 000 for cheap configurations (standard moves). For expensive moves (block Gibbs, weighted SPR), decrease to ~50 for UI responsiveness. Auto-tune based on observed per-iteration wall time from the first batch. Expected: 2–5× wall-clock speedup for standard runs. |
| M-108 | P2 | OPEN | **Cache/patch postorder after NNI.** VTune (M-103) shows `TreeTools::postorder_order` is the #1 C++ hotspot (34% of MkPrime.dll time, 3.2% of total). Called every topology proposal to maintain the preorder invariant. NNI changes only 4 edges; the traversal order can be patched rather than recomputed from scratch. |
| M-104 | P3 | OPEN | **Reduce `clone()` overhead in Gibbs SPR.** Each of the ~100 candidates in `gibbs_spr_impl` allocates 3 new Rcpp vectors. VTune shows clone overhead is only 0.2% in the standard move schedule (M-103) and <0.1% even in Gibbs-specific profiling (mkp-gibbs VTune 2026-03-28). Confirmed low priority across both profiles. |
| M-105 | P2 | OPEN | **Partial likelihood reuse for Gibbs/weighted SPR.** Precomputing the "pruned tree" CL once and only updating the affected path per candidate could reduce per-candidate cost from O(N×C) to O(D×C). Standard-move VTune (M-103) shows pruning is hidden by R↔C++ overhead, but **Gibbs-specific VTune (mkp-gibbs, 2026-03-28) shows `pruning_jc_acrv` is 42% of total CPU** — the dominant hotspot. Each Gibbs SPR candidate triggers a full-tree pruning. Partial reuse is the highest-value Gibbs optimization after OPP-1. |

## Misc / UI improvements

| ID | Priority | Status | Description |
|----|----------|--------|-------------|

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
| S-PROF | dyn | OPEN | **Standing: Performance profiling.** Profile the compiled MCMC hot path using VTune (see `r-package-profiling` skill) or `bench::mark()` microbenchmarks. Identify the current top hotspot after OPP-1–6. Check whether `pruning_jc_flat` / `pruning_jc_acrv_flat` show further vectorisation opportunities, whether chain-swap overhead is visible at scale, or whether R↔C++ boundary crossings dominate for small datasets. File any actionable findings as new `M-nnn` tasks. When completed, record the focus and key finding in Notes and reset to OPEN. Priority: same dynamic rule as S-RED. \| Last run: 2026-03-28 (C). Focus: VTune hotspot profiling (M-103). Key finding: R↔C++ boundary overhead dominates (90%+ CPU); batch size 200 is the bottleneck. C++ pruning functions absent from profile (M-063 effective). Filed M-106 (P1, adaptive batch size) and M-108 (P2, cache postorder). |
| S-COORD | dyn | OPEN | **Standing: Coordination review.** Review `to-do.md` and `completed-tasks.md` for consistency: stale ASSIGNED statuses, tasks that are done but not archived, priorities that need adjusting. Check agent log files (`agent-*.md`) for blocked work or stale context. Scan `u.nnn` issue files and triage any that have accumulated. When completed, record round number and actions taken in Notes and reset to OPEN. Priority: same dynamic rule as S-RED. \| Last run: 2026-03-28 round 3. Actions: (1) pruned completed Phase 10 section from to-do.md (M-094/095/096 all DONE and archived); (2) updated coordination.md project state: Phase 10 "planned" → "complete", added logseries branch note, bumped date; (3) agent logs: all four stale (agent-a: Phase 2 era, agent-b: M-062, agent-c: IDLE with logseries branch ready to merge, agent-e: ACTIVE but all tasks done); (4) no u.nnn files; (5) logseries k' prior on separate branch has no task ID — noted in coordination.md for awareness. 14 OPEN specific tasks → standing tasks remain P3. |
