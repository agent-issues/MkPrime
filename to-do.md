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
| M-122 | P2 | ASSIGNED (D) | **Fix `test-bayesian-module.R` failures (7 tests).** `status()` not found in `testServer()` context — likely a Shiny module export or namespace issue introduced when the module was refactored. All 7 tests fail with the same `could not find function "status"` error. |
| M-123 | P3 | OPEN | **Verify `test-m092-adaptive-scheduler.R:360` passes in full suite.** Appeared to fail with `nEdge` missing in a full-suite run but passes in isolation — likely a test-ordering artifact from stale correction tests. Confirm in next full `R CMD check` or full-suite run; close if clean. |

## Optimization roadmap (Gibbs/weighted moves)

| ID | Priority | Status | Description |
|----|----------|--------|-------------|

| M-121 | P3 | DONE | **Node-level CL dirty flags for standard NNI/SPR.** Extend the partial CL framework (used by Gibbs moves, M-105/M-111) to standard NNI and SPR proposals. NNI invalidates CLs only along the path from modified nodes to root — O(depth) instead of O(n_nodes). Requires careful engineering of CL rollback on rejection. See `research-mcmc-optimizations.md` §4. |
| M-125 | P2 | OPEN | **Block branch-length proposal for better mixing.** The current one-at-a-time beta_simplex on individual relative branch lengths is the main mixing bottleneck (log_likelihood ESS ~15–20 on Vinther2008 fixed-topology). Consider a block Dirichlet or Hamiltonian-within-Gibbs proposal that updates multiple edges jointly. |



## Misc / UI improvements

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-124 | P2 | OPEN | **Submit RevBayes relabelling correction patch.** The same inverted relabelling formula exists in `PhyloCTMCSiteHomogeneousMkPrime.h`. Proof and patch instructions in `relabelling-correction-proof.md`. |

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
| S-RED | dyn | OPEN | **Standing: Red-team review.** Review recent code changes for correctness bugs, edge cases, and safety issues. Focus areas: C++ MCMC hot path (likelihood, proposals, flat-buffer pruning, OPP-1–6 changes), MCMC acceptance logic (rollback correctness for each move type), R orchestration (convergence checking, streaming buffers, checkpoint/resume). File any bugs found as `u.nnn` issue files. When completed, record the focus area and outcome in the Notes column and reset to OPEN. Priority: ≥6 open tasks → P3, 3–5 → P2, <3 → P1. \| Last run: 2026-03-29 round 5. Focus: relabelling correction fix + M-120 (2D joint Bactrian). Found 1 bug: u.123 — NaN correlation propagation in `.EstimateJointRhos()` (constant samples → `cor()` NaN → silently rejected joint moves). Fixed with `is.finite()` guard + `suppressWarnings()`. Also: dead variable cleanup in `corrections.cpp`; set `joint2d = FALSE` in stepping-stone path. Relabelling correction math verified (falling factorial P(k',kObs) correct). M-120 C++ kernel symmetry verified analytically; Hastings ratio correct; rollback correct. Design note: rho not estimated during warmup (no samples saved to R during warmup batches) — joint moves use ρ=0 until Tuning phase. |
| S-PROF | dyn | OPEN | **Standing: Performance profiling.** Profile the compiled MCMC hot path using VTune (see `r-package-profiling` skill) or `bench::mark()` microbenchmarks. Identify the current top hotspot after OPP-1–6. Check whether `pruning_jc_flat` / `pruning_jc_acrv_flat` show further vectorisation opportunities, whether chain-swap overhead is visible at scale, or whether R↔C++ boundary crossings dominate for small datasets. File any actionable findings as new `M-nnn` tasks. When completed, record the focus and key finding in Notes and reset to OPEN. Priority: same dynamic rule as S-RED. \| Last run: 2026-03-28 round 2 (C). Focus: M-106 optimization round — Rprof + bench::mark profiling of real workload (Sun2018, 20k iter). Key finding: VTune "90% R overhead" was artifact of 0% acceptance rates; real bottleneck is C++ Felsenstein pruning (90% of wall time). Implemented 5 optimizations for 1.80× cumulative speedup. Remaining C++ time dominated by pruning traversal (diminishing returns). |
| S-COORD | dyn | OPEN | **Standing: Coordination review.** Review `to-do.md` and `completed-tasks.md` for consistency: stale ASSIGNED statuses, tasks that are done but not archived, priorities that need adjusting. Check agent log files (`agent-*.md`) for blocked work or stale context. Scan `u.nnn` issue files and triage any that have accumulated. When completed, record round number and actions taken in Notes and reset to OPEN. Priority: same dynamic rule as S-RED. \| Last run: 2026-03-28 round 3. Actions: (1) pruned completed Phase 10 section from to-do.md (M-094/095/096 all DONE and archived); (2) updated coordination.md project state: Phase 10 "planned" → "complete", added logseries branch note, bumped date; (3) agent logs: all four stale (agent-a: Phase 2 era, agent-b: M-062, agent-c: IDLE with logseries branch ready to merge, agent-e: ACTIVE but all tasks done); (4) no u.nnn files; (5) logseries k' prior on separate branch has no task ID — noted in coordination.md for awareness. 14 OPEN specific tasks → standing tasks remain P3. |
