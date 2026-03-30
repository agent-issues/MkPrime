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




## Optimization roadmap (Gibbs/weighted moves)

| ID | Priority | Status | Description |
|----|----------|--------|-------------|






## MCMC infrastructure

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-134 | P2 | OPEN | **`treeThin`: differential thinning for tree samples.** New `MkPrimeMCMC(treeThin=)` parameter to store tree topology samples less frequently than scalar parameters. Must be a multiple of `thin`. Motivation: tree samples are ~10× more expensive to store and have higher autocorrelation; 500–1000 posterior trees suffice for CID comparison. *(from u.001)* |
| M-135 | P3 | OPEN | **`thin = "auto"`: adapt thinning to observed autocorrelation.** Estimate ACT from early samples and set `thin ≈ ACT * log(2)`. Benefits from M-134 being done first so tree and scalar thinning can be tuned independently. *(from u.002)* |
| M-131 | P2 | OPEN | **Warmup stabilisation detector: adaptive `windowSize` + validation study.** `.CheckStabilisation()` uses hardcoded `windowSize=10` (Geweke comparison windows), giving a minimum 11,500-iteration warmup floor regardless of tree size. Propose scaling `windowSize` with `nEdge` (`max(5, min(20, nEdge %/% 10))`). Requires empirical validation: run 8 benchmark datasets (20–88 tips) × 4 seeds on Hamilton with warmup disabled (200k iter), then replay the detector offline with a grid of `windowSize` and `nStableRequired` values to verify no false positives. **Briefing:** `.positai/plans/2026-03-29-m131-warmup-stabilisation-validation.md` |

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
| S-RED | dyn | OPEN | **Standing: Red-team review.** Review recent code changes for correctness bugs, edge cases, and safety issues. Focus areas: C++ MCMC hot path (likelihood, proposals, flat-buffer pruning, OPP-1–6 changes), MCMC acceptance logic (rollback correctness for each move type), R orchestration (convergence checking, streaming buffers, checkpoint/resume). File any bugs found as `u.nnn` issue files. When completed, record the focus area and outcome in the Notes column and reset to OPEN. Priority: ≥6 open tasks → P3, 3–5 → P2, <3 → P1. \| Last run: 2026-03-29 round 6. Focus: M-125 (Block Dirichlet) + M-126 (warmup rho). M-125 C++: Hastings ratio math verified (symmetric +1 shift, Jacobian cancels with preserved selSum, floor 0.01 never triggered). Rollback correct at all 3 rejection points (inf Hastings, inf prior, MH reject) — full-vector snapshot restore. `dirEdges` properly flows into `find_dirty_dirichlet()` for partial CL eval; 80% coverage fallback to full eval is sound; cache invalidation correct on acceptance. M-125 R: inverted adaptation (divide by adj) correct for concentration param; weight formula reasonable. M-126: state extraction moved before phase logic; `.AccumulateRhoSnapshot()` accumulates from cold-chain state during warmup; dual estimation (warmup snapshots + tuning buffer). Rho flows to C++ via `.BuildJointRhoMatrix()` rebuilt each batch. 0 bugs found. Also closed M-123 (test-m092 full-suite pass verified, 0 failures). |
| S-PROF | dyn | OPEN | **Standing: Performance profiling.** Profile the compiled MCMC hot path using VTune (see `r-package-profiling` skill) or `bench::mark()` microbenchmarks. Identify the current top hotspot after OPP-1–6. Check whether `pruning_jc_flat` / `pruning_jc_acrv_flat` show further vectorisation opportunities, whether chain-swap overhead is visible at scale, or whether R↔C++ boundary crossings dominate for small datasets. File any actionable findings as new `M-nnn` tasks. When completed, record the focus and key finding in Notes and reset to OPEN. Priority: same dynamic rule as S-RED. \| Last run: 2026-03-28 round 2 (C). Focus: M-106 optimization round — Rprof + bench::mark profiling of real workload (Sun2018, 20k iter). Key finding: VTune "90% R overhead" was artifact of 0% acceptance rates; real bottleneck is C++ Felsenstein pruning (90% of wall time). Implemented 5 optimizations for 1.80× cumulative speedup. Remaining C++ time dominated by pruning traversal (diminishing returns). |
| S-COORD | dyn | OPEN | **Standing: Coordination review.** Review to-do.md and completed-tasks.md for consistency. Check agent logs and u.nnn files. When completed, record round number and actions in Notes and reset to OPEN. Priority: same dynamic rule as S-RED. \| Last run: 2026-03-29 round 4. Actions: (1) deleted M-121/M-122 (DONE) from to-do; (2) archived M-122 in completed-tasks; (3) removed stray EOF line; (4) coordination.md Phase 7 deferred items marked complete; (5) agent logs all stale; (6) no u.nnn; (7) M-080 (C) still ASSIGNED. 4 OPEN + 1 ASSIGNED = standing tasks P2. |
