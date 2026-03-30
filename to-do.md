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

## Display / UI

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-151 | P2 | OPEN | **`plot()` and `print()` broken for streaming multi-run results.** S-RED round 14 finds: **(a) `plot()` crashes:** `plot.MkPosterior()` line 190 accesses `pb$per_run[[run]]$samples[, colIdx]`, but streaming per-run `samples` is NULL (line 1800 of RunMkPrime.R). Crashes with subscript error on any `plot(result)` call for streaming + nRuns > 1. Fix: fall back to combined samples for multi-run overlay, or auto-load per-run from log files. **(b) `print()` displays empty per-run count:** `nrow(x$per_run[[1]]$samples)` → `nrow(NULL)` → NULL → "Runs: 2 ( samples each)". Fix: use `x$per_run[[1]]$saved_idx` for streaming. **(c) `.PostBurninData()` crashes for streaming multi-run with burnin > 0:** `nrow(r$samples)` at burnin.R:169 where `r$samples` is NULL. Low impact (requires explicit `result$burnin <- N`). |


## MCMC infrastructure

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-148 | P2 | OPEN | **Serial multi-run mode ignores `maxRhat` stopping criterion.** *Likely resolved:* `.RunSerialRuns()` two-phase orchestrator now dispatches at line 345 of `RunMkPrime.R` when `nRuns >= 2 && maxRhat` is set. Needs formal verification: confirm tests cover cross-run R-hat convergence, then archive to completed-tasks.md. |
| M-149 | P1 | OPEN | **Checkpoint/resume: warmup/tuning state and move weights not persisted.** Seven related bugs found in S-RED rounds 12–13: **(1) CRITICAL — Tuning-phase resume crashes:** `tuningBuf` is only allocated at Warmup→Tuning transition (line 937). Resume with `r$phase == "Tuning"` leaves `tuningBuf = NULL`; first saved sample hits `nrow(NULL)` → crash. **(2) Move weights not restored:** `checkpoint$moveWeights` is saved (line 1863) but never read back in `ResumeMkPrime()`. Resume always builds fresh default weights, losing all tuning gains. **(3) Warmup state not persisted:** `logPostHistory` and `nStableConsecutive` read from `r` on init (lines 684–685) but never written back → warmup stabilisation restarts from scratch on resume. **(4) Tuning counters not persisted:** `tuningIterUsed` and `tuningRoundsDone` similarly read from `r` (694–695) but never written back → tuning budget tracking lost. **(5) `effectiveTuningBudget` lost:** Computed at Warmup→Tuning transition from `remainingIter`; on Tuning resume, defaults to `mcmc$tuningBudget` (may differ). **(6) `.RunSerialRuns()` orchestrator phase not in checkpoint:** Phase 2 checkpoints (lines 1302, 1368) don't pass `phase` or `moveWeights`. On resume, `ResumeMkPrime()` re-enters `.RunSerialRuns()` from Phase 1, re-running all runs before cross-run R-hat check. Wasteful (not incorrect). **(7) Per-run `startIters` not preserved:** All runs resume from `max(actual_iter) + 1` instead of their individual last iterations. Runs that stopped earlier lose iteration budget. **Fix plan:** (a) Write all 4 local vars back to `r` before every `.SaveCheckpoint()` call; (b) re-allocate `tuningBuf` + `tuningCandidates` on Tuning-phase resume; (c) read `checkpoint$moveWeights` in `ResumeMkPrime()` and pass through to `.RunMkPrimeSingleRun()` (new parameter); (d) persist `effectiveTuningBudget` in `r`; (e) save per-run `startIters` and orchestrator phase in `.RunSerialRuns()` checkpoints; (f) pass `moveWeights` in serial checkpoint calls. |
| M-150 | P3 | OPEN | **Checkpoint/resume: `maxTime` break doesn't save final checkpoint.** Cancel path (lines 1120–1133) flushes streaming buffer and saves checkpoint before break; `maxTime` path (lines 1112–1116) just breaks. Post-loop flush (line 1213) saves samples to log, but no checkpoint is saved for single-run. On resume, `.TruncateLogToN()` discards those flushed-but-uncheckpointed samples. Not corruption — just data loss between last `checkEvery` and `maxTime`. **Fix:** Add flush + checkpoint save in `maxTime` break path, mirroring cancel path. |
| M-131 | P2 | OPEN | **Warmup stabilisation detector: adaptive `windowSize` + validation study.** `.CheckStabilisation()` uses hardcoded `windowSize=10` (Geweke comparison windows), giving a minimum 11,500-iteration warmup floor regardless of tree size. Propose scaling `windowSize` with `nEdge` (`max(5, min(20, nEdge %/% 10))`). Requires empirical validation: run 8 benchmark datasets (20–88 tips) × 4 seeds on Hamilton with warmup disabled (200k iter), then replay the detector offline with a grid of `windowSize` and `nStableRequired` values to verify no false positives. **Briefing:** `.positai/plans/2026-03-29-m131-warmup-stabilisation-validation.md` |


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
| S-RED | dyn | OPEN | **Standing: Red-team review.** Review recent code changes for correctness bugs, edge cases, and safety issues. Focus areas: C++ MCMC hot path (likelihood, proposals, flat-buffer pruning, OPP-1–6 changes), MCMC acceptance logic (rollback correctness for each move type), R orchestration (convergence checking, streaming buffers, checkpoint/resume). File any bugs found as `u.nnn` issue files. When completed, record the focus area and outcome in the Notes column and reset to OPEN. Priority: ≥6 open tasks → P3, 3–5 → P2, <3 → P1. \| Last run: 2026-03-30 round 14 (B). Focus: **`.BuildResult()` / `MkPosterior` output construction**. Audited result assembly, print/summary/plot methods, and `.PostBurninData()` for streaming and multi-run paths. **Filed M-151 (3 bugs):** (a) `plot()` crashes for streaming multi-run — `per_run[[run]]$samples` is NULL, subscript error. (b) `print()` displays empty per-run sample count — `nrow(NULL)` → empty string in cli. (c) `.PostBurninData()` crashes for streaming multi-run with burnin > 0 — `nrow(NULL)` on per-run samples. Verified correct: streaming buffer flush in `.BuildResult()`; `pmax(..., 1L)` prevents division by zero in acceptance rates; tree trimming handles `seq_len(0)` gracefully; ESS/R-hat in summary auto-loads from log file; `MkPosterior()` constructor is clean. Also fixed M-149 #1-5 (tuning crash + state persistence, commit `e56feb5`). Prev: round 13 — `.RunSerialRuns()` audit (M-149 #6-7). Round 12 — checkpoint audit (M-149 #1–5, M-150). Round 11 — cache invalidation audit. Round 10 — M-145. Round 9 — M-144. Round 8 — M-143. Round 7 — M-142. Round 6 — test coverage paths (logHastings -Inf, prior -Inf) and MH rejection. (10) `partLogLik` cache lifecycle consistent: populated at init, cleared by Gibbs/weighted/partial-CL acceptance, correctly rebuilt through per-partition eval when non-empty. (11) Hastings ratios verified: SPR (`log(lRegraft/lMerge)` symmetric candidate counts), TBR (`+log(lSubEdge/lMergeSub)` subtree re-rooting Jacobian), Dirichlet (`logRev - logFwd` Dir density ratio), joint moves (bivariate Jacobian `log(m1) + log(m2)`). (12) Heated acceptance `β × ΔlogLik + ΔlogPrior + logHR` correct (unheated prior, unheated Hastings). (13) Prior computation correct (`Dir(1,...,1)` constant, `lgamma(n)` normalisation). NNI correctly skips prior eval. |
| S-PROF | dyn | OPEN | **Standing: Performance profiling.** Profile the compiled MCMC hot path using VTune (see `r-package-profiling` skill) or `bench::mark()` microbenchmarks. Identify the current top hotspot after OPP-1–6. Check whether `pruning_jc_flat` / `pruning_jc_acrv_flat` show further vectorisation opportunities, whether chain-swap overhead is visible at scale, or whether R↔C++ boundary crossings dominate for small datasets. File any actionable findings as new `M-nnn` tasks. When completed, record the focus and key finding in Notes and reset to OPEN. Priority: same dynamic rule as S-RED. \| Last run: 2026-03-28 round 2 (C). Focus: M-106 optimization round — Rprof + bench::mark profiling of real workload (Sun2018, 20k iter). Key finding: VTune "90% R overhead" was artifact of 0% acceptance rates; real bottleneck is C++ Felsenstein pruning (90% of wall time). Implemented 5 optimizations for 1.80× cumulative speedup. Remaining C++ time dominated by pruning traversal (diminishing returns). |
| S-COORD | dyn | OPEN | **Standing: Coordination review.** Review to-do.md and completed-tasks.md for consistency. Check agent logs and u.nnn files. When completed, record round number and actions in Notes and reset to OPEN. Priority: same dynamic rule as S-RED. \| Last run: 2026-03-30 round 6 (A). Actions: (1) Verified completed-tasks: M-140, M-141, M-142–M-145, M-124, M-135 all archived. 147 completed tasks total. (2) Updated coordination.md: added partial-CL bug fix series (M-142–M-145) and streaming/UI tasks (M-135, M-140, M-141) to recent milestones; corrected standing task priority to P1. (3) Test suite: 3977 pass, 0 fail, 20 skip, 10 warnings (up from 3960). (4) No u.nnn files, no agent logs. (5) Uncommitted: 3 stale man/ files (need `roxygen2::roxygenise()`) + `vignettes/hyoliths_posterior.rds` (generated artifact, consider adding to .gitignore). (6) Vignette PSRF reference (line 233) is correct — explanatory context ("R-hat supersedes the classical PSRF"). (7) M-080 (C) still ASSIGNED. 1 OPEN + 1 ASSIGNED specific = standing tasks P1. |


