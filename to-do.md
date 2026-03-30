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
| M-146 | P3 | OPEN | **Move weights display: silver names, colour-coded values only.** Follow-up to M-140. Current scheme colour-codes entire entries; instead, all parameter names should be silver and only the percentage values should be coloured: green (>10%), orange (5–10%), yellow/white (<5%). In `.FormatMoveWeights()` (`R/RunMkPrime.R`). |
| M-147 | P3 | OPEN | **Fix Unicode escapes in adapted-thin cli message.** `cli::cli_alert_info("... \\u2192 ...")` renders as literal `\u2192` instead of `→`. Same for `\\u2248` (`≈`). In `.AdaptThinning()` around line 1170 of `R/RunMkPrime.R`. Fix: use raw Unicode characters directly or `\u2192` (single backslash) outside the glue string. |

## MCMC infrastructure

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-148 | P2 | OPEN | **Serial multi-run mode ignores `maxRhat` stopping criterion.** In non-parallel mode with `nRuns > 1`, `.CheckConvergence()` receives `list(r)` (single-run list), so R-hat (which requires ≥ 2 runs) is never computed during individual runs. Each run stops at `minEss`/`nIter` independently; between-run R-hat is silently skipped. Options: (a) after run 1 completes, check cross-run R-hat against completed runs at each convergence check of subsequent runs; extend earlier runs if R-hat still high; (b) require `parallel = TRUE` when `nRuns > 1` and `maxRhat` is set; (c) run all runs to `minEss`, then do a post-hoc R-hat check and warn/extend. |
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
| S-RED | dyn | OPEN | **Standing: Red-team review.** Review recent code changes for correctness bugs, edge cases, and safety issues. Focus areas: C++ MCMC hot path (likelihood, proposals, flat-buffer pruning, OPP-1–6 changes), MCMC acceptance logic (rollback correctness for each move type), R orchestration (convergence checking, streaming buffers, checkpoint/resume). File any bugs found as `u.nnn` issue files. When completed, record the focus area and outcome in the Notes column and reset to OPEN. Priority: ≥6 open tasks → P3, 3–5 → P2, <3 → P1. \| Last run: 2026-03-30 round 11 (B). Focus: **systematic cache invalidation audit** — motivated by M-143/M-144/M-145 bug pattern (state changes without cache invalidation). Systematically traced every state mutation path in `mcmc.cpp` (~4150 lines). **0 bugs found.** Verified: (1) All 9 Gibbs/weighted move acceptance paths have `nodeCL.valid = false` (M-143 fix confirmed). (2) Slice sampler has `nodeCL.valid = false` at line 2829 (M-145 fix confirmed). (3) `do_move_impl` general acceptance correctly invalidates via `if (!usedPartialCL && likChanges) state->nodeCL.valid = false`. (4) NNI partial CL: TreeNav update/rollback via `update_topo_nni` correct in both accept and reject paths; `edgeToPar` indices stable across in-place NNI; preorder safety condition (`wRow > edgeRow`) correctly checked. (5) Beta_simplex partial CL: dirty path from 2 modified edge parents to root correct; `absEdgeLen` parameter (not stale `topo.edgeLen`) used for dirty recompute. (6) Dirichlet partial CL: dirty set from K edge parents to root; fullback heuristic (>80% dirty) correctly falls through to full eval with `usedPartialCL = false` → cache invalidated on acceptance. (7) Chain swap `std::swap(*states[i], *states[j])` safe — nodeCL travels with its state. (8) Checkpoint/resume: `init_mcmc_state` creates `nodeCL.valid = false` by default; no stale cache on resume. (9) State rollback correct for all early-rejection paths (logHastings -Inf, prior -Inf) and MH rejection. (10) `partLogLik` cache lifecycle consistent: populated at init, cleared by Gibbs/weighted/partial-CL acceptance, correctly rebuilt through per-partition eval when non-empty. (11) Hastings ratios verified: SPR (`log(lRegraft/lMerge)` symmetric candidate counts), TBR (`+log(lSubEdge/lMergeSub)` subtree re-rooting Jacobian), Dirichlet (`logRev - logFwd` Dir density ratio), joint moves (bivariate Jacobian `log(m1) + log(m2)`). (12) Heated acceptance `β × ΔlogLik + ΔlogPrior + logHR` correct (unheated prior, unheated Hastings). (13) Prior computation correct (`Dir(1,...,1)` constant, `lgamma(n)` normalisation). NNI correctly skips prior eval. |
| S-PROF | dyn | OPEN | **Standing: Performance profiling.** Profile the compiled MCMC hot path using VTune (see `r-package-profiling` skill) or `bench::mark()` microbenchmarks. Identify the current top hotspot after OPP-1–6. Check whether `pruning_jc_flat` / `pruning_jc_acrv_flat` show further vectorisation opportunities, whether chain-swap overhead is visible at scale, or whether R↔C++ boundary crossings dominate for small datasets. File any actionable findings as new `M-nnn` tasks. When completed, record the focus and key finding in Notes and reset to OPEN. Priority: same dynamic rule as S-RED. \| Last run: 2026-03-28 round 2 (C). Focus: M-106 optimization round — Rprof + bench::mark profiling of real workload (Sun2018, 20k iter). Key finding: VTune "90% R overhead" was artifact of 0% acceptance rates; real bottleneck is C++ Felsenstein pruning (90% of wall time). Implemented 5 optimizations for 1.80× cumulative speedup. Remaining C++ time dominated by pruning traversal (diminishing returns). |
| S-COORD | dyn | OPEN | **Standing: Coordination review.** Review to-do.md and completed-tasks.md for consistency. Check agent logs and u.nnn files. When completed, record round number and actions in Notes and reset to OPEN. Priority: same dynamic rule as S-RED. \| Last run: 2026-03-30 round 6 (A). Actions: (1) Verified completed-tasks: M-140, M-141, M-142–M-145, M-124, M-135 all archived. 147 completed tasks total. (2) Updated coordination.md: added partial-CL bug fix series (M-142–M-145) and streaming/UI tasks (M-135, M-140, M-141) to recent milestones; corrected standing task priority to P1. (3) Test suite: 3977 pass, 0 fail, 20 skip, 10 warnings (up from 3960). (4) No u.nnn files, no agent logs. (5) Uncommitted: 3 stale man/ files (need `roxygen2::roxygenise()`) + `vignettes/hyoliths_posterior.rds` (generated artifact, consider adding to .gitignore). (6) Vignette PSRF reference (line 233) is correct — explanatory context ("R-hat supersedes the classical PSRF"). (7) M-080 (C) still ASSIGNED. 1 OPEN + 1 ASSIGNED specific = standing tasks P1. |


