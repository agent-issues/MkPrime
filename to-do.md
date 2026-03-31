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



## MCMC infrastructure

| ID | Priority | Status | Description |
|----|----------|--------|-------------|


| M-155 | P1 | OPEN | **Gibbs kPrime sweep: batch-by-k' optimization.** S-PROF round 3 found the Gibbs sweep is **10× slower** than equivalent int_walk coverage (340 ms vs 33 ms for 99 trans chars on Sun2018). Root cause: `single_char_loglik_jc()` does a **per-character, per-candidate-k' tree traversal** with per-call heap allocation, while int_walk uses batched flat-buffer partition pruning. Inner loop = 99.5% of sweep time. nCat scales linearly (70 ms at nCat=1, 340 ms at nCat=6). **Fix:** restructure the Gibbs inner loop to batch characters by candidate k' value and evaluate them using the existing flat-buffer partition pruning infrastructure. For each distinct candidate k', gather all characters with that candidate into a temporary partition-style buffer and call the partition-level pruning once. Expected: ~K_distinct (5–10) partition traversals × ~0.6 ms ≈ 3–6 ms total — a **50–100× speedup** of the sweep. The three RED-flagged concerns (charToLocalIdx rebuild, CL heap allocation, Het IntegerMatrix copy) are **secondary** to this algorithmic issue and will be resolved as a side effect of batching. |
| M-131 | P2 | BLOCKED (M-155) | **Warmup stabilisation detector: adaptive `windowSize` + validation study.** `.CheckStabilisation()` uses hardcoded `windowSize=10` (Geweke comparison windows), giving a minimum 11,500-iteration warmup floor regardless of tree size. Propose scaling `windowSize` with `nEdge` (`max(5, min(20, nEdge %/% 10))`). Requires empirical validation: run 8 benchmark datasets (20–88 tips) × 4 seeds on Hamilton with warmup disabled (200k iter), then replay the detector offline with a grid of `windowSize` and `nStableRequired` values to verify no false positives. **Briefing:** `.positai/plans/2026-03-29-m131-warmup-stabilisation-validation.md` |


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
| S-RED | dyn | OPEN | **Standing: Red-team review.** Review recent code changes for correctness bugs, edge cases, and safety issues. Focus areas: C++ MCMC hot path (likelihood, proposals, flat-buffer pruning, OPP-1–6 changes), MCMC acceptance logic (rollback correctness for each move type), R orchestration (convergence checking, streaming buffers, checkpoint/resume). File any bugs found as `u.nnn` issue files. When completed, record the focus area and outcome in the Notes column and reset to OPEN. Priority: ≥6 open tasks → P3, 3–5 → P2, <3 → P1. \| Last run: 2026-03-31 round 15 (B). Focus: **M-154 Gibbs kPrime sweep + block shift.** Audited: (1) `jc1_acrv` pruning correctness — tip CL init/reuse across ACRV categories, flg mechanism for internal nodes; (2) `single_char_loglik_jc` — Het path delegates to `pruning_f81_het_acrv_flat` with correct stride/buffer; JC path uses inline pruning; (3) Gibbs full conditional: β×logLik + logPrior (prior untempered, matches cpp_log_prior per-character decomposition for both geometric and logseries); (4) Prior formula cross-check: geometric `logP + u*log1mP` and logseries `k*logC - log(k) - logNorm` both match cpp_log_prior exactly; Beta hyperprior on p correctly excluded (constant in conditional); (5) Block shift: symmetric proposal, feasibility check, partition-level recompute (only type==1 when partLogLik cached); (6) Cache invalidation: both moves rebuild partLogLik, logPrior, and set nodeCL.valid=false; (7) Edge cases: single trans char, all-binary (kObs=2), sparse data, mixed neo+trans. **No bugs found.** All 8 consistency tests pass (logLik exact match with `eval_full_loglik_cpp`, logPrior exact match with R `LogPrior`). **Performance notes for S-PROF:** charToLocalIdx rebuilt every sweep (should cache in McmcData); heavy heap allocation in inner loop (CL buffers per char per candidate k'); IntegerMatrix construction in Het path per candidate. Prev: round 14 (B) — Focus: **`.BuildResult()` / `MkPosterior` output construction**. Audited result assembly, print/summary/plot methods, and `.PostBurninData()` for streaming and multi-run paths. **Filed M-151 (3 bugs):** (a) `plot()` crashes for streaming multi-run — `per_run[[run]]$samples` is NULL, subscript error. (b) `print()` displays empty per-run sample count — `nrow(NULL)` → empty string in cli. (c) `.PostBurninData()` crashes for streaming multi-run with burnin > 0 — `nrow(NULL)` on per-run samples. Verified correct: streaming buffer flush in `.BuildResult()`; `pmax(..., 1L)` prevents division by zero in acceptance rates; tree trimming handles `seq_len(0)` gracefully; ESS/R-hat in summary auto-loads from log file; `MkPosterior()` constructor is clean. Also fixed M-149 #1-5 (tuning crash + state persistence, commit `e56feb5`). Prev: round 13 — `.RunSerialRuns()` audit (M-149 #6-7). Round 12 — checkpoint audit (M-149 #1–5, M-150). Round 11 — cache invalidation audit. Round 10 — M-145. Round 9 — M-144. Round 8 — M-143. Round 7 — M-142. Round 6 — test coverage paths (logHastings -Inf, prior -Inf) and MH rejection. (10) `partLogLik` cache lifecycle consistent: populated at init, cleared by Gibbs/weighted/partial-CL acceptance, correctly rebuilt through per-partition eval when non-empty. (11) Hastings ratios verified: SPR (`log(lRegraft/lMerge)` symmetric candidate counts), TBR (`+log(lSubEdge/lMergeSub)` subtree re-rooting Jacobian), Dirichlet (`logRev - logFwd` Dir density ratio), joint moves (bivariate Jacobian `log(m1) + log(m2)`). (12) Heated acceptance `β × ΔlogLik + ΔlogPrior + logHR` correct (unheated prior, unheated Hastings). (13) Prior computation correct (`Dir(1,...,1)` constant, `lgamma(n)` normalisation). NNI correctly skips prior eval. |
| S-PROF | dyn | OPEN | **Standing: Performance profiling.** Profile the compiled MCMC hot path using VTune (see `r-package-profiling` skill) or `bench::mark()` microbenchmarks. Identify the current top hotspot after OPP-1–6. Check whether `pruning_jc_flat` / `pruning_jc_acrv_flat` show further vectorisation opportunities, whether chain-swap overhead is visible at scale, or whether R↔C++ boundary crossings dominate for small datasets. File any actionable findings as new `M-nnn` tasks. When completed, record the focus and key finding in Notes and reset to OPEN. Priority: same dynamic rule as S-RED. \| Last run: 2026-03-31 round 3 (B). Focus: **M-154 Gibbs kPrime sweep performance.** Benchmarked on Sun2018 (54 taxa, 99 trans, 126 neo, nCat=6). Key finding: **Gibbs sweep 10× slower than equivalent int_walk coverage** (340 ms vs 33 ms). Root cause: `single_char_loglik_jc()` does per-character per-candidate-k' individual tree traversals (99.5% of sweep time), vs batched flat-buffer partition pruning used by int_walk. nCat scales linearly (70ms@1, 124@2, 236@4, 340@6). Final-rebuild negligible (1.8 ms, 0.5%). Block shift is fast (0.66 ms). kObs distribution: 68×kObs=2, 23×kObs=3, 8×kObs=4 (3 partitions). Per-partition int_walk amortized cost: 0.007–0.07 ms/char (batch advantage). Filed M-155 (P1): batch-by-k' optimization for ~50–100× sweep speedup. RED-flagged concerns (charToLocalIdx, CL alloc, Het IntegerMatrix) secondary to algorithmic issue. Prev: round 2 (C) — M-106 optimization, 1.80× cumulative speedup. |
| S-COORD | dyn | OPEN | **Standing: Coordination review.** Review to-do.md and completed-tasks.md for consistency. Check agent logs and u.nnn files. When completed, record round number and actions in Notes and reset to OPEN. Priority: same dynamic rule as S-RED. \| Last run: 2026-03-30 round 7 (A). Actions: (1) Archived M-149 (interrupt-safe checkpointing, 7 sub-bugs) and M-150 (maxTime checkpoint) to completed-tasks.md; removed from to-do.md. 143 completed tasks total. (2) Triaged u.126 — build log with 3 unused-function C++ warnings (minor, not filed as task). Removed from tracking. (3) Untracked accidentally committed artifacts: pcl_diag.txt (3 copies, 7k lines), vignettes/hyoliths.trees, vignettes/hyoliths_{1,2}.log, vignettes/hyoliths_posterior.rds. Added patterns to .gitignore. (4) Updated coordination.md: checkpoint series complete. (5) Test suite: 3972 pass, 0 fail, 15 skip, 0 warnings. (6) Open tasks: M-148 (needs verification), M-151 (display bugs), M-131 (warmup validation). M-080 (C) still ASSIGNED. Standing tasks at P1. |


