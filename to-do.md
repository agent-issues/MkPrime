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



## Performance optimization

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-158 | P2 | OPEN | **Partial CL evaluation for SPR moves.** SPR currently invalidates the entire CL cache, forcing O(nEdge) full rebuild. Extend the dirty-flag CL system (`node_cl_cache.h`) to handle SPR: incrementally update TreeNav at prune + regraft points, identify dirty set (two paths to root + detachment/reattachment nodes), partial recompute. For 54-tip tree: dirty set ~15–20 nodes vs 105 edges → 5–7× speedup for the most heavily weighted topology move. High engineering effort; validate carefully against full-eval path. |
| M-159 | P1 | OPEN | **Cache-aware move scheduling.** When the CL cache is valid, NNI/beta-simplex run ~15× faster via partial CL. Currently move selection is cache-unaware. Add a cache-bonus multiplier to partial-CL-eligible moves (NNI, beta_simplex, Dirichlet) that activates when `nodeCL.valid == true`. Statistically valid: state-dependent mixture kernels preserve detailed balance as long as each component kernel is valid. Can integrate with existing warmup adaptation and tuning bandit. |
| M-160 | P2 | OPEN | **Delayed rejection for topology moves.** When an SPR proposal is rejected, immediately try a cheaper NNI at the regraft edge as a fallback before discarding the iteration. Delayed-rejection MH preserves detailed balance with a modified acceptance ratio (Green & Mira 2001). The full likelihood from the rejected SPR can be partially recycled for the NNI evaluation. Potential for significant tree-mixing improvement at modest cost. |
| M-161 | P2 | OPEN | **Per-partition CL cache validity.** `nodeCL.valid = false` invalidates all CacheUnits on any parameter change, even when only one partition is affected (e.g. `rate_loss` only affects neomorphic partitions). Maintain per-CacheUnit validity flags so NNI proposals can reuse cached CLs for unaffected partitions after a parameter-only move. Benefit scales with partition diversity. |
| M-162 | P2 | OPEN | **Hot-path micro-optimizations: Rcpp bounds checks + vector copies.** (a) In pruning hot paths, replace `Rcpp::NumericVector::operator[]` with raw `double*` pointer access to eliminate `check_index` overhead (1.1% CPU). (b) Pre-allocate scratch `std::vector` buffers in McmcState instead of per-call construction (1.5% CPU). Mechanical changes, low risk, ~2.5% combined savings. |


## MCMC infrastructure

| ID | Priority | Status | Description |
|----|----------|--------|-------------|


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
| S-RED | dyn | OPEN | **Standing: Red-team review.** Review recent code changes for correctness bugs, edge cases, and safety issues. Focus areas: C++ MCMC hot path (likelihood, proposals, flat-buffer pruning, OPP-1–6 changes), MCMC acceptance logic (rollback correctness for each move type), R orchestration (convergence checking, streaming buffers, checkpoint/resume). File any bugs found as `u.nnn` issue files. When completed, record the focus area and outcome in the Notes column and reset to OPEN. Priority: ≥6 open tasks → P3, 3–5 → P2, <3 → P1. \| Last run: 2026-03-31 round 16. Focus: **Adaptive tree ESS convergence path.** Audited: (1) `.CheckConvergence()` adaptive three-tier logic (skip/coarse/fine) — tier selection correct for all combinations of minEss, maxRhat, minTreeEss; (2) `.ComputeTreeEssInLoop()` — **BUG FOUND:** returned `Inf` when all chains had 5–6 trees (C++ `median_pseudo_ess_cpp` needs n≥7, returns NA for shorter chains; `min(NA, na.rm=TRUE)` = `Inf`, falsely satisfying any `minTreeEss` threshold). Fixed: filter non-finite values before `min()`; return `NA_real_` when no finite ESS remains. Same fix applied to `.ComputeTreeEss()` in `Convergence.R`. Regression test added (`97819de`). (3) `.TreeESS()` + C++ `frechet_correlation_ess_cpp` / `median_pseudo_ess_cpp` — correct Geyer IPS, monotone correction, tau clamping; (4) Tree construction at 3 sample-handling sites — confirmed M-153 fix (Preorder + correct Nnode) in place; (5) ETA binding constraint logic — correct when both minEss and minTreeEss set; **noted usability gap:** ETA unavailable when only minTreeEss set (etaTarget=NULL); (6) `.CheckConvergenceFromLogs()` — by design doesn't enforce minTreeEss (trees not in logs); (7) `.MinEssPerSec()` tree ESS integration in tuning bandit — correct. Prev: round 15 (B) — M-154 Gibbs kPrime sweep + block shift. Prev: round 14 (B) — `.BuildResult()` / `MkPosterior` output construction**. Audited result assembly, print/summary/plot methods, and `.PostBurninData()` for streaming and multi-run paths. **Filed M-151 (3 bugs):** (a) `plot()` crashes for streaming multi-run — `per_run[[run]]$samples` is NULL, subscript error. (b) `print()` displays empty per-run sample count — `nrow(NULL)` → empty string in cli. (c) `.PostBurninData()` crashes for streaming multi-run with burnin > 0 — `nrow(NULL)` on per-run samples. Verified correct: streaming buffer flush in `.BuildResult()`; `pmax(..., 1L)` prevents division by zero in acceptance rates; tree trimming handles `seq_len(0)` gracefully; ESS/R-hat in summary auto-loads from log file; `MkPosterior()` constructor is clean. Also fixed M-149 #1-5 (tuning crash + state persistence, commit `e56feb5`). Prev: round 13 — `.RunSerialRuns()` audit (M-149 #6-7). Round 12 — checkpoint audit (M-149 #1–5, M-150). Round 11 — cache invalidation audit. Round 10 — M-145. Round 9 — M-144. Round 8 — M-143. Round 7 — M-142. Round 6 — test coverage paths (logHastings -Inf, prior -Inf) and MH rejection. (10) `partLogLik` cache lifecycle consistent: populated at init, cleared by Gibbs/weighted/partial-CL acceptance, correctly rebuilt through per-partition eval when non-empty. (11) Hastings ratios verified: SPR (`log(lRegraft/lMerge)` symmetric candidate counts), TBR (`+log(lSubEdge/lMergeSub)` subtree re-rooting Jacobian), Dirichlet (`logRev - logFwd` Dir density ratio), joint moves (bivariate Jacobian `log(m1) + log(m2)`). (12) Heated acceptance `β × ΔlogLik + ΔlogPrior + logHR` correct (unheated prior, unheated Hastings). (13) Prior computation correct (`Dir(1,...,1)` constant, `lgamma(n)` normalisation). NNI correctly skips prior eval. |
| S-PROF | dyn | OPEN | **Standing: Performance profiling.** Profile the compiled MCMC hot path using VTune (see `r-package-profiling` skill) or `bench::mark()` microbenchmarks. Identify the current top hotspot after OPP-1–6. Check whether `pruning_jc_flat` / `pruning_jc_acrv_flat` show further vectorisation opportunities, whether chain-swap overhead is visible at scale, or whether R↔C++ boundary crossings dominate for small datasets. File any actionable findings as new `M-nnn` tasks. When completed, record the focus and key finding in Notes and reset to OPEN. Priority: same dynamic rule as S-RED. \| Last run: 2026-03-31 round 4. Focus: **VTune hotspot collection on Sun2018** (54 taxa, 225 all-trans chars, nCat=6). VTune 2025.10 at `C:/Program Files (x86)/Intel/oneAPI/vtune/latest/bin64/`, user-mode sampling (hw PMU needs admin). 90s CPU captured. **Top hotspot: `pruning_jc_acrv_persite` at 59.1%** — the batched per-site Felsenstein pruning used by Gibbs kPrime sweep (M-155). Regular MH proposals (`pruning_jc_acrv_flat`) only 3.9%. **`_expl_internal` (exp()) at 11.7%** — one `exp()` per edge × rate category × traversal; already amortized across characters but unavoidable across candidate k' (eigenvalue depends on k). **`constant_site_prob_jc` (ascertainment) at 6.6%** — separate tree traversals for constant-site correction; could batch like main pruning. Minor: vector copies 1.5%, Rcpp bounds checks 1.1%. Optimization opportunities: (1) fast approx exp could save ~5-8%, (2) batched ascertainment correction, (3) eliminate hot-path vector copies. VTune results saved in `vtune-out/`. Prev: round 3 (B) — M-154 Gibbs sweep 10× overhead, filed M-155. |
| S-COORD | dyn | OPEN | **Standing: Coordination review.** Review to-do.md and completed-tasks.md for consistency. Check agent logs and u.nnn files. When completed, record round number and actions in Notes and reset to OPEN. Priority: same dynamic rule as S-RED. \| Last run: 2026-03-30 round 7 (A). Actions: (1) Archived M-149 (interrupt-safe checkpointing, 7 sub-bugs) and M-150 (maxTime checkpoint) to completed-tasks.md; removed from to-do.md. 143 completed tasks total. (2) Triaged u.126 — build log with 3 unused-function C++ warnings (minor, not filed as task). Removed from tracking. (3) Untracked accidentally committed artifacts: pcl_diag.txt (3 copies, 7k lines), vignettes/hyoliths.trees, vignettes/hyoliths_{1,2}.log, vignettes/hyoliths_posterior.rds. Added patterns to .gitignore. (4) Updated coordination.md: checkpoint series complete. (5) Test suite: 3972 pass, 0 fail, 15 skip, 0 warnings. (6) Open tasks: M-148 (needs verification), M-151 (display bugs), M-131 (warmup validation). M-080 (C) still ASSIGNED. Standing tasks at P1. |


