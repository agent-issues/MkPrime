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

| M-160 | P2 | OPEN | **Delayed rejection for topology moves.** When an SPR proposal is rejected, immediately try a cheaper NNI at the regraft edge as a fallback before discarding the iteration. Delayed-rejection MH preserves detailed balance with a modified acceptance ratio (Green & Mira 2001). The full likelihood from the rejected SPR can be partially recycled for the NNI evaluation. Potential for significant tree-mixing improvement at modest cost. |
| M-161 | P2 | OPEN | **Per-partition CL cache validity.** `nodeCL.valid = false` invalidates all CacheUnits on any parameter change, even when only one partition is affected (e.g. `rate_loss` only affects neomorphic partitions). Maintain per-CacheUnit validity flags so NNI proposals can reuse cached CLs for unaffected partitions after a parameter-only move. Benefit scales with partition diversity. |
| M-167 | P3 | OPEN | **Add BG hyperparameter moves to `moveWeights` validation.** `MkPrimeMCMC()` `validNames` list (line ~464) doesn't include `kprime_alpha`, `kprime_beta`, `slice_kprime_alpha`, `slice_kprime_beta`, or `block_kPrime`, `neo_joint`, `gibbs_kPrime`, `beta_scale`, `slice_rate_loss`, `slice_rate_neo`, `slice_rate_log_sd`, `slice_tree_length`, `slice_beta_scale`. Users can't override their weights. Audit the full move list in `.BuildMoves()` and sync `validNames` to match. Quick fix. |
| M-166 | P3 | OPEN | **VTune: verify M-164 pre-filter reduced Gibbs sweep CPU share.** Re-profile after M-163 (slice sampler) + M-164 (prior-ceiling pre-filter + LOG_CUTOFF tightening) to check whether the Gibbs kPrime sweep (`pruning_jc_acrv_persite` at 59.1%) CPU share decreased. Same config as S-PROF round 4: Sun2018, 54 taxa, 225 all-trans chars, nCat=6, 15k iterations. Compare against saved VTune baseline in `vtune-out/`. If Gibbs share dropped meaningfully, the next bottleneck may shift to ascertainment (`constant_site_prob_jc` at 6.6%) or exp() calls (11.7%). |
| M-162 | P2 | OPEN | **Hot-path micro-optimizations: pre-allocate scratch vectors.** Part A (raw pointers, 1.1% CPU) DONE in `717818c`. Part B remaining: pre-allocate scratch `std::vector` buffers (`site_lik_sum`, ascertainment buffers) in McmcState/McmcData instead of per-call heap allocation (~1.5% CPU savings). Mechanical change, low risk. |



## Bugs

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-168 | P3 | OPEN | **Segfault in Q-het + both Gibbs topology moves.** `test-gibbs-het.R` test 4 crashes (exit code 139) when `qHeterogeneity = TRUE`, `gibbsSpr = TRUE`, and `gibbsSubtreeSwap = TRUE` simultaneously (8-tip, 12-char dataset, seed 9102). Tests 1–3 pass (F81 unit, Q-het + single Gibbs move). Pre-existing: confirmed at `bc1e2ee~1`. Likely buffer overrun or stale CL cache access in dual-Gibbs Q-het path. Filed from u.165 (S-RED round 17). |


## MCMC infrastructure

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-165 | P2 | ASSIGNED (B) | **Fix BG hyperparameter adaptation + mixing.** Investigation (2026-03-31) found three root causes for low kprime_alpha/beta ESS: (1) **Adaptation never fires:** `.AdaptTuning()` requires ≥20 proposals per tuning round, but kprime_alpha/beta get ~18 (1.8% weight × 1000 budget < threshold). `.AdaptSliceWidths()` requires ≥10, but slice moves get ~2 (0.2% weight). Scales stuck at 0.05 default, slice widths stuck at 1.0. Fix: lower `.AdaptTuning()` threshold from 20 to 10 for hyperparameter moves (or globally). (2) **Default scale too small:** scale=0.05 gives ±2.5% log-scale steps → 100% acceptance (should be ~35%). Posterior log-SD is ~0.14 (alpha) and ~0.21 (beta). Increase defaults to 0.3 (alpha) and 0.5 (beta). Manually setting these improved β/α ESS ratio from 0.27 to 0.85. (3) **moveWeights validation missing BG moves:** `MkPrimeMCMC()` `validNames` list doesn't include kprime_alpha, kprime_beta, slice_kprime_alpha, slice_kprime_beta — users can't override their weights. Add them. Secondary: initial move weights (0.05 for scale, 1.0 for slice) are low relative to other moves; consider increasing to ensure ≥20 proposals per tuning round. |

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
| S-RED | dyn | OPEN | **Standing: Red-team review.** Review recent code changes for correctness bugs, edge cases, and safety issues. Focus areas: C++ MCMC hot path (likelihood, proposals, flat-buffer pruning, OPP-1–6 changes), MCMC acceptance logic (rollback correctness for each move type), R orchestration (convergence checking, streaming buffers, checkpoint/resume). File any bugs found as `u.nnn` issue files. When completed, record the focus area and outcome in the Notes column and reset to OPEN. Priority: ≥6 open tasks → P3, 3–5 → P2, <3 → P1. \| Last run: 2026-03-31 round 17 (B). Focus: **M-163 (BG slice sampler) + M-164 (Gibbs kPrime pre-filter).** Audited commit `bc1e2ee`. (1) `slice_kprime_hyper_impl` (moveType 29): log-scale target `logPrior + Jacobian` is correct; evalTarget restore logic verified; no CL cache invalidation needed (prior-only, correct); logPrior cache updated on acceptance; logLik unchanged (correct); logPosterior not cached (computed on-the-fly, no staleness); heated chain handling correct (beta/tempering not needed since prior-only). (2) R wiring: `.kMoveTypes` mapping (both alpha/beta → moveType 29, distinguished by sliceParamCodes); `.BuildMoves()` creates moves with `type = "slice_kprime_hyper"`, `sliceParamIdx = 0/1`; `alwaysAcceptTypes` includes `"slice_kprime_hyper"`; `.AdaptSliceWidths` keys correct; `.AdaptTuning` correctly skips slice moves via NA tuningKeys. (3) Batch runner dispatch (moveType 29) correct. (4) M-164 `LOG_CUTOFF` tightened -57.5 → -25.0: safe (exp(-25) ≈ 1.4e-11, well below RNG resolution). (5) Prior-ceiling pre-filter: uses `charMaxLL[ti]` as LL upper bound; assumption "LL decreases with k" is heuristic but 25-log-unit margin makes incorrect pruning vanishingly unlikely; `charMaxLL` tracking correct (R_FINITE guard); nActive bookkeeping consistent; heated chain interaction safe (beta consistent on both sides). (6) **Pre-existing segfault found**: `test-gibbs-het.R` test 4 (Q-het + both Gibbs SPR + Gibbs swap) crashes with exit code 139 — confirmed pre-existing (also crashes at `bc1e2ee~1`), filed as `u.165`. Not related to M-163/M-164. (7) Targeted probes: BG slice sampler produces 90+ unique values for both alpha/beta; Gibbs sweep kPrime distributions show healthy variation (kPrime up to 14, no truncation artifacts); all log_posterior values finite. **No bugs found in M-163/M-164.** Tests: 23 BG pass, 16 Gibbs kPrime pass, 9 batched Gibbs pass. Prev: round 16 (B). Focus: **Adaptive tree ESS convergence path.** Audited: (1) `.CheckConvergence()` adaptive three-tier logic (skip/coarse/fine) — tier selection correct for all combinations of minEss, maxRhat, minTreeEss; (2) `.ComputeTreeEssInLoop()` — **BUG FOUND:** returned `Inf` when all chains had 5–6 trees (C++ `median_pseudo_ess_cpp` needs n≥7, returns NA for shorter chains; `min(NA, na.rm=TRUE)` = `Inf`, falsely satisfying any `minTreeEss` threshold). Fixed: filter non-finite values before `min()`; return `NA_real_` when no finite ESS remains. Same fix applied to `.ComputeTreeEss()` in `Convergence.R`. Regression test added (`97819de`). (3) `.TreeESS()` + C++ `frechet_correlation_ess_cpp` / `median_pseudo_ess_cpp` — correct Geyer IPS, monotone correction, tau clamping; (4) Tree construction at 3 sample-handling sites — confirmed M-153 fix (Preorder + correct Nnode) in place; (5) ETA binding constraint logic — correct when both minEss and minTreeEss set; **noted usability gap:** ETA unavailable when only minTreeEss set (etaTarget=NULL); (6) `.CheckConvergenceFromLogs()` — by design doesn't enforce minTreeEss (trees not in logs); (7) `.MinEssPerSec()` tree ESS integration in tuning bandit — correct. Prev: round 15 (B) — M-154 Gibbs kPrime sweep + block shift. Prev: round 14 (B) — `.BuildResult()` / `MkPosterior` output construction**. Audited result assembly, print/summary/plot methods, and `.PostBurninData()` for streaming and multi-run paths. **Filed M-151 (3 bugs):** (a) `plot()` crashes for streaming multi-run — `per_run[[run]]$samples` is NULL, subscript error. (b) `print()` displays empty per-run sample count — `nrow(NULL)` → empty string in cli. (c) `.PostBurninData()` crashes for streaming multi-run with burnin > 0 — `nrow(NULL)` on per-run samples. Verified correct: streaming buffer flush in `.BuildResult()`; `pmax(..., 1L)` prevents division by zero in acceptance rates; tree trimming handles `seq_len(0)` gracefully; ESS/R-hat in summary auto-loads from log file; `MkPosterior()` constructor is clean. Also fixed M-149 #1-5 (tuning crash + state persistence, commit `e56feb5`). Prev: round 13 — `.RunSerialRuns()` audit (M-149 #6-7). Round 12 — checkpoint audit (M-149 #1–5, M-150). Round 11 — cache invalidation audit. Round 10 — M-145. Round 9 — M-144. Round 8 — M-143. Round 7 — M-142. Round 6 — test coverage paths (logHastings -Inf, prior -Inf) and MH rejection. (10) `partLogLik` cache lifecycle consistent: populated at init, cleared by Gibbs/weighted/partial-CL acceptance, correctly rebuilt through per-partition eval when non-empty. (11) Hastings ratios verified: SPR (`log(lRegraft/lMerge)` symmetric candidate counts), TBR (`+log(lSubEdge/lMergeSub)` subtree re-rooting Jacobian), Dirichlet (`logRev - logFwd` Dir density ratio), joint moves (bivariate Jacobian `log(m1) + log(m2)`). (12) Heated acceptance `β × ΔlogLik + ΔlogPrior + logHR` correct (unheated prior, unheated Hastings). (13) Prior computation correct (`Dir(1,...,1)` constant, `lgamma(n)` normalisation). NNI correctly skips prior eval. |
| S-PROF | dyn | OPEN | **Standing: Performance profiling.** Profile the compiled MCMC hot path using VTune (see `r-package-profiling` skill) or `bench::mark()` microbenchmarks. Identify the current top hotspot after OPP-1–6. Check whether `pruning_jc_flat` / `pruning_jc_acrv_flat` show further vectorisation opportunities, whether chain-swap overhead is visible at scale, or whether R↔C++ boundary crossings dominate for small datasets. File any actionable findings as new `M-nnn` tasks. When completed, record the focus and key finding in Notes and reset to OPEN. Priority: same dynamic rule as S-RED. \| Last run: 2026-03-31 round 4. Focus: **VTune hotspot collection on Sun2018** (54 taxa, 225 all-trans chars, nCat=6). VTune 2025.10 at `C:/Program Files (x86)/Intel/oneAPI/vtune/latest/bin64/`, user-mode sampling (hw PMU needs admin). 90s CPU captured. **Top hotspot: `pruning_jc_acrv_persite` at 59.1%** — the batched per-site Felsenstein pruning used by Gibbs kPrime sweep (M-155). Regular MH proposals (`pruning_jc_acrv_flat`) only 3.9%. **`_expl_internal` (exp()) at 11.7%** — one `exp()` per edge × rate category × traversal; already amortized across characters but unavoidable across candidate k' (eigenvalue depends on k). **`constant_site_prob_jc` (ascertainment) at 6.6%** — separate tree traversals for constant-site correction; could batch like main pruning. Minor: vector copies 1.5%, Rcpp bounds checks 1.1%. Optimization opportunities: (1) fast approx exp could save ~5-8%, (2) batched ascertainment correction, (3) eliminate hot-path vector copies. VTune results saved in `vtune-out/`. Prev: round 3 (B) — M-154 Gibbs sweep 10× overhead, filed M-155. |
| S-COORD | dyn | OPEN | **Standing: Coordination review.** Review to-do.md and completed-tasks.md for consistency. Check agent logs and u.nnn files. When completed, record round number and actions in Notes and reset to OPEN. Priority: same dynamic rule as S-RED. \| Last run: 2026-03-31 round 8 (B). Actions: (1) Triaged u.165 → M-168 (P3, pre-existing Q-het segfault with dual Gibbs moves); deleted u.165. (2) Verified completed-tasks.md current (162 entries, M-163/M-164 archived in `ea01a51`). (3) No stale Rscript processes. No pending remote jobs. (4) Uncommitted M-165 work in progress (B): R/MkPrimeMCMC.R, R/RunMkPrime.R, src/mcmc.cpp, src/proposals.cpp — BG hyperparameter adaptation fix. (5) u.124 resolved in user commit `c76ea44` (move-weights display restyle). (6) Updated coordination.md with M-163/M-164/M-165 milestones. (7) Task counts: 8 OPEN specific (M-158, M-160, M-161, M-162, M-131, M-167, M-166, M-168) + 2 ASSIGNED (M-165 B, M-080 C). Standing tasks at P3 (≥6 open). Prev: round 7 (A). |


