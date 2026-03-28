# MkPrime — Completed Tasks

| ID | Description | Agent | Date | Notes |
|----|-------------|-------|------|-------|
| M-001 | Package skeleton (DESCRIPTION, Rcpp stub, testthat harness) | A | 2026-03-26 | `693c6ee` |
| M-002 | `MkPrimeData()` — phyDat input and character classification | A | 2026-03-26 | `5afd39e`. 19 tests. |
| M-003 | Character partitioning by type and kObs | A | 2026-03-26 | `efc1ba0`. 10 new tests. |
| M-004 | Data validation and edge cases | A | 2026-03-26 | `46386b6`. Invariant drop, index remapping. |
| M-005 | Comprehensive unit tests for data layer | A | 2026-03-26 | `46386b6`. 70 total tests. |
| M-006 | JC(k') rate matrix and analytical P(t) | A | 2026-03-26 | `da4492f`. O(1) eigendecomposition. |
| M-007 | MkN asymmetric binary rate matrix and P(t) | A | 2026-03-26 | `da4492f`. rate_loss parameterization. |
| M-008 | Felsenstein pruning (JC + MkN) | A | 2026-03-26 | `8e0fd20`. Post-order traversal, hand-verified. |
| M-009 | ACRV discretized lognormal rate categories | A | 2026-03-26 | `f2d5ee3`. 6-category midpoint quantiles. |
| M-010 | Ascertainment bias correction (variable coding) | A | 2026-03-26 | `ffb5769`. Constant site probability. |
| M-011 | Mk' relabelling correction | A | 2026-03-26 | `e001c74`. lgamma-based, scalar + batch. |
| M-012 | Likelihood validation + mkp_loglikelihood() | A | 2026-03-26 | `25e08ad`. Matches phangorn to machine epsilon. 194 total tests. |
| M-013 | Fix per-character k' in likelihood | A | 2026-03-26 | `41105eb`. kPrime > kObs now uses correct JC(k') matrix. |
| M-014 | Prior specification & log-prior (MkPrimeModel) | A | 2026-03-26 | `41105eb`. Gamma/LogNormal/Geometric priors with boundary checks. |
| M-015 | Proposal functions (Scale, BetaSimplex, BoundedIntegerWalk) | A | 2026-03-26 | `41105eb`. Correct Hastings ratios, reversibility tested. |
| M-016 | MCMC state object & initialization | A | 2026-03-26 | `41105eb`. tree_length + rel_br_lengths parameterization. |
| M-017 | Single-chain MH engine (RunMkPrime) | A | 2026-03-26 | `41105eb`. R-side loop, weighted move schedule. |
| M-018 | Sample accumulation & progress bar | A | 2026-03-26 | `41105eb`. In-memory matrix, cli progress. |
| M-019 | Proposal adaptation | A | 2026-03-26 | `41105eb`. Multiplicative tuning during warmup. |
| M-020 | MkPosterior result object | A | 2026-03-26 | `41105eb`. print/summary/plot methods, ESS via coda. |
| M-021 | Integration validation | A | 2026-03-26 | `41105eb`. Posterior recovery on simulated data. 466 total tests. |
| M-022 | Topology as mutable MCMC state | A | 2026-03-26 | `7dcd780`. state$tree stores topology; .do_move() updated. |
| M-023 | NNI proposal on unrooted binary trees | A | 2026-03-26 | `7dcd780`. Symmetric (Hastings=0), swaps across internal edges. |
| M-024 | SPR proposal on unrooted binary trees | A | 2026-03-26 | `7dcd780`. Node u moves; Jacobian correction log(l_regraft/l_merge). |
| M-025 | Tree move integration into MCMC loop | A | 2026-03-26 | `7dcd780`. NNI wt=nEdge/2, SPR wt=nEdge/4, fix_topology param. |
| M-026 | exp_steps default from Fitch parsimony | A | 2026-03-26 | `266dc07`. Matches phangorn exactly. |
| M-027 | Tree logging to Newick file | A | 2026-03-26 | `266dc07`. tree_file param in MkPrimeMCMC, append-mode. |
| M-028 | Tree topology recovery test | A | 2026-03-26 | `fd38bd2`. Recovers true 5-tip tree from simulated data. ~2059 tests. |
| M-029 | Parallel tempering core (multi-chain + temp ladder) | A | 2026-03-26 | `4fe7fe8`. Geometric beta ladder, heated MH acceptance, per-chain tuning. |
| M-030 | Chain swap proposals | A | 2026-03-26 | `4fe7fe8`. Adjacent pair swaps with correct Metropolis criterion. |
| M-031 | Adaptive temperature tuning | A | 2026-03-26 | `77042ae`. Multiplicative heat adjustment targeting 25% swap rate. |
| M-032 | Independent runs (nRuns) | A | 2026-03-26 | `911464b`. Per-run state, .perturb_start(), .combine_runs(). |
| M-033 | Convergence monitoring (ESS + PSRF) | A | 2026-03-26 | `2e05242`. convergence_diagnostics() via coda. |
| M-034 | Stopping rules + interleaved runs | A | 2026-03-26 | `266f571`. max_time, min_ess, max_psrf, check_every. Lockstep run architecture. |
| M-035 | Checkpointing | A | 2026-03-26 | `79ccb6a`. checkpoint_file, resume_mkprime(). |
| M-036 | MkPosterior multi-run updates | A | 2026-03-26 | `623b9a9`. ESS/PSRF in summary, per-run traces in plot. 2898 total tests. |
| M-037 | Progress callback infrastructure | A | 2026-03-26 | plot_every + progress_fn params in MkPrimeMCMC; .build_progress_info() in RunMkPrime + resume. |
| M-038 | PlotDuringMCMC — console trace plots | A | 2026-03-26 | mkp_trace_plot(): multi-panel base-R traces, per-run colors, warmup handling. |
| M-039 | PNG progress writer | A | 2026-03-26 | mkp_png_progress(): factory fn, atomic PNG rename, hand-rolled JSON status. 46 new tests. |
| M-040 | Standalone Shiny app (EasyMkPrime) | A | 2026-03-26 | bslib page_sidebar, callr background MCMC, PNG progress polling, traces/summary/consensus tabs. |
| M-042 | _pkgdown.yml for MkPrime | A | 2026-03-26 | Created by human; confirmed present. |
| M-043 | Auto-detect neomorphic characters | A | 2026-03-26 | auto_detect_neomorphic(): checks level labels {0,1}; wired into EasyMkPrime data loading. 14 new tests. |
| M-041 | Migrate issues protocol to u.nnn files | A | 2026-03-26 | Updated parent, TreeSearch, MkPrime AGENTS.md. Deprecated issues.md files. |
| M-044 | C++ singleton_site_prob_jc() | A | 2026-03-26 | n pseudo-chars × k*(k-1) JC symmetry. ACRV averaging. |
| M-045 | C++ singleton_site_prob_mkn() | A | 2026-03-26 | 2n pseudo-chars (asymmetric states). |
| M-046 | R-level coding="informative" | A | 2026-03-26 | MkpLogLikelihood + MkPrimeModel support. P(uninf) = P(const) + P(singleton). |
| M-047 | Tests for coding="informative" | A | 2026-03-26 | 20 new tests in test-ascertainment.R. |
| M-048 | Partition rate scalar (rate_neo) | A | 2026-03-26 | Multiplies neomorphic branch lengths. LogNormal prior, scale proposal, adaptive tuning. |
| M-049 | Tests for rate_neo | A | 2026-03-26 | Likelihood, prior, MCMC smoke test, backward compat. |
| M-050 | Stepping-stone marginal likelihood | A | 2026-03-26 | mkp_stepping_stone(): power posterior with Beta(α,1) schedule. |
| M-051 | Tests for stepping-stone | A | 2026-03-26 | 11 tests: finite output, consistency, neomorphic, fixed topology. |
| M-055 | Eliminate redundant ape::reorder.phylo from likelihood hot path | C | 2026-03-26 | .MkpLogLikelihood() internal fast-path; .DoMove() and .InitState() use it. |
| M-056 | Pre-compute move weight vector | (prev) | 2026-03-26 | Already implemented — weights pre-computed before main loop. |
| M-057 | Replace ape::reorder.phylo with TreeTools::Postorder | C | 2026-03-26 | proposals.R, RunMkPrime.R, MkPrimeModel.R. Exported MkpLogLikelihood retains ape for user-facing use. |
| M-069 | Fix stepping-stone postorder invariant + NaN SE | B | 2026-03-26 | mkp_stepping_stone() now reorders tree before .InitState(). .EssVector() handles NaN sd. SE loop handles degenerate stones (all -Inf logLiks). 3199 tests pass. |
| M-058 | C++ NNI proposal (nni_proposal in tree_moves.cpp) | B | 2026-03-26 | TreeTools::postorder_order() for reordering. R wrapper reconstructs phylo. Fixes heredoc contamination in proposals.cpp/.R. 3197 tests pass. |
| M-059 | C++ SPR proposal (spr_proposal in proposals.cpp) | C | 2026-03-26 | BFS descendant exclusion, Jacobian Hastings ratio. |
| M-060 | C++ BetaSimplex proposal (beta_simplex_proposal in proposals.cpp) | C | 2026-03-26 | In-place via R::rbeta/dbeta. |
| M-061 | C++ MCMC inner loop (do_move_cpp, McmcState XPtr) | B/C | 2026-03-26 | mcmc.cpp: McmcState XPtr, do_move_cpp() (all proposals + prior + lik in C++). mcmc_likelihood.cpp: prepare_mcmc_data(), cpp_log_likelihood(). mcmc_state.h: McmcData/McmcState. R: .InitMcmcData, .InitMcmcChain, .DoMove (XPtr path), .SaveCheckpoint serializes XPtrs. SteppingStone updated to use new engine. 3197 tests pass. |
| M-066 | Automatic burnin selection for MkPosterior | C | 2026-03-26 | SetBurnin(), AutoBurnin(), .PostBurninData(). summary/print/plot/ConvergenceDiagnostics all respect burnin. 25 new tests. |
| M-067 | Vignette and package doc references via Rdpack | C | 2026-03-26 | inst/REFERENCES.bib (9 entries), Rdpack in Imports/RdMacros, \insertCite in MkpLogLikelihood, bibliography in hyoliths.qmd. |
| M-068 | Rogue taxon suppression in hyoliths.qmd | C | 2026-03-26 | Rogue::QuickRogue() in consensus tree section, both production and demo blocks. Rogue added to Suggests. |
| M-062 | C++ outer iteration loop (run_mcmc_batch_cpp) | B | 2026-03-26 | `002f603`. Refactored `do_move_cpp` into `do_move_impl` (raw ptrs) + thin wrapper. `run_mcmc_batch_cpp()` runs 200-iter batches: weighted move selection, chain swaps (std::swap), sample collection. R loop calls it per-batch, handles adaptation/convergence/progress. 335+ tests pass. |
| M-064 | Partial likelihood recalculation for parameter-only moves | C | 2026-03-26 | `de0d7c3`. Per-partition caching: p moves skip lik entirely, rate_loss/rate_neo recompute only neomorphic partitions, kPrime recompute only affected partition. cpp_partition_log_likelihood() + fill_partition_cache(). MCMC engine tests 3.4s→2.2s (~35% speedup). |
| M-065 | Eliminate edge matrix round-trips in C++ MCMC | C | 2026-03-26 | `a34d674`. Vector-based cpp_log_likelihood, nni/spr_proposal_impl; -71 lines net. |
| M-063 | C++ state struct with pre-allocated CL workspace | B | 2026-03-27 | `ClWorkspace` struct (flat `double[]` buf + `uint8_t[]` init, nNodeMax × strideMax) added to `mcmc_state.h` and `McmcState`. Four static flat-buffer pruning helpers (`pruning_jc_flat`, `pruning_jc_acrv_flat`, `pruning_mkn_flat`, `pruning_mkn_acrv_flat`) in `mcmc_likelihood.cpp`; used automatically when workspace fits, fall back to allocating versions when not. `allocate_cl_workspace()` Rcpp export sizes workspace from partition info + kPrimeMax+4 headroom. All three entry points (RunMkPrime, ResumeMkPrime, SteppingStone) call it after fill_partition_cache. Also fixed pre-existing resume bug: fill_partition_cache now runs after init_mcmc_state (not before). 3197 tests pass. |
| M-070 | RunMkPrime auto-resume from checkpoint | C | 2026-03-27 | `c01df8d`. `overwrite = FALSE` param; delegates to `ResumeMkPrime()` when checkpointFile exists. 2 new tests. |
| M-072 | hyoliths.qmd tree-summary-demo fixes | C | 2026-03-27 | `7b4b1bc`. `nzchar(keepNA=FALSE)` NA filter for rogueNames; `par(mar=rep(0,4))` before plot; conditional "rogues excluded" title. Production chunk updated to match. |
| M-071 | Progress bar redesign + nIter=Inf default | C | 2026-03-27 | `e3919a7`, `f3ada92`. MkPrimeMCMC nIter=Inf default (warmup=5000L); for→repeat loop in RunMkPrime+ResumeMkPrime; new format "iter N \| logP: N \| ESS: N \| PSRF: N"; dynamic buffer doubling; convergence check added to ResumeMkPrime. Fixed latent `.CheckConvergence` bug: `converged=TRUE` when no criteria → now `FALSE`. 3211 tests pass. |
| M-073 | Tree ESS in ConvergenceDiagnostics | C | 2026-03-27 | `67892fd`. `.ComputeTreeEss()`: treess+RF, subsample to 1000/run, interactive message. `print.MkpDiagnostics` always shows topology rows. TreeDist added to Suggests. |
| M-075 | Streaming log output (logFile / bufferSize / ReadMkLog) | C | 2026-03-27 | `5b48a79`. Tracer-compatible TSV log during MCMC. Flush-buffer design; per-run log files for nRuns>1; checkpoint v2 includes buffer state; convergence uses rolling window in streaming mode. ReadMkLog() to reload. |
| M-074 | Per-parameter progress table during MCMC | C | 2026-03-27 | `fc5f391`. `.PrintProgressTable()` prints ESS/PSRF per parameter at each checkEvery interval (mirrors print.MkpDiagnostics). `.CheckConvergence()` extended: nRuns=1 support, returns ess/psrf vectors. 4 new tests. |
| M-076 | `MkpWatchLog()` — live trace plot from disk log files | C | 2026-03-27 | `2c9f2a3`. Polling watcher + `.WatchLogPlot()` internal; `MkLogPaths()` path expander; 10 tests. |
| M-077 | EasyMkPrime: auto-detect neomorphic characters | C | 2026-03-27 | Fixed `auto_detect_neomorphic` typo → `AutoDetectNeomorphic()` in app.R; added detection notification. Standalone helper + 5 tests already in place (MkPrimeData.R / test-MkPrimeData.R). |
| M-081 | Cancel-file support in RunMkPrime/ResumeMkPrime | C | 2026-03-27 | `bffd844`. `cancelFile` in MkPrimeMCMC(); checked every batch; flush+checkpoint+break with stop_reason='cancelled'. .FlushAndSaveCheckpoint() refactors doCheck duplication. MkCancelPath() exported. 5 tests. |
| M-093 | PID liveness check on Reconnect after OS-kill | E | 2026-03-27 | `.PidIsAlive(pid)`: tries `ps` package, falls back to Windows `tasklist` or POSIX signal-0. Integrated into Reconnect handler. `ps` added to Suggests. 3 new unit tests. |
| M-096 | Tests and docs for parallel mode | E | 2026-03-27 | `c2527ce` on `feature/parallel-runs`. test-parallel.R: MkPrimeMCMC validation (2 tests), nRuns=1 fallback, sequential-plan orchestration test (exercises full .RunParallelRuns() path), skippable multisession test (guarded via pkgload::is_dev_package()). @section HPC usage added to RunMkPrime() roxygen. man/ not regenerated (stale installed pkg — do after merge). |
| M-095 | Parallel orchestration via `future` | E | 2026-03-27 | `329c61d` on `feature/parallel-runs`. `.RunParallelRuns()`: validates future, forces streaming, per-run cancel files, future::future() workers, polling loop (sleep/user-cancel/maxTime/.CheckConvergenceFromLogs/resolved), collects via future::value(). `.CheckConvergenceFromLogs()`: disk-based ESS+PSRF from log files. `MkPrimeMCMC()` gains `parallel`+`pollInterval`. `future` in Suggests. Sequential path unchanged. |
| M-094 | Extract `.RunMkPrimeSingleRun()`; fix multi-run checkpoint regression | E | 2026-03-27 | `01255d4` on `feature/parallel-runs`. Pure refactor: both RunMkPrime() and ResumeMkPrime() delegate to new internal. Regression fix: .SaveCheckpoint() guards against already-serialized chainStates; RunMkPrime() saves combined checkpoint after all sequential runs for nRuns > 1. Checkpoint, stopping, streaming tests all pass. |
| M-082 | Checkpoint integration in EasyMkPrime (detached-process design) | E | 2026-03-27 | `.RelaunchFromCheckpoint(job, rv)`: clears stale signals, relaunches `mkp_run.R` (RunMkPrime auto-resumes via overwrite=FALSE + checkpointFile in mcmc config), updates job.rds with new PID. Reconnect handler now has full decision tree: done→show; error/cancel/OS-kill + checkpoint→relaunch; error/cancel/OS-kill + no checkpoint→terminal status. 4 new testServer tests for signal-file paths. |
| M-079 | Refactor app.R to use MkBayesianServer | E | 2026-03-27 | Removed callr/PNG inline logic (−200 lines). Data + Starting tree accordions kept app-owned. MkBayesianUi("bayes") below sidebar accordion. Progress tab removed (progress in sidebar). Traces/Summary/Consensus fed from result.rds read on module status=="done". EasyMkPrime() updated: requires processx not callr. |
| M-078 | MkBayesianUi / MkBayesianServer — detached-process Shiny module | E | 2026-03-27 | New R/BayesianModule.R. MkBayesianUi(): accordion (MCMC config, logDir, neomorphic), Run/Stop/Reconnect buttons, live trace plot + ESS table. MkBayesianServer(): launches detached Rscript via processx, writes job.rds, polls log files every 5s, cancel via cancelFile, Reconnect reads job.rds from logDir. Returns jobFile/trees/status reactives. processx added to Suggests. 12 tests. |
| M-083 | Extract `compute_full_loglik()` / `compute_full_loglik_at()` C++ helpers | B | 2026-03-27 | `dad1fcc` on `feature/gibbs-weighted-moves`. Pure C++ helpers for evaluating full posterior likelihood from within proposal code. Blocker for M-085–M-089. |
| M-084 | Subtree swap primitives (`swap_subtrees`, `get_valid_swap_partners`) | B | 2026-03-27 | `fc37ccd` on `feature/gibbs-weighted-moves`. Round-trip and hand-verified tests. Used by M-086 and M-089. |
| M-085 | GibbsSPR (moveType 10) | B | 2026-03-27 | `974a4cf` on `feature/gibbs-weighted-moves`. O(N) Gibbs sampling of SPR reattachment candidates. |
| M-086 | GibbsSubtreeSwap (moveType 11) | B | 2026-03-27 | `974a4cf` on `feature/gibbs-weighted-moves`. O(N) Gibbs sampling with subtree swap operation. |
| M-087 | WeightedBranchLengthScale (moveType 12) | B | 2026-03-27 | `cfd8097` on `feature/gibbs-weighted-moves`. O(B) discretised branch fractions with Beta(0.25,0.25) bins. |
| M-088 | WeightedSPR (moveType 13) | B | 2026-03-27 | `1d49f69` on `feature/gibbs-weighted-moves`. O(N×B) combined Gibbs topology + marginalised branch fractions. |
| M-089 | WeightedSubtreeSwap (moveType 14) | B | 2026-03-28 | `b894fe6` on `feature/gibbs-weighted-moves`. O(N×B) combined Gibbs subtree swap + marginalised branch fractions. |
| M-090 | Wire Gibbs/Weighted moves into MkPrimeMCMC API and dispatch | E | 2026-03-28 | `5533749` on `feature/gibbs-weighted-moves`. 6 new MkPrimeMCMC() params (gibbsSpr, gibbsSubtreeSwap, weightedBranchScale, weightedSpr, weightedSubtreeSwap, nBranchBins). nBranchBins stored in McmcData via set_branch_bins() setter. .BuildMoves() conditionally includes 5 new moves. .kMoveTypes updated (codes 10,11). .AdaptTuning() skips Gibbs/Weighted (NA keys). 49 new tests. |
| M-092 | Adaptive move scheduler | B | 2026-03-28 | `1d1c275` on `feature/gibbs-weighted-moves`. C++ per-move wall-time tracking; MkPrimeMCMC() moveWeights param for pinned weights; log-space softmax with temperature annealing (2.0→0.5) during warmup; per-move floor enforcement; weights frozen post-warmup. Checkpoint includes moveWeights for resume. 45 tests. |
| M-091 | Mixing validation: Gibbs/Weighted vs standard moves | B | 2026-03-28 | Vinther 2008 (23 taxa). Gibbs SPR 32% acc vs standard SPR 2%, but 7.5x slower → baseline wins ESS/s on small trees. WeightedSPR 24x slower, no net gain. Found wMin floor bug (148d31c). Production-build + larger dataset comparison recommended. |
| M-054 | Block Gibbs branch-length sweep (reframed from HMC) | A | 2026-03-28 | `2696b5d` on `feature/gibbs-weighted-moves`. moveType 15: random-permutation-scan MH-within-Gibbs over all edge pairs using bin-based approximate conditionals. Dimension-adjusted adaptive scheduler (score = accept_rate × dim / cost). Default OFF. 24 tests. |
| M-097 | Progress display fixes: warmup label + table overwrite | F,C | 2026-03-28 | F: `3cff5f6` on `feature/tree-ess` (cli pause/recreate, trailing rule, 15 tests). C: `fca2f44` on `main` (ANSI cursor-up overwrite, `prevLines` return value, 8 tests). Merge of tree-ess needs conflict resolution — take union of both approaches. |
| M-098 | Exclude kPrime from convergence criteria | C | 2026-03-28 | `fca2f44`. `.CheckConvergence()`, `.CheckConvergenceFromLogs()`, and `ConvergenceDiagnostics()` exclude kPrime_ from minEss/maxPsrf. kPrime are discrete nuisance parameters being marginalized over — remain in ESS output for display but no longer gate convergence. 7 tests. |
