# Red-team focus areas — mkp

Rotation list. The /red-team skill picks the next area each invocation.

**Rescoped 2026-09-18.** Areas 1 and 2 were written for the empirical-geometric
prior phase and were still tagged `[EG]` "while that arm runs on Hamilton". That
arm has concluded — EG-001/002 resolved as not-a-bug (a prior is pre-data, so
`kObs` must not enter it), EG-003 recorded as a model property — and the released
scope is MkNT. Both rows now cover the seams that actually moved since: the whole
k′ prior family rather than one member of it, and the marginal-k evaluator
alongside the p-samplers. Scopes were widened, not replaced, so the `area:1` and
`area:2` labels already carried by issues #1–#5 remain correct.

**Known uncovered seam:** the partition API and the hierarchical hyperprior
(`R/partition-api.R`, `R/partition.R`, `hyper_tau` / `class_rate_log_sd_z`,
`unlink=`, the `neoScale`/`transScale` normalisation, case 34 and
`dirichlet_simplex_class_w`) belong to no area. It needs its own row rather than
being folded into a neighbour — see area 12.

## Rotation

**The rule is not recorded here.** It lives in the `/red-team` skill (Normal
run, steps 1-2): the next area is the one whose most recent Discussion has the
oldest `createdAt`, an area with no Discussion taking precedence. The skill
carries the query too. Restating either here would create a second copy that
can drift from it -- which is what `last_focus:` was.

Two local facts the skill cannot know:

- The queue was normalised on 2026-09-19 by re-posting all twelve areas in
  sequence, so `createdAt` order is now area order: 7, 8, 9, 10, 11, 12, 1, 2,
  3, 4, 5, 6.
- **Discussions #118-#123 are ordering markers, not round records.** Each says
  so and names the real record. Delete one as soon as a genuine round record
  supersedes it.

`dev/red-team/log.md` and its `last_focus:` were retired the same day; the
pre-Discussions history is archived at discussion #124.

| # | Area | Files | Key questions |
|---|------|-------|---------------|
| 1 | k′ prior families & normalisation | `R/MkPrimeModel.R` (`LogPrior`, `.LogPriorEmpiricalGeometric`, `LogPemp`, `MkPrimeEmpiricalPrior`, all four `kPrimePrior` branches); `R/data.R`; `data-raw/empirical_n_obs.R`; `data/empiricalNObs.rda`; `src/mcmc.cpp` (`cpp_log_prior`); `tests/testthat/test-priors.R`, `test-empirical-geometric-prior.R`, `test-logseries-prior.R` | Covers **all four** priors — `geometric`, `empirical_geometric`, `beta_geometric`, `logseries` — not just EG, and their shared `LogPrior` spine. **The governing principle is that a prior is pre-data**: `kObs_i` is an observation and must not enter it. That resolved EG-001/EG-002 and retired LS-001. Does any surviving code path still condition on `kObs` — `priorVariant="conditional"` (the retained opt-in, called *Model B* in the older records; the pre-data default is *Model A*), the truncation normaliser, the relabelling correction? Is each prior normalised over its own support, and does the empirical body→tail join stay consistent (#2)? Do R and C++ agree to ~1e-9 on every arm, including `beta_geometric`, whose two hyperparameter columns have already caused an off-by-one elsewhere (STREAM-004)? Numerical stability in log space at the tails. Does `empiricalNObs` cover the range real datasets hit? Edge cases: `n_obs = 0`, `n_obs` above the table maximum, `NA` in `kObs` (#1), dead plumbing (#3). |
| 2 | k′ marginalisation & the p-samplers | `src/mcmc.cpp` (`compute_per_kprime_log_lik`, marginal_k branches, case 30 `mh_logit_p`, case 35 `gibbs_p_marginal`, the move-gating switch); `src/mcmc_likelihood.cpp`; `src/node_cl_cache.h`; `R/MkPrimeModel.R` (`likelihoodMode`); `R/RunMkPrime.R` (EG + marginal paths, `.BuildMoves`); `tests/testthat/test-marginal-k-*.R`; `dev/notes/2026-05-28-marginal-k-plan.md`; `dev/red-team/MARGINAL-K-*.md` | The `marginal_k` evaluator is the newest and highest-risk surface on `main`, and **`fixTopology = TRUE` SBC has already failed to catch two topology defects here** (FREEZE-003 Bug B: the evaluator read `state->parent/child` instead of its passed arguments) — so every question below must be asked under a **free** topology. **Two moves remain gated off** — `gibbs_spr` and `gibbs_subtree_swap` (`R/RunMkPrime.R:3570,3576`); the weighted/block family was re-enabled via scratch-eval in `c68f3dc`. Both remaining gates are necessary for a deeper reason than cache coherence: they select candidates at fixed kPrime and write `state->logLik` directly, so no cache fix makes them target-correct (verified 2026-09-18). Is the scratch-eval re-enable sound — `compute_full_loglik_at` does not force cold, and the `useCache` gate carries no topology or branch fingerprint. Is the `u_max(p)` truncation and its normaliser correct, and is the tail mass it discards bounded and stated? Is case 35 (`gibbs_p_marginal`, opt-in) invariant — derive it, do not assume; the claimed ~3.2× p-ESS win says nothing about correctness. Does case 30 still hold up (Jacobian, detailed balance, gating, #4, #5)? Do the v1 unsupported combinations (`empirical_geometric`/`beta_geometric`/`logseries`, `qHeterogeneity`, `usePartitioned`) actually abort at the boundary rather than silently mis-evaluating? |
| 3 | Logseries prior + `LogPrior` `treeLengthRate=NULL` fix | `R/MkPrimeModel.R` (logseries branch); `dev/plans/2026-05-12-logprior-treelengthrate-null.md`; `tests/testthat/test-logseries-prior.R` | Recent fix returns numeric(0) when treeLengthRate is NULL — does the fix cover all call sites? Is treeLengthRate NULL semantics consistent between R/cpp? |
| 4 | Hamilton harness: combine mode, recovery, ETAs | `data-raw/hamilton/run_one.R` (+97 lines in EG phase); `data-raw/hamilton/{mk,mkp,mkp_eg,combine}_array.slurm`; `data-raw/hamilton/mkp_check.slurm` | Heat 0.1 swap-cold logic correct? maxTime 10h with checkpoint resume — does state survive? Combine-mode CID computation against true tree handles unrooted/reordered correctly? Race conditions if two array tasks share scratch? |
| 5 | Move scheduler & adaptive weights | `R/MkPrimeMCMC.R`; `tests/testthat/test-m092-adaptive-scheduler.R` | After EG additions, are weights renormalised? Adaptive scheduler converges? m092 regression tests cover EG case? |
| 6 | Streaming, checkpoint, interrupt recovery | `R/streaming.R`; `R/RunMkPrime.R` (`ResumeMkPrime`, m175 interrupt handler, `.RunParallelRuns`); `R/burnin.R` | Stop-button handler (m175) leaves valid state? Resume from checkpoint reproduces continuation correctly? Trace file appends don't double-write the header on resume? Warmup ETA anchor (m177) doesn't slide back? **Known limitation (PAR-002, deferred):** parallel runs (`nCore > 1`) cannot checkpoint worker state — `ResumeMkPrime` restarts from scratch; documented in interrupt handler and `@param nCore`. Revisit if/when a user actually needs resumable parallel runs. |
| 7 | Tree-move proposals & Rcpp surface | `R/proposals.R`; `src/mcmc.cpp` (tree-move cases); `src/mcmc_state.h` | NNI/SPR Hastings ratios correct? Branch-length proposals respect bounds? Rcpp ↔ R state round-trip after EG additions? |
| 8 | Likelihood & partial CLs | `R/likelihood.R`; `src/mcmc_likelihood.cpp`; `src/mcmc.cpp` | Ascertainment correction under EG prior on k'? Partial-CL cache invalidation correct when k' changes per-character? Constant-site handling? |
| 9 | Convergence diagnostics & ESS | `R/Convergence.R`; `R/ess.R`; `R/treeESS.R`; `R/acrv.R` | Tree ESS computation correct for variable k'? Multi-chain convergence (R-hat) handles within-chain heat changes? Diagnostics emit on EG arm? |
| 10 | Test-suite health | `tests/testthat/test-empirical-geometric-prior.R`; `tests/testthat/test-priors.R`; `tests/testthat/test-gibbs.R`; `tests/testthat/test-m092-adaptive-scheduler.R` | New EG tests cover the math and the p-sampler? Snapshot tests stable? Tests that broke under EG default — were they fixed by re-baselining or by hiding a real regression? |
| 11 | **RB-oracle equivalence** | `dev/rb-equivalence/*`; `R/MkPrimeModel.R`; `R/MkPrimeMCMC.R`; `R/RunMkPrime.R`; `src/*` | Do MkPrime posteriors agree statistically with RevBayes on the 5 frozen matrices × 2 models (by_nt_9v, by_nt_kv)? Cross-sampler R-hat < 1.025 and pooled ESS > 128 on every scalar (tree_length, rate_log_sd, rate_loss, rate_neo) and on CID-to-pooled-median? Wall-seconds to target trend stable across commits? Tree ESS comparable (reported, not gating, since RB lacks the diagnostic)? Re-run on Hamilton (same node per cell) when any in-scope file changes — cache invalidation by `mtime(out/<rds>) < git log -1 --format=%cd R/ src/`. Pid 950 cid_to_median is the documented exception (near-prior posterior, degenerate tree CID). |
| 12 | **Red-team process meta-review** | `dev/red-team/focus-areas.md`, `dev/red-team/README.md`, `dev/red-team/discussion-categories.md`, `dev/red-team/findings-archive.md`, `dev/red-team/migration-map.tsv`; `gh issue list --label red-team` | Reviews the rotation itself, not the package. Are any areas **too broad** — spanning distinct seams, so a finder concentrating on one file family systematically misses another? Too **narrow** — a single-feature scope better merged into a neighbour? Do any **overlap**, auditing the same source under two headings? Has any gone persistently **dry** (three consecutive rounds, zero confirmed sev:med+) and should be retired, merged or down-tiered — remembering that dry is model-scoped, so check the model stamped in the round titles before calling it? Are there **new seams** not covered by any area: recently merged features, new `src/` or `R/` files, `dev/rb-equivalence/`, `dev/profiling/` (scope question — profiling has its own rotation; does red-team own its harnesses?), the Hamilton harness under `data-raw/hamilton/`? Are **tier assignments calibrated to recorded yield** in the round Discussions and in closed issues? Are areas 1–2 (tagged `[EG]`) still live, given the EG arm concluded and the released scope is MkNT? **Check Discussion-category headroom before proposing a new area** — GitHub caps a repo at 25, and each new area needs both an `area:N` label and a hand-created category (`discussion-categories.md`). Output is concrete restructuring actions — split, merge, retire, add, re-tier — each tied to yield evidence, not taste. |
