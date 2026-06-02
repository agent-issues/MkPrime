#!/usr/bin/env Rscript
# marginal-k-freeze-repro.R
# =========================================================================
# Reproduce + LOCALIZE the marginal_k FREE-TOPOLOGY whole-chain freeze.
#
# FINDING (2026-06-02, T-OVL): under FREE topology, likelihoodMode="marginal_k"
# whole-chain FREEZES in ~7/16 datasets (sd=0 on tree_length / rate_log_sd / p,
# stuck at init, 0% MH+slice acceptance); sampled_k never freezes (0/16). The
# Stage-2 sampled_k SBC PASSED, but it runs fixTopology=TRUE, so it structurally
# cannot exercise the topology moves that trigger this. => the freeze is the v1
# marginal-k GATING bug, and a fixed-topology SBC alone cannot catch it.
#
# CODE-GROUNDED HYPOTHESIS (src/mcmc.cpp):
#   The Gibbs topology moves -- case 10 (gibbs_spr_impl) and case 11
#   (gibbs_subtree_swap_impl), both default-ON (gibbsSpr/gibbsSubtreeSwap=TRUE)
#   -- compute their candidate likelihoods by grouping transformational chars at
#   the SINGLE current state->kPrime (lines 1254-1282) through the sampled-k
#   partial-CL machinery: a FIXED-kPrime likelihood, with NO marginal-over-k'
#   logsumexp and NO geometric P(u|p) weight. They then write
#       state->logLik = candLL[chosen]            (line 1410 / 2132)
#   and RETURN early from the dispatch (4976 / 4979), BEFORE any marginal
#   re-evaluation. do_move_impl invalidates the per-(char,k') cache at its top
#   (line 4694) but NOTHING restores the correct MARGINAL baseline. The inflated
#   state->logLik (the same ~+10 nat INIT-001 inflation, now RECURRING on every
#   accepted topology move) becomes the logY0 / MH-baseline for the next
#   continuous-parameter move => every proposal looks worse => 0% acceptance =>
#   the continuous params freeze at (near) init.
#
# This driver tests that hypothesis directly:
#   PART B  (fast, decisive) -- a DISCRIMINATING move-type coherence-gap sweep:
#     fire each topology move until-accept on a frozen cell and measure
#       baseline_gap = state->logLik (what the move LEFT) - cold marginal LL.
#     PREDICT: gibbs moves (10,11) leave a large +ve gap (~nats); the MH moves
#     (5 nni, 6 spr, 20 pspr) and 30 (mh_p) leave gap ~= 0 because they route
#     through compute_full_loglik_at (the marginal evaluator).
#   PART B.2 -- dataset dependence: the gap is data-dependent (= fixed-k LL -
#     marginal LL), which is WHY only ~7/16 datasets freeze. Compare case-10 gap
#     on a frozen cell (r02) vs a non-frozen cell (r03).
#   PART B.3 -- causal close: corrupt state->logLik via one case-10 move, then
#     show mh_p (case 30) acceptance collapses to ~0 (vs >0 uncorrupted) on the
#     SAME data -- "wrong baseline => freeze", isolated to the single corruption.
#   PART A  (slower, end-to-end) -- RunMkPrime free topology, marginal_k vs
#     sampled_k on r02 + r03: sd(tl/rls/p) and per-move accept rates in situ.
#     Set MK_FREEZE_PARTA=0 to skip.
#
# Run:  Rscript dev/red-team/numerical/marginal-k-freeze-repro.R
# Scope: REPRO + LOCALIZE only. The fix has real forks (disable Gibbs under
#        marginal_k / make the Gibbs candidate eval marginal-aware /
#        recompute-the-baseline-after) -- reported, NOT chosen here.
# =========================================================================

suppressPackageStartupMessages({
  pkgload::load_all(".", quiet = TRUE)
  library("ape")
  library("TreeTools")
})

OUT_DIR <- "dev/red-team/numerical/marginal-k-freeze-results"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
sink_con <- file(file.path(OUT_DIR, "freeze-repro.log"), open = "wt")
sink(sink_con, split = TRUE)
on.exit({ sink(); close(sink_con) }, add = TRUE)

cat(sprintf("marginal-k-freeze-repro.R  |  %s\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

# ---- T-OVL dataset construction (verbatim, so cells match the frozen run) ----
simulate_dataset <- function(nTip, nChar, seed) {
  set.seed(seed)
  tree <- ape::rtree(nTip, br = function(n) runif(n, 0.02, 0.15))
  mat <- matrix(sample(0:1, nTip * nChar, replace = TRUE),
                nrow = nTip, ncol = nChar,
                dimnames = list(tree$tip.label, NULL))
  list(tree = tree, mkd = MkPrimeData(MatrixToPhyDat(mat)))
}
cell_seed <- function(nTip, nChar, rep) as.integer(1e3 * rep + 10 * nTip + nChar)

# ---- low-level marginal_k state builder (mirrors test-marginal-k-cache-*.R) --
build_ptrs <- function(sim, p = 0.5, coding = "variable") {
  tree  <- TreeTools::Preorder(sim$tree)
  mkd   <- sim$mkd
  model <- suppressMessages(MkPrimeModel(kPrimePrior    = "geometric",
                                         likelihoodMode = "marginal_k",
                                         coding         = coding))
  model  <- MkPrime:::.FinalizeModel(model, tree, mkd)
  state0 <- MkPrime:::.InitState(tree, mkd, model)
  state0$p              <- p
  state0$tree_length    <- sum(tree$edge.length)
  state0$rel_br_lengths <- tree$edge.length / state0$tree_length
  dataPtr  <- MkPrime:::.InitMcmcData(mkd, model)
  statePtr <- MkPrime:::.InitMcmcChain(state0)
  fill_partition_cache(dataPtr, statePtr)
  list(dataPtr = dataPtr, statePtr = statePtr, tree = tree, mkd = mkd)
}

# Coherence probe. Reads state->logLik (what the LAST move left as the MH/slice
# baseline) BEFORE touching the cache, then a warm eval, then a COLD eval (the
# true marginal LL). baseline_gap is the freeze signal; warm_gap (warm != cold)
# is the separately-tracked cache-option-a anomaly (root unpinned: cache-
# invalidation gap vs proposal/commit edge bookkeeping).
probe <- function(fx) {
  Lstate <- get_state_log_lik(fx$statePtr)                  # read FIRST
  Lwarm  <- eval_full_loglik_cpp(fx$dataPtr, fx$statePtr)   # warm (cache if ready)
  invalidate_marginal_cache(fx$statePtr)
  Lcold  <- eval_full_loglik_cpp(fx$dataPtr, fx$statePtr)   # cold = true marginal
  c(Lstate = Lstate, Lwarm = Lwarm, Lcold = Lcold,
    baseline_gap = Lstate - Lcold, warm_gap = Lwarm - Lcold)
}

fire_until_accept <- function(fx, moveType, maxTries = 2000L,
                              scaleTuning = 0.5, bsTuning = 0.5) {
  for (i in seq_len(maxTries)) {
    ok <- do_move_cpp(fx$dataPtr, fx$statePtr,
                      moveType          = as.integer(moveType),
                      charIdx           = 0L,
                      scaleTuning       = scaleTuning,
                      betaSimplexTuning = bsTuning,
                      intWalkWindow     = 1L,
                      beta              = 1.0)
    if (isTRUE(ok)) return(i)
  }
  NA_integer_
}

MOVE_NAMES <- c(`5`  = "nni (MH)",        `6`  = "spr (MH)",
                `17` = "tbr (MH)",        `20` = "pspr (MH)",
                `30` = "mh_p (MH)",
                `10` = "gibbs_spr",       `11` = "gibbs_subtree_swap",
                `12` = "weighted_branch_scale", `15` = "block_gibbs_branch",
                `13` = "weighted_spr",    `14` = "weighted_subtree_swap",
                `25` = "gibbs_kprime_sweep", `26` = "block_kprime_shift")
# FREEZE-003 cache-coherence audit: sweep EVERY move type that can corrupt the
# marginal charLLCache. The multi-eval moves that call compute_full_loglik_at in
# a topology/branch loop WITHOUT forcing the cache cold per-eval read a stale
# per-(char,k') cache after the first fill (no topology/branch fingerprint in the
# useCache gate) -> state==warm!=cold. Beyond the four gated topology moves, this
# adds 12 (weighted_branch_scale) + 15 (block_gibbs_branch) -- both UNGATED under
# marginal_k pre-this-audit -- plus the k'-moves 25/26 (k' is integrated out under
# marginal_k; expect coherent-or-no-op). MH moves (5/6/17/20) eval once -> gap 0.
SWEEP_ORDER <- c(30L, 5L, 6L, 20L, 17L,           # MH (single eval -> expect gap 0)
                 12L, 15L,                          # multi-eval BRANCH moves (audit)
                 10L, 11L, 13L, 14L,                # gibbs/weighted topology (gated)
                 25L, 26L)                          # k'-moves (integrated under marginal_k)

# =========================================================================
# PART B -- discriminating move-type coherence-gap sweep (the smoking gun)
# =========================================================================
cat("\n========================================================\n")
cat("PART B: move-type coherence-gap sweep (frozen cell n16_c16_r02)\n")
cat("========================================================\n")
seed_r02 <- cell_seed(16L, 16L, 2L)
sim_r02  <- simulate_dataset(16L, 16L, seed_r02)

# gotcha 2: what does state->kPrime hold under marginal_k (k' is integrated out)?
fx0 <- build_ptrs(sim_r02)
st0 <- get_mcmc_state(fx0$statePtr)
cat(sprintf("\nstate$kPrime under marginal_k (n=%d chars): %s ...\n",
            length(st0$kPrime), paste(head(st0$kPrime, 12), collapse = " ")))
cat(sprintf("  (range %d..%d)\n", min(st0$kPrime), max(st0$kPrime)))
b0 <- probe(fx0)
cat(sprintf("INIT coherence (post fill_partition_cache, INIT-001 fix in effect):\n"))
cat(sprintf("  state->logLik=%.4f  cold=%.4f  baseline_gap=%.4f  warm_gap=%.4g\n",
            b0["Lstate"], b0["Lcold"], b0["baseline_gap"], b0["warm_gap"]))

na_row <- function(mt) data.frame(
  moveType = mt, name = unname(MOVE_NAMES[as.character(mt)]),
  accepted_on = NA_integer_, gap_before = NA_real_, gap_after = NA_real_,
  warm_after = NA_real_, stringsAsFactors = FALSE)
sweep_rows <- lapply(SWEEP_ORDER, function(mt) {
  tryCatch({
    fx <- build_ptrs(sim_r02)              # FRESH state per move type (isolate)
    gB <- probe(fx)                        # baseline (must be ~0 after INIT-001)
    tries <- fire_until_accept(fx, mt)     # fire that move until it accepts
    gA <- if (is.na(tries)) c(baseline_gap = NA, warm_gap = NA) else probe(fx)
    data.frame(moveType    = mt,
               name        = unname(MOVE_NAMES[as.character(mt)]),
               accepted_on = tries,
               gap_before  = unname(gB["baseline_gap"]),
               gap_after   = unname(gA["baseline_gap"]),
               warm_after  = unname(gA["warm_gap"]),
               stringsAsFactors = FALSE)
  }, error = function(e) { message(sprintf("  [sweep] moveType %d threw: %s",
                                           mt, conditionMessage(e))); na_row(mt) })
})
sweep <- do.call(rbind, sweep_rows)
cat("\nbaseline_gap = state->logLik (LEFT by move) - cold marginal LL\n")
cat("warm_after   = warm eval - cold eval (warm!=cold = cache-option-a anomaly)\n\n")
print(sweep, row.names = FALSE, digits = 4)
saveRDS(sweep, file.path(OUT_DIR, "sweep-r02.rds"))

bad <- sweep[!is.na(sweep$gap_after) & abs(sweep$gap_after) > 0.05, ]
cat(sprintf("\n>> moves leaving |baseline_gap| > 0.05 nat: %s\n",
            if (nrow(bad)) paste(sprintf("%s(%.3f)", bad$name, bad$gap_after),
                                 collapse = ", ") else "NONE"))

# =========================================================================
# PART B.2 -- dataset dependence: case-10 gap, frozen (r02) vs non-frozen (r03)
# =========================================================================
cat("\n========================================================\n")
cat("PART B.2: case-10 (gibbs_spr) baseline_gap across cells\n")
cat("========================================================\n")
dep_rows <- lapply(1:4, function(rep) {
  seed <- cell_seed(16L, 16L, rep)
  sim  <- simulate_dataset(16L, 16L, seed)
  fx   <- build_ptrs(sim)
  tries <- fire_until_accept(fx, 10L)
  g <- if (is.na(tries)) c(baseline_gap = NA) else probe(fx)
  data.frame(cell = sprintf("n16_c16_r%02d", rep), seed = seed,
             gibbs_spr_accepted_on = tries,
             baseline_gap = unname(g["baseline_gap"]),
             stringsAsFactors = FALSE)
})
dep <- do.call(rbind, dep_rows)
print(dep, row.names = FALSE, digits = 4)
saveRDS(dep, file.path(OUT_DIR, "dataset-dependence.rds"))

# =========================================================================
# PART B.3 -- causal close: case-10 corruption => mh_p (case 30) freezes
# =========================================================================
cat("\n========================================================\n")
cat("PART B.3: causal close -- mh_p accept rate, uncorrupted vs case-10-corrupted\n")
cat("========================================================\n")
count_mh_p <- function(fx, n = 300L) {
  acc <- 0L
  for (i in seq_len(n))
    if (isTRUE(do_move_cpp(fx$dataPtr, fx$statePtr, 30L, 0L, 0.5, 0.5, 1L, 1.0)))
      acc <- acc + 1L
  acc / n
}
fxU <- build_ptrs(sim_r02)                         # uncorrupted
rateU <- count_mh_p(fxU)
fxC <- build_ptrs(sim_r02)                         # corrupt via ONE gibbs_spr
tC  <- fire_until_accept(fxC, 10L)
gC  <- probe(fxC)["baseline_gap"]
# probe() invalidated the cache + did not change state->logLik; mh_p path rebuilds
rateC <- count_mh_p(fxC)
cat(sprintf("\nuncorrupted state->logLik baseline -> mh_p accept rate = %.3f\n", rateU))
cat(sprintf("after ONE gibbs_spr (accepted on try %s, baseline_gap=%.3f nat)\n",
            as.character(tC), gC))
cat(sprintf("   corrupted state->logLik baseline -> mh_p accept rate = %.3f\n", rateC))
cat(sprintf(">> %s\n", if (rateU > 0.02 && rateC < 0.01)
            "CAUSAL CHAIN CONFIRMED: one fixed-kPrime topology move freezes mh_p"
            else "INCONCLUSIVE -- inspect rates above"))

# =========================================================================
# PART B.4 -- p-gating: is the Bug-A gap the omitted geometric weight,
#            gap(p) ~= -n_char * log(p)?  (advisor's discriminating test)
# If yes, the freeze is gated by p (gap large at low p, ->0 as p->1), which
# explains the seed/data-dependent 7/16 incidence WITHOUT the (circular)
# "a slice rescues the baseline" story.
# =========================================================================
cat("\n========================================================\n")
cat("PART B.4: p-gating -- gibbs_spr gap vs -n_char*log(p)\n")
cat("========================================================\n")
nChar <- length(st0$kPrime)
pgate_rows <- lapply(c(0.3, 0.5, 0.7, 0.9, 0.97), function(pp) {
  fxC <- build_ptrs(sim_r02, p = pp)
  tC  <- fire_until_accept(fxC, 10L)
  gap <- if (is.na(tC)) NA else unname(probe(fxC)["baseline_gap"])
  rateC <- { fxD <- build_ptrs(sim_r02, p = pp); fire_until_accept(fxD, 10L)
             count_mh_p(fxD) }
  data.frame(p = pp,
             gibbs_gap        = gap,
             predicted_neglogp = -nChar * log(pp),
             mh_p_accept_after_corrupt = rateC,
             stringsAsFactors = FALSE)
})
pgate <- do.call(rbind, pgate_rows)
cat(sprintf("(n_char = %d; predicted gap = -n_char*log(p) = missing geometric weight)\n\n", nChar))
print(pgate, row.names = FALSE, digits = 4)
saveRDS(pgate, file.path(OUT_DIR, "p-gating.rds"))
trk <- with(pgate, cor(gibbs_gap, predicted_neglogp))
cat(sprintf("\n>> gap-vs-(-n_char*log p) correlation = %.4f %s\n", trk,
            if (!is.na(trk) && trk > 0.99)
              "=> Bug-A gap IS the omitted geometric weight; freeze is p-gated"
            else "=> relationship weaker than expected; inspect"))

# =========================================================================
# PART A -- end-to-end RunMkPrime (free topology): marginal_k vs sampled_k
# =========================================================================
if (!identical(Sys.getenv("MK_FREEZE_PARTA", unset = "1"), "0")) {
  cat("\n========================================================\n")
  cat("PART A: end-to-end RunMkPrime free topology (sd + accept in situ)\n")
  cat("========================================================\n")
  N_ITER <- as.integer(Sys.getenv("MK_FREEZE_NITER", unset = "4000"))
  N_WARM <- as.integer(Sys.getenv("MK_FREEZE_NWARM", unset = "1000"))

  run_one <- function(sim, mode, seed) {
    model <- suppressMessages(MkPrimeModel(kPrimePrior = "geometric",
                                           likelihoodMode = mode,
                                           coding = "variable"))
    mcmc <- MkPrimeMCMC(nIter = N_ITER, thin = "auto",
                        minWarmup = N_WARM, maxWarmup = N_WARM,
                        autoTune = FALSE, nRuns = 1L, nChains = 1L)
    set.seed(seed)
    suppressMessages(suppressWarnings(
      RunMkPrime(sim$mkd, sim$tree, model = model, mcmc = mcmc, overwrite = TRUE)))
  }
  sd3 <- function(fit) c(
    tl  = sd(fit$samples[, "tree_length"]),
    rls = sd(fit$samples[, "rate_log_sd"]),
    p   = sd(fit$samples[, "p"]))

  partA <- list()
  for (rep in c(2L, 3L)) {
    seed <- cell_seed(16L, 16L, rep)
    sim  <- simulate_dataset(16L, 16L, seed)
    fitM <- run_one(sim, "marginal_k", seed)
    fitS <- run_one(sim, "sampled_k",  seed)
    cell <- sprintf("n16_c16_r%02d", rep)
    sm <- sd3(fitM); ss <- sd3(fitS)
    cat(sprintf("\n[%s] marginal_k sd(tl/rls/p) = %.4f / %.4f / %.4f %s\n",
                cell, sm["tl"], sm["rls"], sm["p"],
                if (sm["rls"] < 1e-3 || sm["tl"] < 1e-3) "<-FROZEN" else ""))
    cat(sprintf("[%s] sampled_k  sd(tl/rls/p) = %.4f / %.4f / %.4f %s\n",
                cell, ss["tl"], ss["rls"], ss["p"],
                if (ss["rls"] < 1e-3 || ss["tl"] < 1e-3) "<-FROZEN" else ""))
    acc <- tryCatch(fitM$acceptance, error = function(e) NULL)
    if (!is.null(acc)) {
      cat(sprintf("[%s] marginal_k per-move accept:\n", cell))
      print(round(unlist(acc), 4))
    }
    partA[[cell]] <- list(marginal_sd = sm, sampled_sd = ss,
                          marginal_accept = acc,
                          marginal_fit = fitM, sampled_fit = fitS)
  }
  saveRDS(partA, file.path(OUT_DIR, "partA-runmkprime.rds"))
} else {
  cat("\n[PART A skipped: MK_FREEZE_PARTA=0]\n")
}

cat("\nDONE.\n")
