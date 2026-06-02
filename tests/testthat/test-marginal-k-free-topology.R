# Regression guard for MARGINAL-K-FREEZE-003 (free-topology freeze).
# Diagnosis: dev/red-team/MARGINAL-K-FREEZE-003-diagnosis.md.
#
# Bug B (FIXED, src/mcmc.cpp): under marginal_k a topology-changing MH move
# evaluated the OLD topology, because cpp_log_likelihood_marginal /
# compute_per_kprime_log_lik read state->parent/child instead of the passed
# parent/child. The proposal is evaluated BEFORE the topology is committed, so
# the committed state->logLik was wrong -> the continuous samplers froze (or
# sampled the wrong posterior) under FREE topology. The fixed-topology SBC
# could not catch this. This file asserts that after an accepted nni/spr move
# under marginal_k, the committed state->logLik equals a fresh COLD marginal
# recompute of the committed tree (gap ~ 0).
#
# Bug A (MITIGATED, R/RunMkPrime.R): gibbs_spr / gibbs_subtree_swap /
# weighted_spr / weighted_subtree_swap write a fixed-kPrime state->logLik (no
# marginal-over-k' sum) and are disabled under marginal_k in .BuildMoves. This
# file asserts they are absent from the marginal_k schedule but present (when
# enabled) under sampled_k.

library("ape")
library("TreeTools")

# --- shared 8-tip / 4-trans-char fixture (mirrors test-marginal-k-cache-*.R) --
.ft_make_tree <- function() {
  read.tree(text = paste0(
    "(((t1:0.05,t2:0.07):0.04,(t3:0.06,t4:0.05):0.03):0.05,",
    "((t5:0.04,t6:0.06):0.05,(t7:0.05,t8:0.04):0.06):0.04);"
  ))
}
.ft_make_mkd <- function() {
  set.seed(42L)
  tips <- paste0("t", 1:8)
  mat <- matrix(
    c(0, 0, 1, 1, 0, 1, 0, 1,
      0, 1, 1, 0, 1, 0, 0, 1,
      0, 1, 0, 0, 1, 1, 1, 0,
      1, 1, 0, 1, 0, 0, 1, 0),
    nrow = 8, ncol = 4, dimnames = list(tips, NULL)
  )
  MatrixToPhyDat(mat)
}
.ft_build <- function(p = 0.5) {
  tree  <- TreeTools::Preorder(.ft_make_tree())
  mkd   <- MkPrimeData(.ft_make_mkd())
  model <- MkPrime:::.FinalizeModel(
    MkPrimeModel(kPrimePrior = "geometric", likelihoodMode = "marginal_k",
                 coding = "variable"),
    tree, mkd)
  s0 <- MkPrime:::.InitState(tree, mkd, model)
  s0$p              <- p
  s0$tree_length    <- sum(tree$edge.length)
  s0$rel_br_lengths <- tree$edge.length / s0$tree_length
  dataPtr  <- MkPrime:::.InitMcmcData(mkd, model)
  statePtr <- MkPrime:::.InitMcmcChain(s0)
  fill_partition_cache(dataPtr, statePtr)
  list(dataPtr = dataPtr, statePtr = statePtr)
}
.ft_cold <- function(fx) {            # true marginal LL of the current state tree
  invalidate_marginal_cache(fx$statePtr)
  eval_full_loglik_cpp(fx$dataPtr, fx$statePtr)
}

test_that("marginal-k: accepted topology moves leave state->logLik coherent (FREEZE-003 Bug B)", {
  # The MH topology moves enabled under marginal_k -- nni(5), spr(6), tbr(17),
  # pspr(20) -- all change topology and route through the marginal evaluator.
  # Pre-fix these committed an OLD-topology likelihood (gap 1.3-4.8 nat). The
  # gap-sweep (marginal-k-freeze-repro.R) confirms all four read 0.000 post-fix;
  # this locks every enabled topology move into the always-on guard (tbr/pspr
  # were gap-swept manually but not previously regression-tested).
  for (mt in c(5L, 6L, 17L, 20L)) {
    fx <- .ft_build()
    accepted <- FALSE
    for (i in seq_len(2000L)) {
      if (isTRUE(do_move_cpp(fx$dataPtr, fx$statePtr, moveType = mt, charIdx = 0L,
                             scaleTuning = 0.5, betaSimplexTuning = 0.5,
                             intWalkWindow = 1L, beta = 1.0))) {
        accepted <- TRUE
        break
      }
    }
    expect_true(accepted, info = paste("moveType", mt, "should accept within 2000 tries"))
    committed <- get_state_log_lik(fx$statePtr)   # read BEFORE the cold recompute
    cold      <- .ft_cold(fx)
    expect_equal(committed, cold, tolerance = 1e-8,
                 info = paste("moveType", mt,
                              "committed state->logLik must equal a cold marginal recompute"))
  }
})

test_that("marginal-k: Gibbs/weighted topology moves are disabled in the schedule (FREEZE-003 Bug A)", {
  nEdge <- 13L   # 8-tip unrooted binary -> 2*8-3 edges
  nTrans <- 4L
  # suppressWarnings: tiny nIter clamps minWarmup->maxWarmup (warmup config,
  # irrelevant to the move-list assertions this test makes).
  mcmc <- suppressWarnings(MkPrimeMCMC(nIter = 10L, weightedSpr = TRUE,
                                       weightedSubtreeSwap = TRUE))
  banned <- c("gibbs_spr", "gibbs_subtree_swap", "weighted_spr", "weighted_subtree_swap")

  mv_marg <- MkPrime:::.BuildMoves(nEdge, nTrans, hasNeo = FALSE, mcmc,
                                   fixTopology = FALSE, likelihoodMode = "marginal_k")
  nm_marg <- vapply(mv_marg, function(m) m$name, character(1))
  for (g in banned)
    expect_false(g %in% nm_marg, info = paste(g, "must be absent under marginal_k"))
  # Topology is still searched by the marginal-correct MH moves (all four enabled).
  for (keep in c("nni", "spr", "tbr", "pspr"))
    expect_true(keep %in% nm_marg, info = paste(keep, "must remain under marginal_k"))

  # sampled_k keeps them (no regression to the legacy path).
  mv_samp <- MkPrime:::.BuildMoves(nEdge, nTrans, hasNeo = FALSE, mcmc,
                                   fixTopology = FALSE, likelihoodMode = "sampled_k")
  nm_samp <- vapply(mv_samp, function(m) m$name, character(1))
  for (g in banned)
    expect_true(g %in% nm_samp, info = paste(g, "must be present under sampled_k"))
})
