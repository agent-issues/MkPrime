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
# Bug A + cache-coherence audit (MITIGATED, R/RunMkPrime.R): six moves write an
# incoherent committed state->logLik under marginal_k and are disabled in
# .BuildMoves -- the gibbs pair (gibbs_spr / gibbs_subtree_swap, fixed-k'
# likelihood) and the four multi-eval moves (weighted_spr / weighted_subtree_swap
# / weighted_branch_scale / block_gibbs_branch, which read a stale marginal
# charLLCache because compute_full_loglik_at does not force the cache cold between
# topology/branch evals). The k'-sampling moves (kPrime / gibbs_kPrime /
# block_kPrime) are also gated out -- k' is integrated out under marginal_k.
# Test 2 asserts all are absent from the marginal_k schedule but present under
# sampled_k.

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

test_that("marginal-k: every cache-incoherent / k'-sampling move is gated out of the schedule (FREEZE-003)", {
  nEdge <- 13L   # 8-tip unrooted binary -> 2*8-3 edges
  nTrans <- 4L
  # Enable EVERY move that must be gated off under marginal_k, so this test would
  # FAIL if any were left ungated. suppressWarnings: tiny nIter clamps minWarmup.
  mcmc <- suppressWarnings(MkPrimeMCMC(
    nIter = 10L, gibbsSpr = TRUE, gibbsSubtreeSwap = TRUE,
    weightedSpr = TRUE, weightedSubtreeSwap = TRUE,
    weightedBranchScale = TRUE, blockGibbsBranch = TRUE))

  # (a) Multi-eval moves that write an incoherent committed state->logLik under
  #     marginal_k. They call compute_full_loglik_at across topology/branch
  #     configs without forcing the marginal charLLCache cold per-eval, so they
  #     read a stale per-(char,k') cache (state==warm!=cold) -- or, for the gibbs
  #     pair, write a fixed-k' likelihood (state!=warm==cold). Coherence sweep
  #     (marginal-k-freeze-repro.R) gaps: gibbs_spr +10.7, gibbs_subtree +10.5,
  #     weighted_branch_scale -1.26, block_gibbs_branch -1.60, weighted_spr -0.90,
  #     weighted_subtree_swap -2.26 nat.
  corruptors <- c("gibbs_spr", "gibbs_subtree_swap", "weighted_spr",
                  "weighted_subtree_swap", "weighted_branch_lengths",
                  "block_gibbs_branch")
  # (b) k'-sampling moves: k' is analytically integrated out under marginal_k, so
  #     sampling it is semantically void AND fires the fixed-k' Bug-A path (sweep:
  #     gibbs_kprime_sweep gap = +6.85 nat if fired). .BuildMoves replaces them
  #     with mh_logit_p under marginal_k; this asserts that exclusion is live.
  kprime_moves <- c("kPrime", "gibbs_kPrime", "block_kPrime")
  banned <- c(corruptors, kprime_moves)

  expect_message(
    MkPrime:::.BuildMoves(nEdge, nTrans, hasNeo = FALSE, mcmc,
                          fixTopology = FALSE, likelihoodMode = "marginal_k"),
    regexp = "not used under"
  )
  mv_marg <- suppressMessages(
    MkPrime:::.BuildMoves(nEdge, nTrans, hasNeo = FALSE, mcmc,
                          fixTopology = FALSE, likelihoodMode = "marginal_k")
  )
  nm_marg <- vapply(mv_marg, function(m) m$name, character(1))
  for (g in banned)
    expect_false(g %in% nm_marg, info = paste(g, "must be absent under marginal_k"))
  # Topology is still searched by the marginal-correct MH moves; p by mh_logit_p.
  for (keep in c("nni", "spr", "tbr", "pspr", "mh_logit_p"))
    expect_true(keep %in% nm_marg, info = paste(keep, "must remain under marginal_k"))

  # sampled_k keeps them all (no regression to the legacy path).
  mv_samp <- MkPrime:::.BuildMoves(nEdge, nTrans, hasNeo = FALSE, mcmc,
                                   fixTopology = FALSE, likelihoodMode = "sampled_k")
  nm_samp <- vapply(mv_samp, function(m) m$name, character(1))
  for (g in banned)
    expect_true(g %in% nm_samp, info = paste(g, "must be present under sampled_k"))
})
