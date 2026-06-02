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
# Bug A + cache-coherence audit (R/RunMkPrime.R + src/mcmc.cpp):
#  * The four MULTI-EVAL moves (weighted_branch_scale / weighted_spr /
#    weighted_subtree_swap / block_gibbs_branch) select candidates via the
#    marginal evaluator (compute_full_loglik_at) with a full MH accept; their only
#    defect was a stale per-(char,k') charLLCache reused across the many configs
#    each evaluates. FIXED (FREEZE-003 follow-up) by the scratch-eval flag
#    (fillCharLLCache=false) on every intra-move eval, and RE-ENABLED under
#    marginal_k. Test 3 asserts each leaves a coherent committed state->logLik AND
#    a cache coherent with the committed state (warm==cold), on accept and reject.
#  * The GIBBS pair (gibbs_spr / gibbs_subtree_swap) select candidates via the
#    fixed-k' partial-CL path (no marginal sum, no P(u|p) weight), so a cache fix
#    alone does NOT make them target-correct: they stay DISABLED under marginal_k
#    (marginal-aware candidate eval deferred).
#  * The k'-sampling moves (kPrime / gibbs_kPrime / block_kPrime) stay gated out --
#    k' is integrated out under marginal_k.
# Test 2 asserts the gibbs pair + k'-moves are absent from the marginal_k
# schedule while the four weighted/block moves + nni/spr/tbr/pspr/mh_logit_p are
# present; all are present under sampled_k.

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

test_that("marginal-k: gibbs pair + k'-moves gated out; weighted/block moves re-enabled (FREEZE-003)", {
  nEdge <- 13L   # 8-tip unrooted binary -> 2*8-3 edges
  nTrans <- 4L
  # Enable every gateable move so the test would FAIL if the gating drifted.
  # suppressWarnings: tiny nIter clamps minWarmup.
  mcmc <- suppressWarnings(MkPrimeMCMC(
    nIter = 10L, gibbsSpr = TRUE, gibbsSubtreeSwap = TRUE,
    weightedSpr = TRUE, weightedSubtreeSwap = TRUE,
    weightedBranchScale = TRUE, blockGibbsBranch = TRUE))

  # (a) BANNED under marginal_k:
  #   gibbs pair -- select candidates via the fixed-k' partial-CL path (no
  #     marginal sum, no P(u|p) weight): a cache fix cannot make them
  #     target-correct, so they remain disabled (marginal-aware eval deferred).
  #   k'-moves   -- k' is analytically integrated out under marginal_k, so
  #     sampling it is semantically void; .BuildMoves substitutes mh_logit_p.
  banned <- c("gibbs_spr", "gibbs_subtree_swap",
              "kPrime", "gibbs_kPrime", "block_kPrime")
  # (b) PRESENT under marginal_k: the four weighted/block moves are now
  #     cache-coherent (scratch-eval fix, Test 3) and re-enabled; topology is
  #     also searched by nni/spr/tbr/pspr; p by mh_logit_p.
  present <- c("weighted_spr", "weighted_subtree_swap", "weighted_branch_lengths",
               "block_gibbs_branch", "nni", "spr", "tbr", "pspr", "mh_logit_p")

  # The dropped-move message still fires (gibbs pair was requested).
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
  for (keep in present)
    expect_true(keep %in% nm_marg, info = paste(keep, "must be present under marginal_k"))

  # sampled_k keeps every move (no regression to the legacy path).
  mv_samp <- MkPrime:::.BuildMoves(nEdge, nTrans, hasNeo = FALSE, mcmc,
                                   fixTopology = FALSE, likelihoodMode = "sampled_k")
  nm_samp <- vapply(mv_samp, function(m) m$name, character(1))
  for (g in c("gibbs_spr", "gibbs_subtree_swap", present[1:4]))
    expect_true(g %in% nm_samp, info = paste(g, "must be present under sampled_k"))
})

test_that("marginal-k: re-enabled weighted/block moves stay cache-coherent (FREEZE-003 scratch-eval)", {
  # The four moves re-enabled under marginal_k -- weighted_branch_scale(12),
  # weighted_spr(13), weighted_subtree_swap(14), block_gibbs_branch(15) -- each
  # evaluate many topology/branch configs per call. The scratch-eval fix
  # (compute_full_loglik_at(..., fillCharLLCache=false)) recomputes every config's
  # marginal LL cold. We fire each move many times (mix of accept + reject) and
  # after EACH fire assert two invariants, on BOTH paths:
  #   (i)  committed state->logLik == cold marginal recompute  (baseline_gap ~ 0)
  #   (ii) a warm eval == cold eval                            (warm_gap ~ 0)
  # (ii) is the cross-move post-condition: eval_full_loglik_cpp uses the
  # charLLCache exactly as the next mh_logit_p (case 30, no entry-invalidation)
  # would, so warm==cold proves no stale cache can poison it -- whether the move
  # left the cache cold (13/14/15) or ready-and-coherent (12, whose generic-MH
  # commit re-fills it for the accepted config). Pre-fix these gapped by
  # -1.26/-0.90/-2.26/-1.60 nat respectively.
  for (mt in c(12L, 13L, 14L, 15L)) {
    fx <- .ft_build()
    accepts <- 0L
    for (i in seq_len(300L)) {
      ok <- do_move_cpp(fx$dataPtr, fx$statePtr, moveType = mt, charIdx = 0L,
                        scaleTuning = 0.5, betaSimplexTuning = 0.5,
                        intWalkWindow = 1L, beta = 1.0)
      if (isTRUE(ok)) accepts <- accepts + 1L
      committed <- get_state_log_lik(fx$statePtr)                 # read FIRST
      warm      <- eval_full_loglik_cpp(fx$dataPtr, fx$statePtr)  # uses cache if ready
      invalidate_marginal_cache(fx$statePtr)
      cold      <- eval_full_loglik_cpp(fx$dataPtr, fx$statePtr)  # true marginal LL
      expect_equal(committed, cold, tolerance = 1e-8,
                   info = paste("moveType", mt, "committed==cold, iter", i))
      expect_equal(warm, cold, tolerance = 1e-8,
                   info = paste("moveType", mt, "warm==cold (no stale cache), iter", i))
    }
    expect_gt(accepts, 0L)   # ensure the accept path was actually exercised
  }
})

test_that("marginal-k: mh_logit_p after a weighted move reads a coherent cache (FREEZE-003)", {
  # The pointed cross-move hazard: case 30 (mh_logit_p) does NOT invalidate the
  # charLLCache at entry (it is the one move that reuses it). If a preceding
  # weighted move left a stale-ready cache, mh_p would consume it and commit a
  # wrong marginal LL. Fire [weighted move -> mh_logit_p] and assert mh_p's
  # committed LL equals a cold recompute at the resulting (topology, p).
  for (mt in c(12L, 13L, 14L, 15L)) {
    fx <- .ft_build()
    for (i in seq_len(60L)) {
      do_move_cpp(fx$dataPtr, fx$statePtr, moveType = mt, charIdx = 0L,
                  scaleTuning = 0.5, betaSimplexTuning = 0.5,
                  intWalkWindow = 1L, beta = 1.0)
      do_move_cpp(fx$dataPtr, fx$statePtr, moveType = 30L, charIdx = 0L,
                  scaleTuning = 0.5, betaSimplexTuning = 0.5,
                  intWalkWindow = 1L, beta = 1.0)
      committed <- get_state_log_lik(fx$statePtr)
      cold      <- .ft_cold(fx)
      expect_equal(committed, cold, tolerance = 1e-8,
                   info = paste("mh_logit_p after moveType", mt,
                                "must commit a coherent marginal LL, iter", i))
    }
  }
})

test_that("marginal-k: candidate (preorder_into) and commit (preorder_weighted) paths give identical LL (FREEZE-003 Phase 2 proposal-selection check)", {
  # The deterministic settler for proposal-SELECTION correctness (the dimension
  # committed==cold / warm==cold do NOT test). The weighted topology moves score
  # candidates via preorder_into -> compute_full_loglik_at (the selection weight)
  # but COMMIT via preorder_weighted_impl. If those two canonicalisations gave
  # different likelihoods for the same tree, candidate selection would be skewed
  # invisibly. eval_preorder_paths_cpp computes the marginal LL of the SAME tree
  # through BOTH paths (scratch evals); assert bit-equality across many evolving
  # topologies AND arbitrary input edge orders (the realistic mid-edit input).
  set.seed(11L)
  fx <- .ft_build()
  for (step in seq_len(90L)) {
    mt <- c(6L, 5L, 17L)[(step %% 3L) + 1L]   # spr / nni / tbr -> diverse topologies
    do_move_cpp(fx$dataPtr, fx$statePtr, moveType = mt, charIdx = 0L,
                scaleTuning = 0.5, betaSimplexTuning = 0.5,
                intWalkWindow = 1L, beta = 1.0)
    st  <- get_mcmc_state(fx$statePtr)
    par <- st$edge[, 1]; ch <- st$edge[, 2]
    el  <- st$treeLength * st$relBrLengths
    base <- eval_preorder_paths_cpp(fx$dataPtr, fx$statePtr, par, ch, el)
    expect_equal(unname(base[["preorder_into"]]), unname(base[["preorder_weighted"]]),
                 tolerance = 1e-9,
                 info = paste("step", step, "moveType", mt,
                              ": candidate-eval (preorder_into) and commit (preorder_weighted) LL must agree"))
    # Arbitrary input edge order (preorder_into rebuilds adjacency from scratch):
    # both canonicalisers must recover the same tree -> same LL.
    perm <- sample.int(length(par))
    pp <- eval_preorder_paths_cpp(fx$dataPtr, fx$statePtr, par[perm], ch[perm], el[perm])
    expect_equal(unname(pp[["preorder_into"]]), unname(base[["preorder_into"]]),
                 tolerance = 1e-9, info = paste("step", step, ": preorder_into edge-order invariance"))
    expect_equal(unname(pp[["preorder_weighted"]]), unname(base[["preorder_weighted"]]),
                 tolerance = 1e-9, info = paste("step", step, ": preorder_weighted edge-order invariance"))
  }
})
