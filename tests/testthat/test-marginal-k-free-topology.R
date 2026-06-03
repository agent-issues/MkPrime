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

# ---------------------------------------------------------------------------
# gibbs_p_marginal (case 35): data-augmentation Metropolis-within-Gibbs p-update
# for marginal_k. Derivation + math-prover verification + numerical invariance:
#   dev/red-team/proofs/marginal-k-gibbs-p.md
#   dev/red-team/numerical/gibbs-p-identity-check.R
# Impute u_i from the cached per-(char,k') weights, propose p* from the
# untruncated conjugate Beta(a+n, b+Sigma u + c_A), accept with the truncation-
# normaliser ratio Sigma_i[logZ_i(p) - logZ_i(p*)]. p-only move: preserves the
# charLL cache (like case 30). Tested for BOTH Model A (unconditional) and Model
# B (conditional).
# ---------------------------------------------------------------------------

# Multistate fixture: two kObs=3 characters (states {0,1,2}) + two binary, so
# Model A's c_A = Sum(kObs_i - 2) = 2 > 0. REQUIRED to exercise the c_A branch of
# case 35 -- the shared binary fixture has all kObs=2 => c_A=0 => Model A == Model
# B, leaving the most error-prone branch (the one disabled case-9 got wrong)
# untested. Mirrors the 8-tip layout; do NOT fold into .ft_make_mkd (other tests
# depend on its exact binary content).
.gp_make_mkd_ms <- function() {
  set.seed(43L); tips <- paste0("t", 1:8)
  m <- matrix(c(
    0, 0, 1, 1, 2, 2, 0, 1,    # kObs = 3
    0, 1, 2, 0, 1, 2, 0, 1,    # kObs = 3
    0, 0, 1, 1, 0, 1, 0, 1,    # kObs = 2
    1, 1, 0, 1, 0, 0, 1, 0),   # kObs = 2
    nrow = 8, ncol = 4, dimnames = list(tips, NULL))
  TreeTools::MatrixToPhyDat(m)
}
.gp_build <- function(p = 0.5, variant = "conditional", mkdFn = .ft_make_mkd) {
  tree  <- TreeTools::Preorder(.ft_make_tree())
  mkd   <- MkPrimeData(mkdFn())
  model <- MkPrime:::.FinalizeModel(
    MkPrimeModel(kPrimePrior = "geometric", likelihoodMode = "marginal_k",
                 priorVariant = variant, coding = "variable"),
    tree, mkd)
  s0 <- MkPrime:::.InitState(tree, mkd, model)
  s0$p              <- p
  s0$tree_length    <- sum(tree$edge.length)
  s0$rel_br_lengths <- tree$edge.length / s0$tree_length
  dataPtr  <- MkPrime:::.InitMcmcData(mkd, model)
  statePtr <- MkPrime:::.InitMcmcChain(s0)
  fill_partition_cache(dataPtr, statePtr)
  # Prime the per-(char,k') charLL cache COLD at this p: case 35 requires a warm
  # cache and (unlike case 30) never fills it itself, so an unprimed cache would
  # make it a no-op. The cold prime also gives the correct marginal LL at p.
  invisible(eval_full_loglik_cpp(dataPtr, statePtr))
  list(dataPtr = dataPtr, statePtr = statePtr,
       a = model$kprimeHyperA, b = model$kprimeHyperB)
}

# Fire one p-only move repeatedly on a FIXED tree/rates fixture; record p (thinned)
# and the acceptance fraction. The cache stays warm throughout (p-only moves
# preserve it), so the chain explores pi(p | theta, tree) exactly.
.gp_run_pmove <- function(fx, moveType, nIter, scaleTuning = 1.0, thin = 5L) {
  nRec <- nIter %/% thin
  ps <- numeric(nRec); j <- 0L; acc <- 0L
  for (i in seq_len(nIter)) {
    ok <- do_move_cpp(fx$dataPtr, fx$statePtr, moveType = moveType, charIdx = 0L,
                      scaleTuning = scaleTuning, betaSimplexTuning = 0.5,
                      intWalkWindow = 1L, beta = 1.0)
    if (isTRUE(ok)) acc <- acc + 1L
    if (i %% thin == 0L) { j <- j + 1L; ps[j] <- get_mcmc_state(fx$statePtr)$p }
  }
  list(p = ps, accept = acc / nIter)
}

test_that("marginal-k gibbs_p_marginal (case 35): committed==cold & warm==cold (cache coherence)", {
  # p-only move: on accept it re-marginalises at p* via the cache fast-path and
  # leaves the cache ready; on reject nothing changes. After EACH fire the
  # committed state->logLik must equal a cold recompute, and a warm eval must
  # equal a cold eval (no stale cache for a subsequent mh_logit_p to consume).
  for (variant in c("conditional", "unconditional")) {
    fx <- .gp_build(p = 0.5, variant = variant)
    accepts <- 0L
    for (i in seq_len(300L)) {
      ok <- do_move_cpp(fx$dataPtr, fx$statePtr, moveType = 35L, charIdx = 0L,
                        scaleTuning = 0.5, betaSimplexTuning = 0.5,
                        intWalkWindow = 1L, beta = 1.0)
      if (isTRUE(ok)) accepts <- accepts + 1L
      committed <- get_state_log_lik(fx$statePtr)
      warm      <- eval_full_loglik_cpp(fx$dataPtr, fx$statePtr)
      invalidate_marginal_cache(fx$statePtr)
      cold      <- eval_full_loglik_cpp(fx$dataPtr, fx$statePtr)
      # Re-prime the cache for the next iteration (the invalidate above cleared it).
      expect_equal(committed, cold, tolerance = 1e-8,
                   info = paste(variant, "committed==cold, iter", i))
      expect_equal(warm, cold, tolerance = 1e-8,
                   info = paste(variant, "warm==cold, iter", i))
    }
    # K=200 default: Z ~ 1, so the move is near-pure Gibbs and accepts readily.
    expect_gt(accepts, 150L)
  }
})

# Ground-truth oracle: the marginal_k posterior on p at the FIXED fixture tree is
# pi(p|theta,tree) ∝ exp(marginal_LL_cold(p)) · Beta(p; a, b). Tabulate it by COLD
# marginal evals on a p-grid (each .gp_build cold-primes the marginal at its p),
# add the Beta prior, normalise, take E[p]. (A pure mh_logit_p chain is NOT a valid
# oracle in isolation: with no cache-refilling move interleaved its warm cache
# sticks at the prime-p support and undercounts as p wanders -- the very drift
# case 35 fixes via force-cold. The grid oracle sidesteps both moves.)
.gp_analytic_Ep <- function(variant, mkdFn) {
  grid <- seq(0.02, 0.98, by = 0.02)
  a <- NULL; b <- NULL
  logpost <- vapply(grid, function(pg) {
    fxg <- .gp_build(p = pg, variant = variant, mkdFn = mkdFn)
    a <<- fxg$a; b <<- fxg$b
    eval_full_loglik_cpp(fxg$dataPtr, fxg$statePtr) + stats::dbeta(pg, a, b, log = TRUE)
  }, numeric(1L))
  w <- exp(logpost - max(logpost)); w <- w / sum(w)   # uniform grid -> spacing cancels
  sum(grid * w)
}

test_that("marginal-k gibbs_p_marginal (case 35): p-target matches analytic (binary fixture, c_A=0)", {
  # Binary fixture -> all kObs=2 -> c_A=0 -> Model A == Model B. The case-35 chain
  # (force-cold-on-commit -> correct target, even under large p-jumps) must match
  # the grid oracle within MC error. A flipped Z-sign biases E[p] by ~0.2
  # (gibbs-p-identity-check.R controls), far outside the tolerance.
  for (variant in c("conditional", "unconditional")) {
    Ep_analytic <- .gp_analytic_Ep(variant, .ft_make_mkd)
    set.seed(2026L)
    fx <- .gp_build(p = 0.5, variant = variant)
    r  <- .gp_run_pmove(fx, moveType = 35L, nIter = 40000L)
    Ep_emp <- mean(r$p[-seq_len(length(r$p) %/% 5L)])
    expect_equal(Ep_emp, Ep_analytic, tolerance = 0.02,
                 info = paste("binary", variant, ": case-35 E[p]=", round(Ep_emp, 4),
                              "vs analytic E[p]=", round(Ep_analytic, 4)))
    expect_gt(r$accept, 0.5)
  }
})

test_that("marginal-k gibbs_p_marginal (case 35): Model-A c_A shift is exercised (multistate fixture)", {
  # Multistate fixture -> c_A = Sum(kObs_i - 2) = 2 > 0, so Model A's proposal
  # shape2 = b + Sum(u) + c_A differs from Model B's b + Sum(u). This is the ONLY
  # test that executes the c_A branch of case 35 (the disabled case-9 bug).
  Ep_A <- .gp_analytic_Ep("unconditional", .gp_make_mkd_ms)   # Model A (c_A>0)
  Ep_B <- .gp_analytic_Ep("conditional",   .gp_make_mkd_ms)   # Model B (c_A=0)
  # Sanity: c_A actually shifts the posterior, so the test is SENSITIVE to it.
  # (If A==B here the comparison below could not detect a dropped c_A.)
  expect_gt(abs(Ep_A - Ep_B), 0.01,
            label = sprintf("Model-A vs Model-B E[p] gap (A=%.4f B=%.4f)", Ep_A, Ep_B))
  # case-35 chain under Model A must match the Model-A analytic (NOT the Model-B
  # one): a dropped c_A would make the chain track Ep_B instead -> caught.
  set.seed(2027L)
  fxA <- .gp_build(p = 0.5, variant = "unconditional", mkdFn = .gp_make_mkd_ms)
  rA  <- .gp_run_pmove(fxA, moveType = 35L, nIter = 40000L)
  EpA_emp <- mean(rA$p[-seq_len(length(rA$p) %/% 5L)])
  expect_equal(EpA_emp, Ep_A, tolerance = 0.02,
               info = sprintf("Model A: case-35 E[p]=%.4f vs analytic A=%.4f (B=%.4f)",
                              EpA_emp, Ep_A, Ep_B))
  # And Model B matches the Model-B analytic.
  set.seed(2028L)
  fxB <- .gp_build(p = 0.5, variant = "conditional", mkdFn = .gp_make_mkd_ms)
  rB  <- .gp_run_pmove(fxB, moveType = 35L, nIter = 40000L)
  EpB_emp <- mean(rB$p[-seq_len(length(rB$p) %/% 5L)])
  expect_equal(EpB_emp, Ep_B, tolerance = 0.02,
               info = sprintf("Model B: case-35 E[p]=%.4f vs analytic B=%.4f", EpB_emp, Ep_B))
})

test_that("marginal-k gibbs_p_marginal is opt-in (default off; present only when enabled)", {
  # Conservative scheduling: the move is correctness-verified but unproven on
  # ESS/second, so production is UNCHANGED by default (mh_logit_p only). The
  # `gibbsPMarginal=TRUE` flag adds it (low pinned weight) for measurement.
  off <- suppressWarnings(MkPrimeMCMC(nIter = 10L))
  nm_off <- vapply(suppressMessages(MkPrime:::.BuildMoves(
    13L, 4L, hasNeo = FALSE, off, fixTopology = FALSE, likelihoodMode = "marginal_k")),
    function(m) m$name, character(1))
  expect_false("gibbs_p_marginal" %in% nm_off)
  expect_true("mh_logit_p" %in% nm_off)   # primary p-move always present

  on <- suppressWarnings(MkPrimeMCMC(nIter = 10L, gibbsPMarginal = TRUE))
  nm_on <- vapply(suppressMessages(MkPrime:::.BuildMoves(
    13L, 4L, hasNeo = FALSE, on, fixTopology = FALSE, likelihoodMode = "marginal_k")),
    function(m) m$name, character(1))
  expect_true("gibbs_p_marginal" %in% nm_on)
  expect_true("mh_logit_p" %in% nm_on)    # kept alongside, as primary

  # No-op under sampled_k (marginal evaluator + case 35 are geometric marginal_k only).
  nm_samp <- vapply(MkPrime:::.BuildMoves(
    13L, 4L, hasNeo = FALSE, on, fixTopology = FALSE, likelihoodMode = "sampled_k"),
    function(m) m$name, character(1))
  expect_false("gibbs_p_marginal" %in% nm_samp)
})
