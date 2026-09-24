# Tests for Layer 1 per-class MCMC moves:
#   case 31: scale_class_rate_log_sd (Bactrian MH on classRateLogSd[c])
#   case 32: dirichlet_simplex_class_w (Dirichlet simplex on class_w)
#
# IMPORTANT: cpp_log_prior does not yet include the per-class Gamma (shape)
# or Dirichlet (ratemultiplier) prior contributions — that is the parallel
# agent's work. Acceptance ratios are therefore driven by LL-ratio only.
# These tests assert:
#   (a) chains run to completion with all-finite LLs (moves fire without error),
#   (b) per-class columns are present in result$samples,
#   (c) per-class values evolve across samples (moves are accepted sometimes),
#   (d) simplex constraint: classW columns sum to ~1 every sample.
#
# Precise posterior correctness is NOT asserted here; the integration commit
# (third agent) will add the prior terms and a numeric-accuracy test.

# ---- shared data setup -------------------------------------------------------

.setup_2class_data <- function(seed = 77L, nChar = 10L, nTip = 6L) {
  set.seed(seed)
  mat <- matrix(
    sample(0:1, nTip * nChar, replace = TRUE),
    nrow = nTip, ncol = nChar,
    dimnames = list(paste0("t", seq_len(nTip)), NULL)
  )
  # Ensure every character is variable
  for (j in seq_len(nChar)) {
    if (length(unique(mat[, j])) < 2L) mat[1L, j] <- 1L - mat[1L, j]
  }
  pd  <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)
  tree <- Preorder(
    NJTree(pd, edgeLengths = TRUE)
  )
  if (is.null(tree$edge.length) || any(tree$edge.length <= 0)) {
    tree$edge.length <- pmax(
      tree$edge.length %||% rep(0.1, nrow(tree$edge)), 1e-8
    )
  }
  # 2-class partition: first half class 1, second half class 2
  n  <- mkd$nChar
  part <- c(rep(1L, floor(n / 2)), rep(2L, n - floor(n / 2)))
  list(mkd = mkd, tree = tree, part = part)
}

.short_mcmc <- function() {
  MkPrimeMCMC(
    nIter       = 100L,
    maxWarmup   = 50L,
    minWarmup   = 50L,
    nChains     = 1L,
    thin        = 1L,
    autoTune    = FALSE,
    gibbsSubtreeSwap = FALSE   # avoid non-determinism from extra topology moves
  )
}

# ==============================================================================
# 1. unlink = "shape"
# ==============================================================================

test_that("unlink=shape chain runs, per-class shape columns present and vary", {
  d <- .setup_2class_data()
  set.seed(42L)
  result <- allow_warning(
    RunMkPrime(
      data      = d$mkd,
      tree      = d$tree,
      mcmc      = .short_mcmc(),
      partition = d$part,
      unlink    = "shape"
    ),
    "without stabilisation"
  )

  samp <- result$samples

  # (a) All LLs finite
  expect_true(all(is.finite(samp[, "log_likelihood"])))
  expect_gt(result$nSamples, 0L)

  # (b) Per-class columns present
  expect_true("class1_rate_log_sd" %in% colnames(samp))
  expect_true("class2_rate_log_sd" %in% colnames(samp))

  # (c) Values are positive (rateLogSd is always positive)
  expect_true(all(samp[, "class1_rate_log_sd"] > 0))
  expect_true(all(samp[, "class2_rate_log_sd"] > 0))

  # (d) Values evolve: not all identical (moves were accepted at least sometimes)
  # Use a relaxed test — at least one unique value beyond the starting value.
  expect_gt(length(unique(round(samp[, "class1_rate_log_sd"], 8))), 1L)
  expect_gt(length(unique(round(samp[, "class2_rate_log_sd"], 8))), 1L)
})

# ==============================================================================
# 2. unlink = "ratemultiplier"
# ==============================================================================

test_that("unlink=ratemultiplier chain runs, w columns present on simplex", {
  d <- .setup_2class_data()
  set.seed(43L)
  result <- allow_warning(
    RunMkPrime(
      data      = d$mkd,
      tree      = d$tree,
      mcmc      = .short_mcmc(),
      partition = d$part,
      unlink    = "ratemultiplier"
    ),
    "without stabilisation"
  )

  samp <- result$samples

  # (a) All LLs finite
  expect_true(all(is.finite(samp[, "log_likelihood"])))
  expect_gt(result$nSamples, 0L)

  # (b) Per-class columns present
  expect_true("w_1" %in% colnames(samp))
  expect_true("w_2" %in% colnames(samp))

  # (c) Simplex constraint: w_1 + w_2 ~= 1 every sample
  w_sum <- samp[, "w_1"] + samp[, "w_2"]
  expect_true(all(abs(w_sum - 1.0) < 1e-8))

  # (d) Weights are in (0, 1)
  expect_true(all(samp[, "w_1"] > 0 & samp[, "w_1"] < 1))
  expect_true(all(samp[, "w_2"] > 0 & samp[, "w_2"] < 1))

  # (e) Values evolve (at least some acceptance)
  expect_gt(length(unique(round(samp[, "w_1"], 8))), 1L)
})

# ==============================================================================
# 3. unlink = c("shape", "ratemultiplier")
# ==============================================================================

test_that("unlink=c(shape,ratemultiplier) chain runs, all per-class columns present", {
  d <- .setup_2class_data()
  set.seed(44L)
  result <- allow_warning(
    RunMkPrime(
      data      = d$mkd,
      tree      = d$tree,
      mcmc      = .short_mcmc(),
      partition = d$part,
      unlink    = c("shape", "ratemultiplier")
    ),
    "without stabilisation"
  )

  samp <- result$samples

  # (a) All LLs finite
  expect_true(all(is.finite(samp[, "log_likelihood"])))

  # (b) All per-class columns present
  for (col in c("class1_rate_log_sd", "class2_rate_log_sd", "w_1", "w_2")) {
    expect_true(col %in% colnames(samp), label = paste("column", col, "present"))
  }

  # (c) Simplex constraint holds
  w_sum <- samp[, "w_1"] + samp[, "w_2"]
  expect_true(all(abs(w_sum - 1.0) < 1e-8))

  # (d) Per-class shape values are positive
  expect_true(all(samp[, "class1_rate_log_sd"] > 0))
  expect_true(all(samp[, "class2_rate_log_sd"] > 0))
})

# ==============================================================================
# 4. class1_rate_log_sd random-walk regression
# ==============================================================================

# Guards the lockstep invariant `classRateLogSd[0] == rateLogSd`. The
# partitioned prior in mcmc.cpp explicitly skips the c==0 Gamma term on the
# assumption this holds (init_mcmc_state sets it; the move sites in
# do_move_impl maintain it). If that invariant is dropped — as it was when
# Layer 1 landed — the c==0 slot becomes unconstrained, a symmetric
# log-scale Bactrian proposal random-walks it upward unboundedly while the
# likelihood saturates under ACRV (rates → 0/∞), leaving log_likelihood
# bounded and the corruption invisible until the trace explodes to ~1e+46
# at production scale.
test_that("class1_rate_log_sd stays in lockstep with rate_log_sd (no slot-0 drift)", {
  d <- .setup_2class_data(seed = 91L, nChar = 12L, nTip = 8L)
  set.seed(91L)
  result <- allow_warning(
    RunMkPrime(
      data      = d$mkd,
      tree      = d$tree,
      mcmc      = MkPrimeMCMC(
        nIter            = 5000L,
        maxWarmup        = 1000L,
        minWarmup        = 1000L,
        nChains          = 1L,
        thin             = 1L,
        autoTune         = FALSE,
        gibbsSubtreeSwap = FALSE
      ),
      partition = d$part,
      unlink    = c("shape", "ratemultiplier")
    ),
    "without stabilisation"
  )

  samp <- result$samples

  # All class shapes must stay finite and modest in magnitude. The slot-0
  # drift bug pushed class1 to 1e+25..1e+46 in production runs; even at
  # 5000 iter the unbounded RW will leave class1 well outside [1e-3, 1e3].
  for (col in c("class1_rate_log_sd", "class2_rate_log_sd")) {
    v <- samp[, col]
    expect_true(all(is.finite(v)), label = paste(col, "finite"))
    expect_lt(max(v), 1e3, label = paste(col, "max bounded"))
    expect_gt(min(v), 1e-3, label = paste(col, "min bounded"))
  }

  # Lockstep: class1_rate_log_sd and the legacy rate_log_sd column are the
  # same underlying state vector slot, so every sample must agree exactly.
  expect_equal(
    unname(samp[, "class1_rate_log_sd"]),
    unname(samp[, "rate_log_sd"]),
    tolerance = 0,
    label = "class1 lockstep with legacy rate_log_sd"
  )
})

# ==============================================================================
# 5. Every evaluator scores the partitioned model
# ==============================================================================

# The partitioned likelihood gives each class its own shape and rate
# multiplier. Every path that scores or commits a likelihood must use them,
# or an accepted move leaves state$logLik holding the likelihood of a
# different model (#153). Each move type runs on a fresh chain, so a partial-
# CL accept that clears the partition cache cannot divert the moves after it
# onto a different evaluator.
.Ternary12 <- function(seed = 153L, nTip = 8L) {
  set.seed(seed)
  tree <- Preorder(ape::rtree(nTip, rooted = FALSE))
  mat <- matrix(sample(0:2, nTip * 12L, replace = TRUE), nrow = nTip,
                dimnames = list(tree$tip.label, NULL))
  list(tree = tree, mkd = suppressWarnings(MkPrimeData(MatrixToPhyDat(mat))))
}

test_that("every move commits the partitioned likelihood", {
  d <- .Ternary12()
  tree <- d$tree
  mkd <- d$mkd

  chain <- .PartitionedChain(tree, mkd)
  expect_lt(abs(get_state_log_lik(chain$statePtr) - .PartitionedLogLik(chain)),
            1e-9)

  moves <- c(scale_tree_length = 0L, beta_simplex = 4L, nni = 5L, spr = 6L,
             int_walk = 7L, gibbs_spr = 10L, gibbs_subtree_swap = 11L,
             weighted_branch_scale = 12L, weighted_spr = 13L,
             weighted_subtree_swap = 14L, block_gibbs_branch = 15L, tbr = 17L,
             slice_tree_length = 19L, pspr = 20L, dirichlet_branch = 23L,
             local_dirichlet = 24L, gibbs_kprime_sweep = 25L,
             block_kprime_shift = 26L, scale_class_rate_log_sd = 31L,
             dirichlet_simplex_class_w = 32L)
  for (move in names(moves)) {
    chain <- .PartitionedChain(tree, mkd)
    set.seed(1L)
    accepted <- 0L
    drift <- 0
    for (i in seq_len(60L)) {
      charIdx <- switch(move,
                        int_walk = sample(chain$nChar, 1L) - 1L,
                        scale_class_rate_log_sd = sample(2L, 1L),
                        0L)
      if (do_move_cpp(chain$dataPtr, chain$statePtr, moves[[move]], charIdx,
                      0.5, 0.5, 3L, 1.0)) {
        accepted <- accepted + 1L
        drift <- max(drift, abs(get_state_log_lik(chain$statePtr) -
                                  .PartitionedLogLik(chain)))
      }
    }
    expect_gt(accepted, 0L, label = paste(move, "acceptances"))
    expect_lt(drift, 1e-9, label = paste(move, "logLik drift"))
  }
})

# The Gibbs k' sweep draws each k'_i from a conditional it computes for
# itself, which must also be the partitioned one (#153). With the tree and
# every other parameter held fixed, successive sweeps are independent draws
# from it, so the draws must be likelier under the partitioned conditional
# than under the class-blind one.
test_that("gibbs_kprime_sweep samples the partitioned conditional", {
  d <- .Ternary12()
  chain <- .PartitionedChain(d$tree, d$mkd)
  st <- get_mcmc_state(chain$statePtr)
  edgeLen <- st$treeLength * st$relBrLengths
  kObs <- d$mkd$kObs
  shifts <- 0:8

  # The k' prior does not involve the class parameters, so a legacy state
  # supplies it.
  model <- MkPrime:::.FinalizeModel(MkPrimeModel(), d$tree, d$mkd)
  legacyData <- MkPrime:::.InitMcmcData(d$mkd, model)
  legacyState <- MkPrime:::.InitState(d$tree, d$mkd, model)
  LogConditional <- function(i, LogLik) {
    lw <- vapply(kObs[[i]] + shifts, function(k) {
      kPrime <- st$kPrime
      kPrime[[i]] <- k
      priorState <- legacyState
      priorState$kPrime[[i]] <- k
      LogLik(kPrime) + eval_log_prior_cpp(
        legacyData, MkPrime:::.InitMcmcChain(priorState))
    }, double(1))
    # Return:
    lw - max(lw) - log(sum(exp(lw - max(lw))))
  }
  Partitioned <- function(kPrime) {
    cpp_log_likelihood_partitioned_xptr(
      chain$dataPtr, st$edge[, 1], st$edge[, 2], edgeLen, kPrime,
      st$rateLoss, st$classRateLogSd, st$classRate, st$etaNeo, st$betaScale)
  }
  ClassBlind <- function(kPrime) {
    cpp_log_likelihood_xptr(
      chain$dataPtr, st$edge[, 1], st$edge[, 2], edgeLen, kPrime,
      st$rateLoss, st$rateLogSd, st$rateNeo, st$betaScale)
  }

  set.seed(25L)
  draws <- vapply(seq_len(1000L), function(j) {
    do_move_cpp(chain$dataPtr, chain$statePtr, 25L, 0L, 0.5, 0.5, 3L, 1.0)
    get_mcmc_state(chain$statePtr)$kPrime
  }, integer(d$mkd$nChar))

  logRatio <- sum(vapply(seq_len(d$mkd$nChar), function(i) {
    counts <- tabulate(draws[i, ] - kObs[[i]] + 1L, length(shifts))
    sum(counts * (LogConditional(i, Partitioned) -
                    LogConditional(i, ClassBlind)))
  }, double(1)))
  expect_gt(logRatio, 0)
})
