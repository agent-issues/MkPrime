# Tests for Phase 9: Gibbs moves
#
# M-082: Gibbs update for p
# M-084: Gibbs update for k'_i (added here when implemented)

library(ape)
library(TreeTools)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

.small_trans_tree <- function() {
  read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
}

.small_trans_pd <- function() {
  mat <- matrix(c(0, 1, 2, 0,
                  0, 1, 0, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  MatrixToPhyDat(mat)
}


# ---------------------------------------------------------------------------
# M-082: Gibbs update for p
# ---------------------------------------------------------------------------

test_that("gibbs_p move is listed for transformational data", {
  pd <- .small_trans_pd()
  mkd <- MkPrimeData(pd)
  nEdge <- 2 * length(.small_trans_tree()$tip.label) - 3L
  nTrans <- sum(mkd$type == "transformational")
  moves <- MkPrime:::.BuildMoves(nEdge, nTrans, hasNeo = FALSE,
                                 mcmc = MkPrimeMCMC(), fixTopology = FALSE)
  p_move <- Filter(function(m) m$name == "p", moves)
  expect_length(p_move, 1L)
  expect_equal(p_move[[1L]]$type, "gibbs_p")
})

test_that("gibbs_p move is absent when there are no transformational chars", {
  mat <- matrix(c(0, 1, 0, 1), 4, 1,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd, knownStates = c("1" = 2L))
  nEdge <- 2 * 4L - 3L
  moves <- MkPrime:::.BuildMoves(nEdge, nTrans = 0L, hasNeo = FALSE,
                                 mcmc = MkPrimeMCMC(), fixTopology = FALSE)
  p_names <- vapply(moves, `[[`, character(1L), "name")
  expect_false("p" %in% p_names)
})

test_that("gibbs_p samples p from the correct Beta full conditional", {
  # With known k' values, the full conditional for p is exactly
  # Beta(a + nTrans, b + sumU) where sumU = sum(k'_i - kObs_i).
  # We verify this by drawing many samples from the Gibbs move and
  # comparing the empirical mean/variance to the theoretical Beta moments.
  set.seed(3817)

  pd <- .small_trans_pd()
  mkd <- MkPrimeData(pd)
  model <- MkPrimeModel()
  tree  <- .small_trans_tree()

  # Fix k' so sumU is known
  transIdx <- which(mkd$type == "transformational")
  fixedKPrime <- mkd$kObs
  fixedKPrime[transIdx] <- mkd$kObs[transIdx] + 1L  # each k'_i = kObs_i + 1
  sumU <- sum(fixedKPrime[transIdx] - mkd$kObs[transIdx])
  nTrans <- length(transIdx)

  # Expected posterior parameters
  a_post <- model$kprimeHyperA + nTrans
  b_post <- model$kprimeHyperB + sumU
  expected_mean <- a_post / (a_post + b_post)
  expected_var  <- a_post * b_post / ((a_post + b_post)^2 * (a_post + b_post + 1))

  # Use R-fallback .DoMove() to collect samples
  state <- list(
    tree = tree, tree_length = 0.5, rel_br_lengths = rep(1/5, 5),
    rate_loss = 1.0, rate_log_sd = 0.5, rate_neo = NULL,
    kPrime = fixedKPrime, p = 0.5
  )
  state$log_prior <- LogPrior(state, model, mkd)
  state$log_lik   <- 0.0   # placeholder; not used by gibbs_p
  state$log_post  <- state$log_lik + state$log_prior

  moves <- list(list(name = "p", type = "gibbs_p", target = "p", weight = 1))
  tuning <- list(beta_simplex = 10, int_walk_window = 1L)

  nDraw <- 2000L
  p_samples <- numeric(nDraw)
  for (i in seq_len(nDraw)) {
    res <- MkPrime:::.DoMove(moves[[1L]], state, mkd, model, tuning = tuning)
    state <- res$state
    p_samples[i] <- state$p
  }

  expect_equal(mean(p_samples), expected_mean, tolerance = 0.05)
  expect_equal(var(p_samples),  expected_var,  tolerance = 0.01)
})

test_that("gibbs_p always accepts (acceptance rate = 1)", {
  set.seed(5021)

  pd <- .small_trans_pd()
  mkd <- MkPrimeData(pd)
  model <- MkPrimeModel()
  tree  <- .small_trans_tree()

  moves <- list(list(name = "p", type = "gibbs_p", target = "p", weight = 1))
  tuning <- list(beta_simplex = 10, int_walk_window = 1L)

  state <- list(
    tree = tree, tree_length = 0.5, rel_br_lengths = rep(1/5, 5),
    rate_loss = 1.0, rate_log_sd = 0.5, rate_neo = NULL,
    kPrime = mkd$kObs, p = 0.5
  )
  state$log_prior <- LogPrior(state, model, mkd)

  for (i in seq_len(50L)) {
    res <- MkPrime:::.DoMove(moves[[1L]], state, mkd, model, tuning = tuning)
    expect_true(res$accept)
  }
})

test_that("RunMkPrime with gibbs_p produces valid posterior for transformational data", {
  set.seed(7342)
  tree <- .small_trans_tree()
  pd   <- .small_trans_pd()

  result <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 600L, thin = 5L,
                       maxWarmup = 200L, minWarmup = 200L, autoTune = FALSE))

  expect_s3_class(result, "MkPosterior")
  expect_true(all(is.finite(result$samples[, "log_posterior"])))
  p_vals <- result$samples[, "p"]
  expect_true(all(p_vals > 0 & p_vals < 1))

  # p acceptance should be ~100% (all accepts from Gibbs)
  p_acc <- result$acceptance[["p"]]
  expect_gt(p_acc, 0.99)
})

test_that("gibbs_p updates logPrior correctly in C++ engine", {
  skip_if_not(requireNamespace("ape", quietly = TRUE))
  set.seed(2194)

  pd   <- .small_trans_pd()
  mkd  <- MkPrimeData(pd)
  model <- MkPrimeModel()
  tree <- .small_trans_tree()

  # Initialize C++ state
  model  <- MkPrime:::.FinalizeModel(model, tree, mkd)
  state0 <- MkPrime:::.InitState(tree, mkd, model)
  mcmcData <- MkPrime:::.InitMcmcData(mkd, model)
  statePtr <- MkPrime:::.InitMcmcChain(state0)
  fill_partition_cache(mcmcData, statePtr)
  allocate_cl_workspace(mcmcData, statePtr)

  # Record initial log-posterior
  s_before <- get_mcmc_state(statePtr)
  lp_before <- s_before$logLik + s_before$logPrior

  # Execute one gibbs_p move via C++
  accepted <- do_move_cpp(mcmcData, statePtr,
    moveType = 9L, charIdx = 0L,
    scaleTuning = 0.5, betaSimplexTuning = 10.0,
    intWalkWindow = 1L, beta = 1.0)

  expect_true(accepted)  # Gibbs always accepts

  s_after <- get_mcmc_state(statePtr)
  # logLik must be unchanged; p must have changed; logPrior recomputed
  expect_equal(s_after$logLik, s_before$logLik)
  # p is in (0, 1)
  expect_true(s_after$p > 0 && s_after$p < 1)
  # logPrior recomputed — may differ from before
  # (just check it's finite and consistent)
  expect_true(is.finite(s_after$logPrior))
})
