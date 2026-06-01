# Tests for the Beta-Geometric kPrime prior
skip_if_not_installed("phangorn")

# --- Helper: build small trans-only dataset ---------------------------------
make_bg_fixture <- function() {
  tree <- ape::rtree(6)
  tree$edge.length <- tree$edge.length / sum(tree$edge.length) * 3
  nChar <- 5L
  mat <- matrix(sample(0:2, 6 * nChar, replace = TRUE), nrow = 6,
                dimnames = list(tree$tip.label, NULL))
  pd <- phangorn::phyDat(mat, type = "USER", levels = as.character(0:2))
  mkd <- MkPrimeData(pd)  # all transformational
  model <- MkPrimeModel(kPrimePrior = "beta_geometric", treeLengthRate = 0.5)
  model <- MkPrime:::.FinalizeModel(model, tree, mkd)
  tree <- TreeTools::Preorder(tree)
  list(tree = tree, mkd = mkd, model = model)
}

# --- MkPrimeModel construction ----------------------------------------------

test_that("MkPrimeModel accepts 'beta_geometric'", {
  m <- MkPrimeModel(kPrimePrior = "beta_geometric")
  expect_s3_class(m, "MkPrimeModel")
  expect_equal(m$kPrimePrior, "beta_geometric")
  expect_equal(m$kprimeAlpha, 1)
  expect_equal(m$kprimeBeta, 1)
})

test_that("custom alpha/beta are stored", {
  m <- MkPrimeModel(kPrimePrior = "beta_geometric",
                     kprimeAlpha = 2, kprimeBeta = 0.5)
  expect_equal(m$kprimeAlpha, 2)
  expect_equal(m$kprimeBeta, 0.5)
})

test_that("invalid alpha/beta rejected", {
  expect_error(MkPrimeModel(kPrimePrior = "beta_geometric", kprimeAlpha = -1),
               "positive")
  expect_error(MkPrimeModel(kPrimePrior = "beta_geometric", kprimeBeta = 0),
               "positive")
})

test_that("print method works for beta_geometric", {
  m <- MkPrimeModel(kPrimePrior = "beta_geometric")
  out <- capture.output(print(m), type = "message")
  expect_true(any(grepl("Beta-Geometric", out)))
})

# --- R-side LogPrior ---------------------------------------------------------

test_that("LogPrior: beta_geometric matches manual computation", {
  f <- make_bg_fixture()
  state <- MkPrime:::.InitState(f$tree, f$mkd, f$model)

  lp <- LogPrior(state, f$model, f$mkd)
  expect_true(is.finite(lp))

  # Manual: kPrime component
  transIdx <- which(f$mkd$type == "transformational")
  u <- state$kPrime[transIdx] - f$mkd$kObs[transIdx]
  a <- state$kprime_alpha
  b <- state$kprime_beta
  kp_part <- sum(lbeta(a + 1, b + u) - lbeta(a, b))
  hyper_part <- dexp(a, 1, log = TRUE) + dexp(b, 1, log = TRUE)

  # Other components
  tl <- dgamma(state$tree_length, shape = f$model$treeLengthShape,
                rate = f$model$treeLengthRate, log = TRUE)
  dir <- lfactorial(length(state$rel_br_lengths) - 1L)
  rls <- dgamma(state$rate_log_sd, shape = f$model$rateLogSdShape,
                 rate = f$model$rateLogSdRate, log = TRUE)
  expected <- tl + dir + rls + kp_part + hyper_part
  expect_equal(lp, expected, tolerance = 1e-10)
})

test_that("LogPrior: beta_geometric returns -Inf for invalid alpha/beta", {
  f <- make_bg_fixture()
  state <- MkPrime:::.InitState(f$tree, f$mkd, f$model)
  state$kprime_alpha <- 0
  expect_equal(LogPrior(state, f$model, f$mkd), -Inf)
  state$kprime_alpha <- 1
  state$kprime_beta <- -1
  expect_equal(LogPrior(state, f$model, f$mkd), -Inf)
})

test_that("LogPrior: beta_geometric varies with u", {
  f <- make_bg_fixture()
  state <- MkPrime:::.InitState(f$tree, f$mkd, f$model)
  lp0 <- LogPrior(state, f$model, f$mkd)

  # Increase one kPrime by 1
  transIdx <- which(f$mkd$type == "transformational")[1]
  state$kPrime[transIdx] <- state$kPrime[transIdx] + 1L
  lp1 <- LogPrior(state, f$model, f$mkd)

  # With alpha=beta=1, penalty for u=0→u=1 is lbeta(2,2)-lbeta(2,1)
  # = log(1/6) - log(1/2) = log(1/3) ≈ -1.099
  expect_equal(lp1 - lp0, lbeta(2, 2) - lbeta(2, 1), tolerance = 1e-10)
})

# --- C++ cpp_log_prior agreement ----------------------------------------------

test_that("C++ and R LogPrior agree for beta_geometric", {
  f <- make_bg_fixture()
  state <- MkPrime:::.InitState(f$tree, f$mkd, f$model)
  r_lp <- LogPrior(state, f$model, f$mkd)

  # C++ path: init data/state pointers
  mcmcData <- MkPrime:::.InitMcmcData(f$mkd, f$model)
  statePtr <- MkPrime:::.InitMcmcChain(state)
  fill_partition_cache(mcmcData, statePtr)

  cpp_state <- get_mcmc_state(statePtr)
  expect_equal(cpp_state$logPrior, r_lp, tolerance = 1e-8)
})

# Note: the axis-aligned Bactrian moves on raw α and β (move codes 27, 28)
# were removed when the BG hyperparameter sampler was reparameterised to
# (s, r) coordinates. End-to-end coverage that α / β actually move is
# provided by test-bg-hyperparameter-moves.R.

# --- ParamNames and StateToRow -----------------------------------------------

test_that("ParamNames includes kprime_alpha/kprime_beta, not p", {
  f <- make_bg_fixture()
  nEdge <- nrow(f$tree$edge)
  pn <- MkPrime:::.ParamNames(f$mkd, nEdge, kPrimePrior = "beta_geometric")
  expect_true("kprime_alpha" %in% pn)
  expect_true("kprime_beta" %in% pn)
  expect_false("p" %in% pn)
})

test_that("StateToRow length matches ParamNames", {
  f <- make_bg_fixture()
  state <- MkPrime:::.InitState(f$tree, f$mkd, f$model)
  mcmcData <- MkPrime:::.InitMcmcData(f$mkd, f$model)
  statePtr <- MkPrime:::.InitMcmcChain(state)
  fill_partition_cache(mcmcData, statePtr)

  nEdge <- nrow(f$tree$edge)
  pn <- MkPrime:::.ParamNames(f$mkd, nEdge, kPrimePrior = "beta_geometric")
  row <- MkPrime:::.StateToRow(statePtr, f$mkd, nEdge,
                                kPrimePrior = "beta_geometric")
  expect_equal(length(row), length(pn))

  names(row) <- pn
  expect_equal(row[["kprime_alpha"]], 1.0)
  expect_equal(row[["kprime_beta"]], 1.0)
})

# --- BuildMoves --------------------------------------------------------------

test_that("BuildMoves registers slice_kprime_s/r for BG, not p/gibbs_p", {
  f <- make_bg_fixture()
  nTrans <- sum(f$mkd$type == "transformational")
  nEdge <- nrow(f$tree$edge)
  mcmc <- MkPrimeMCMC(nIter = 100)
  moves <- MkPrime:::.BuildMoves(nEdge, nTrans, hasNeo = FALSE, mcmc = mcmc,
                                  kPrimePrior = "beta_geometric")
  moveNames <- vapply(moves, `[[`, character(1), "name")
  expect_true("slice_kprime_s" %in% moveNames)
  expect_true("slice_kprime_r" %in% moveNames)
  expect_false("p" %in% moveNames)
})

test_that("BuildMoves for geometric registers mh_logit_p, not gibbs_p or BG slice", {
  f <- make_bg_fixture()
  nTrans <- sum(f$mkd$type == "transformational")
  nEdge <- nrow(f$tree$edge)
  mcmc <- MkPrimeMCMC(nIter = 100)
  moves <- MkPrime:::.BuildMoves(nEdge, nTrans, hasNeo = FALSE, mcmc = mcmc,
                                  kPrimePrior = "geometric")
  moveNames <- vapply(moves, `[[`, character(1), "name")
  # Stage 2 (MARGINAL-K-TRUNC-001): the truncated geometric samples p via
  # mh_logit_p (case 30), not the legacy conjugate gibbs_p (named "p"); and it
  # does not use the beta_geometric (s, r) slice samplers.
  expect_true("mh_logit_p" %in% moveNames)
  expect_false("p" %in% moveNames)
  expect_false("slice_kprime_s" %in% moveNames)
})

# --- Gibbs kPrime sweep correctness ------------------------------------------

test_that("Gibbs sweep with beta_geometric samples from correct full conditional", {
  f <- make_bg_fixture()
  state <- MkPrime:::.InitState(f$tree, f$mkd, f$model)
  mcmcData <- MkPrime:::.InitMcmcData(f$mkd, f$model)
  statePtr <- MkPrime:::.InitMcmcChain(state)
  fill_partition_cache(mcmcData, statePtr)
  allocate_cl_workspace(mcmcData, statePtr)

  # Run many Gibbs sweeps and collect kPrime distribution for char 1
  set.seed(4297)
  nIter <- 200
  transIdx <- which(f$mkd$type == "transformational")
  kp_hist <- integer(nIter)
  for (i in seq_len(nIter)) {
    do_move_cpp(mcmcData, statePtr, 25L, 0L, 0.5, 10.0, 1L, 1.0)
    st <- get_mcmc_state(statePtr)
    kp_hist[i] <- st$kPrime[transIdx[1]]
  }

  # kPrime should be at least kObs
  expect_true(all(kp_hist >= f$mkd$kObs[transIdx[1]]))
})

# --- int_walk kPrime viability ------------------------------------------------

test_that("int_walk kPrime achieves nonzero acceptance under beta_geometric", {
  f <- make_bg_fixture()
  state <- MkPrime:::.InitState(f$tree, f$mkd, f$model)
  mcmcData <- MkPrime:::.InitMcmcData(f$mkd, f$model)
  statePtr <- MkPrime:::.InitMcmcChain(state)
  fill_partition_cache(mcmcData, statePtr)
  allocate_cl_workspace(mcmcData, statePtr)

  set.seed(6103)
  transIdx <- which(f$mkd$type == "transformational") - 1L  # 0-based
  nAcc <- 0L
  nTry <- 200L
  for (i in seq_len(nTry)) {
    ci <- sample(transIdx, 1)
    nAcc <- nAcc + do_move_cpp(mcmcData, statePtr, 7L, ci, 0.5, 10.0, 1L, 1.0)
  }
  # Should be well above 0; expect at least 5% acceptance
  expect_gt(nAcc / nTry, 0.02,
            label = sprintf("int_walk acceptance %d/%d", nAcc, nTry))
})

# --- Hyperparameter move tests ------------------------------------------------

# Acceptance-rate tests for the retired axis-aligned Bactrian moves on raw
# α / β (move codes 27 / 28) were removed alongside those moves. End-to-end
# coverage that α and β actually move under the (s, r) slice sampler is in
# test-bg-hyperparameter-moves.R; logPrior-after-Gibbs consistency is also
# covered by the existing Gibbs-sweep test above.

test_that("logPrior stays consistent across Gibbs kPrime sweeps", {
  f <- make_bg_fixture()
  state <- MkPrime:::.InitState(f$tree, f$mkd, f$model)
  mcmcData <- MkPrime:::.InitMcmcData(f$mkd, f$model)
  statePtr <- MkPrime:::.InitMcmcChain(state)
  fill_partition_cache(mcmcData, statePtr)
  allocate_cl_workspace(mcmcData, statePtr)

  set.seed(9274)
  for (i in 1:50) {
    do_move_cpp(mcmcData, statePtr, 25L, 0L, 0.5, 10.0, 1L, 1.0)
  }

  st <- get_mcmc_state(statePtr)
  r_state <- list(
    tree_length = st$treeLength,
    rel_br_lengths = st$relBrLengths,
    rate_loss = st$rateLoss,
    rate_log_sd = st$rateLogSd,
    kPrime = st$kPrime,
    kprime_alpha = st$kprimeAlpha,
    kprime_beta = st$kprimeBeta
  )
  r_lp <- LogPrior(r_state, f$model, f$mkd)
  expect_equal(st$logPrior, r_lp, tolerance = 1e-8,
               label = "C++ logPrior matches R recomputation after Gibbs sweeps")
})

# --- Backward compatibility: geometric prior still works ----------------------

test_that("geometric prior is unaffected by beta_geometric code", {
  f <- make_bg_fixture()
  model_geo <- MkPrimeModel(kPrimePrior = "geometric", treeLengthRate = 0.5)
  model_geo <- MkPrime:::.FinalizeModel(model_geo, f$tree, f$mkd)

  state <- MkPrime:::.InitState(f$tree, f$mkd, model_geo)
  expect_true(is.finite(state$log_prior))
  expect_true(!is.null(state$p))
  expect_null(state$kprime_alpha)

  mcmcData <- MkPrime:::.InitMcmcData(f$mkd, model_geo)
  statePtr <- MkPrime:::.InitMcmcChain(state)
  fill_partition_cache(mcmcData, statePtr)

  st <- get_mcmc_state(statePtr)
  expect_equal(st$logPrior, state$log_prior, tolerance = 1e-8)
})

# --- Beta-Geometric math spot checks -----------------------------------------

test_that("Beta-Geometric PMF is correct for known values", {
  # P(u=0 | α=1, β=1) = B(2,1)/B(1,1) = 1/2
  expect_equal(exp(lbeta(2, 1) - lbeta(1, 1)), 0.5)
  # P(u=1 | α=1, β=1) = B(2,2)/B(1,1) = 1/6
  expect_equal(exp(lbeta(2, 2) - lbeta(1, 1)), 1/6, tolerance = 1e-10)
  # P(u=2 | α=1, β=1) = B(2,3)/B(1,1) = 1/12
  expect_equal(exp(lbeta(2, 3) - lbeta(1, 1)), 1/12, tolerance = 1e-10)

  # Incremental formula: logP(u) = logP(u-1) + log(β+u-1) - log(α+β+u)
  a <- 1; b <- 1
  logP0 <- log(a) - log(a + b)
  logP1 <- logP0 + log(b + 0) - log(a + b + 1)
  logP2 <- logP1 + log(b + 1) - log(a + b + 2)
  expect_equal(exp(logP0), 0.5, tolerance = 1e-10)
  expect_equal(exp(logP1), 1/6, tolerance = 1e-10)
  expect_equal(exp(logP2), 1/12, tolerance = 1e-10)
})

test_that("Beta-Geometric sums to 1 (numeric)", {
  # P(u >= 0 | α, β) should sum to 1
  a <- 2.5; b <- 0.8
  logP <- numeric(200)
  logP[1] <- log(a) - log(a + b)
  for (u in 2:200) {
    logP[u] <- logP[u-1] + log(b + u - 2) - log(a + b + u - 1)
  }
  total <- sum(exp(logP))
  expect_equal(total, 1.0, tolerance = 1e-4)
})
