test_that("fused JC ascertainment matches standalone", {
  set.seed(7284)
  tree <- ape::rtree(8, br = runif)
  tree$edge.length <- abs(tree$edge.length) + 0.01
  tree <- ape::unroot(tree)
  nTip <- ape::Ntip(tree)
  parent <- tree$edge[, 1]
  child  <- tree$edge[, 2]

  kStates <- 3L
  nChar <- 5
  tipStates <- matrix(sample(0:(kStates - 1), nTip * nChar, replace = TRUE),
                      nrow = nTip, ncol = nChar,
                      dimnames = list(tree$tip.label, NULL))
  pd <- TreeTools::MatrixToPhyDat(tipStates)

  mkd <- MkPrimeData(pd, knownStates = rep(kStates, nChar))

  ll_none <- MkpLogLikelihood(tree, mkd, coding = "none", relabel = FALSE)
  ll_var  <- MkpLogLikelihood(tree, mkd, coding = "variable", relabel = FALSE)

  rootFreqs <- rep(1 / kStates, kStates)
  p_standalone <- constant_site_prob_jc(parent, child, tree$edge.length,
                                         nTip, kStates, rootFreqs, 1.0)

  expected_correction <- -nChar * log(1 - p_standalone)
  expect_equal(ll_var - ll_none, expected_correction, tolerance = 1e-10)
})


test_that("fused MkN ascertainment matches standalone", {
  set.seed(3419)
  tree <- ape::rtree(10, br = runif)
  tree$edge.length <- abs(tree$edge.length) + 0.01
  tree <- ape::unroot(tree)
  nTip <- ape::Ntip(tree)
  parent <- tree$edge[, 1]
  child  <- tree$edge[, 2]

  nChar <- 8
  tipStates <- matrix(sample(0:1, nTip * nChar, replace = TRUE),
                      nrow = nTip, ncol = nChar,
                      dimnames = list(tree$tip.label, NULL))
  pd <- TreeTools::MatrixToPhyDat(tipStates)

  mkd <- MkPrimeData(pd, neomorphic = seq_len(nChar))

  rateLoss <- 1.5

  ll_none <- MkpLogLikelihood(tree, mkd, coding = "none", relabel = FALSE,
                               rate_loss = rateLoss)
  ll_var  <- MkpLogLikelihood(tree, mkd, coding = "variable", relabel = FALSE,
                               rate_loss = rateLoss)

  rootFreqs <- c(rateLoss / (1 + rateLoss), 1 / (1 + rateLoss))
  p_standalone <- constant_site_prob_mkn(parent, child, tree$edge.length,
                                          nTip, rateLoss, rootFreqs, 1.0)

  expected_correction <- -nChar * log(1 - p_standalone)
  expect_equal(ll_var - ll_none, expected_correction, tolerance = 1e-10)
})


test_that("fused JC ACRV ascertainment matches standalone", {
  set.seed(8176)
  tree <- ape::rtree(12, br = runif)
  tree$edge.length <- abs(tree$edge.length) + 0.01
  tree <- ape::unroot(tree)
  nTip <- ape::Ntip(tree)
  parent <- tree$edge[, 1]
  child  <- tree$edge[, 2]

  kStates <- 2L
  nChar <- 6
  tipStates <- matrix(sample(0:1, nTip * nChar, replace = TRUE),
                      nrow = nTip, ncol = nChar,
                      dimnames = list(tree$tip.label, NULL))
  pd <- TreeTools::MatrixToPhyDat(tipStates)

  mkd <- MkPrimeData(pd, knownStates = rep(kStates, nChar))

  rateLogSd <- 0.5
  nCat <- 4L

  ll_none <- MkpLogLikelihood(tree, mkd, coding = "none", relabel = FALSE,
                               nCat = nCat, rate_log_sd = rateLogSd)
  ll_var  <- MkpLogLikelihood(tree, mkd, coding = "variable", relabel = FALSE,
                               nCat = nCat, rate_log_sd = rateLogSd)

  rootFreqs <- rep(1 / kStates, kStates)
  acrv_rates <- MkPrime:::DiscreteLognormalRates(rateLogSd, nCat)
  p_standalone <- constant_site_prob_jc(parent, child, tree$edge.length,
                                         nTip, kStates, rootFreqs,
                                         acrv_rates)

  expected_correction <- -nChar * log(1 - p_standalone)
  expect_equal(ll_var - ll_none, expected_correction, tolerance = 1e-10)
})


test_that("fused Het ascertainment gives valid variable correction", {
  # The fused ascertainment in pruning_f81_het_acrv_flat replaces the
  # standalone het_constant_site_prob (buggy per-edge-product shortcut)
  # with proper Felsenstein pruning of k pseudo-characters.
  # Existing test-het.R (28+ tests) validates Het likelihood numerics
  # thoroughly; here we verify the sanity property that the
  # variable-coding correction is positive and finite.
  set.seed(5982)
  tree <- ape::rtree(6, br = runif)
  tree$edge.length <- abs(tree$edge.length) + 0.01
  tree <- ape::unroot(tree)
  nTip <- ape::Ntip(tree)

  nChar <- 4
  kStates <- 3L
  tipStates <- matrix(sample(0:(kStates - 1), nTip * nChar, replace = TRUE),
                      nrow = nTip, ncol = nChar,
                      dimnames = list(tree$tip.label, NULL))
  pd <- TreeTools::MatrixToPhyDat(tipStates)

  mkd <- MkPrimeData(pd, knownStates = rep(kStates, nChar))

  # Build model + C++ state via internals (same pattern as test-het.R)
  model_var  <- MkPrimeModel(coding = "variable", nCat = 2L,
                              qHeterogeneity = TRUE, nBetaCat = 3L)
  model_none <- MkPrimeModel(coding = "none", nCat = 2L,
                              qHeterogeneity = TRUE, nBetaCat = 3L)

  .build <- function(model) {
    model <- MkPrime:::.FinalizeModel(model, tree, mkd)
    tree0 <- TreeTools::Preorder(tree)
    state0 <- MkPrime:::.InitState(tree0, mkd, model)
    state0$beta_scale <- 2.0
    state0$rate_log_sd <- 0.3
    state0$log_prior <- MkPrime:::LogPrior(state0, model, mkd)
    dataPtr  <- MkPrime:::.InitMcmcData(mkd, model)
    statePtr <- MkPrime:::.InitMcmcChain(state0)
    fill_partition_cache(dataPtr, statePtr)
    allocate_cl_workspace(dataPtr, statePtr)
    eval_full_loglik_cpp(dataPtr, statePtr)
  }

  ll_none <- .build(model_none)
  ll_var  <- .build(model_var)

  expect_true(is.finite(ll_none))
  expect_true(is.finite(ll_var))
  # Variable correction increases log-likelihood (denominator < 1)
  expect_gt(ll_var, ll_none)
})
