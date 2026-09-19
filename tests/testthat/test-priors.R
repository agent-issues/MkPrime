# Tests for MkPrimeModel, LogPrior (M-014), and Fitch parsimony (M-026)

test_that("MkPrimeModel creates valid object with defaults", {
  model <- MkPrimeModel()
  expect_s3_class(model, "MkPrimeModel")
  expect_equal(model$coding, "variable")
  expect_equal(model$nCat, 6L)
  expect_true(model$relabel)
  expect_equal(model$treeLengthShape, 2)
  expect_null(model$treeLengthRate)
  expect_null(model$expSteps)
})


test_that("MkPrimeModel accepts custom parameters", {
  model <- MkPrimeModel(
    coding = "none", nCat = 4L,
    treeLengthShape = 3, expSteps = 50,
    rateLossMeanlog = 1, rateLossSdlog = 0.5
  )
  expect_equal(model$coding, "none")
  expect_equal(model$nCat, 4L)
  expect_equal(model$expSteps, 50)
  expect_equal(model$treeLengthRate, 2 / 50)
})


test_that(".FinalizeModel computes defaults", {
  model <- MkPrimeModel()
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,t3:0.3);")

  mat <- matrix(c(0, 1, 0, 0, 1, 2), 3, 2,
                dimnames = list(c("t1", "t2", "t3"), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)

  model <- MkPrime:::.FinalizeModel(model, tree, mkd)
  expect_true(!is.null(model$expSteps))
  expect_true(model$expSteps > 0)
  expect_true(!is.null(model$treeLengthRate))
  expect_true(model$treeLengthRate > 0)
})


test_that("LogPrior matches manual density calculations", {
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1), 4, 1,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)
  # Explicit `geometric` prior because the manual calculation below assumes
  # the unconvoluted geometric density on k'.  The default
  # (`empirical_geometric`) uses a convolution and would give a different
  # value.
  model <- MkPrimeModel(expSteps = 10, kPrimePrior = "geometric")

  state <- list(
    tree_length = 0.5,
    rel_br_lengths = tree$edge.length / sum(tree$edge.length),
    rate_loss = 1.5,
    rate_log_sd = 0.3,
    kPrime = 3L,
    p = 0.4
  )

  lp <- MkPrime:::LogPrior(state, model, mkd)

  # Manual calculation
  expected <- dgamma(0.5, shape = 2, rate = 2 / 10, log = TRUE) +
    lfactorial(length(state$rel_br_lengths) - 1L) +
    dgamma(0.3, shape = 1, rate = 1, log = TRUE) +
    # k' geometric: u = 3 - 2 = 1, log(p) + u*log(1-p)
    log(0.4) + 1 * log1p(-0.4) +
    dbeta(0.4, 1, 1, log = TRUE)
  # No rate_loss prior since no neomorphic chars

  expect_equal(lp, expected, tolerance = 1e-12)
})


test_that("LogPrior includes rate_loss for neomorphic chars", {
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1), 4, 1,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd, neomorphic = 1L)
  model <- MkPrimeModel(expSteps = 10)

  state <- list(
    tree_length = 1.0,
    rel_br_lengths = tree$edge.length / sum(tree$edge.length),
    rate_loss = 2.0,
    rate_log_sd = 0.5,
    kPrime = 2L,
    p = 0.5
  )

  lp <- MkPrime:::LogPrior(state, model, mkd)

  # Should include rate_loss LogNormal(0, 2) density
  lp_rate_loss <- dlnorm(2.0, meanlog = 0, sdlog = 2, log = TRUE)
  # rate_loss component is included since there are neomorphic chars
  expect_true(is.finite(lp))

  # Remove rate_loss from expected and check separately
  state2 <- state
  state2$rate_loss <- 1.0
  lp2 <- MkPrime:::LogPrior(state2, model, mkd)
  # Different rate_loss should give different prior

  expect_false(isTRUE(all.equal(lp, lp2)))
})


test_that("LogPrior returns -Inf for invalid parameter values", {
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,t3:0.3);")
  mat <- matrix(c(0, 1, 0), 3, 1,
                dimnames = list(c("t1", "t2", "t3"), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)
  model <- MkPrimeModel(expSteps = 10)

  base_state <- list(
    tree_length = 0.5,
    rel_br_lengths = tree$edge.length / sum(tree$edge.length),
    rate_loss = 1.0,
    rate_log_sd = 0.3,
    kPrime = 2L,
    p = 0.5
  )

  # Negative tree_length
  bad <- base_state
  bad$tree_length <- -1
  expect_equal(MkPrime:::LogPrior(bad, model, mkd), -Inf)

  # Negative rate_log_sd
  bad <- base_state
  bad$rate_log_sd <- -0.1
  expect_equal(MkPrime:::LogPrior(bad, model, mkd), -Inf)
})


test_that("LogPrior: k' prior favors kObs when p is high", {
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,t3:0.3);")
  mat <- matrix(c(0, 1, 0), 3, 1,
                dimnames = list(c("t1", "t2", "t3"), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)
  model <- MkPrimeModel(expSteps = 10)

  state_low_k <- list(
    tree_length = 0.5,
    rel_br_lengths = tree$edge.length / sum(tree$edge.length),
    rate_loss = 1.0,
    rate_log_sd = 0.3,
    kPrime = 2L,
    p = 0.9
  )
  state_high_k <- state_low_k
  state_high_k$kPrime <- 10L

  lp_low <- MkPrime:::LogPrior(state_low_k, model, mkd)
  lp_high <- MkPrime:::LogPrior(state_high_k, model, mkd)

  # High p means geometric is concentrated near 0 (kObs)
  expect_gt(lp_low, lp_high)
})


test_that("Fitch parsimony matches hand calculation", {
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  # Char 1: (0,0,1,1) -> 1 change; Char 2: (0,1,0,1) -> 2 changes
  mat <- matrix(c(0, 0, 1, 1, 0, 1, 0, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)

  expect_equal(MkPrime:::.FitchScore(tree, mkd), 3L)
})


test_that("Fitch score used for expSteps default", {
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 0, 1, 1, 0, 1, 0, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)

  model <- MkPrimeModel()
  finalized <- MkPrime:::.FinalizeModel(model, tree, mkd)
  expect_equal(finalized$expSteps, 3)
  expect_equal(finalized$treeLengthRate, 2 / 3)
})
