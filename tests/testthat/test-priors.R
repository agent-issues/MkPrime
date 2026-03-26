# Tests for MkPrimeModel, log_prior (M-014), and Fitch parsimony (M-026)

test_that("MkPrimeModel creates valid object with defaults", {
  model <- MkPrimeModel()
  expect_s3_class(model, "MkPrimeModel")
  expect_equal(model$coding, "variable")
  expect_equal(model$nCat, 6L)
  expect_true(model$relabel)
  expect_equal(model$tree_length_shape, 2)
  expect_null(model$tree_length_rate)
  expect_null(model$exp_steps)
})


test_that("MkPrimeModel accepts custom parameters", {
  model <- MkPrimeModel(
    coding = "none", nCat = 4L,
    tree_length_shape = 3, exp_steps = 50,
    rate_loss_meanlog = 1, rate_loss_sdlog = 0.5
  )
  expect_equal(model$coding, "none")
  expect_equal(model$nCat, 4L)
  expect_equal(model$exp_steps, 50)
  expect_equal(model$tree_length_rate, 2 / 50)
})


test_that(".finalize_model computes defaults", {
  library(ape)
  model <- MkPrimeModel()
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,t3:0.3);")

  mat <- matrix(c(0, 1, 0, 0, 1, 2), 3, 2,
                dimnames = list(c("t1", "t2", "t3"), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)

  model <- MkPrime:::.finalize_model(model, tree, mkd)
  expect_true(!is.null(model$exp_steps))
  expect_true(model$exp_steps > 0)
  expect_true(!is.null(model$tree_length_rate))
  expect_true(model$tree_length_rate > 0)
})


test_that("log_prior matches manual density calculations", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1), 4, 1,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)
  model <- MkPrimeModel(exp_steps = 10)

  state <- list(
    tree_length = 0.5,
    rel_br_lengths = tree$edge.length / sum(tree$edge.length),
    rate_loss = 1.5,
    rate_log_sd = 0.3,
    kPrime = 3L,
    p = 0.4
  )

  lp <- MkPrime:::log_prior(state, model, mkd)

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


test_that("log_prior includes rate_loss for neomorphic chars", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1), 4, 1,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd, neomorphic = 1L)
  model <- MkPrimeModel(exp_steps = 10)

  state <- list(
    tree_length = 1.0,
    rel_br_lengths = tree$edge.length / sum(tree$edge.length),
    rate_loss = 2.0,
    rate_log_sd = 0.5,
    kPrime = 2L,
    p = 0.5
  )

  lp <- MkPrime:::log_prior(state, model, mkd)

  # Should include rate_loss LogNormal(0, 2) density
  lp_rate_loss <- dlnorm(2.0, meanlog = 0, sdlog = 2, log = TRUE)
  # rate_loss component is included since there are neomorphic chars
  expect_true(is.finite(lp))

  # Remove rate_loss from expected and check separately
  state2 <- state
  state2$rate_loss <- 1.0
  lp2 <- MkPrime:::log_prior(state2, model, mkd)
  # Different rate_loss should give different prior

  expect_false(isTRUE(all.equal(lp, lp2)))
})


test_that("log_prior returns -Inf for invalid parameter values", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,t3:0.3);")
  mat <- matrix(c(0, 1, 0), 3, 1,
                dimnames = list(c("t1", "t2", "t3"), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)
  model <- MkPrimeModel(exp_steps = 10)

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
  expect_equal(MkPrime:::log_prior(bad, model, mkd), -Inf)

  # Negative rate_log_sd
  bad <- base_state
  bad$rate_log_sd <- -0.1
  expect_equal(MkPrime:::log_prior(bad, model, mkd), -Inf)
})


test_that("log_prior: k' prior favors kObs when p is high", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,t3:0.3);")
  mat <- matrix(c(0, 1, 0), 3, 1,
                dimnames = list(c("t1", "t2", "t3"), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)
  model <- MkPrimeModel(exp_steps = 10)

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

  lp_low <- MkPrime:::log_prior(state_low_k, model, mkd)
  lp_high <- MkPrime:::log_prior(state_high_k, model, mkd)

  # High p means geometric is concentrated near 0 (kObs)
  expect_gt(lp_low, lp_high)
})


test_that("Fitch parsimony matches hand calculation", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  # Char 1: (0,0,1,1) -> 1 change; Char 2: (0,1,0,1) -> 2 changes
  mat <- matrix(c(0, 0, 1, 1, 0, 1, 0, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)

  expect_equal(MkPrime:::.fitch_score(tree, mkd), 3L)
})


test_that("Fitch score used for exp_steps default", {
  library(ape)
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 0, 1, 1, 0, 1, 0, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  pd <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)

  model <- MkPrimeModel()
  finalized <- MkPrime:::.finalize_model(model, tree, mkd)
  expect_equal(finalized$exp_steps, 3)
  expect_equal(finalized$tree_length_rate, 2 / 3)
})
