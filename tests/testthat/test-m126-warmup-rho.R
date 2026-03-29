# Tests for M-126: Warmup rho estimation for 2D joint Bactrian moves

# --- .AccumulateRhoSnapshot unit tests ---

test_that(".AccumulateRhoSnapshot builds a matrix with correct columns", {
  paramNames <- c("log_posterior", "log_likelihood", "tree_length",
                  "rate_log_sd", "rel_br_1", "rel_br_2")
  state <- list(treeLength = 5.0, rateLogSd = 0.8, rateLoss = 0.3)

  buf <- MkPrime:::.AccumulateRhoSnapshot(NULL, state, hasNeo = FALSE,
                                          paramNames)
  expect_equal(nrow(buf), 1L)
  expect_equal(colnames(buf), paramNames)
  expect_equal(unname(buf[1, "tree_length"]), 5.0)
  expect_equal(unname(buf[1, "rate_log_sd"]), 0.8)
})

test_that(".AccumulateRhoSnapshot accumulates rows", {
  paramNames <- c("log_posterior", "log_likelihood", "tree_length",
                  "rate_loss", "rate_log_sd")
  state <- list(treeLength = 5.0, rateLogSd = 0.8, rateLoss = 0.3)

  buf <- NULL
  for (i in 1:10) {
    state$treeLength <- 5.0 + i * 0.1
    buf <- MkPrime:::.AccumulateRhoSnapshot(buf, state, hasNeo = TRUE,
                                            paramNames)
  }
  expect_equal(nrow(buf), 10L)
  expect_equal(unname(buf[10, "tree_length"]), 6.0)
  expect_equal(unname(buf[10, "rate_loss"]), 0.3)
})

test_that(".AccumulateRhoSnapshot caps at 500 rows", {
  paramNames <- c("log_posterior", "log_likelihood", "tree_length",
                  "rate_log_sd")
  state <- list(treeLength = 5.0, rateLogSd = 0.8, rateLoss = 0.3)

  buf <- NULL
  for (i in 1:550) {
    state$treeLength <- i * 0.01
    buf <- MkPrime:::.AccumulateRhoSnapshot(buf, state, hasNeo = FALSE,
                                            paramNames)
  }
  expect_equal(nrow(buf), 500L)
  # Should retain the LAST 500 rows (rows 51–550)
  expect_equal(unname(buf[500, "tree_length"]), 550 * 0.01)
  expect_equal(unname(buf[1, "tree_length"]), 51 * 0.01)
})

test_that(".AccumulateRhoSnapshot includes rate_loss when hasNeo", {
  paramNames <- c("log_posterior", "log_likelihood", "tree_length",
                  "rate_loss", "rate_log_sd")
  state <- list(treeLength = 5.0, rateLogSd = 0.8, rateLoss = 0.3)

  buf <- MkPrime:::.AccumulateRhoSnapshot(NULL, state, hasNeo = TRUE,
                                          paramNames)
  expect_equal(unname(buf[1, "rate_loss"]), 0.3)
})

# --- Integration: rho estimated from correlated warmup snapshots ---

test_that(".EstimateJointRhos detects positive correlation", {
  set.seed(8147)
  n <- 100
  tl <- exp(rnorm(n, 1.5, 0.3))
  # rate_log_sd positively correlated with tree_length
  rls <- exp(0.5 * log(tl) + rnorm(n, 0, 0.2))
  paramNames <- c("log_posterior", "log_likelihood", "tree_length",
                  "rate_log_sd")
  buf <- matrix(0, nrow = n, ncol = length(paramNames),
                dimnames = list(NULL, paramNames))
  buf[, "tree_length"] <- tl
  buf[, "rate_log_sd"] <- rls
  rhos <- MkPrime:::.EstimateJointRhos(buf, hasNeo = FALSE)
  expect_gt(rhos$rho_tl_rls, 0.3)
})

test_that(".EstimateJointRhos detects negative correlation", {
  set.seed(2739)
  n <- 100
  tl <- exp(rnorm(n, 1.5, 0.3))
  # rate_log_sd negatively correlated with tree_length
  rls <- exp(-0.8 * log(tl) + rnorm(n, 2, 0.2))
  paramNames <- c("log_posterior", "log_likelihood", "tree_length",
                  "rate_log_sd")
  buf <- matrix(0, nrow = n, ncol = length(paramNames),
                dimnames = list(NULL, paramNames))
  buf[, "tree_length"] <- tl
  buf[, "rate_log_sd"] <- rls
  rhos <- MkPrime:::.EstimateJointRhos(buf, hasNeo = FALSE)
  expect_lt(rhos$rho_tl_rls, -0.3)
})

test_that(".EstimateJointRhos returns 0 with insufficient samples", {
  paramNames <- c("log_posterior", "log_likelihood", "tree_length",
                  "rate_log_sd")
  buf <- matrix(1, nrow = 10, ncol = length(paramNames),
                dimnames = list(NULL, paramNames))
  rhos <- MkPrime:::.EstimateJointRhos(buf, hasNeo = FALSE)
  expect_equal(rhos$rho_tl_rls, 0.0)
})
