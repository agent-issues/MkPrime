# Unit tests for .ConvergenceTier() and its consumers

test_that(".ConvergenceTier classifies every sampled column", {
  cols <- c("log_posterior", "log_likelihood", "tree_length", "rate_loss",
            "rate_log_sd", "p", "kprime_alpha", "kprime_beta", "rate_neo",
            "beta_scale", "swap_cold", "topo_hash", "kPrime_1", "br_1")
  tier <- MkPrime:::.ConvergenceTier(cols)

  expect_equal(unname(tier[c("log_posterior", "tree_length", "rate_loss",
                             "rate_log_sd", "p", "rate_neo")]),
               rep("gate", 6))
  expect_equal(unname(tier[c("log_likelihood", "kPrime_1")]),
               rep("nuisance", 2))
  expect_equal(unname(tier["br_1"]), "branch")
  expect_equal(unname(tier[c("swap_cold", "topo_hash")]),
               rep("bookkeeping", 2))
})

test_that(".ConvergenceTier gates the beta-geometric and Het hyperparameters", {
  # These were invisible to the old allow-list regexp, so neither ESS nor
  # R-hat was ever computed for them.
  tier <- MkPrime:::.ConvergenceTier(c("kprime_alpha", "kprime_beta",
                                       "beta_scale"))
  expect_true(all(tier == "gate"))
})

test_that(".ConvergenceTier gates an unrecognised parameter family", {
  expect_equal(unname(MkPrime:::.ConvergenceTier("hyper_tau")), "gate")
})

test_that(".GateCols and .ReportCols select the expected columns", {
  cols <- c("log_posterior", "log_likelihood", "swap_cold", "kPrime_1", "br_1")
  expect_equal(MkPrime:::.GateCols(cols), 1L)
  expect_equal(MkPrime:::.ReportCols(cols), c(1L, 2L, 4L))
})

test_that(".KeyParamCols keeps the hash and branch columns out", {
  cols <- c("log_posterior", "topo_hash", "beta_scale", "br_1", "kPrime_2")
  mat <- matrix(0, 1, 5, dimnames = list(NULL, cols))
  expect_equal(colnames(mat)[MkPrime:::.KeyParamCols(mat)],
               c("log_posterior", "beta_scale", "kPrime_2"))
  expect_equal(colnames(mat)[MkPrime:::.PlotParamCols(mat)],
               c("log_posterior", "beta_scale"))
})

test_that(".MinEssRate ignores the topology hash", {
  set.seed(90210)
  n <- 200L
  mat <- cbind(log_posterior = rnorm(n), tree_length = rnorm(n),
               topo_hash = cumsum(rnorm(n)))
  # A random walk has an ESS of order 1; were topo_hash scored, it would set
  # the minimum and so drive both the bandit and the thinning interval.
  withHash <- MkPrime:::.MinEssRate(mat, 1)
  noHash   <- MkPrime:::.MinEssRate(mat[, 1:2], 1)
  expect_equal(withHash, noHash)
})
