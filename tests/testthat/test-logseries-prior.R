# Tests for Logseries k' prior (opt-in alternative to hierarchical geometric)

# ---------------------------------------------------------------------------
# 1. MkPrimeModel stores logseries settings correctly
# ---------------------------------------------------------------------------
test_that("MkPrimeModel stores logseries prior settings", {
  model <- MkPrimeModel(kPrimePrior = "logseries", kprimeLogseriesC = 0.7)
  expect_s3_class(model, "MkPrimeModel")
  expect_equal(model$kPrimePrior, "logseries")
  expect_equal(model$kprimeLogseriesC, 0.7)
})


test_that("MkPrimeModel kPrimePrior defaults to empirical_geometric", {
  model <- MkPrimeModel()
  expect_equal(model$kPrimePrior, "empirical_geometric")
})


test_that("MkPrimeModel rejects invalid kPrimePrior", {
  expect_error(MkPrimeModel(kPrimePrior = "uniform"), class = "error")
})


# ---------------------------------------------------------------------------
# 2. LogPrior logseries: matches manual density
# ---------------------------------------------------------------------------
test_that("LogPrior logseries matches manual density calculation", {
  library("ape")
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1), 4, 1,
                dimnames = list(paste0("t", 1:4), NULL))
  pd  <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)

  model <- MkPrimeModel(kPrimePrior = "logseries",
                        kprimeLogseriesC = 0.7,
                        expSteps = 10)

  state <- list(
    tree_length    = 0.5,
    rel_br_lengths = tree$edge.length / sum(tree$edge.length),
    rate_loss      = 1.5,
    rate_log_sd    = 0.3,
    kPrime         = 3L   # kObs = 2, so k' = 3 is valid
    # no p — logseries
  )

  lp <- MkPrime:::LogPrior(state, model, mkd)

  # Manual:
  # log P(k' = 3; c = 0.7) = 3*log(0.7) - log(3) - log(-log(1-0.7))
  c_ls  <- 0.7
  kp    <- 3L
  expected_kprime <- kp * log(c_ls) - log(kp) - log(-log1p(-c_ls))

  expected <- dgamma(0.5, shape = 2, rate = 2 / 10, log = TRUE) +
    lfactorial(length(state$rel_br_lengths) - 1L) +
    dgamma(0.3, shape = 1, rate = 1, log = TRUE) +
    expected_kprime
  # No rate_loss (no neomorphic), no Beta(p) term

  expect_equal(lp, expected, tolerance = 1e-12)
})


# ---------------------------------------------------------------------------
# 3. LogPrior logseries: higher k' has lower prior density (for c < 1)
# ---------------------------------------------------------------------------
test_that("LogPrior logseries: higher k' has lower prior density", {
  library("ape")
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,t3:0.3);")
  mat  <- matrix(c(0, 1, 0), 3, 1,
                 dimnames = list(c("t1", "t2", "t3"), NULL))
  pd   <- TreeTools::MatrixToPhyDat(mat)
  mkd  <- MkPrimeData(pd)
  model <- MkPrimeModel(kPrimePrior = "logseries",
                        kprimeLogseriesC = 0.7,
                        expSteps = 10)

  base <- list(
    tree_length    = 0.5,
    rel_br_lengths = tree$edge.length / sum(tree$edge.length),
    rate_loss      = 1.0,
    rate_log_sd    = 0.3,
    kPrime         = 2L   # minimum valid (kObs = 2)
  )

  state_high_k    <- base
  state_high_k$kPrime <- 10L

  lp_low  <- MkPrime:::LogPrior(base,         model, mkd)
  lp_high <- MkPrime:::LogPrior(state_high_k, model, mkd)

  expect_gt(lp_low, lp_high)
})


# ---------------------------------------------------------------------------
# 4. LogPrior logseries: c out of bounds returns -Inf
# ---------------------------------------------------------------------------
test_that("LogPrior logseries: c out of bounds returns -Inf", {
  library("ape")
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,t3:0.3);")
  mat  <- matrix(c(0, 1, 0), 3, 1,
                 dimnames = list(c("t1", "t2", "t3"), NULL))
  pd   <- TreeTools::MatrixToPhyDat(mat)
  mkd  <- MkPrimeData(pd)

  base_state <- list(
    tree_length    = 0.5,
    rel_br_lengths = tree$edge.length / sum(tree$edge.length),
    rate_loss      = 1.0,
    rate_log_sd    = 0.3,
    kPrime         = 2L
  )

  bad_c_zero <- MkPrimeModel(kPrimePrior = "logseries",
                              kprimeLogseriesC = 0.0, expSteps = 10)
  bad_c_one  <- MkPrimeModel(kPrimePrior = "logseries",
                              kprimeLogseriesC = 1.0, expSteps = 10)

  expect_equal(MkPrime:::LogPrior(base_state, bad_c_zero, mkd), -Inf)
  expect_equal(MkPrime:::LogPrior(base_state, bad_c_one,  mkd), -Inf)
})


# ---------------------------------------------------------------------------
# 5. LogPrior logseries: k' < kObs returns -Inf
# ---------------------------------------------------------------------------
test_that("LogPrior logseries: k' < kObs returns -Inf", {
  library("ape")
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,t3:0.3);")
  mat  <- matrix(c(0, 1, 0), 3, 1,
                 dimnames = list(c("t1", "t2", "t3"), NULL))
  pd   <- TreeTools::MatrixToPhyDat(mat)
  mkd  <- MkPrimeData(pd)
  model <- MkPrimeModel(kPrimePrior = "logseries",
                        kprimeLogseriesC = 0.7,
                        expSteps = 10)

  bad_state <- list(
    tree_length    = 0.5,
    rel_br_lengths = tree$edge.length / sum(tree$edge.length),
    rate_loss      = 1.0,
    rate_log_sd    = 0.3,
    kPrime         = 1L   # kObs = 2, so this is invalid
  )
  expect_equal(MkPrime:::LogPrior(bad_state, model, mkd), -Inf)
})


# ---------------------------------------------------------------------------
# 6. LogPrior logseries: works when state has no p field
# ---------------------------------------------------------------------------
test_that("LogPrior logseries: works correctly with no p in state", {
  library("ape")
  tree  <- read.tree(text = "((t1:0.1,t2:0.2):0.15,t3:0.3);")
  mat   <- matrix(c(0, 1, 0), 3, 1,
                  dimnames = list(c("t1", "t2", "t3"), NULL))
  pd    <- TreeTools::MatrixToPhyDat(mat)
  mkd   <- MkPrimeData(pd)
  model <- MkPrimeModel(kPrimePrior = "logseries",
                        kprimeLogseriesC = 0.7, expSteps = 10)

  # State deliberately has no p field
  state <- list(
    tree_length    = 0.5,
    rel_br_lengths = tree$edge.length / sum(tree$edge.length),
    rate_loss      = 1.0,
    rate_log_sd    = 0.3,
    kPrime         = 2L
  )
  expect_null(state$p)
  lp <- MkPrime:::LogPrior(state, model, mkd)
  expect_true(is.finite(lp))
})


# ---------------------------------------------------------------------------
# 7. print.MkPrimeModel runs without error for both prior types
# ---------------------------------------------------------------------------
test_that("print.MkPrimeModel runs without error for logseries model", {
  model <- MkPrimeModel(kPrimePrior = "logseries", kprimeLogseriesC = 0.7)
  # cli output goes to the console directly; just verify no error is thrown
  expect_invisible(print(model))
  # And verify the field values that print would show
  expect_equal(model$kPrimePrior, "logseries")
  expect_equal(model$kprimeLogseriesC, 0.7)
})


test_that("print.MkPrimeModel runs without error for default model", {
  model <- MkPrimeModel()
  expect_invisible(print(model))
  expect_equal(model$kPrimePrior, "empirical_geometric")
})


# ---------------------------------------------------------------------------
# 8. .ParamNames omits p for logseries
# ---------------------------------------------------------------------------
test_that(".ParamNames omits p column for logseries prior", {
  library("ape")
  mat <- matrix(c(0, 1, 0, 1), 4, 1,
                dimnames = list(paste0("t", 1:4), NULL))
  pd  <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)

  nms_geo <- MkPrime:::.ParamNames(mkd, nEdge = 6L, kPrimePrior = "geometric")
  nms_ls  <- MkPrime:::.ParamNames(mkd, nEdge = 6L, kPrimePrior = "logseries")

  expect_true("p" %in% nms_geo)
  expect_false("p" %in% nms_ls)
  # Logseries should have exactly one fewer column
  expect_equal(length(nms_ls), length(nms_geo) - 1L)
})


# ---------------------------------------------------------------------------
# 9. RunMkPrime smoke test with logseries prior
#    Short run; verify no error and samples contain no p column.
# ---------------------------------------------------------------------------
test_that("RunMkPrime smoke test: logseries prior runs and produces correct columns", {
  skip_on_cran()
  library("ape")

  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat  <- matrix(c(0, 1, 0, 1,
                   0, 0, 1, 1), 4, 2,
                 dimnames = list(paste0("t", 1:4), NULL))
  pd   <- TreeTools::MatrixToPhyDat(mat)
  mkd  <- MkPrimeData(pd)

  model <- MkPrimeModel(
    kPrimePrior = "logseries",
    kprimeLogseriesC = 0.7,
    expSteps = 5
  )

  mcmc <- MkPrimeMCMC(nIter = 100L, maxWarmup = 50L, minWarmup = 50L,
                      thin = 5L, autoTune = FALSE)

  set.seed(4217)
  result <- RunMkPrime(mkd, tree = tree, model = model, mcmc = mcmc)

  expect_s3_class(result, "MkPosterior")
  # result$samples is a matrix; colnames gives the parameter names
  cols <- colnames(result$samples)

  # No p column
  expect_false("p" %in% cols)

  # kPrime columns present for transformational characters
  expect_true(any(grepl("^kPrime_", cols)))

  # Branch length columns present
  expect_true(any(grepl("^br_", cols)))
})


# ---------------------------------------------------------------------------
# 10. R and C++ agree on the logseries prior
# ---------------------------------------------------------------------------
test_that("LogPrior logseries agrees with cpp_log_prior", {
  tree <- TreeTools::Preorder(
    ape::read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  )
  # Two transformational characters with different kObs (2 and 3), so a
  # per-character error would not cancel across the sum.
  mat <- matrix(c(0, 1, 0, 1,
                  0, 1, 2, 0), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  mkd <- MkPrimeData(TreeTools::MatrixToPhyDat(mat))

  for (logseriesC in c(0.05, 0.3, 0.7, 0.95, 0.999)) {
    model <- MkPrimeModel(kPrimePrior = "logseries",
                          kprimeLogseriesC = logseriesC, expSteps = 10)
    model <- MkPrime:::.FinalizeModel(model, tree, mkd)
    dataPtr <- MkPrime:::.InitMcmcData(mkd, model)

    for (kp in list(c(2L, 3L), c(3L, 5L), c(4L, 8L), c(2L, 25L), c(20L, 30L))) {
      state <- list(
        tree           = tree,
        tree_length    = 0.5,
        rel_br_lengths = tree$edge.length / sum(tree$edge.length),
        rate_loss      = 1.0,
        rate_log_sd    = 0.2,
        rate_neo       = 1.0,
        kPrime         = kp,
        log_lik        = 0.0,
        log_prior      = 0.0
      )
      lpR <- MkPrime:::LogPrior(state, model, mkd)
      lpCpp <- eval_log_prior_cpp(dataPtr, MkPrime:::.InitMcmcChain(state))
      expect_equal(lpCpp, lpR, tolerance = 1e-10,
                   info = sprintf("c=%g, kPrime=%s", logseriesC,
                                  paste(kp, collapse = ",")))
    }
  }
})


test_that("LogPrior logseries and cpp_log_prior both reject k' < kObs", {
  tree <- TreeTools::Preorder(
    ape::read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  )
  mat <- matrix(c(0, 1, 0, 1,
                  0, 1, 2, 0), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  mkd <- MkPrimeData(TreeTools::MatrixToPhyDat(mat))
  model <- MkPrimeModel(kPrimePrior = "logseries", kprimeLogseriesC = 0.7,
                        expSteps = 10)
  model <- MkPrime:::.FinalizeModel(model, tree, mkd)
  dataPtr <- MkPrime:::.InitMcmcData(mkd, model)

  state <- list(
    tree           = tree,
    tree_length    = 0.5,
    rel_br_lengths = tree$edge.length / sum(tree$edge.length),
    rate_loss      = 1.0,
    rate_log_sd    = 0.2,
    rate_neo       = 1.0,
    kPrime         = c(2L, 2L),   # second character has kObs = 3
    log_lik        = 0.0,
    log_prior      = 0.0
  )
  expect_equal(MkPrime:::LogPrior(state, model, mkd), -Inf)
  expect_equal(eval_log_prior_cpp(dataPtr, MkPrime:::.InitMcmcChain(state)),
               -Inf)
})
