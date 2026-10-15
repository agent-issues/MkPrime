.InitValuesFixture <- function(neo = integer(0)) {
  set.seed(4120)
  tree <- Preorder(RandomTree(8, root = FALSE))
  tree$edge.length <- rexp(nrow(tree$edge), 5)
  mat <- matrix(sample(0:2, 8 * 12, replace = TRUE), 8,
                dimnames = list(tree$tip.label, NULL))
  mat[1:2, ] <- 0L
  mat[3, ] <- 1L
  mkd <- suppressWarnings(MkPrimeData(MatrixToPhyDat(mat), neomorphic = neo))
  model <- MkPrime:::.FinalizeModel(MkPrimeModel(), tree, mkd)
  list(tree = tree, mkd = mkd, model = model)
}

test_that("initValues = NULL leaves the initial state unchanged", {
  fx <- .InitValuesFixture()
  expect_identical(
    MkPrime:::.InitState(fx$tree, fx$mkd, fx$model, initValues = NULL),
    MkPrime:::.InitState(fx$tree, fx$mkd, fx$model)
  )
})

test_that("initValues set the state and its likelihood", {
  fx <- .InitValuesFixture()
  nEdge <- nrow(fx$tree$edge)
  rel <- seq_len(nEdge)
  kp <- fx$mkd$kObs + 1L
  st <- MkPrime:::.InitState(fx$tree, fx$mkd, fx$model, initValues = list(
    tree_length = 2.5, rel_br_lengths = rel, rate_loss = 1.7,
    rate_log_sd = 0.2, kPrime = kp
  ))
  expect_equal(st$tree_length, 2.5)
  expect_equal(st$rel_br_lengths, rel / sum(rel))
  expect_equal(st$tree$edge.length, 2.5 * rel / sum(rel))
  expect_identical(st$kPrime, as.integer(kp))
  expected <- MkPrime:::.MkpLogLikelihood(
    st$tree, fx$mkd, kPrime = kp, rate_loss = 1.7, rate_log_sd = 0.2,
    nCat = fx$model$nCat, coding = fx$model$coding, rate_neo = 1,
    relabel = fx$model$relabel
  )
  expect_identical(st$log_lik, expected)
  expect_identical(st$log_prior, LogPrior(st, fx$model, fx$mkd))
})

test_that("invalid initValues are rejected", {
  fx <- .InitValuesFixture()
  init <- function(...) {
    MkPrime:::.InitState(fx$tree, fx$mkd, fx$model, initValues = list(...))
  }
  expect_error(init(rate_los = 2), "Unknown")
  expect_error(init(tree_length = -1), "positive")
  expect_error(init(rate_log_sd = -0.1), "non-negative")
  expect_error(init(rel_br_lengths = 1:3), "one positive value per edge")
  expect_error(init(kPrime = fx$mkd$kObs - 1L), "valid k'")
  expect_error(init(rate_neo = 2), "no character is neomorphic")
  expect_error(MkPrime:::.InitState(fx$tree, fx$mkd, fx$model,
                                    initValues = list(2)), "named list")
})

test_that("RunMkPrime starts every chain from initValues", {
  fx <- .InitValuesFixture()
  mcmc <- MkPrimeMCMC(nIter = 20L, thin = 10L, minWarmup = 10L,
                      maxWarmup = 10L, autoTune = FALSE, nRuns = 1L,
                      nChains = 2L, logFile = tempfile(fileext = ".log"),
                      checkpointFile = tempfile(fileext = ".ckp"))
  captured <- NULL
  local_mocked_bindings(
    .RunMkPrimeSingleRun = function(mkd, model, mcmc, initialState, ...) {
      captured <<- initialState
      stop("captured")
    }
  )
  expect_error(
    RunMkPrime(fx$mkd, fx$tree, mcmc = mcmc,
               initValues = list(rate_loss = 3, rate_log_sd = 0.1)),
    "captured"
  )
  expect_identical(vapply(captured$chains, `[[`, 0, "rate_loss"), c(3, 3))
  expect_identical(vapply(captured$chains, `[[`, 0, "rate_log_sd"), c(0.1, 0.1))
})
