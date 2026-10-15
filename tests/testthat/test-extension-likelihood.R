# Extension likelihood terms and custom moves (#441 S2b, extension ABI v2).
# toy_lik_extension_cpp(): N = 3 Poisson events at rate toy_lambda along the
# tree, toy_lambda ~ Gamma(2, 4), so toy_lambda | T ~ Gamma(5, 4 + T) for
# tree length T. Moves: "gibbs" (draw from that conditional) and "jump".

ToyLik <- function() MkPrimeExtension(function() toy_lik_extension_cpp())

ToyLogLik <- function(lambda, treeLength) {
  stats::dpois(3L, lambda * treeLength, log = TRUE)
}

ToyLogPrior <- function(lambda) stats::dgamma(lambda, 2, 4, log = TRUE)

test_that("a likelihood extension describes its term and moves", {
  info <- extension_info_cpp(toy_lik_extension_cpp())
  expect_identical(info$abi_version, 2L)
  expect_true(info$log_lik)
  expect_identical(info$moves, c("gibbs", "jump"))
  expect_identical(info$move_weights, c(1, 0.5))
  ext <- ToyLik()
  expect_true(ext$log_lik)
  expect_identical(ext$moves, c(gibbs = 1, jump = 0.5))
})

test_that("extension moves join, and core-likelihood moves leave", {
  model <- MkPrimeModel(extensions = ToyLik())
  core <- list(
    list(name = "nni", type = "nni", weight = 3),
    list(name = "gibbs_spr", type = "gibbs_spr", weight = 2),
    list(name = "weighted_spr", type = "weighted_spr", weight = 1),
    list(name = "slice_tree_length", type = "slice", weight = 1),
    list(name = "slice_rate_loss", type = "slice", weight = 1)
  )
  moves <- MkPrime:::.WithExtensionMoves(core, model)
  names <- vapply(moves, `[[`, "", "name")
  expect_identical(names, c("nni", "slice_rate_loss",
                            "scale_param_toy_lambda", "toylik_gibbs",
                            "toylik_jump"))
  own <- moves[4:5]
  expect_identical(vapply(own, MkPrime:::.MoveTypeCode, 0L), c(41L, 41L))
  expect_identical(vapply(own, `[[`, 0L, "sliceParamIdx"), 0:1)
  expect_identical(vapply(own, `[[`, 0, "weight"), c(1, 0.5))
  # A prior-only extension leaves the core moves alone.
  priorOnly <- MkPrimeModel(extensions = MkPrimeExtension(
    function() toy_extension_cpp()))
  expect_length(MkPrime:::.WithExtensionMoves(core, priorOnly), 7L)
})

test_that("the extension likelihood is tracked through every kind of move", {
  fx <- ExtFixture()
  attach_extensions_cpp(fx$dataPtr, list(toy_lik_extension_cpp()))
  statePtr <- ExtState(fx, 0.5)
  T0 <- sum(fx$tree$edge.length)
  expect_equal(refresh_extension_log_lik_cpp(fx$dataPtr, statePtr),
               ToyLogLik(0.5, T0), tolerance = 1e-12)
  start <- get_mcmc_state(statePtr)
  expect_equal(start$extLogLik, ToyLogLik(0.5, T0), tolerance = 1e-12)
  expect_equal(start$logPost, start$logLik + start$extLogLik + start$logPrior,
               tolerance = 1e-12)

  Check <- function(label) {
    s <- get_mcmc_state(statePtr)
    expect_equal(s$extLogLik, ToyLogLik(s$extParams, s$treeLength),
                 tolerance = 1e-10, label = label)
    expect_equal(s$logLik, eval_full_loglik_cpp(fx$dataPtr, statePtr),
                 tolerance = 1e-8, label = label)
    expect_equal(s$logPrior, eval_log_prior_cpp(fx$dataPtr, statePtr),
                 tolerance = 1e-10, label = label)
    s
  }
  Move <- function(type, idx = 0L, beta = 1.0) {
    do_move_cpp(fx$dataPtr, statePtr, type, idx, 0.6, 1.0, 1L, beta)
  }

  set.seed(4416L)
  # Gibbs draws from the exact conditional: always accepted in the cold chain.
  for (i in 1:20) {
    before <- get_mcmc_state(statePtr)
    expect_true(Move(41L, 0L))
    after <- Check("gibbs")
    expect_false(identical(after$extParams, before$extParams))
    expect_identical(after$logLik, before$logLik)
    expect_identical(after$treeLength, before$treeLength)
  }
  # The random walk is accepted sometimes, and rejection restores the block.
  acc <- logical(40)
  for (i in seq_along(acc)) {
    before <- get_mcmc_state(statePtr)
    acc[i] <- Move(41L, 1L)
    after <- Check("jump")
    if (!acc[i]) {
      expect_identical(after$extParams, before$extParams)
      expect_identical(after$extLogLik, before$extLogLik)
    }
  }
  expect_true(any(acc) && !all(acc))
  # Tree length (and branch, topology) moves rescore the extension term.
  tl <- logical(30)
  for (i in seq_along(tl)) {
    tl[i] <- Move(0L)
    Check("tree_length")
    Move(4L)
    Check("beta_simplex")
    Move(5L)
    Check("nni")
    Move(40L, length(MkPrime:::.kScalarParams))
    Check("scale_param")
  }
  expect_true(any(tl))
  # Tempered chains run the same moves.
  for (i in 1:10) {
    Move(41L, 0L, beta = 0.4)
    Check("tempered gibbs")
  }

  # Refused: moves scoring the core likelihood alone, and unknown moves.
  before <- get_mcmc_state(statePtr)
  for (type in 10:15) expect_false(Move(type))
  expect_false(Move(41L, 2L))
  expect_false(Move(41L, -1L))
  expect_identical(get_mcmc_state(statePtr), before)
})

test_that("exact Gibbs draws of an extension parameter hit the target", {
  skip_on_cran()
  fx <- ExtFixture()
  attach_extensions_cpp(fx$dataPtr, list(toy_lik_extension_cpp()))
  statePtr <- ExtState(fx, 0.5)
  # toy_lambda | T ~ Gamma(5, 4 + T): after each draw toy_lambda * (4 + T) is
  # a fresh Gamma(5, 1) variate, whatever tree-length moves did in between.
  set.seed(4426L)
  n <- 4000L
  z <- numeric(n)
  for (i in seq_len(n)) {
    do_move_cpp(fx$dataPtr, statePtr, 0L, 0L, 0.6, 1.0, 1L, 1.0)
    expect_true(do_move_cpp(fx$dataPtr, statePtr, 41L, 0L, 0.6, 1.0, 1L, 1.0))
    s <- get_mcmc_state(statePtr)
    z[i] <- s$extParams * (4 + s$treeLength)
  }
  expect_lt(abs(mean(z) - 5), 0.15)
  expect_lt(abs(var(z) - 5), 0.8)
})

test_that("RunMkPrime runs an extension likelihood", {
  skip_on_cran()
  skip_under_memcheck()
  d <- ExtRunData()
  model <- MkPrimeModel(extensions = list(ToyLik()))
  set.seed(4417L)
  res <- allow_warning(RunMkPrime(
    d$pd, d$tree, model = model,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 6000L, thin = 5L,
                       minWarmup = 1000L, maxWarmup = 1000L,
                       autoTune = FALSE, maxTime = 120)
  ), "maxWarmup|stabilis")
  s <- res$samples
  expect_true("toy_lambda" %in% colnames(s))
  expect_true(all(is.finite(s[, "log_likelihood"])))
  expect_false(any(c("gibbs_spr", "gibbs_subtree_swap") %in%
                     names(res$acceptance)))
  expect_gt(res$acceptance[["toylik_gibbs"]], 0.9)

  expect_error(mkp_stepping_stone(d$pd, d$tree, model = model),
               "extension likelihoods")
})
