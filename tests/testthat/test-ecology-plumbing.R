library("TreeTools")


# Build an ecology-enabled MkPrimeData + tree fixture used by the plumbing tests.
.MakeEcologyFixture <- function() {
  set.seed(2026)
  tips <- paste0("t", 1:6)
  mat <- matrix(c(
    0, 1, 2, 0, 1, 2,   # transformational
    1, 0, 1, 2, 2, 0,
    2, 1, 0, 1, 0, 2,
    0, 1, 0, 1, 0, 1,   # neomorphic (binary)
    1, 0, 1, 0, 1, 0,
    0, 0, 1, 1, 2, 2    # ecology (3 states)
  ), nrow = 6, ncol = 6, byrow = FALSE,
     dimnames = list(tips, NULL))
  pd <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd, neomorphic = c(4L, 5L), ecology = 6L)
  tree <- Preorder(ape::rtree(6, tip.label = tips))
  list(mkd = mkd, tree = tree)
}


test_that(".InitState initialises phi, pi0, z for ecology-aware models", {
  f <- .MakeEcologyFixture()
  model <- MkPrimeModel(ecologyAware = TRUE, expSteps = 10,
                        kPrimePrior = "geometric")
  state <- MkPrime:::.InitState(f$tree, f$mkd, model)

  expect_equal(state$phi, 1.0)
  expect_equal(state$pi0, model$rho0Alpha / (model$rho0Alpha + model$rho0Beta))
  expect_true(is.matrix(state$z))
  expect_equal(dim(state$z), c(f$mkd$nChar, f$mkd$kEcology))
  expect_true(all(state$z == 0L))
  expect_true(is.finite(state$log_lik))
})


test_that(".InitState initialises phi as a kEcology vector in per_ecology mode", {
  f <- .MakeEcologyFixture()
  model <- MkPrimeModel(ecologyAware = TRUE, magnitudeMode = "per_ecology",
                        expSteps = 10, kPrimePrior = "geometric")
  state <- MkPrime:::.InitState(f$tree, f$mkd, model)

  expect_length(state$phi, f$mkd$kEcology)
  expect_true(all(state$phi == 1.0))
})


test_that(".InitState log_lik uses ecology-aware orchestrator when ecologyAware", {
  f <- .MakeEcologyFixture()
  model <- MkPrimeModel(ecologyAware = TRUE, expSteps = 10,
                        kPrimePrior = "geometric")
  state_eco <- MkPrime:::.InitState(f$tree, f$mkd, model)

  # With z = 0 the ecology likelihood equals the standard one
  # (modulo coding correction — orchestrator currently runs without it).
  ll_ref <- MkPrime:::.MkpEcologyLogLikelihood(
    f$tree, f$mkd,
    kPrime = state_eco$kPrime,
    rate_loss = state_eco$rate_loss,
    rate_log_sd = state_eco$rate_log_sd,
    nCat = model$nCat,
    rate_neo = state_eco$rate_neo %||% 1.0,
    relabel = model$relabel,
    phi = state_eco$phi, zMat = state_eco$z,
    magnitudeMode = model$magnitudeMode
  )
  expect_equal(state_eco$log_lik, ll_ref, tolerance = 1e-12)
})


test_that(".InitState skips ecology fields when ecologyAware = FALSE", {
  f <- .MakeEcologyFixture()
  model <- MkPrimeModel(kPrimePrior = "geometric", expSteps = 10)
  state <- MkPrime:::.InitState(f$tree, f$mkd, model)

  expect_null(state$phi)
  expect_null(state$pi0)
  expect_null(state$z)
})


test_that("prepare_mcmc_data + init_mcmc_state round-trip ecology values", {
  f <- .MakeEcologyFixture()
  model <- MkPrimeModel(ecologyAware = TRUE, expSteps = 10,
                        kPrimePrior = "geometric")

  # Build partitions list matching prepare_mcmc_data's expected format.
  parts <- lapply(f$mkd$partitions, function(p) {
    list(type = p$type, k = p$k, kObs = p$kObs,
         char_indices = p$char_indices,
         tip_states = p$tip_states,
         unique_tip_states = p$unique_tip_states,
         pattern_index = p$pattern_index)
  })

  dataPtr <- prepare_mcmc_data(
    parts, as.integer(f$mkd$kObs), f$mkd$type,
    any(f$mkd$type == "neomorphic"),
    model$nCat, model$coding, model$relabel,
    model$treeLengthShape, model$treeLengthRate,
    model$rateLossMeanlog, model$rateLossSdlog,
    model$rateLogSdShape, model$rateLogSdRate,
    model$rateNeoMeanlog, model$rateNeoSdlog,
    model$kprimeHyperA, model$kprimeHyperB,
    identical(model$kPrimePrior, "logseries"),
    model$kprimeLogseriesC %||% 0.7,
    FALSE, FALSE, 4L, 1.0, 1.0,
    FALSE, numeric(0), 1L, 0L, 0, -1e308,
    TRUE,                                  # ecologyAware
    as.integer(f$mkd$ecology),
    as.integer(f$mkd$kEcology),
    "global",
    model$rho0Alpha, model$rho0Beta,
    model$sigmaPhi, as.integer(model$gibbsZEvery)
  )
  expect_true(inherits(dataPtr, "externalptr"))

  state <- MkPrime:::.InitState(f$tree, f$mkd, model)
  statePtr <- MkPrime:::.InitMcmcChain(state)
  cppState <- get_mcmc_state(statePtr)

  expect_equal(cppState$phi, state$phi)
  expect_equal(cppState$pi0, state$pi0)
  expect_equal(cppState$zMatrix, state$z)
})


test_that("init_mcmc_state allows empty phi (non-ecology mode)", {
  # When ecologyAware is FALSE, init_mcmc_state should accept empty phi.
  f <- .MakeEcologyFixture()
  model <- MkPrimeModel(kPrimePrior = "geometric", expSteps = 10)
  state <- MkPrime:::.InitState(f$tree, f$mkd, model)
  statePtr <- MkPrime:::.InitMcmcChain(state)
  cppState <- get_mcmc_state(statePtr)
  expect_length(cppState$phi, 0L)
})
