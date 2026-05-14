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
  expect_equal(dim(state$z), c(f$mkd$nChar, f$mkd$kEcology - 1L))
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
    magnitudeMode = model$magnitudeMode,
    coding = model$coding
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
    model$sigmaPhi, as.integer(model$gibbsZEvery),
    model$thetaAlpha, model$thetaBeta
  )
  expect_true(inherits(dataPtr, "externalptr"))

  state <- MkPrime:::.InitState(f$tree, f$mkd, model)
  statePtr <- MkPrime:::.InitMcmcChain(state)
  cppState <- get_mcmc_state(statePtr)

  expect_equal(cppState$phi, state$phi)
  expect_equal(cppState$pi0, state$pi0)
  expect_equal(cppState$zMatrix, state$z)
})


# ===== C++ vs R-level ecology orchestrator =====


# Build the data pointer needed for cpp_log_likelihood_ecology tests.
.MakeEcoDataPtr <- function(mkd, model) {
  parts <- lapply(mkd$partitions, function(p) {
    list(type = p$type, k = p$k, kObs = p$kObs,
         char_indices = p$char_indices,
         tip_states = p$tip_states,
         unique_tip_states = p$unique_tip_states,
         pattern_index = p$pattern_index)
  })
  ecoTip <- as.integer(mkd$ecology)
  ecoTip[is.na(ecoTip)] <- -1L
  prepare_mcmc_data(
    parts, as.integer(mkd$kObs), mkd$type,
    any(mkd$type == "neomorphic"),
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
    TRUE,
    ecoTip, as.integer(mkd$kEcology),
    model$magnitudeMode %||% "global",
    model$rho0Alpha %||% 75.0, model$rho0Beta %||% 25.0,
    model$sigmaPhi %||% 1.0, as.integer(model$gibbsZEvery %||% 50L),
    model$thetaAlpha %||% 1.0, model$thetaBeta %||% 1.0
  )
}


test_that("cpp_log_likelihood_ecology matches R orchestrator (z = 0)", {
  f <- .MakeEcologyFixture()
  model <- MkPrimeModel(ecologyAware = TRUE, expSteps = 10,
                        kPrimePrior = "geometric", coding = "none")
  dataPtr <- .MakeEcoDataPtr(f$mkd, model)
  zMat <- matrix(0L, nrow = f$mkd$nChar, ncol = f$mkd$kEcology - 1L)

  parent <- f$tree$edge[, 1]
  child  <- f$tree$edge[, 2]
  edgeLen <- f$tree$edge.length
  kPrime <- as.integer(f$mkd$kObs)

  ll_cpp <- MkPrime:::.CppLogLikelihoodEcology(
    dataPtr, parent, child, edgeLen, kPrime,
    rateLoss = 1.0, rateLogSd = 0, rateNeo = 1.0,
    phi = 1.0, zMatrix = zMat)
  ll_r <- MkPrime:::.MkpEcologyLogLikelihood(
    f$tree, f$mkd, kPrime,
    rate_loss = 1.0, rate_log_sd = 0, nCat = model$nCat,
    rate_neo = 1.0, relabel = model$relabel,
    phi = 1.0, zMat = zMat, magnitudeMode = "global")
  expect_equal(ll_cpp, ll_r, tolerance = 1e-10)
})


test_that("cpp_log_likelihood_ecology matches R orchestrator (mixed z, global phi)", {
  f <- .MakeEcologyFixture()
  model <- MkPrimeModel(ecologyAware = TRUE, expSteps = 10,
                        kPrimePrior = "geometric", coding = "none")
  dataPtr <- .MakeEcoDataPtr(f$mkd, model)
  set.seed(99)
  zMat <- matrix(as.integer(sample(0:2, f$mkd$nChar * (f$mkd$kEcology - 1L),
                                    replace = TRUE)),
                 nrow = f$mkd$nChar, ncol = f$mkd$kEcology - 1L)

  parent <- f$tree$edge[, 1]
  child  <- f$tree$edge[, 2]
  edgeLen <- f$tree$edge.length

  # v2: pi0 default differs between C++ wrapper (0.0) and R orchestrator
  # (0.5); pass explicitly to align.
  ll_cpp <- MkPrime:::.CppLogLikelihoodEcology(
    dataPtr, parent, child, edgeLen,
    kPrime = as.integer(f$mkd$kObs),
    rateLoss = 1.2, rateLogSd = 0, rateNeo = 0.9,
    phi = 2.0, zMatrix = zMat, pi0 = 0.5)
  ll_r <- MkPrime:::.MkpEcologyLogLikelihood(
    f$tree, f$mkd, kPrime = as.integer(f$mkd$kObs),
    rate_loss = 1.2, rate_log_sd = 0, nCat = model$nCat,
    rate_neo = 0.9, relabel = model$relabel,
    phi = 2.0, zMat = zMat, magnitudeMode = "global")
  expect_equal(ll_cpp, ll_r, tolerance = 1e-10)
})


test_that("cpp_log_likelihood_ecology matches R orchestrator with ACRV", {
  f <- .MakeEcologyFixture()
  model <- MkPrimeModel(ecologyAware = TRUE, expSteps = 10,
                        kPrimePrior = "geometric", coding = "none",
                        nCat = 4L)
  dataPtr <- .MakeEcoDataPtr(f$mkd, model)
  set.seed(13)
  zMat <- matrix(as.integer(sample(0:2, f$mkd$nChar * (f$mkd$kEcology - 1L),
                                    replace = TRUE)),
                 nrow = f$mkd$nChar, ncol = f$mkd$kEcology - 1L)

  parent <- f$tree$edge[, 1]
  child  <- f$tree$edge[, 2]
  edgeLen <- f$tree$edge.length

  ll_cpp <- MkPrime:::.CppLogLikelihoodEcology(
    dataPtr, parent, child, edgeLen,
    kPrime = as.integer(f$mkd$kObs),
    rateLoss = 1.0, rateLogSd = 0.5, rateNeo = 1.0,
    phi = 1.8, zMatrix = zMat, pi0 = 0.5)
  ll_r <- MkPrime:::.MkpEcologyLogLikelihood(
    f$tree, f$mkd, kPrime = as.integer(f$mkd$kObs),
    rate_loss = 1.0, rate_log_sd = 0.5, nCat = 4L,
    rate_neo = 1.0, relabel = model$relabel,
    phi = 1.8, zMat = zMat, magnitudeMode = "global")
  expect_equal(ll_cpp, ll_r, tolerance = 1e-10)
})


test_that("cpp_log_likelihood_ecology matches R orchestrator with per_ecology phi", {
  f <- .MakeEcologyFixture()
  model <- MkPrimeModel(ecologyAware = TRUE, expSteps = 10,
                        kPrimePrior = "geometric", coding = "none",
                        magnitudeMode = "per_ecology")
  dataPtr <- .MakeEcoDataPtr(f$mkd, model)
  set.seed(31)
  zMat <- matrix(as.integer(sample(0:2, f$mkd$nChar * (f$mkd$kEcology - 1L),
                                    replace = TRUE)),
                 nrow = f$mkd$nChar, ncol = f$mkd$kEcology - 1L)
  phi <- c(1.4, 0.7, 2.1)
  expect_equal(length(phi), f$mkd$kEcology)

  parent <- f$tree$edge[, 1]
  child  <- f$tree$edge[, 2]
  edgeLen <- f$tree$edge.length

  ll_cpp <- MkPrime:::.CppLogLikelihoodEcology(
    dataPtr, parent, child, edgeLen,
    kPrime = as.integer(f$mkd$kObs),
    rateLoss = 1.0, rateLogSd = 0, rateNeo = 1.0,
    phi = phi, zMatrix = zMat, pi0 = 0.5)
  ll_r <- MkPrime:::.MkpEcologyLogLikelihood(
    f$tree, f$mkd, kPrime = as.integer(f$mkd$kObs),
    rate_loss = 1.0, rate_log_sd = 0, nCat = model$nCat,
    rate_neo = 1.0, relabel = model$relabel,
    phi = phi, zMat = zMat, magnitudeMode = "per_ecology")
  expect_equal(ll_cpp, ll_r, tolerance = 1e-10)
})


test_that("cpp_log_likelihood_ecology with variable coding matches R orchestrator", {
  f <- .MakeEcologyFixture()
  model <- MkPrimeModel(ecologyAware = TRUE, expSteps = 10,
                        kPrimePrior = "geometric", coding = "variable")
  dataPtr <- .MakeEcoDataPtr(f$mkd, model)
  set.seed(401)
  zMat <- matrix(as.integer(sample(0:2, f$mkd$nChar * (f$mkd$kEcology - 1L),
                                    replace = TRUE)),
                 nrow = f$mkd$nChar, ncol = f$mkd$kEcology - 1L)

  parent <- f$tree$edge[, 1]
  child  <- f$tree$edge[, 2]
  edgeLen <- f$tree$edge.length

  ll_cpp <- MkPrime:::.CppLogLikelihoodEcology(
    dataPtr, parent, child, edgeLen,
    kPrime = as.integer(f$mkd$kObs),
    rateLoss = 1.0, rateLogSd = 0, rateNeo = 1.0,
    phi = 1.8, zMatrix = zMat, pi0 = 0.5)
  ll_r <- MkPrime:::.MkpEcologyLogLikelihood(
    f$tree, f$mkd, kPrime = as.integer(f$mkd$kObs),
    rate_loss = 1.0, rate_log_sd = 0, nCat = model$nCat,
    rate_neo = 1.0, relabel = model$relabel,
    phi = 1.8, zMat = zMat, magnitudeMode = "global",
    coding = "variable")
  expect_equal(ll_cpp, ll_r, tolerance = 1e-10)
})


test_that("cpp_log_likelihood_ecology variable coding + z = 0 matches standard", {
  f <- .MakeEcologyFixture()
  model <- MkPrimeModel(ecologyAware = TRUE, expSteps = 10,
                        kPrimePrior = "geometric", coding = "variable")
  dataPtr <- .MakeEcoDataPtr(f$mkd, model)
  zMat <- matrix(0L, nrow = f$mkd$nChar, ncol = f$mkd$kEcology - 1L)

  parent <- f$tree$edge[, 1]
  child  <- f$tree$edge[, 2]
  edgeLen <- f$tree$edge.length

  # v2: z = 0 only matches the non-ecology baseline when phi = 1 (otherwise
  # the non-ref ecology rate factor 1/gamma_e drifts away from 1).
  ll_cpp <- MkPrime:::.CppLogLikelihoodEcology(
    dataPtr, parent, child, edgeLen,
    kPrime = as.integer(f$mkd$kObs),
    rateLoss = 1.0, rateLogSd = 0, rateNeo = 1.0,
    phi = 1.0, zMatrix = zMat)
  ll_baseline <- MkPrime:::.MkpLogLikelihood(
    f$tree, f$mkd, kPrime = as.integer(f$mkd$kObs),
    rate_loss = 1.0, rate_log_sd = 0, nCat = model$nCat,
    coding = "variable", rate_neo = 1.0, relabel = model$relabel)
  expect_equal(ll_cpp, ll_baseline, tolerance = 1e-10)
})


test_that("cpp_log_likelihood_ecology errors on informative coding", {
  f <- .MakeEcologyFixture()
  model <- MkPrimeModel(ecologyAware = TRUE, expSteps = 10,
                        kPrimePrior = "geometric", coding = "informative")
  dataPtr <- .MakeEcoDataPtr(f$mkd, model)
  zMat <- matrix(0L, nrow = f$mkd$nChar, ncol = f$mkd$kEcology - 1L)
  expect_error(
    MkPrime:::.CppLogLikelihoodEcology(
      dataPtr,
      f$tree$edge[, 1], f$tree$edge[, 2], f$tree$edge.length,
      as.integer(f$mkd$kObs),
      1.0, 0, 1.0, 1.0, zMat),
    "informative.*not yet"
  )
})


test_that("cpp_log_likelihood_ecology errors when ecology not enabled", {
  f <- .MakeEcologyFixture()
  model <- MkPrimeModel(ecologyAware = FALSE, expSteps = 10,
                        kPrimePrior = "geometric")
  # Build a non-ecology-aware data pointer
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
    FALSE, integer(0), 0L, "global",
    7.0, 3.0, 0.5, 50L
  )
  zMat <- matrix(0L, nrow = f$mkd$nChar, ncol = 1L)
  expect_error(
    MkPrime:::.CppLogLikelihoodEcology(
      dataPtr,
      f$tree$edge[, 1], f$tree$edge[, 2], f$tree$edge.length,
      as.integer(f$mkd$kObs),
      1.0, 0, 1.0,
      phi = 1.0, zMatrix = zMat),
    "ecologyAware = TRUE"
  )
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


# ===== Phase 3-phi: Bactrian on log(phi) =====


# Run a single phi Bactrian move (moveType = 34) and return the post-move
# C++ state plus the proposal seed details.  Used by the tests below.
.PhiStateAfterMove <- function(f, model, seed = 7L,
                               scaleTuning = 0.5, beta = 1.0) {
  dataPtr  <- .MakeEcoDataPtr(f$mkd, model)
  state    <- MkPrime:::.InitState(f$tree, f$mkd, model)
  statePtr <- MkPrime:::.InitMcmcChain(state)
  pre      <- get_mcmc_state(statePtr)
  set.seed(seed)
  accepted <- do_move_cpp(dataPtr, statePtr,
                          moveType = 34L, charIdx = 0L,
                          scaleTuning = scaleTuning,
                          betaSimplexTuning = 1.0,
                          intWalkWindow = 1L, beta = beta)
  post <- get_mcmc_state(statePtr)
  list(pre = pre, post = post, accepted = accepted,
       dataPtr = dataPtr, statePtr = statePtr, refState = state)
}


test_that("scale_phi move keeps logLik in sync with ecology likelihood (global)", {
  f <- .MakeEcologyFixture()
  model <- MkPrimeModel(ecologyAware = TRUE, expSteps = 10,
                        kPrimePrior = "geometric", coding = "none")
  res <- .PhiStateAfterMove(f, model, seed = 13L, scaleTuning = 0.4)

  expect_length(res$post$phi, 1L)
  # Bactrian draws are non-zero, so phi must move when accepted.
  if (res$accepted) {
    expect_false(isTRUE(all.equal(res$post$phi, res$pre$phi)))
  } else {
    expect_equal(res$post$phi, res$pre$phi)
    expect_equal(res$post$logLik, res$pre$logLik)
    expect_equal(res$post$logPrior, res$pre$logPrior)
  }

  # state->logLik must match a fresh ecology evaluation at the post-move phi.
  ll_fresh <- MkPrime:::.CppLogLikelihoodEcology(
    res$dataPtr,
    f$tree$edge[, 1], f$tree$edge[, 2], f$tree$edge.length,
    as.integer(f$mkd$kObs),
    rateLoss = res$post$rateLoss,
    rateLogSd = res$post$rateLogSd,
    rateNeo = res$post$rateNeo,
    phi = res$post$phi,
    zMatrix = res$post$zMatrix,
    pi0 = res$post$pi0, theta = res$post$theta)
  expect_equal(res$post$logLik, ll_fresh, tolerance = 1e-10)
})


test_that("scale_phi prior delta equals dlnorm difference (only phi moves)", {
  f <- .MakeEcologyFixture()
  model <- MkPrimeModel(ecologyAware = TRUE, expSteps = 10,
                        kPrimePrior = "geometric", coding = "none")
  res <- .PhiStateAfterMove(f, model, seed = 9L, scaleTuning = 0.4)
  if (!res$accepted) skip("Phi move was rejected at this seed; rerun with another.")

  # All other params unchanged.
  expect_equal(res$post$treeLength, res$pre$treeLength)
  expect_equal(res$post$rateLoss, res$pre$rateLoss)
  expect_equal(res$post$zMatrix, res$pre$zMatrix)
  expect_equal(res$post$pi0, res$pre$pi0)

  # Prior delta should be entirely the LogNormal(0, sigmaPhi) contribution.
  d_lp_phi <- dlnorm(res$post$phi, 0, model$sigmaPhi, log = TRUE) -
              dlnorm(res$pre$phi,  0, model$sigmaPhi, log = TRUE)
  expect_equal(res$post$logPrior - res$pre$logPrior, d_lp_phi,
               tolerance = 1e-10)
})


test_that("scale_phi in per_ecology mode mutates exactly one phi entry", {
  f <- .MakeEcologyFixture()
  model <- MkPrimeModel(ecologyAware = TRUE, magnitudeMode = "per_ecology",
                        expSteps = 10, kPrimePrior = "geometric",
                        coding = "none")
  res <- .PhiStateAfterMove(f, model, seed = 31L, scaleTuning = 0.4)
  expect_length(res$post$phi, f$mkd$kEcology)

  if (res$accepted) {
    nChanged <- sum(abs(res$post$phi - res$pre$phi) > 1e-12)
    expect_equal(nChanged, 1L)
    ll_fresh <- MkPrime:::.CppLogLikelihoodEcology(
      res$dataPtr,
      f$tree$edge[, 1], f$tree$edge[, 2], f$tree$edge.length,
      as.integer(f$mkd$kObs),
      rateLoss = res$post$rateLoss,
      rateLogSd = res$post$rateLogSd,
      rateNeo = res$post$rateNeo,
      phi = res$post$phi,
      zMatrix = res$post$zMatrix,
      pi0 = res$post$pi0, theta = res$post$theta)
    expect_equal(res$post$logLik, ll_fresh, tolerance = 1e-10)
  } else {
    expect_equal(res$post$phi, res$pre$phi)
  }
})


test_that("scale_phi declines outside ecology mode", {
  f <- .MakeEcologyFixture()
  model <- MkPrimeModel(kPrimePrior = "geometric", expSteps = 10)
  state    <- MkPrime:::.InitState(f$tree, f$mkd, model)
  statePtr <- MkPrime:::.InitMcmcChain(state)
  # Non-ecology data pointer.
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
    FALSE, integer(0), 0L, "global", 7.0, 3.0, 0.5, 50L
  )
  accepted <- do_move_cpp(dataPtr, statePtr,
                          moveType = 34L, charIdx = 0L,
                          scaleTuning = 0.5,
                          betaSimplexTuning = 1.0,
                          intWalkWindow = 1L, beta = 1.0)
  expect_false(accepted)
})


test_that("MkPrimeModel rejects ecologyAware + qHeterogeneity", {
  expect_error(
    MkPrimeModel(ecologyAware = TRUE, qHeterogeneity = TRUE,
                 expSteps = 10, kPrimePrior = "geometric"),
    "ecologyAware.*qHeterogeneity"
  )
})


# ===== Phase 3-pi0: logit-Bactrian on pi0 =====


# Run a single pi0 logit-Bactrian move (moveType = 35) on a chain with at
# least one non-zero z cell, so the prior actually depends on pi0.
.Pi0StateAfterMove <- function(f, model, zInit, seed = 7L,
                               scaleTuning = 0.5, beta = 1.0) {
  dataPtr  <- .MakeEcoDataPtr(f$mkd, model)
  state    <- MkPrime:::.InitState(f$tree, f$mkd, model)
  state$z  <- zInit
  storage.mode(state$z) <- "integer"
  state$log_prior <- MkPrime:::LogPrior(state, model, f$mkd)
  statePtr <- MkPrime:::.InitMcmcChain(state)
  pre <- get_mcmc_state(statePtr)
  set.seed(seed)
  accepted <- do_move_cpp(dataPtr, statePtr,
                          moveType = 35L, charIdx = 0L,
                          scaleTuning = scaleTuning,
                          betaSimplexTuning = 1.0,
                          intWalkWindow = 1L, beta = beta)
  post <- get_mcmc_state(statePtr)
  list(pre = pre, post = post, accepted = accepted,
       dataPtr = dataPtr, statePtr = statePtr)
}


test_that("scale_pi0 move stays in (0, 1) and leaves logLik untouched", {
  f <- .MakeEcologyFixture()
  model <- MkPrimeModel(ecologyAware = TRUE, expSteps = 10,
                        kPrimePrior = "geometric", coding = "none")
  zInit <- matrix(0L, f$mkd$nChar, f$mkd$kEcology - 1L)
  zInit[1, 1] <- 1L  # one slab cell so pi0 matters in the prior
  zInit[2, 2] <- 2L
  res <- .Pi0StateAfterMove(f, model, zInit, seed = 17L, scaleTuning = 0.5)

  expect_true(res$post$pi0 > 0 && res$post$pi0 < 1)
  expect_equal(res$post$logLik, res$pre$logLik)        # prior-only move
  if (res$accepted) {
    expect_false(isTRUE(all.equal(res$post$pi0, res$pre$pi0)))
  } else {
    expect_equal(res$post$pi0, res$pre$pi0)
    expect_equal(res$post$logPrior, res$pre$logPrior)
  }
})


test_that("scale_pi0 prior delta matches dbeta + spike-and-slab contribution", {
  f <- .MakeEcologyFixture()
  model <- MkPrimeModel(ecologyAware = TRUE, expSteps = 10,
                        kPrimePrior = "geometric", coding = "none")
  zInit <- matrix(0L, f$mkd$nChar, f$mkd$kEcology - 1L)
  zInit[1, 1] <- 1L
  zInit[3, 2] <- 2L
  zInit[4, 1] <- 1L
  res <- .Pi0StateAfterMove(f, model, zInit, seed = 22L, scaleTuning = 0.5)
  if (!res$accepted) skip("Pi0 move rejected at this seed; rerun with another.")

  nNone <- sum(zInit == 0L)
  nSlab <- length(zInit) - nNone

  # Analytical prior delta: dbeta(pi0_new) - dbeta(pi0_old)
  # plus the spike-and-slab change: nNone * (log(p_new) - log(p_old))
  #                                + nSlab * (log(1 - p_new) - log(1 - p_old))
  # The slab includes a -log(2) per cell, which cancels in the delta.
  pOld <- res$pre$pi0
  pNew <- res$post$pi0
  d_beta  <- dbeta(pNew, model$rho0Alpha, model$rho0Beta, log = TRUE) -
             dbeta(pOld, model$rho0Alpha, model$rho0Beta, log = TRUE)
  d_spike <- nNone * (log(pNew) - log(pOld)) +
             nSlab * (log1p(-pNew) - log1p(-pOld))
  expect_equal(res$post$logPrior - res$pre$logPrior,
               d_beta + d_spike, tolerance = 1e-10)
})


# ===== Phase 3f scaffolding: per-character ecology log-likelihoods =====


# Invariant: sum of .CppLogLikelihoodEcologyPerChar equals
# .CppLogLikelihoodEcology under all combinations of z, ACRV, phi mode,
# and coding.  This guards the per-character helper that the z Gibbs
# sweep is built on.
.PerCharCases <- list(
  list(name = "z = 0, coding = none, phi global",
       coding = "none", magMode = "global", rateLogSd = 0,
       z = function(nC, kE) matrix(0L, nC, kE),
       phi = function(kE, mode)
         if (mode == "per_ecology") rep(1.0, kE) else 1.0),
  list(name = "z mixed, coding = variable, phi global",
       coding = "variable", magMode = "global", rateLogSd = 0,
       z = function(nC, kE) {
         set.seed(101)
         matrix(as.integer(sample(0:2, nC * kE, replace = TRUE)),
                nrow = nC, ncol = kE)
       },
       phi = function(kE, mode)
         if (mode == "per_ecology") rep(1.0, kE) else 1.6),
  list(name = "z mixed, ACRV, phi per_ecology, coding = none",
       coding = "none", magMode = "per_ecology", rateLogSd = 0.5,
       z = function(nC, kE) {
         set.seed(202)
         matrix(as.integer(sample(0:2, nC * kE, replace = TRUE)),
                nrow = nC, ncol = kE)
       },
       phi = function(kE, mode) {
         set.seed(303)
         exp(rnorm(kE, 0, 0.4))
       }),
  list(name = "z mixed, ACRV, variable coding, per_ecology phi",
       coding = "variable", magMode = "per_ecology", rateLogSd = 0.5,
       z = function(nC, kE) {
         set.seed(404)
         matrix(as.integer(sample(0:2, nC * kE, replace = TRUE)),
                nrow = nC, ncol = kE)
       },
       phi = function(kE, mode) {
         set.seed(505)
         exp(rnorm(kE, 0, 0.3))
       })
)

for (case in .PerCharCases) {
  local({
    cs <- case
    test_that(paste0("per-char ecology log-liks sum to total (", cs$name, ")"), {
      f <- .MakeEcologyFixture()
      model <- MkPrimeModel(ecologyAware = TRUE, expSteps = 10,
                            kPrimePrior = "geometric",
                            coding = cs$coding,
                            magnitudeMode = cs$magMode)
      dataPtr <- .MakeEcoDataPtr(f$mkd, model)
      # v2: zMatrix is nChar x (kEcology - 1); phi length matches kEcology.
      zMat <- cs$z(f$mkd$nChar, f$mkd$kEcology - 1L)
      phi  <- cs$phi(f$mkd$kEcology, cs$magMode)

      parent <- f$tree$edge[, 1]
      child  <- f$tree$edge[, 2]
      edgeLen <- f$tree$edge.length
      kPrime <- as.integer(f$mkd$kObs)

      ll_total <- MkPrime:::.CppLogLikelihoodEcology(
        dataPtr, parent, child, edgeLen, kPrime,
        rateLoss = 1.0, rateLogSd = cs$rateLogSd, rateNeo = 1.0,
        phi = phi, zMatrix = zMat)
      ll_per <- MkPrime:::.CppLogLikelihoodEcologyPerChar(
        dataPtr, parent, child, edgeLen, kPrime,
        rateLoss = 1.0, rateLogSd = cs$rateLogSd, rateNeo = 1.0,
        phi = phi, zMatrix = zMat)
      expect_length(ll_per, f$mkd$nChar)
      expect_equal(sum(ll_per), ll_total, tolerance = 1e-10)
    })
  })
}


# ===== Phase 3f: z Gibbs sweep =====


test_that("gibbs_z sweep updates z and syncs logLik/logPrior to fresh eval", {
  f <- .MakeEcologyFixture()
  model <- MkPrimeModel(ecologyAware = TRUE, expSteps = 10,
                        kPrimePrior = "geometric", coding = "none")
  dataPtr <- .MakeEcoDataPtr(f$mkd, model)
  state    <- MkPrime:::.InitState(f$tree, f$mkd, model)
  state$z[1, 1] <- 1L
  state$z[2, 2] <- 2L
  storage.mode(state$z) <- "integer"
  state$log_prior <- MkPrime:::LogPrior(state, model, f$mkd)
  statePtr <- MkPrime:::.InitMcmcChain(state)
  pre <- get_mcmc_state(statePtr)

  set.seed(7L)
  accepted <- do_move_cpp(dataPtr, statePtr,
                          moveType = 36L, charIdx = 0L,
                          scaleTuning = 0.5,
                          betaSimplexTuning = 1.0,
                          intWalkWindow = 1L, beta = 1.0)
  expect_true(accepted)
  post <- get_mcmc_state(statePtr)

  # All z entries are valid spike-and-slab values.
  expect_true(all(post$zMatrix %in% 0:2))

  # logLik / logPrior match a fresh evaluation at the post-sweep state.
  ll_fresh <- MkPrime:::.CppLogLikelihoodEcology(
    dataPtr,
    f$tree$edge[, 1], f$tree$edge[, 2], f$tree$edge.length,
    as.integer(post$kPrime),
    rateLoss = post$rateLoss, rateLogSd = post$rateLogSd,
    rateNeo = post$rateNeo,
    phi = post$phi, zMatrix = post$zMatrix,
    pi0 = post$pi0, theta = post$theta)
  expect_equal(post$logLik, ll_fresh, tolerance = 1e-10)

  # Some cells likely changed (sanity — not strictly required, but the
  # sweep should not leave z untouched given the spike-and-slab prior
  # and finite data on the small fixture).
  expect_true(!identical(post$zMatrix, pre$zMatrix))
})


test_that("gibbs_z sweep collapses z to the spike at large pi0 + no data weight", {
  # Test: with pi0 ~ 1 and beta = 0 (so the likelihood is ignored), the
  # Gibbs full conditional is dominated by the spike and almost every
  # cell should land in z = 0.  This isolates the prior-side logic.
  f <- .MakeEcologyFixture()
  model <- MkPrimeModel(ecologyAware = TRUE, expSteps = 10,
                        kPrimePrior = "geometric", coding = "none")
  dataPtr <- .MakeEcoDataPtr(f$mkd, model)
  state    <- MkPrime:::.InitState(f$tree, f$mkd, model)
  state$pi0 <- 0.999  # spike dominates
  state$z[] <- 1L     # everything starts in encouraged
  storage.mode(state$z) <- "integer"
  state$log_prior <- MkPrime:::LogPrior(state, model, f$mkd)
  statePtr <- MkPrime:::.InitMcmcChain(state)

  set.seed(123L)
  accepted <- do_move_cpp(dataPtr, statePtr,
                          moveType = 36L, charIdx = 0L,
                          scaleTuning = 0.5,
                          betaSimplexTuning = 1.0,
                          intWalkWindow = 1L, beta = 0.0)
  expect_true(accepted)
  post <- get_mcmc_state(statePtr)

  # With beta = 0, P(z = 0) = pi0 = 0.999.  Out of nChar*kEco = 18 cells,
  # essentially all should be 0; allow at most 2 slabs (≈ 2% tail).
  nSlab <- sum(post$zMatrix != 0L)
  expect_lt(nSlab, 3L)
})


test_that("gibbs_z declines outside ecology mode", {
  f <- .MakeEcologyFixture()
  model <- MkPrimeModel(kPrimePrior = "geometric", expSteps = 10)
  state    <- MkPrime:::.InitState(f$tree, f$mkd, model)
  statePtr <- MkPrime:::.InitMcmcChain(state)
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
    FALSE, integer(0), 0L, "global", 7.0, 3.0, 0.5, 50L
  )
  accepted <- do_move_cpp(dataPtr, statePtr,
                          moveType = 36L, charIdx = 0L,
                          scaleTuning = 0.5,
                          betaSimplexTuning = 1.0,
                          intWalkWindow = 1L, beta = 1.0)
  expect_false(accepted)
})


test_that("scale_pi0 declines outside ecology mode", {
  f <- .MakeEcologyFixture()
  model <- MkPrimeModel(kPrimePrior = "geometric", expSteps = 10)
  state    <- MkPrime:::.InitState(f$tree, f$mkd, model)
  statePtr <- MkPrime:::.InitMcmcChain(state)
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
    FALSE, integer(0), 0L, "global", 7.0, 3.0, 0.5, 50L
  )
  accepted <- do_move_cpp(dataPtr, statePtr,
                          moveType = 35L, charIdx = 0L,
                          scaleTuning = 0.5,
                          betaSimplexTuning = 1.0,
                          intWalkWindow = 1L, beta = 1.0)
  expect_false(accepted)
})
