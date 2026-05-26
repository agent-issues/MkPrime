# Plan v4 §5.1 + §5.4: per-class prior contributions.
#
# Four test groups (mirroring the §7b structure for the likelihood):
#
#   (A) Trivial spec equivalence: partitioned prior == legacy prior to ~1e-10.
#   (B) Multi-class finite prior at the .InitStatePartitioned starting point.
#   (C) Sensitivity: perturbing one class_rate_log_sd[c] changes the prior by
#       exactly the Gamma log-density difference (sanity check for §5.4).
#   (D) Cross-validation: R-side LogPrior() agrees with C++
#       eval_log_prior_partitioned_cpp() to 1e-10.
#
# All tests use the same helper .setup_prior() which builds a minimal mkd,
# tree, model, and McmcData XPtr in one call.

library("TreeTools")


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

.setup_prior <- function(partition = NULL, nChar = 8L, nTip = 6L, seed = 42L) {
  set.seed(seed)
  mat <- matrix(sample(0:1, nTip * nChar, replace = TRUE),
                nrow = nTip, ncol = nChar,
                dimnames = list(paste0("t", seq_len(nTip)), NULL))
  # Ensure every character is variable (at least two distinct states)
  for (j in seq_len(nChar)) {
    if (length(unique(mat[, j])) < 2L) mat[1L, j] <- 1L - mat[1L, j]
  }
  pd  <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)

  if (!is.null(partition)) {
    mkd$partitions <- .BuildPartitions(mkd, partition = partition)
  }

  tree <- TreeTools::Preorder(
    TreeTools::NJTree(pd, edgeLengths = TRUE) %||%
    TreeTools::RandomTree(pd, root = TRUE)
  )
  if (is.null(tree$edge.length) || any(tree$edge.length <= 0)) {
    tree$edge.length <- rep(0.1, nrow(tree$edge))
  }

  # Use "geometric" kPrime prior so that both the R LogPrior() and the C++
  # cpp_log_prior use the same prior path. The empirical_geometric default
  # requires passing the full empirical body to prepare_mcmc_data, which adds
  # complexity orthogonal to what this test file validates.
  model <- MkPrimeModel(kPrimePrior = "geometric")
  # Finalise so treeLengthRate is set (needed by LogPrior)
  model <- .FinalizeModel(model, tree, mkd)

  dataPtr <- prepare_mcmc_data(
    partitions_r             = mkd$partitions,
    kObs_r                   = mkd$kObs,
    charTypes_r              = mkd$type,
    hasNeo                   = any(mkd$type == "neomorphic"),
    nCat                     = model$nCat,
    codingStr                = model$coding,
    relabelFlag              = isTRUE(model$relabel),
    treeLengthShape          = model$treeLengthShape,
    treeLengthRate           = model$treeLengthRate %||% 1,
    rateLossMeanlog          = model$rateLossMeanlog,
    rateLossSdlog            = model$rateLossSdlog,
    rateLogSdShape           = model$rateLogSdShape,
    rateLogSdRate            = model$rateLogSdRate,
    rateNeoMeanlog           = model$rateNeoMeanlog,
    rateNeoSdlog             = model$rateNeoSdlog,
    kprimeHyperA             = model$kprimeHyperA,
    kprimeHyperB             = model$kprimeHyperB,
    kPriorLogseries          = identical(model$kPrimePrior, "logseries"),
    kprimeLogseriesC         = model$kprimeLogseriesC,
    kPriorBetaGeometric      = identical(model$kPrimePrior, "beta_geometric"),
    qHeterogeneity           = isTRUE(model$qHeterogeneity),
    nBetaCat                 = model$nBetaCat,
    betaScaleShape           = model$betaScaleShape,
    betaScaleRate            = model$betaScaleRate,
    kPriorEmpiricalGeometric = FALSE,
    empLogBody               = numeric(0)
  )

  # Build a minimal state list (mirrors what .InitState produces)
  state <- .InitState(tree, mkd, model)

  list(
    mkd     = mkd,
    tree    = tree,
    model   = model,
    dataPtr = dataPtr,
    state   = state,
    kPrime  = as.integer(ifelse(mkd$type == "known", mkd$known_k, mkd$kObs))
  )
}


# Helper: build an XPtr<McmcState> from a state list for eval_log_prior_cpp
.make_state_ptr <- function(s_list, dataPtr) {
  state <- s_list$state
  init_mcmc_state(
    state$tree$edge[, 1], state$tree$edge[, 2],
    state$rel_br_lengths, state$tree_length,
    state$rate_loss,
    state$rate_log_sd,
    state$rate_neo %||% 1.0,
    state$p %||% 0.5,
    as.integer(state$kPrime),
    state$log_lik, state$log_prior,
    state$beta_scale %||% 1.0,
    state$kprime_alpha %||% 1.0,
    state$kprime_beta  %||% 1.0
  )
}


# ---------------------------------------------------------------------------
# (A) Trivial spec equivalence
# ---------------------------------------------------------------------------

test_that("(A) trivial spec: eval_log_prior_partitioned_cpp == eval_log_prior_cpp to 1e-10", {
  s <- .setup_prior()
  statePtr <- .make_state_ptr(s, s$dataPtr)

  lp_legacy <- eval_log_prior_cpp(s$dataPtr, statePtr)
  lp_part   <- eval_log_prior_partitioned_cpp(
    s$dataPtr, statePtr,
    classRateLogSd = s$state$rate_log_sd,  # length-1, == rate_log_sd
    classW         = 1.0,                   # length-1 trivial simplex
    etaNeo         = 1.0
  )

  expect_true(is.finite(lp_legacy))
  expect_equal(lp_part, lp_legacy, tolerance = 1e-10)
})


test_that("(A) trivial spec: R-side LogPrior(state, model, mkd) matches legacy eval_log_prior_cpp", {
  s <- .setup_prior()
  statePtr <- .make_state_ptr(s, s$dataPtr)

  lp_cpp <- eval_log_prior_cpp(s$dataPtr, statePtr)
  lp_r   <- LogPrior(s$state, s$model, s$mkd)

  expect_true(is.finite(lp_cpp))
  expect_equal(lp_r, lp_cpp, tolerance = 1e-10)
})


# ---------------------------------------------------------------------------
# (B) Multi-class finite prior at .InitStatePartitioned starting point
# ---------------------------------------------------------------------------

test_that("(B) multi-class partition: partitioned prior is finite at w = nChar_c / nChar", {
  partition <- c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L)
  s <- .setup_prior(partition = partition)

  # Mimic .InitStatePartitioned per-class initialisation
  nChar   <- s$mkd$nChar
  nChar_c <- tabulate(partition, nbins = 2L)
  class_w <- nChar_c / nChar   # simplex, class_rate == 1 by construction
  class_rls <- rep(s$state$rate_log_sd, 2L)   # unlinked: equal at start

  statePtr <- .make_state_ptr(s, s$dataPtr)

  lp_part <- eval_log_prior_partitioned_cpp(
    s$dataPtr, statePtr,
    classRateLogSd = class_rls,
    classW         = class_w,
    etaNeo         = 1.0
  )

  expect_true(is.finite(lp_part))
})


test_that("(B) multi-class: R-side LogPrior with per-class fields is finite at initial point", {
  partition <- c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L)
  s <- .setup_prior(partition = partition)

  nChar   <- s$mkd$nChar
  nChar_c <- tabulate(partition, nbins = 2L)
  class_w <- nChar_c / nChar

  state2 <- s$state
  state2$class_w          <- class_w
  state2$class_rate_log_sd <- rep(state2$rate_log_sd, 2L)
  state2$class_rate        <- as.numeric(class_w) * nChar / as.numeric(nChar_c)
  state2$eta_neo           <- 1.0

  lp_r <- LogPrior(state2, s$model, s$mkd)
  expect_true(is.finite(lp_r))
})


# ---------------------------------------------------------------------------
# (C) Sensitivity: perturbing class_rate_log_sd[c] changes prior by the
#     expected Gamma log-density difference
# ---------------------------------------------------------------------------

test_that("(C) perturbing class_rate_log_sd[2] changes prior by Gamma density diff", {
  partition <- c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L)
  s <- .setup_prior(partition = partition)

  nChar   <- s$mkd$nChar
  nChar_c <- tabulate(partition, nbins = 2L)
  class_w <- nChar_c / nChar

  sd_base <- s$state$rate_log_sd   # e.g. 0.5
  sd_pert <- sd_base * 2.0          # perturbed class 2

  statePtr <- .make_state_ptr(s, s$dataPtr)

  lp_base <- eval_log_prior_partitioned_cpp(
    s$dataPtr, statePtr,
    classRateLogSd = c(sd_base, sd_base),
    classW         = class_w,
    etaNeo         = 1.0
  )
  lp_pert <- eval_log_prior_partitioned_cpp(
    s$dataPtr, statePtr,
    classRateLogSd = c(sd_base, sd_pert),
    classW         = class_w,
    etaNeo         = 1.0
  )

  # The difference should equal the Gamma log-density difference for class 2
  shape <- s$model$rateLogSdShape
  rate  <- s$model$rateLogSdRate
  expected_diff <- dgamma(sd_pert, shape = shape, rate = rate, log = TRUE) -
                   dgamma(sd_base, shape = shape, rate = rate, log = TRUE)

  expect_equal(lp_pert - lp_base, expected_diff, tolerance = 1e-10)
})


# ---------------------------------------------------------------------------
# (D) Cross-validation: R-side LogPrior() agrees with C++ to 1e-10
# ---------------------------------------------------------------------------

test_that("(D) R-side LogPrior agrees with eval_log_prior_partitioned_cpp to 1e-10", {
  partition <- c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L)
  s <- .setup_prior(partition = partition)

  nChar   <- s$mkd$nChar
  nChar_c <- tabulate(partition, nbins = 2L)
  class_w <- nChar_c / nChar
  class_rls <- c(s$state$rate_log_sd, s$state$rate_log_sd * 1.5)

  state2 <- s$state
  state2$class_w           <- class_w
  state2$class_rate_log_sd <- class_rls
  state2$class_rate        <- as.numeric(class_w) * nChar / as.numeric(nChar_c)
  state2$eta_neo           <- 1.0

  statePtr <- .make_state_ptr(s, s$dataPtr)

  lp_r   <- LogPrior(state2, s$model, s$mkd)
  lp_cpp <- eval_log_prior_partitioned_cpp(
    s$dataPtr, statePtr,
    classRateLogSd = class_rls,
    classW         = class_w,
    etaNeo         = 1.0
  )

  expect_true(is.finite(lp_r))
  expect_equal(lp_r, lp_cpp, tolerance = 1e-10)
})


# ---------------------------------------------------------------------------
# Regression: §7a bit-identity test must stay green (not repeated here —
# covered by test-partition-bitcompat-null.R which is run at every step).
# ---------------------------------------------------------------------------
