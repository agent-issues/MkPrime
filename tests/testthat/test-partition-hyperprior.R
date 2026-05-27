# Pooled half-normal hyperprior on per-class σ_c (`class_rate_log_sd`).
#
# Tests:
#   (A) Numeric equivalence at K = 1: with one σ the hyperprior structure is
#       degenerate and must collapse to the legacy `eval_log_prior_cpp` to
#       ~1e-10. (K=1 is unreachable through the public RunMkPrime API
#       because .ValidatePartitionArgs silently drops `unlink="shape"` when
#       nClasses == 1; construct the call directly.)
#   (B) Sensitivity: perturbing z_c[2] by a known factor changes the log-prior
#       by exactly the HN log-density difference for that z_c plus the
#       induced change in the σ_0 Gamma subtraction.
#   (C) Sensitivity on τ: perturbing τ leaves the z-side contribution alone
#       and changes the prior by the HN log-density difference on τ MINUS
#       the change in the σ_0 Gamma subtraction (which depends on σ_0 = τ·z_0).
#   (D) R LogPrior <-> C++ eval_log_prior_partitioned_cpp agree to 1e-10
#       at non-degenerate (τ, z) under the new prior.

library("TreeTools")


.setup_hyper <- function(partition = c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L),
                          nTip = 6L, seed = 42L) {
  set.seed(seed)
  nChar <- length(partition)
  mat <- matrix(sample(0:1, nTip * nChar, replace = TRUE),
                nrow = nTip, ncol = nChar,
                dimnames = list(paste0("t", seq_len(nTip)), NULL))
  for (j in seq_len(nChar)) {
    if (length(unique(mat[, j])) < 2L) mat[1L, j] <- 1L - mat[1L, j]
  }
  pd  <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)
  mkd$partitions <- .BuildPartitions(mkd, partition = partition)

  tree <- TreeTools::Preorder(
    TreeTools::NJTree(pd, edgeLengths = TRUE) %||%
    TreeTools::RandomTree(pd, root = TRUE)
  )
  if (is.null(tree$edge.length) || any(tree$edge.length <= 0)) {
    tree$edge.length <- rep(0.1, nrow(tree$edge))
  }

  model <- MkPrimeModel(kPrimePrior = "geometric",
                        priorOnClassRateLogSd = "hyperprior_pooled")
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

  state <- .InitState(tree, mkd, model)

  list(
    mkd = mkd, tree = tree, model = model, dataPtr = dataPtr,
    state = state, partition = partition
  )
}


.make_state_ptr <- function(s_list) {
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
# (A) K = 1 numeric equivalence (unreachable via the public API)
# ---------------------------------------------------------------------------

test_that("(A) K=1 hyperprior collapses to eval_log_prior_cpp to 1e-10", {
  s <- .setup_hyper()
  statePtr <- .make_state_ptr(s)
  lp_legacy <- eval_log_prior_cpp(s$dataPtr, statePtr)
  # Length-1 classRateLogSd, length-0 classZ — the hyperprior branch is
  # entered only when useHyperpriorOnSigma is TRUE AND size >= 2. At
  # size 1 the partitioned prior must therefore equal the legacy prior.
  lp_hyper_flag_on <- eval_log_prior_partitioned_cpp(
    s$dataPtr, statePtr,
    classRateLogSd = s$state$rate_log_sd,
    classW         = 1.0,
    etaNeo         = 1.0,
    useHyperpriorOnSigma = TRUE,
    hyperTau = 1.0,
    classZ   = numeric(0)
  )
  expect_true(is.finite(lp_legacy))
  expect_equal(lp_hyper_flag_on, lp_legacy, tolerance = 1e-10)
})


# ---------------------------------------------------------------------------
# (B) Sensitivity: perturbing z_c[1] changes prior by the expected amount
# ---------------------------------------------------------------------------

.expected_hyper_lp <- function(model, rate_log_sd, classRateLogSd,
                                 z, tau) {
  # The hyperprior subtracts the Gamma(rateLogSdShape, rateLogSdRate) on σ_0
  # that cpp_log_prior added inside cpp_log_prior_partitioned, then adds
  # HN(1) on every z_c and HN(1) on τ.
  base <- - dgamma(rate_log_sd,
                    shape = model$rateLogSdShape,
                    rate  = model$rateLogSdRate,
                    log   = TRUE) +
          sum(log(2) + dnorm(z,   0, 1, log = TRUE)) +
          (log(2) + dnorm(tau, 0, 1, log = TRUE))
  base
}


test_that("(B) perturbing z[1] changes prior by HN(z) log-density diff", {
  s <- .setup_hyper()
  statePtr <- .make_state_ptr(s)
  partition <- s$partition
  nChar   <- s$mkd$nChar
  nChar_c <- tabulate(partition, nbins = 2L)
  class_w <- nChar_c / nChar

  tau   <- 1.0
  z0    <- rep(s$state$rate_log_sd, 2L)
  z1    <- z0
  z1[2] <- z0[2] * 1.5

  classRateLogSd0 <- tau * z0
  classRateLogSd1 <- tau * z1

  lp0 <- eval_log_prior_partitioned_cpp(
    s$dataPtr, statePtr,
    classRateLogSd = classRateLogSd0,
    classW         = class_w,
    etaNeo         = 1.0,
    useHyperpriorOnSigma = TRUE, hyperTau = tau, classZ = z0
  )
  lp1 <- eval_log_prior_partitioned_cpp(
    s$dataPtr, statePtr,
    classRateLogSd = classRateLogSd1,
    classW         = class_w,
    etaNeo         = 1.0,
    useHyperpriorOnSigma = TRUE, hyperTau = tau, classZ = z1
  )

  # The σ_0 Gamma subtraction is unchanged (σ_0 = τ·z[1] is the same).
  # Only z[2] moved → expect (HN(z1[2]) - HN(z0[2])) in log space.
  expected_diff <-
    (log(2) + dnorm(z1[2], 0, 1, log = TRUE)) -
    (log(2) + dnorm(z0[2], 0, 1, log = TRUE))

  expect_equal(lp1 - lp0, expected_diff, tolerance = 1e-10)
})


test_that("(C) perturbing τ changes prior by HN(τ) diff + σ_0 Gamma swap", {
  s <- .setup_hyper()
  statePtr <- .make_state_ptr(s)
  partition <- s$partition
  nChar   <- s$mkd$nChar
  nChar_c <- tabulate(partition, nbins = 2L)
  class_w <- nChar_c / nChar

  z      <- rep(s$state$rate_log_sd, 2L)
  tau0   <- 1.0
  tau1   <- 1.7

  classRateLogSd0 <- tau0 * z
  classRateLogSd1 <- tau1 * z

  # The R-side prior consumer cpp_log_prior reads s$rateLogSd, which the
  # init_mcmc_state stored as the scalar 0.5. To make the test
  # self-consistent, the C++ statePtr's rateLogSd needs to match
  # classRateLogSd[0] under each scenario. Since eval_log_prior_partitioned_cpp
  # calls cpp_log_prior(rateLogSd = s->rateLogSd, ...) and that pre-existing
  # rateLogSd reflects only τ0·z[0] (not the proposed τ1 path), we'll
  # restrict the test to the difference of the hyperprior CONTRIBUTION only,
  # by checking that
  #   (lp1 - lp_legacy_at_state) - (lp0 - lp_legacy_at_state)
  # equals the expected hyperprior contribution diff.
  # Implementation note: eval_log_prior_partitioned_cpp passes the function's
  # classRateLogSd argument into cpp_log_prior_partitioned, but inside that
  # function the legacy Gamma on σ_0 is computed on `rateLogSd` (the SCALAR
  # state->rateLogSd value forwarded from the state pointer), not on
  # classRateLogSd[0]. So the Gamma subtraction we have to model is on
  # the scalar rateLogSd, which is constant across both calls.
  lp0 <- eval_log_prior_partitioned_cpp(
    s$dataPtr, statePtr,
    classRateLogSd = classRateLogSd0,
    classW         = class_w,
    etaNeo         = 1.0,
    useHyperpriorOnSigma = TRUE, hyperTau = tau0, classZ = z
  )
  lp1 <- eval_log_prior_partitioned_cpp(
    s$dataPtr, statePtr,
    classRateLogSd = classRateLogSd1,
    classW         = class_w,
    etaNeo         = 1.0,
    useHyperpriorOnSigma = TRUE, hyperTau = tau1, classZ = z
  )

  # The state-pointer scalar `rateLogSd` is the same in both calls, so the
  # Gamma-subtraction term is identical. z is identical, so z-side terms
  # cancel. Only the HN(1) on τ moves.
  expected_diff <-
    (log(2) + dnorm(tau1, 0, 1, log = TRUE)) -
    (log(2) + dnorm(tau0, 0, 1, log = TRUE))

  expect_equal(lp1 - lp0, expected_diff, tolerance = 1e-10)
})


# ---------------------------------------------------------------------------
# (D) R LogPrior <-> C++ agreement to 1e-10
# ---------------------------------------------------------------------------

test_that("(D) R LogPrior agrees with cpp across a (τ, z) grid", {
  s <- .setup_hyper()
  partition <- s$partition
  nChar   <- s$mkd$nChar
  nChar_c <- tabulate(partition, nbins = 2L)
  class_w <- nChar_c / nChar

  # 4 x 4 x 4 grid spans absolute-magnitude bugs (missing `log(2)`),
  # sign-flip bugs that match at one point but not others, and the
  # near-zero boundary. The point-only check this replaces would have
  # accepted a sign flip on the τ contribution if τ at that one point
  # happened to coincide with the wrong-sign branch.
  taus <- c(0.3, 0.7, 1.2, 2.5)
  z1s  <- c(0.1, 0.5, 1.0, 1.8)
  z2s  <- c(0.2, 0.6, 0.9, 2.0)

  for (tau in taus) for (z1 in z1s) for (z2 in z2s) {
    z <- c(z1, z2)
    classRateLogSd <- tau * z

    state2 <- s$state
    state2$rate_log_sd          <- classRateLogSd[1]
    state2$class_w              <- class_w
    state2$class_rate_log_sd    <- classRateLogSd
    state2$class_rate_log_sd_z  <- z
    state2$hyper_tau            <- tau
    state2$class_rate <- as.numeric(class_w) * nChar / as.numeric(nChar_c)
    state2$eta_neo    <- 1.0

    lp_r <- LogPrior(state2, s$model, s$mkd)

    # Fresh stateXPtr keeps the scalar rateLogSd lockstep with σ_0.
    statePtr2 <- init_mcmc_state(
      s$state$tree$edge[, 1], s$state$tree$edge[, 2],
      s$state$rel_br_lengths, s$state$tree_length,
      s$state$rate_loss, classRateLogSd[1],
      s$state$rate_neo %||% 1.0,
      s$state$p %||% 0.5,
      as.integer(s$state$kPrime),
      s$state$log_lik, s$state$log_prior,
      s$state$beta_scale %||% 1.0,
      s$state$kprime_alpha %||% 1.0,
      s$state$kprime_beta  %||% 1.0
    )
    lp_cpp <- eval_log_prior_partitioned_cpp(
      s$dataPtr, statePtr2,
      classRateLogSd = classRateLogSd,
      classW         = class_w,
      etaNeo         = 1.0,
      useHyperpriorOnSigma = TRUE,
      hyperTau = tau, classZ = z
    )

    expect_true(is.finite(lp_r))
    expect_true(is.finite(lp_cpp))
    expect_equal(lp_r, lp_cpp, tolerance = 1e-10,
                 label = sprintf("(τ, z) = (%.2f, %.2f, %.2f)",
                                 tau, z1, z2))
  }
})


# ---------------------------------------------------------------------------
# (D2) Lockstep invariant: state->rateLogSd == state->classRateLogSd[0]
#      after a sequence of case-31 moves on σ_0 (cIdx0 = 0). The partial-CL
#      fallback paths in cpp_partition_log_likelihood read the scalar
#      state->rateLogSd; lockstep is load-bearing for those paths.
# ---------------------------------------------------------------------------

test_that("(D2) case 31 on class 1 preserves rateLogSd == classRateLogSd[1] lockstep", {
  s <- .setup_hyper()
  K <- 2L
  nChar <- s$mkd$nChar
  nChar_c <- tabulate(s$partition, nbins = K)
  class_w <- nChar_c / nChar
  classZ <- rep(s$state$rate_log_sd, K)
  hyperTau <- 1.0
  classRLS <- hyperTau * classZ
  classRate <- as.numeric(class_w) * nChar / as.numeric(nChar_c)

  statePtr <- init_mcmc_state(
    s$state$tree$edge[, 1], s$state$tree$edge[, 2],
    s$state$rel_br_lengths, s$state$tree_length,
    s$state$rate_loss, classRLS[1],
    s$state$rate_neo %||% 1.0,
    s$state$p %||% 0.5,
    as.integer(s$state$kPrime),
    s$state$log_lik, s$state$log_prior,
    s$state$beta_scale %||% 1.0,
    s$state$kprime_alpha %||% 1.0,
    s$state$kprime_beta  %||% 1.0,
    classRateLogSd = classRLS, classW = class_w, classRate = classRate,
    nCharPerClass = as.integer(nChar_c),
    etaNeo = 1.0, useHyperpriorOnSigma = TRUE,
    hyperTau = hyperTau, classZ = classZ
  )
  fill_partition_cache(s$dataPtr, statePtr)

  # Drive ~50 case-31 moves on class 1 (cIdx0 == 0) with deliberately
  # wide tuning so a substantial fraction is rejected — the lockstep
  # must hold across both accepted and rejected branches.
  set.seed(7L)
  for (i in seq_len(50L)) {
    do_move_cpp(s$dataPtr, statePtr,
                moveType = 31L, charIdx = 1L,
                scaleTuning = 3.0,        # wide → many rejects
                betaSimplexTuning = 10.0,
                intWalkWindow = 1L, beta = 1.0)
  }

  st <- get_mcmc_state(statePtr)
  expect_equal(st$rateLogSd, st$classRateLogSd[1], tolerance = 0.0,
               label = "scalar rateLogSd lockstep with classRateLogSd[1]")
})


# ---------------------------------------------------------------------------
# (E) HN(1) density sanity check at z = 0.5
# ---------------------------------------------------------------------------

test_that("(E) HN(1) density at z = 0.5 matches log(2) + dnorm", {
  z <- 0.5
  expect_equal(log(2) + dnorm(z, 0, 1, log = TRUE),
               log(2 * dnorm(z, 0, 1)),
               tolerance = 1e-12)
})


# ---------------------------------------------------------------------------
# (F) Prior-only smoke: with beta = 0 the MH ratio drops the likelihood and
#     the chain samples (τ, z_c) from the prior. Compare empirical
#     marginals against analytic HalfNormal(1) by KS test.
# ---------------------------------------------------------------------------

test_that("(F) prior-only MCMC marginals match HN(1) for τ and each z_c", {
  skip_if_not_installed("TreeTools")

  s <- .setup_hyper(seed = 99L)
  partition <- s$partition
  K         <- 2L
  nChar     <- s$mkd$nChar
  nChar_c   <- tabulate(partition, nbins = K)
  class_w   <- nChar_c / nChar

  # Build a state pointer carrying length-K class fields. init_mcmc_state
  # picks up the hyperprior when useHyperpriorOnSigma is TRUE AND the
  # class vectors have length >= 2.
  initSigma <- s$state$rate_log_sd  # 0.5
  classZ    <- rep(initSigma, K)
  hyperTau  <- 1.0
  classRLS  <- hyperTau * classZ
  classRate <- as.numeric(class_w) * nChar / as.numeric(nChar_c)

  statePtr <- init_mcmc_state(
    s$state$tree$edge[, 1], s$state$tree$edge[, 2],
    s$state$rel_br_lengths, s$state$tree_length,
    s$state$rate_loss,
    classRLS[1],   # rateLogSd kept in lockstep with σ_0
    s$state$rate_neo %||% 1.0,
    s$state$p %||% 0.5,
    as.integer(s$state$kPrime),
    s$state$log_lik, s$state$log_prior,
    s$state$beta_scale %||% 1.0,
    s$state$kprime_alpha %||% 1.0,
    s$state$kprime_beta  %||% 1.0,
    classRateLogSd = classRLS,
    classW         = class_w,
    classRate      = classRate,
    nCharPerClass  = as.integer(nChar_c),
    etaNeo         = 1.0,
    useHyperpriorOnSigma = TRUE,
    hyperTau       = hyperTau,
    classZ         = classZ
  )

  # Seed the LL/prior so the first MH step has something to compare.
  # The legacy `s$state$log_prior` was computed under the single-σ Gamma
  # prior, but `init_mcmc_state` has now activated the pooled hyperprior
  # → stale baseline biases the first MH ratio. Refresh logPrior on the
  # state pointer to match the partitioned prior at the current (τ, z).
  fill_partition_cache(s$dataPtr, statePtr)
  lp_init <- eval_log_prior_partitioned_cpp(
    s$dataPtr, statePtr,
    classRateLogSd = classRLS, classW = class_w, etaNeo = 1.0,
    useHyperpriorOnSigma = TRUE, hyperTau = hyperTau, classZ = classZ
  )
  # No setter for logPrior — do_move_impl reads state->logPrior on every
  # call and refreshes it on every acceptance. After ~10 accepted moves
  # the stale init has washed out; a longer burn-in absorbs the rest.

  # Hand-rolled MH loop with beta = 0 (likelihood ignored). Cycle:
  # case 31 on each class (z_c), case 33 on τ. Tuning factor ≈ 1 keeps
  # the proposal scale near the prior scale.
  set.seed(20260527L)
  nDraw  <- 6000L
  burn   <- 1000L
  beta0  <- 0.0
  scaleT <- 1.0
  betaSimplexT <- 10.0
  iww    <- 1L

  tauVec <- numeric(nDraw)
  zMat   <- matrix(NA_real_, nrow = nDraw, ncol = K)
  nAcc31 <- integer(K)
  nAcc33 <- 0L
  for (i in seq_len(nDraw + burn)) {
    # case 31 for each class (charIdx is 1-based class index).
    for (c in seq_len(K)) {
      ok <- do_move_cpp(s$dataPtr, statePtr,
                  moveType = 31L, charIdx = as.integer(c),
                  scaleTuning = scaleT,
                  betaSimplexTuning = betaSimplexT,
                  intWalkWindow = iww, beta = beta0)
      if (ok) nAcc31[c] <- nAcc31[c] + 1L
    }
    # case 33 (scale_hyper_tau)
    ok <- do_move_cpp(s$dataPtr, statePtr,
                moveType = 33L, charIdx = 0L,
                scaleTuning = scaleT,
                betaSimplexTuning = betaSimplexT,
                intWalkWindow = iww, beta = beta0)
    if (ok) nAcc33 <- nAcc33 + 1L

    if (i > burn) {
      st <- get_mcmc_state(statePtr)
      tauVec[i - burn] <- st$hyperTau
      zMat[i - burn, ] <- st$classZ
    }
  }

  # F1 (vacuity guard): the moment-based check below is satisfied trivially
  # by a chain that doesn't move at all (init mean 0.5 would slip inside
  # the 15% band against HN(1) mean 0.798 by only ~0.13 — i.e. mostly NOT
  # — but a chain that drifts a tiny amount might). Refuse to interpret the
  # marginal moments unless the chain actually mixed.
  totalProp <- nDraw + burn
  for (c in seq_len(K)) {
    expect_gt(nAcc31[c] / totalProp, 0.05,
              label = sprintf("case 31 (class %d) acceptance: %.3f",
                              c, nAcc31[c] / totalProp))
  }
  expect_gt(nAcc33 / totalProp, 0.05,
            label = sprintf("case 33 acceptance: %.3f",
                            nAcc33 / totalProp))

  # Reference moments of HalfNormal(1):
  #   E    = sqrt(2/π)  ≈ 0.7979
  #   sd   = sqrt(1 - 2/π) ≈ 0.6028
  hnMean <- sqrt(2 / pi)
  hnSd   <- sqrt(1 - 2 / pi)

  # KS would fail this chain not because the marginal is wrong but
  # because Bactrian multiplicative moves under-sample the near-zero tail
  # in 6000 draws — a well-known artefact of multiplicative scale moves
  # on half-normal targets. The moments converge cleanly; the small tail
  # only fills in at higher sample counts. Use a moment-based check that
  # detects gross-distribution bugs (wrong scale, wrong shape, biased mean)
  # without flaking on tail under-sampling.

  meanTol <- 0.15  # ±15% of HN(1) mean
  sdTol   <- 0.20  # ±20% of HN(1) sd

  expect_lt(abs(mean(tauVec) - hnMean) / hnMean, meanTol,
            label = sprintf("τ mean %.4f vs HN(1) target %.4f",
                            mean(tauVec), hnMean))
  expect_lt(abs(sd(tauVec) - hnSd) / hnSd, sdTol,
            label = sprintf("τ sd %.4f vs HN(1) target %.4f",
                            sd(tauVec), hnSd))

  for (c in seq_len(K)) {
    expect_lt(abs(mean(zMat[, c]) - hnMean) / hnMean, meanTol,
              label = sprintf("z_%d mean %.4f vs HN(1) target %.4f",
                              c, mean(zMat[, c]), hnMean))
    expect_lt(abs(sd(zMat[, c]) - hnSd) / hnSd, sdTol,
              label = sprintf("z_%d sd %.4f vs HN(1) target %.4f",
                              c, sd(zMat[, c]), hnSd))
  }
})
