# The continuous MH moves hold the prior at beta = 0 (#450): the Bactrian
# scales (moveTypes 0-3), the joint 2D scales (21, 22, 33), and the k' shifts
# (7, 26), whose proposal must be symmetric.  Dropping any of their Hastings
# terms leaves every topology and branch-length statistic unchanged.
#
# Each replicate starts from an exact prior draw, so E[x_t] = E[x_0] is the
# prior mean at every step t whatever the mixing.  The replicate statistic,
# mean_t(x_t) less the prior mean -- or less x_0 for the sum of k', which the
# k' moves shift only by a common delta -- is zero in expectation, and as
# replicates are independent its z-score needs no ESS.  With 100 replicates,
# |z| < 5 false-fails with probability 2 pt(-5, 99) = 2.4e-6 per coordinate;
# 3.4e-5 over the 14 coordinates checked.

.ScalarFixture <- function() {
  tree <- Preorder(read.tree(
    text = "(((t1:0.1,t2:0.2):0.1,(t3:0.1,t4:0.3):0.2):0.1,(t5:0.2,t6:0.1):0.1);"
  ))
  mat <- matrix(c(0, 0, 1, 1, 0, 1,
                  1, 0, 0, 1, 1, 0,
                  0, 1, 2, 0, 1, 2,
                  0, 1, 0, 1, 1, 0,
                  2, 0, 1, 1, 0, 2), 6,
                dimnames = list(paste0("t", 1:6), NULL))
  mkd <- MkPrimeData(MatrixToPhyDat(mat), neomorphic = 1:2)
  model <- .FinalizeModel(MkPrimeModel(kPrimePrior = "geometric"), tree, mkd)
  list(mkd = mkd, model = model, state = .InitState(tree, mkd, model),
       dataPtr = .InitMcmcData(mkd, model),
       trans = which(mkd$type == "transformational"))
}

# Each k'_i has a geometric prior given p, truncated to [kObs_i, K].
.KPrimeSupport <- function(fx, p) {
  lapply(fx$mkd$kObs[fx$trans], function(kObs) {
    k <- kObs:fx$model$kprimeTruncK
    list(k = k, prob = (1 - p)^(k - kObs) / sum((1 - p)^(k - kObs)))
  })
}

# Prior means of the logged scalars.
.ScalarPriorMean <- function(m) {
  c(tl = digamma(m$treeLengthShape) - log(m$treeLengthRate),
    rl = m$rateLossMeanlog,
    rls = digamma(m$rateLogSdShape) - log(m$rateLogSdRate),
    rn = m$rateNeoMeanlog)
}

.ScalarPriorDraw <- function(fx, p) {
  m <- fx$model
  st <- fx$state
  st$tree_length <- rgamma(1, m$treeLengthShape, m$treeLengthRate)
  st$rate_loss <- rlnorm(1, m$rateLossMeanlog, m$rateLossSdlog)
  st$rate_log_sd <- rgamma(1, m$rateLogSdShape, m$rateLogSdRate)
  st$rate_neo <- rlnorm(1, m$rateNeoMeanlog, m$rateNeoSdlog)
  st$p <- p
  st$kPrime[fx$trans] <- vapply(.KPrimeSupport(fx, p), function(s) {
    s$k[sample.int(length(s$k), 1L, prob = s$prob)]
  }, integer(1))
  st$log_prior <- LogPrior(st, m, fx$mkd)
  st
}

# log tree_length, rate_loss, rate_log_sd and rate_neo, and the sum of the
# transformational k', from run_mcmc_batch_cpp's scalar_samples (columns:
# log_post, log_lik, tree_length, rate_loss, rate_log_sd, p, rate_neo,
# swap_cold, topo_hash, then each k').
.ScalarTrack <- function(samples, nTrans) {
  logs <- log(samples[, c(3, 4, 5, 7), drop = FALSE])
  colnames(logs) <- c("tl", "rl", "rls", "rn")
  cbind(logs, k = rowSums(samples[, 9L + seq_len(nTrans), drop = FALSE]))
}

# `nMove` applications of `moveType` at beta = 0 through the production
# dispatcher.
.ScalarPath <- function(fx, statePtr, moveType, nMove, rho = 0) {
  run_mcmc_batch_cpp(
    fx$dataPtr, list(statePtr), 0, moveType, fx$trans - 1L,
    0L, 1, matrix(1), 10, 1L,
    0L, matrix(1), matrix(rho),
    nMove, 1L, 0L, 1L, TRUE, length(fx$state$rel_br_lengths)
  )$scalar_samples
}

# z-scores of the replicate statistic over `nReps` replicates.
.ScalarDriftZ <- function(fx, moveType, p = 0.5, rho = 0, nMove = 200L,
                          nReps = 100L) {
  priorMean <- .ScalarPriorMean(fx$model)
  drift <- t(vapply(seq_len(nReps), function(r) {
    start <- .ScalarPriorDraw(fx, p)
    statePtr <- .InitMcmcChain(start)
    fill_partition_cache(fx$dataPtr, statePtr)
    allocate_cl_workspace(fx$dataPtr, statePtr)
    path <- .ScalarPath(fx, statePtr, moveType, nMove, rho)
    colMeans(.ScalarTrack(path, length(fx$trans))) -
      c(priorMean, k = sum(start$kPrime[fx$trans]))
  }, numeric(5)))
  colMeans(drift) / apply(drift, 2, sd) * sqrt(nReps)
}

.ExpectNoDrift <- function(z, coords) {
  expect_true(all(abs(z[coords]) < 5),
              info = paste(names(z), signif(z, 2), collapse = ", "))
}

test_that("prior draws start from the prior the moves target", {
  fx <- .ScalarFixture()
  set.seed(4500)
  st <- .ScalarPriorDraw(fx, p = 0.2)
  statePtr <- .InitMcmcChain(st)
  expect_equal(st$log_prior, eval_log_prior_cpp(fx$dataPtr, statePtr))
  fill_partition_cache(fx$dataPtr, statePtr)
  allocate_cl_workspace(fx$dataPtr, statePtr)
  track <- .ScalarTrack(.ScalarPath(fx, statePtr, 0L, 3L), length(fx$trans))
  now <- get_mcmc_state(statePtr)
  expect_equal(track[3, ], c(
    tl = log(now$treeLength), rl = log(now$rateLoss),
    rls = log(now$rateLogSd), rn = log(now$rateNeo),
    k = sum(now$kPrime[fx$trans])
  ))
})

test_that("Bactrian scale moves hold the prior at beta = 0", {
  skip_under_memcheck()
  fx <- .ScalarFixture()
  set.seed(4501)
  coords <- c("tl", "rl", "rls", "rn")
  for (move in 0:3) {
    .ExpectNoDrift(.ScalarDriftZ(fx, move), coords[move + 1L])
  }
})

test_that("joint 2D moves hold the prior at beta = 0", {
  skip_under_memcheck()
  fx <- .ScalarFixture()
  set.seed(4502)
  pairs <- list(`21` = c("tl", "rls"), `22` = c("tl", "rl"),
                `33` = c("tl", "rn"))
  for (move in names(pairs)) {
    .ExpectNoDrift(.ScalarDriftZ(fx, as.integer(move)), pairs[[move]])
  }
  .ExpectNoDrift(.ScalarDriftZ(fx, 21L, rho = 0.7), pairs[["21"]])
})

test_that("k' shifts are symmetric: int_walk and block shift hold the prior", {
  skip_under_memcheck()
  fx <- .ScalarFixture()
  set.seed(4504)
  .ExpectNoDrift(.ScalarDriftZ(fx, 7L, p = 0.2), "k")
  .ExpectNoDrift(.ScalarDriftZ(fx, 26L, p = 0.05), "k")
})
