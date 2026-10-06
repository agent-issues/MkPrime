# The Gibbs k' sweep (moveType 25) must draw each k'_i from its full
# conditional.  With the tree and every other parameter held fixed, successive
# sweeps are independent exact draws from it, so a chi-squared test against
# the enumerated conditional is valid with no burn-in or mixing assumption.
#
# The conditional is enumerated over the candidates the sweep itself considers
# (kprime_sweep_candidates), which stop where the weight falls below a
# negligible fraction of the maximum: far short of any fixed cap on the k'
# range, so the test does not depend on where that cap sits.
#
# Sampling from the square root of the conditional (a flattened draw) gives
# chi-squared ~ 14000 on ~60 df.

test_that("gibbs_kprime_sweep draws from the k' full conditional", {
  set.seed(153L)
  nTip <- 8L
  tree <- Preorder(ape::rtree(nTip, rooted = FALSE))
  mat <- matrix(sample(0:2, nTip * 12L, replace = TRUE), nrow = nTip,
                dimnames = list(tree$tip.label, NULL))
  mkd <- suppressWarnings(MkPrimeData(MatrixToPhyDat(mat)))
  chain <- .PartitionedChain(tree, mkd)
  st <- get_mcmc_state(chain$statePtr)
  edgeLen <- st$treeLength * st$relBrLengths
  kObs <- mkd$kObs
  nCand <- kprime_sweep_candidates(chain$dataPtr, chain$statePtr)

  # The k' prior does not involve the class parameters, so a legacy state
  # supplies it.
  model <- MkPrime:::.FinalizeModel(MkPrimeModel(), tree, mkd)
  legacyData <- MkPrime:::.InitMcmcData(mkd, model)
  legacyState <- MkPrime:::.InitState(tree, mkd, model)
  Conditional <- function(i) {
    lw <- vapply(kObs[[i]] + seq_len(nCand[[i]]) - 1L, function(k) {
      kPrime <- st$kPrime
      kPrime[[i]] <- k
      priorState <- legacyState
      priorState$kPrime[[i]] <- k
      cpp_log_likelihood_partitioned_xptr(
        chain$dataPtr, st$edge[, 1], st$edge[, 2], edgeLen, kPrime,
        st$rateLoss, st$classRateLogSd, st$classRate, st$etaNeo,
        st$betaScale) +
        eval_log_prior_cpp(legacyData, MkPrime:::.InitMcmcChain(priorState))
    }, double(1))
    w <- exp(lw - max(lw))
    # Return:
    w / sum(w)
  }

  nDraw <- 1000L
  draws <- vapply(seq_len(nDraw), function(j) {
    do_move_cpp(chain$dataPtr, chain$statePtr, 25L, 0L, 0.5, 0.5, 3L, 1.0)
    get_mcmc_state(chain$statePtr)$kPrime
  }, integer(mkd$nChar))

  chiSq <- 0
  df <- 0
  for (i in seq_len(mkd$nChar)) {
    expect_true(all(draws[i, ] >= kObs[[i]] &
                      draws[i, ] < kObs[[i]] + nCand[[i]]))
    observed <- tabulate(draws[i, ] - kObs[[i]] + 1L, nCand[[i]])
    expected <- Conditional(i) * nDraw
    # Pool the sparse tail into one bin of expectation >= 5.
    tailSum <- rev(cumsum(rev(expected)))
    m <- max(which(tailSum >= 5))
    obsBins <- c(observed[seq_len(m - 1L)], sum(observed[m:nCand[[i]]]))
    expBins <- c(expected[seq_len(m - 1L)], tailSum[[m]])
    chiSq <- chiSq + sum((obsBins - expBins)^2 / expBins)
    df <- df + length(expBins) - 1L
  }
  expect_gt(pchisq(chiSq, df, lower.tail = FALSE), 1e-4,
            label = sprintf("chi-squared %.1f on %d df", chiSq, df))
})
