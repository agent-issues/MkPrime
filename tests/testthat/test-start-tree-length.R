# Saturated JC data under four-category ACRV: 16 tips, 100 two- or
# three-state characters, tree length 30, rate_log_sd 1.5. Parsimony
# undercounts changes here, so expSteps (about 4.6) sits well below the
# posterior tree length (about 16).
SaturatedData <- function() {
  withr::local_seed(4741L)
  nTip <- 16L
  tree <- ape::rtree(nTip, tip.label = paste0("t", seq_len(nTip)))
  weight <- rexp(nrow(tree$edge))
  tree$edge.length <- 30 * weight / sum(weight)
  tree <- Preorder(tree)
  rates <- DiscreteLognormalRates(1.5, 4L)
  k <- sample(2:3, 100L, replace = TRUE)
  mat <- vapply(k, function(kj) {
    repeat {
      states <- integer(2L * nTip - 1L)
      states[nTip + 1L] <- sample.int(kj, 1L) - 1L
      rate <- rates[sample.int(4L, 1L)]
      for (e in seq_len(nrow(tree$edge))) {
        parent <- tree$edge[e, 1]
        t <- tree$edge.length[e] * rate
        pSame <- 1 / kj + (1 - 1 / kj) * exp(-kj * t / (kj - 1))
        states[tree$edge[e, 2]] <- if (runif(1L) < pSame) {
          states[parent]
        } else {
          sample(setdiff(seq_len(kj) - 1L, states[parent]), 1L)
        }
      }
      tips <- states[seq_len(nTip)]
      if (length(unique(tips)) > 1L) break
    }
    match(tips, unique(tips)) - 1L
  }, integer(nTip))
  rownames(mat) <- tree$tip.label
  MkPrimeData(MatrixToPhyDat(mat), knownStates = setNames(k, seq_along(k)))
}

test_that("edge-less start trees total three times expSteps", {
  expect_equal(.StartEdgeLengths(7L, 2.5), rep(7.5 / 7, 7))
  tree <- BalancedTree(8)
  tree$edge.length <- NULL
  mkd <- MkPrimeData(MatrixToPhyDat(matrix(
    c(0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1), 8,
    dimnames = list(tree$tip.label, NULL)
  )))
  out <- .PrepareStartTree(tree, mkd, expSteps = 2)
  expect_equal(sum(out$edge.length), 6)
  expect_true(all(out$edge.length == out$edge.length[1]))
})

test_that("default start avoids the short-tree ACRV basin (#474)", {
  skip_if_not_installed("TreeSearch")
  mkd <- SaturatedData()
  expSteps <- .FinalizeModel(MkPrimeModel(), NULL, mkd)$expSteps
  expect_gt(expSteps, 4)
  expect_lt(expSteps, 6)

  withr::local_seed(3L)
  mcmc <- MkPrimeMCMC(
    nIter = 1000L, thin = 10L, minWarmup = 500L, maxWarmup = 500L,
    autoTune = FALSE, nRuns = 1L, maxTime = 60
  )
  post <- allow_warning(RunMkPrime(mkd, mcmc = mcmc), "without stabilisation")
  # From 0.1 per edge (tree length 3), chains sat at rate_log_sd ~ 6 and
  # tree length ~ 2 for 20 000+ iterations; the posterior is near 1 and 16.
  expect_lt(median(post$samples[, "rate_log_sd"]), 3)
  expect_gt(median(post$samples[, "tree_length"]), 8)
})
