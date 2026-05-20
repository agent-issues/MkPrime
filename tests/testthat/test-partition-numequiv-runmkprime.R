# Plan v4 §7b integration test: RunMkPrime accepts partition = rep(1L, nChar)
# with unlink = character(0) (trivial spec) and routes the MCMC loop through
# cpp_log_likelihood_partitioned. The chain must run without error and produce
# all-finite log-likelihood values.
#
# The §7b numeric-equivalence guarantee (initial-state LL matches legacy to
# ~1e-10) is verified at the function level by test-partition-numequiv-singleclass.R
# and at the InitState level by test-partition-initstate.R. At the full-run
# level, the sample matrix is NOT required to match bit-for-bit — RNG ordering
# can differ once the partitioned fill_partition_cache recomputes logLik with
# the partitioned function (which changes MH acceptance on the first step).
#
# This test checks:
#   (a) The gate opens for trivial partition.
#   (b) The chain runs to completion and produces all-finite LLs.
#   (c) The run-start logLik from .InitStatePartitioned equals .InitState to ~1e-10.

library("TreeTools")


# Build a small mkd (8 chars, 6 tips, no neomorphic for Casali parity).
.setup_runmkprime_data <- function(seed = 42L, nChar = 8L, nTip = 6L) {
  set.seed(seed)
  mat <- matrix(sample(0:1, nTip * nChar, replace = TRUE),
                nrow = nTip, ncol = nChar,
                dimnames = list(paste0("t", seq_len(nTip)), NULL))
  # Ensure every character is variable (avoids invariant-char drop confounds)
  for (j in seq_len(nChar)) {
    if (length(unique(mat[, j])) < 2L) mat[1L, j] <- 1L - mat[1L, j]
  }
  pd  <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)  # hasNeo = FALSE by construction
  tree <- TreeTools::Preorder(
    TreeTools::NJTree(pd, edgeLengths = TRUE)
  )
  if (is.null(tree$edge.length) || any(tree$edge.length <= 0)) {
    tree$edge.length <- pmax(tree$edge.length %||% rep(0.1, nrow(tree$edge)), 1e-8)
  }
  list(mkd = mkd, tree = tree, pd = pd)
}


# (a) Gate opens: RunMkPrime does not error for trivial partition spec.
test_that("trivial partition gate opens (no error)", {
  d <- .setup_runmkprime_data()
  expect_no_error(
    RunMkPrime(
      data       = d$mkd,
      tree       = d$tree,
      mcmc       = MkPrimeMCMC(
        nIter    = 10L,
        maxWarmup = 5L,
        minWarmup = 5L,
        nChains  = 1L,
        thin     = 1L
      ),
      partition  = rep(1L, d$mkd$nChar),
      unlink     = character(0)
    )
  )
})


# (b) Chain produces all-finite LL values and a non-zero sample count.
test_that("trivial partition chain produces all-finite LL", {
  d <- .setup_runmkprime_data()
  set.seed(17L)
  result <- RunMkPrime(
    data       = d$mkd,
    tree       = d$tree,
    mcmc       = MkPrimeMCMC(
      nIter    = 30L,
      maxWarmup = 15L,
      minWarmup = 15L,
      nChains  = 1L,
      thin     = 1L
    ),
    partition  = rep(1L, d$mkd$nChar),
    unlink     = character(0)
  )
  ll <- result$samples[, "log_likelihood"]
  expect_true(all(is.finite(ll)))
  expect_gt(result$nSamples, 0L)
})


# (c) Initial-state LL from the partitioned initializer matches legacy to ~1e-10.
# This is the §7b contract at the state-init level, verified here end-to-end
# by calling .InitStatePartitioned and .InitState on the same inputs.
test_that("§7b: trivial partition initial LL matches legacy .InitState to ~1e-10", {
  d    <- .setup_runmkprime_data()
  model <- MkPrimeModel()
  model <- MkPrime:::.FinalizeModel(model, d$tree, d$mkd)

  spec <- .ValidatePartitionArgs(partition = rep(1L, d$mkd$nChar),
                                  unlink = character(0), d$mkd)
  # Rebuild partitions with classIdx for the partitioned path
  d$mkd$partitions <- .BuildPartitions(d$mkd, partition = spec$partition)

  s_legacy <- .InitState(d$tree, d$mkd, model)
  s_part   <- .InitStatePartitioned(d$tree, d$mkd, model, spec)

  expect_equal(s_part$log_lik, s_legacy$log_lik, tolerance = 1e-10)
})
