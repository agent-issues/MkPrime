# Regression test: gibbs_subtree_swap must be reproducible across RunMkPrime
# calls in the same R session under the same set.seed().
#
# Bug: two static-local DIAG counters (diagFileCount, preCheckCount) in
# do_move_impl() persisted across RunMkPrime calls.  preCheckCount triggered
# a full cpp_log_likelihood() call that inflated the measured wall-clock time
# for whichever move type was selected at that iteration; because the counter
# carried over from run 1 the inflation landed on a different move type in
# run 2, shifting the timing-based move-weight adaptation and ultimately
# diverging the RNG stream.  Fix: remove both DIAG static blocks.
skip_slow_tests()

library("ape")
library("TreeTools")

# ---------------------------------------------------------------------------
# Shared fixture: 6-tip transformational data + starting tree
# ---------------------------------------------------------------------------

.det_tree <- function() {
  Preorder(read.tree(text =
    "((t1:0.1,t2:0.2):0.15,(t3:0.1,(t4:0.1,t5:0.2):0.12):0.18,t6:0.3);"))
}

.det_mat <- function() {
  matrix(
    c(0, 1, 0, 1, 1, 0,
      0, 0, 1, 1, 0, 1,
      0, 1, 2, 0, 1, 0),
    nrow = 6, ncol = 3,
    dimnames = list(paste0("t", 1:6), NULL)
  )
}

.det_pd <- function() TreeTools::MatrixToPhyDat(.det_mat())

.det_mcmc <- function() {
  MkPrimeMCMC(
    nRuns         = 1L,
    nIter         = 600L,
    thin          = 5L,
    maxWarmup     = 300L,
    minWarmup     = 300L,
    autoTune      = FALSE,
    gibbsSubtreeSwap = TRUE
  )
}


# ---------------------------------------------------------------------------
# Test: identical samples across two RunMkPrime calls (same seed, same session)
# ---------------------------------------------------------------------------

test_that("gibbsSubtreeSwap=TRUE produces identical samples across calls", {
  tree <- .det_tree()
  pd   <- .det_pd()
  mcmc <- .det_mcmc()

  set.seed(424242L)
  r1 <- RunMkPrime(pd, tree, mcmc = mcmc)

  set.seed(424242L)
  r2 <- RunMkPrime(pd, tree, mcmc = mcmc)

  expect_identical(r1$samples, r2$samples,
    info = paste(
      "gibbsSubtreeSwap non-determinism: samples differ between two RunMkPrime",
      "calls with the same set.seed() in the same R session.",
      "Root cause: static-local DIAG counters in do_move_impl() persisted",
      "across calls and altered timing-based move-weight adaptation."
    )
  )
})


# ---------------------------------------------------------------------------
# Test: gibbsSubtreeSwap=FALSE baseline (must also be deterministic)
# ---------------------------------------------------------------------------

test_that("gibbsSubtreeSwap=FALSE also produces identical samples across calls", {
  tree <- .det_tree()
  pd   <- .det_pd()
  mcmc <- MkPrimeMCMC(
    nRuns         = 1L,
    nIter         = 600L,
    thin          = 5L,
    maxWarmup     = 300L,
    minWarmup     = 300L,
    autoTune      = FALSE,
    gibbsSubtreeSwap = FALSE
  )

  set.seed(424242L)
  r1 <- RunMkPrime(pd, tree, mcmc = mcmc)

  set.seed(424242L)
  r2 <- RunMkPrime(pd, tree, mcmc = mcmc)

  expect_identical(r1$samples, r2$samples,
    info = "gibbsSubtreeSwap=FALSE baseline: samples should be deterministic"
  )
})
