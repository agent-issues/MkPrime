# marginal_k must refuse what it does not implement (issue #26).
#
# The marginal evaluator sums over k' for transformational (type 1) partitions
# only: known-state (type 2) partitions are evaluated by neither loop, and the
# partition API never reaches it (compute_full_loglik_at returns from the
# marginal branch before testing usePartitioned, while the prior still routes
# through cpp_log_prior_partitioned). Both are deferred, not supported
# (dev/notes/2026-05-28-marginal-k-plan.md sections 11 and 13) -- so both must
# abort rather than silently produce a posterior that ignores them.

library("TreeTools", quietly = TRUE)

.mku_data <- function(knownStates = integer(0)) {
  m <- rbind(t1 = c("0", "0", "0"), t2 = c("0", "1", "1"),
             t3 = c("1", "1", "2"), t4 = c("1", "0", "2"),
             t5 = c("0", "1", "0"), t6 = c("1", "0", "1"))
  colnames(m) <- c("c1", "c2", "c3")
  MkPrimeData(MatrixToPhyDat(m), knownStates = knownStates)
}

.mku_tree <- function(mkd) {
  tr <- PectinateTree(rownames(mkd$matrix))
  tr$edge.length <- rep(0.1, nrow(tr$edge))
  tr
}

.mku_model <- function(mode) {
  suppressMessages(MkPrimeModel(
    coding = "variable", kPrimePrior = "geometric", likelihoodMode = mode,
    kprimeTruncK = 20L, expSteps = 1.4))
}

# Tiny budget: the guard must fire long before any iteration runs. The warning
# is MkPrimeMCMC clamping minWarmup to the (absurd) maxWarmup this implies.
.mku_mcmc <- function() {
  suppressWarnings(MkPrimeMCMC(nIter = 4L, thin = 1L, nChains = 1L,
                               nRuns = 1L, autoTune = FALSE))
}

test_that("marginal_k aborts on known-state characters", {
  mkd <- .mku_data(knownStates = c("1" = 3L))
  expect_true("known" %in% mkd$type)
  expect_error(
    MkPrime:::.InitMcmcData(mkd, .mku_model("marginal_k")),
    "known-state characters"
  )
  # The same data under sampled_k is fine: known-k partitions are supported there.
  expect_no_error(MkPrime:::.InitMcmcData(mkd, .mku_model("sampled_k")))
})

test_that("RunMkPrime aborts on marginal_k + knownStates", {
  mkd <- .mku_data(knownStates = c("2" = 4L))
  expect_error(
    RunMkPrime(mkd, tree = .mku_tree(mkd), model = .mku_model("marginal_k"),
               mcmc = .mku_mcmc()),
    "known-state characters"
  )
})

test_that("RunMkPrime aborts on marginal_k + the partition API", {
  mkd <- .mku_data()
  tr <- .mku_tree(mkd)
  # Trivial single-class spec: still routes the state through the partitioned
  # initializer (usePartitioned), which the marginal evaluator cannot see.
  expect_error(
    RunMkPrime(mkd, tree = tr, model = .mku_model("marginal_k"),
               mcmc = .mku_mcmc(), partition = rep(1L, mkd$nChar)),
    "partition API"
  )
  # Multi-class, with an unlink token Layer 1 supports.
  expect_error(
    RunMkPrime(mkd, tree = tr, model = .mku_model("marginal_k"),
               mcmc = .mku_mcmc(), partition = c(1L, 1L, 2L),
               unlink = "shape"),
    "partition API"
  )
})

test_that("marginal_k without knownStates or a partition is accepted", {
  mkd <- .mku_data()
  expect_false("known" %in% mkd$type)
  expect_no_error(MkPrime:::.InitMcmcData(mkd, .mku_model("marginal_k")))
  expect_no_error(MkPrime:::.RequireMarginalKSupported(
    .mku_model("marginal_k"), mkd,
    list(partition = NULL, unlink = character(0), nClasses = 1L)))
})
