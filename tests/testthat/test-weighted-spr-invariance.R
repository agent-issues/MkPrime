# GSPR-003: weighted_spr (moveType 13) leaves pi invariant.
#
# The move samples a branch-fraction bin, then draws the committed fraction
# from a Beta centred on that bin's midpoint.  The bin is auxiliary and is not
# carried in the state, so the fraction's proposal density is the Beta MIXTURE
# over bins; the components overlap far too heavily for one of them to stand
# in.  Derivation: dev/red-team/proofs/weighted-spr-hastings.md.
#
# At beta = 0 the target is the prior -- Uniform over labelled unrooted
# topologies, Dirichlet(1, ..., 1) over the edge fractions -- which is exactly
# i.i.d.-sampleable.  Starting each replicate from an exact draw therefore
# tests pi K^B = pi with no burn-in or mixing assumption.

library("ape")
library("TreeTools")

test_that("weighted_spr holds the prior at beta = 0", {
  nTip <- 6L
  nEdge <- 2L * nTip - 3L
  tipLabel <- paste0("t", seq_len(nTip))
  nReps <- 1500L
  nSweeps <- 25L

  set.seed(41L)
  mat <- matrix(sample(0:1, nTip * 4L, replace = TRUE), nrow = nTip,
                dimnames = list(tipLabel, NULL))
  mkd <- suppressWarnings(MkPrimeData(MatrixToPhyDat(mat)))
  model <- MkPrimeModel()

  # An exact draw from the target: uniform topology, flat Dirichlet fractions.
  PiDraw <- function() {
    tree <- TreeTools::Preorder(
      ape::rtree(nTip, rooted = FALSE, tip.label = tipLabel))
    f <- rexp(nEdge)
    tree$edge.length <- f / sum(f)
    tree
  }
  Stats <- function(edge, f) {
    c(min(f), max(f), sum(f[edge[, 2] > nTip]))
  }

  sampled <- t(vapply(seq_len(nReps), function(r) {
    tree <- PiDraw()
    finalModel <- MkPrime:::.FinalizeModel(model, tree, mkd)
    dataPtr  <- MkPrime:::.InitMcmcData(mkd, finalModel)
    statePtr <- MkPrime:::.InitMcmcChain(
      MkPrime:::.InitState(tree, mkd, finalModel))
    fill_partition_cache(dataPtr, statePtr)
    allocate_cl_workspace(dataPtr, statePtr)
    for (i in seq_len(nSweeps))
      do_move_cpp(dataPtr, statePtr, 13L, 0L, 0.5, 0.5, 1L, 0.0)
    st <- get_mcmc_state(statePtr)
    Stats(st$edge, st$relBrLengths)
  }, numeric(3L)))

  reference <- t(vapply(seq_len(nReps * 10L), function(i) {
    tree <- PiDraw()
    Stats(tree$edge, tree$edge.length)
  }, numeric(3L)))

  # A single-component ratio drives fractions to the extremes, putting the
  # smallest 25% below its target mean at p ~ 1e-22.
  for (j in seq_len(3L)) {
    p <- suppressWarnings(ks.test(sampled[, j], reference[, j])$p.value)
    expect_gt(p, 1e-4)
  }
})


# The beta = 0 check above leaves every bin weight at 1, so it says nothing
# about how the weights enter the mixture.  Pin that separately against an
# independent R-side sum, and against the single component it replaces.
test_that("bin_mixture_log_density integrates over the bins", {
  nBins <- 10L
  conc <- 2 * nBins
  breaks <- c(0, qbeta(seq_len(nBins - 1L) / nBins, 0.25, 0.25), 1)
  mids <- 0.5 * (breaks[-1L] + breaks[-(nBins + 1L)])
  Component <- function(f, b) {
    dbeta(f, mids[b] * conc + 1, (1 - mids[b]) * conc + 1)
  }

  set.seed(20L)
  for (trial in seq_len(20L)) {
    weights <- exp(rnorm(nBins, sd = 2))
    f <- runif(1L)
    expected <- log(sum(weights * vapply(seq_len(nBins),
                                         function(b) Component(f, b),
                                         numeric(1L))))
    expect_equal(bin_mixture_log_density(weights, f, nBins), expected,
                 tolerance = 1e-12)
  }

  # Order matters: the weights are not interchangeable.
  weights <- c(5, rep(1, nBins - 1L))
  expect_false(isTRUE(all.equal(
    bin_mixture_log_density(weights, 0.9, nBins),
    bin_mixture_log_density(rev(weights), 0.9, nBins))))

  # A zero weight drops its component rather than poisoning the sum.
  weights <- c(0, rep(1, nBins - 1L))
  expect_equal(bin_mixture_log_density(weights, 0.5, nBins),
               log(sum(weights * vapply(seq_len(nBins),
                                        function(b) Component(0.5, b),
                                        numeric(1L)))),
               tolerance = 1e-12)
  expect_identical(bin_mixture_log_density(rep(0, nBins), 0.5, nBins), -Inf)
})
