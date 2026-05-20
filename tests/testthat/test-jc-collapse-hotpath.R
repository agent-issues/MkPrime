# Direct equivalence of the flat-CL collapsed persite kernel (used in the
# MCMC Gibbs k' sweep) with its uncollapsed counterpart. The PR that added
# the nested-vector pruning_jc_collapsed proved the math via tests/test-jc-collapse.R;
# this file proves the flat-CL port. Chain-trace equivalence is NOT a valid
# oracle because FP reordering in siteLL flips knife-edge multinomial draws
# in the Gibbs sweep — direct per-site comparison is.

test_that("flat-CL persite collapsed kernel matches uncollapsed at <1e-12", {
  library("ape")

  one_case <- function(seed, nTip, nChar, kFull, kObs, useAcrv) {
    set.seed(seed)
    tree <- ape::rtree(nTip, br = function(n) runif(n, 0.05, 0.5))
    tree$tip.label <- paste0("t", seq_len(nTip))
    parent <- as.integer(tree$edge[, 1])
    child  <- as.integer(tree$edge[, 2])
    el     <- tree$edge.length
    # Tip states drawn from [0, kObs-1]; some missing
    mat <- matrix(sample(0:(kObs - 1L), nTip * nChar, replace = TRUE),
                  nTip, nChar)
    if (length(mat) > 10) {
      miss <- sample.int(length(mat), max(1L, length(mat) %/% 10L))
      mat[miss] <- -1L
    }
    storage.mode(mat) <- "integer"
    rates <- if (useAcrv) {
      nCat <- 4L
      z <- qnorm((seq_len(nCat) - 0.5) / nCat)
      r <- exp(z * 0.5); r / mean(r)
    } else 1.0
    ll_unc <- test_persite_uncollapsed(parent, child, el, mat, kFull, rates)
    ll_col <- test_persite_collapsed(parent, child, el, mat, kFull, kObs, rates)
    max(abs(ll_unc - ll_col))
  }

  # Sweep: kObs ∈ {1..4}, lift {1, 2, 5, 10}, tip counts {5, 12, 25},
  # both ACRV settings. Mirrors test-jc-collapse.R's coverage for the
  # nested-vector kernel.
  cases <- expand.grid(
    seed   = c(1L, 7L, 19L, 42L),
    nTip   = c(5L, 12L, 25L),
    kObs   = c(1L, 2L, 3L, 4L),
    uplift = c(1L, 2L, 5L, 10L),
    useAcrv = c(TRUE, FALSE)
  )
  cases$kFull <- cases$kObs + cases$uplift
  cases <- cases[cases$kFull <= 14L, ]
  diffs <- vapply(seq_len(nrow(cases)), function(i) {
    one_case(cases$seed[i], cases$nTip[i], 20L,
             cases$kFull[i], cases$kObs[i], cases$useAcrv[i])
  }, numeric(1))
  expect_true(max(diffs) < 1e-12,
              info = sprintf("max diff = %.3e across %d cases",
                             max(diffs), length(diffs)))
})

test_that("collapsed kernel handles kObs=1 (singleton observed) cleanly", {
  set.seed(2026)
  library("ape")
  tree <- ape::rtree(8, br = function(n) runif(n, 0.05, 0.5))
  tree$tip.label <- paste0("t", seq_len(8))
  parent <- as.integer(tree$edge[, 1])
  child  <- as.integer(tree$edge[, 2])
  el     <- tree$edge.length
  # All tips in state 0 (kObs=1 in practice) plus one missing
  mat <- matrix(0L, 8, 5)
  mat[1, 1] <- -1L
  ll_unc <- test_persite_uncollapsed(parent, child, el, mat, 6L, 1.0)
  ll_col <- test_persite_collapsed(parent, child, el, mat, 6L, 1L, 1.0)
  expect_true(max(abs(ll_unc - ll_col)) < 1e-12)
})
