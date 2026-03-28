# Detailed balance validation for TBR proposals (M-053).
#
# Under flat likelihood, a correctly implemented MH chain using TBR
# proposals should sample uniformly over unrooted tree topologies.
# We verify this with a chi-squared test for uniformity over the 15
# distinct unrooted binary topologies on 5 tips.
#
# This test is slow (~30s) and only runs with MKPRIME_SLOW_TESTS=true.

# Canonical topology key: sorted non-trivial splits, identified by tip names.
# Invariant under root placement and node renumbering.
.topo_key <- function(tr) {
  splits <- TreeTools::as.Splits(tr)
  mat <- as.logical(splits)
  tip_names <- colnames(mat)
  row_keys <- apply(mat, 1, function(r) {
    side_a <- sort(tip_names[r])
    side_b <- sort(tip_names[!r])
    if (side_a[1] < side_b[1]) paste(side_a, collapse = ",")
    else paste(side_b, collapse = ",")
  })
  paste(sort(row_keys), collapse = "|")
}


test_that("TBR satisfies detailed balance (flat-posterior topology chain)", {
  skip_slow_tests()
  library(ape)

  n_tip <- 5L
  n_expected_topos <- 15L  # (2n-5)!! = 15
  n_iter <- 1000000L
  thin <- 200L
  n_samples <- n_iter %/% thin

  set.seed(4729)
  tree <- rtree(n_tip)
  tree <- unroot(tree)
  tree <- TreeTools::Preorder(tree)
  tl <- sum(tree$edge.length)
  rel_br <- tree$edge.length / tl

  current <- tree
  current_rel <- rel_br
  keys <- character(n_samples)

  for (i in seq_len(n_iter)) {
    prop <- MkPrime:::ProposeTbr(current, tl, current_rel)
    # MH accept with flat likelihood: alpha = min(1, exp(logHR))
    if (is.finite(prop$logHastings) &&
        log(runif(1)) < prop$logHastings) {
      current <- prop$tree
      current_rel <- prop$rel_br_lengths
    }
    if (i %% thin == 0L) {
      keys[i %/% thin] <- .topo_key(current)
    }
  }

  freq <- table(keys)

  # All 15 topologies should be visited (irreducibility)
  expect_equal(length(freq), n_expected_topos)

  # Chi-squared test for uniformity (should NOT reject at alpha = 0.005)
  chi2 <- chisq.test(as.numeric(freq))
  expect_gt(chi2$p.value, 0.005)

  # No single topology should dominate (max/min ratio < 2)
  expect_lt(max(freq) / min(freq), 2)
})
