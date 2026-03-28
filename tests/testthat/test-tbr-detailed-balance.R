# Detailed balance validation for topology proposals.
#
# Under flat likelihood, a correct MH chain using any topology proposal
# should sample uniformly over unrooted tree topologies.  We verify this
# with a chi-squared test for uniformity over the 15 distinct unrooted
# binary topologies on 5 tips.
#
# These tests are slow (~30s each) and only run with MKPRIME_SLOW_TESTS=true.
# Intended for CI / Hamilton, not interactive development.

# --- Shared helpers --------------------------------------------------------

# Canonical topology key: sorted non-trivial splits identified by tip names.
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

# Run a flat-posterior MH chain using a single proposal function.
# Returns a frequency table of topology keys.
.flat_posterior_chain <- function(propose_fn, n_iter, thin, seed) {
  n_tip <- 5L
  n_samples <- n_iter %/% thin

  set.seed(seed)
  tree <- ape::rtree(n_tip)
  tree <- ape::unroot(tree)
  tree <- TreeTools::Preorder(tree)
  tl <- sum(tree$edge.length)

  current <- tree
  current_rel <- tree$edge.length / tl
  keys <- character(n_samples)

  for (i in seq_len(n_iter)) {
    prop <- propose_fn(current, tl, current_rel)
    if (is.finite(prop$logHastings) &&
        log(runif(1)) < prop$logHastings) {
      current <- prop$tree
      current_rel <- prop$rel_br_lengths
    }
    if (i %% thin == 0L) {
      keys[i %/% thin] <- .topo_key(current)
    }
  }

  table(keys)
}

# Assert that a topology frequency table is consistent with uniformity
# over the 15 unrooted binary topologies on 5 tips.
.expect_uniform <- function(freq) {
  expect_equal(length(freq), 15L)
  chi2 <- chisq.test(as.numeric(freq))
  expect_gt(chi2$p.value, 0.005)
  expect_lt(max(freq) / min(freq), 2)
}


# --- Tests -----------------------------------------------------------------

test_that("TBR satisfies detailed balance (flat-posterior topology chain)", {
  skip_slow_tests()
  freq <- .flat_posterior_chain(MkPrime:::ProposeTbr,
                                n_iter = 1000000L, thin = 200L,
                                seed = 4729)
  .expect_uniform(freq)
})


test_that("SPR satisfies detailed balance (flat-posterior topology chain)", {
  skip_slow_tests()
  freq <- .flat_posterior_chain(MkPrime:::ProposeSpr,
                                n_iter = 1000000L, thin = 200L,
                                seed = 7831)
  .expect_uniform(freq)
})
