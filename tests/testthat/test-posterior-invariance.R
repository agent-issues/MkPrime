# beta = 1 invariance gates for the weighted moves and the Q-heterogeneity
# paths of the Gibbs topology moves (#481).  The beta = 0 gates leave every
# weight at 1, so a kernel that scores the reverse move with the forward
# weights, or drops the weights from its bin mixture, passes them; these do
# not.  Harness: helper-posterior-invariance.R.
#
# Each test makes 7 or 15 comparisons.  Batch-means z with 50 batches is
# roughly t on 49 df, for which P(|z| > 4.5) ~ 4e-5; over at most 15 cells a
# correct kernel fails less than once in a thousand runs, whatever the seed.

.qHetModel <- function() {
  MkPrimeModel(qHeterogeneity = TRUE, nBetaCat = 4L, nCat = 1L,
               coding = "none")
}

test_that("weighted_spr samples its topology posterior at beta = 1", {
  skip_slow_tests()
  set.seed(1301L)
  z <- .TopologyZ(.PosteriorContext(), moveType = 13L, nIter = 20000L)
  expect_lt(max(abs(z)), 4.5)
})

test_that("weighted_branch_scale samples its edge posterior at beta = 1", {
  skip_slow_tests()
  set.seed(1201L)
  z <- .EdgeFractionZ(.PosteriorContext(), moveType = 12L, nIter = 20000L)
  expect_lt(max(abs(z)), 4.5)
})

test_that("block_gibbs_branch samples its edge posterior at beta = 1", {
  skip_slow_tests()
  set.seed(1501L)
  z <- .EdgeFractionZ(.PosteriorContext(), moveType = 15L, nIter = 5000L)
  expect_lt(max(abs(z)), 4.5)
})

# gibbs_spr_enumerate_cpp refuses qHet, so the deterministic forward/reverse
# weight check in test-gibbs-spr-candidates.R never reaches
# gibbs_spr_impl_het.  Three-state characters exercise its rotation loop.
test_that("gibbs_spr under Q-heterogeneity samples its posterior at beta = 1", {
  skip_slow_tests()
  set.seed(1001L)
  ctx <- .PosteriorContext(.qHetModel(), states = 0:2)
  z <- .TopologyZ(ctx, moveType = 10L, nIter = 20000L)
  expect_lt(max(abs(z)), 4.5)
})

test_that("weighted_subtree_swap under Q-heterogeneity samples at beta = 1", {
  skip_slow_tests()
  set.seed(1401L)
  ctx <- .PosteriorContext(.qHetModel(), states = 0:2)
  z <- .TopologyZ(ctx, moveType = 14L, nIter = 20000L)
  expect_lt(max(abs(z)), 4.5)
})
