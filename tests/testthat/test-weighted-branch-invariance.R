# WBS-001: weighted_branch_scale (moveType 12) and block_gibbs_branch
# (moveType 15) leave pi invariant.
#
# Both draw a pair's new branch fraction by selecting a bin, then sampling a
# Beta centred on its midpoint.  The bin is not carried in the state, so the
# fraction's proposal density is the Beta mixture over bins; the bin weights
# depend only on edges the move leaves alone, so the mixture is the same in
# both directions.  Derivation: dev/red-team/proofs/weighted-branch-hastings.md.
#
# A single-component ratio drives fractions to the extremes: over 25 moves the
# smallest edge fraction lands 38% below its target mean, KS p ~ 1e-49.

test_that("weighted_branch_scale holds the prior at beta = 0", {
  skip_under_memcheck()
  set.seed(1204L)
  p <- .PriorInvarianceP(moveType = 12L, nMove = 25L, nReps = 1000L)
  expect_true(all(p > 1e-4), info = paste(signif(p, 2), collapse = ", "))
})

test_that("block_gibbs_branch holds the prior at beta = 0", {
  skip_under_memcheck()
  set.seed(1504L)
  p <- .PriorInvarianceP(moveType = 15L, nMove = 5L, nReps = 1000L)
  expect_true(all(p > 1e-4), info = paste(signif(p, 2), collapse = ", "))
})
