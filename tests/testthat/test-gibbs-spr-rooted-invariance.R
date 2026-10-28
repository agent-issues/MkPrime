# gibbs_spr (moveType 10) leaves pi invariant on rooted trees (#338).
#
# The default start trees have a degree-2 root, 2n - 2 edges.  gibbs_spr's
# eligible-prune count is then 2n - 4 rather than 2n - 6, but still the same
# for every tree, so the argument of dev/red-team/proofs/hastings-tree-moves.md
# section 5 carries over; this pins it empirically.

test_that("gibbs_spr holds the rooted prior at beta = 0", {
  skip_under_memcheck()
  set.seed(1006L)
  p <- .RootedInvarianceP(moveType = 10L, nMove = 25L, nReps = 1000L)
  expect_true(all(p > 1e-4), info = paste(signif(p, 2), collapse = ", "))
})
