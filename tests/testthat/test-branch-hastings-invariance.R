# The branch-length Hastings terms of the default-on SPR (moveType 6, the
# TreeNav path whenever the node-CL cache is ready), dirichlet_branch (23) and
# local_dirichlet (24) hold the prior at beta = 0 (#143, #324), starting from
# the bifurcating rooted trees a production run stores.
#
# Topology frequencies alone cannot see these terms (the slow-tier
# test-topology-detailed-balance.R passes with both dropped); the edge-fraction
# marginals can: zeroing SPR's Hastings term gives KS p = 0 on move 6, and
# zeroing dirichlet_branch's gives p < 1e-7 on move 23.

test_that("spr holds the prior at beta = 0", {
  skip_under_memcheck()
  set.seed(3246L)
  p <- .PriorInvarianceP(moveType = 6L, nMove = 25L, nReps = 1000L,
                         rooted = TRUE)
  expect_true(all(p > 1e-4), info = paste(signif(p, 2), collapse = ", "))
})

test_that("dirichlet_branch holds the prior at beta = 0", {
  skip_under_memcheck()
  set.seed(3243L)
  p <- .PriorInvarianceP(moveType = 23L, nMove = 25L, nReps = 1000L,
                         rooted = TRUE)
  expect_true(all(p > 1e-4), info = paste(signif(p, 2), collapse = ", "))
})

test_that("local_dirichlet holds the prior at beta = 0", {
  skip_under_memcheck()
  set.seed(3244L)
  p <- .PriorInvarianceP(moveType = 24L, nMove = 25L, nReps = 1000L,
                         rooted = TRUE)
  expect_true(all(p > 1e-4), info = paste(signif(p, 2), collapse = ", "))
})
