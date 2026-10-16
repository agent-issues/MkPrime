# The branch-length Hastings terms of the default-on beta_simplex (moveType 4),
# SPR (6, the TreeNav path whenever the node-CL cache is ready),
# dirichlet_branch (23) and local_dirichlet (24) hold the prior at beta = 0
# (#143, #324, #450), starting from the bifurcating rooted trees a production
# run stores.
#
# Topology frequencies alone cannot see these terms (the slow-tier
# test-topology-detailed-balance.R passes with both dropped); the edge-fraction
# marginals can: zeroing SPR's Hastings term gives KS p = 0 on move 6, and
# zeroing dirichlet_branch's gives p < 1e-7 on move 23 with a block of K = 2.
# At the production block sizes (5 and 6 here) the same mutant gives only
# p ~ 1e-4, so those runs guard the block selection rather than the term.

test_that("beta_simplex holds the prior at beta = 0", {
  skip_under_memcheck()
  set.seed(3248L)
  p <- .PriorInvarianceP(moveType = 4L, nMove = 25L, nReps = 500L,
                         rooted = TRUE)
  expect_true(all(p > 1e-4), info = paste(signif(p, 2), collapse = ", "))
})

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
  p <- c(
    .PriorInvarianceP(moveType = 23L, nMove = 25L, nReps = 1000L,
                      rooted = TRUE, intWalkWindow = 2L),
    .PriorInvarianceP(moveType = 23L, nMove = 25L, nReps = 500L,
                      rooted = TRUE, intWalkWindow = 5L)
  )
  expect_true(all(p > 1e-4), info = paste(signif(p, 2), collapse = ", "))
})

test_that("local_dirichlet holds the prior at beta = 0", {
  skip_under_memcheck()
  set.seed(3244L)
  p <- c(
    .PriorInvarianceP(moveType = 24L, nMove = 25L, nReps = 1000L,
                      rooted = TRUE, intWalkWindow = 2L),
    .PriorInvarianceP(moveType = 24L, nMove = 25L, nReps = 500L,
                      rooted = TRUE, intWalkWindow = 6L)
  )
  expect_true(all(p > 1e-4), info = paste(signif(p, 2), collapse = ", "))
})
