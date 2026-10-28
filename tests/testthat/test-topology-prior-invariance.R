# NNI (moveType 5), TBR (17) and pSPR (20) hold the prior at beta = 0 on the
# bifurcating rooted trees a production run stores: topology and edge
# fractions jointly, which the topology-only detailed-balance tests cannot
# check (#324).

test_that("nni holds the prior at beta = 0", {
  skip_under_memcheck()
  set.seed(3245L)
  p <- .PriorInvarianceP(moveType = 5L, nMove = 25L, nReps = 1000L,
                         rooted = TRUE)
  expect_true(all(p > 1e-4), info = paste(signif(p, 2), collapse = ", "))
})

test_that("tbr holds the prior at beta = 0", {
  skip_under_memcheck()
  set.seed(3247L)
  p <- .PriorInvarianceP(moveType = 17L, nMove = 25L, nReps = 1000L,
                         rooted = TRUE)
  expect_true(all(p > 1e-4), info = paste(signif(p, 2), collapse = ", "))
})

test_that("pspr holds the prior at beta = 0", {
  skip_under_memcheck()
  set.seed(3240L)
  p <- .PriorInvarianceP(moveType = 20L, nMove = 25L, nReps = 1000L,
                         rooted = TRUE)
  expect_true(all(p > 1e-4), info = paste(signif(p, 2), collapse = ", "))
})
