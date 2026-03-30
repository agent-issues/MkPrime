test_that("RunMkPrime forwards ... to MkPrimeMCMC", {
  pd   <- .mkp_test_pd()
  tree <- .mkp_test_tree()
  res <- RunMkPrime(pd, tree,
                    nIter = 500L, maxWarmup = 100L, nRuns = 1L)
  expect_s3_class(res, "MkPosterior")
  expect_gt(nrow(res$samples), 0L)
})

test_that("RunMkPrime errors when both mcmc and ... are supplied", {
  pd   <- .mkp_test_pd()
  tree <- .mkp_test_tree()
  expect_error(
    RunMkPrime(pd, tree, mcmc = MkPrimeMCMC(), nIter = 500L),
    "not both"
  )
})

test_that("RunMkPrime with no mcmc and no dots uses defaults", {
  pd   <- .mkp_test_pd()
  tree <- .mkp_test_tree()
  res <- RunMkPrime(pd, tree,
                    nIter = 300L, maxWarmup = 100L, nRuns = 1L)
  expect_s3_class(res, "MkPosterior")
})

test_that("RunMkPrime dots pass logFile correctly", {
  pd   <- .mkp_test_pd()
  tree <- .mkp_test_tree()
  tmpLog <- tempfile(fileext = ".log")
  on.exit(unlink(tmpLog), add = TRUE)
  res <- RunMkPrime(pd, tree,
                    nIter = 500L, maxWarmup = 100L, nRuns = 1L,
                    logFile = tmpLog)
  expect_true(file.exists(tmpLog))
  expect_equal(res$logFile, tmpLog)
})
