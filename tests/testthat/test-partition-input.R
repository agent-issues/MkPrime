# Partition input handling under coding = "informative", and ignored
# partition/unlink requests (#360).

.InformativePartitionFixture <- function() {
  set.seed(3)
  nTip <- 8L
  mat <- matrix(sample(0:2, nTip * 12L, TRUE), nTip,
                dimnames = list(paste0("t", seq_len(nTip)), NULL))
  mat[, 11:12] <- c(1L, rep(0L, nTip - 1L))
  list(pd = MatrixToPhyDat(mat), model = MkPrimeModel(coding = "informative"))
}

.QuickMcmc <- function(ckp = NULL, nIter = 50L, logFile = tempfile()) {
  MkPrimeMCMC(nRuns = 1L, nChains = 1L, nIter = nIter, thin = 5L,
              maxWarmup = 20L, minWarmup = 20L, autoTune = FALSE,
              checkEvery = 20L, nCore = 1L, maxTime = 60,
              checkpointFile = ckp, logFile = logFile)
}

test_that("informative coding checks partition length against the data (#360)", {
  skip_on_cran()
  fx <- .InformativePartitionFixture()
  setTimeLimit(elapsed = 120, transient = TRUE)
  expect_error(
    suppressWarnings(suppressMessages(RunMkPrime(
      fx$pd, model = fx$model, partition = rep(1:2, 5), mcmc = .QuickMcmc(),
      overwrite = TRUE, verbosity = 0))),
    "length 10.*12 character")
})

test_that("informative coding aborts when a class would be emptied (#360)", {
  skip_on_cran()
  fx <- .InformativePartitionFixture()
  setTimeLimit(elapsed = 120, transient = TRUE)
  expect_error(
    suppressWarnings(suppressMessages(RunMkPrime(
      fx$pd, model = fx$model, partition = c(rep(1:2, 5), 3L, 3L),
      mcmc = .QuickMcmc(), overwrite = TRUE, verbosity = 0))),
    "class 3 empty")
})

test_that("informative coding accepts a correct partition (#360)", {
  skip_on_cran()
  skip_under_memcheck()
  fx <- .InformativePartitionFixture()
  setTimeLimit(elapsed = 120, transient = TRUE)
  expect_no_error(suppressWarnings(suppressMessages(RunMkPrime(
    fx$pd, model = fx$model, partition = c(rep(1:2, 5), 1L, 2L),
    mcmc = .QuickMcmc(), overwrite = TRUE, verbosity = 0))))
})

test_that(".BuildPartitions refuses a partition of the wrong length (#360)", {
  mkd <- MkPrimeData(.InformativePartitionFixture()$pd)
  expect_error(.BuildPartitions(mkd, rep(1:2, 2)), "has length 4")
})

test_that("an ignored unlink warns (#360)", {
  mkd <- MkPrimeData(.InformativePartitionFixture()$pd)
  expect_warning(.ValidatePartitionArgs(NULL, "shape", mkd), "ignored")
  expect_warning(
    .ValidatePartitionArgs(rep(1L, mkd$nChar), "shape", mkd), "ignored")
})

test_that("auto-resume warns only when a request differs from the checkpoint (#360)", {
  skip_on_cran()
  skip_under_memcheck()
  fx <- .InformativePartitionFixture()
  td <- tempfile("mkp_ignored_")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  ckp <- file.path(td, "a.ckp")
  lg <- file.path(td, "a.log")
  part <- c(rep(1:2, 5), 1L, 2L)
  setTimeLimit(elapsed = 240, transient = TRUE)
  suppressWarnings(suppressMessages(RunMkPrime(
    fx$pd, model = fx$model, partition = part,
    mcmc = .QuickMcmc(ckp, logFile = lg), verbosity = 0)))

  Resume <- function(...) {
    suppressMessages(RunMkPrime(
      fx$pd, model = fx$model, mcmc = .QuickMcmc(ckp, 60L, lg), verbosity = 0,
      ...))
  }
  expect_warning(Resume(partition = rev(part)), "partition")
  expect_warning(Resume(partition = part, fixTopology = TRUE), "fixTopology")
  expect_no_warning(Resume(partition = part))
})
