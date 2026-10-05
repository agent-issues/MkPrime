# #287: a checkpoint carries per-character chain state, so a resume on data
# whose characters differ must abort before any chain state is rebuilt.

.FingerprintTree <- function() {
  ape::read.tree(
    text = "(((t1:0.1,t2:0.2):0.1,t3:0.2):0.1,(t4:0.1,(t5:0.1,t6:0.2):0.1):0.1);")
}

.FingerprintMatrix <- function() {
  matrix(c(0, 0, 1, 1, 1, 0,
           0, 1, 0, 1, 0, 1,
           0, 0, 0, 1, 1, 1),
         6, 3, dimnames = list(paste0("t", 1:6), NULL))
}

.FingerprintCheckpoint <- function() {
  tree <- .FingerprintTree()
  cpFile <- tempfile(fileext = ".rds")
  logFile <- tempfile(fileext = ".log")
  setTimeLimit(elapsed = 120, transient = TRUE)
  on.exit(setTimeLimit(), add = TRUE)
  set.seed(287)
  suppressWarnings(RunMkPrime(MatrixToPhyDat(.FingerprintMatrix()), tree,
             mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 5000L, thin = 5L,
                                maxWarmup = 200L, minWarmup = 200L,
                                autoTune = FALSE, checkEvery = 300L,
                                checkpointFile = cpFile, logFile = logFile,
                                maxTime = 0.5)))
  list(cpFile = cpFile, tree = tree,
       files = c(cpFile, logFile, sub("\\.[^.]+$", "_1.log", logFile),
                 sub("\\.[^.]+$", "_trees.nwk", logFile)))
}

test_that("ResumeMkPrime refuses permuted characters (#287)", {
  ck <- .FingerprintCheckpoint()
  on.exit(unlink(ck$files), add = TRUE)
  permuted <- MatrixToPhyDat(.FingerprintMatrix()[, c(2, 3, 1)])
  expect_error(ResumeMkPrime(ck$cpFile, permuted, ck$tree),
               "differs from the data this chain was run on")
})

test_that("ResumeMkPrime refuses a different character count (#287)", {
  ck <- .FingerprintCheckpoint()
  on.exit(unlink(ck$files), add = TRUE)
  fewer <- MatrixToPhyDat(.FingerprintMatrix()[, 1:2])
  expect_error(ResumeMkPrime(ck$cpFile, fewer, ck$tree), "differs from the data this chain was run on")
})

test_that("ResumeMkPrime refuses changed kObs (#287)", {
  ck <- .FingerprintCheckpoint()
  on.exit(unlink(ck$files), add = TRUE)
  mat <- .FingerprintMatrix()
  mat[1:2, 3] <- 2
  expect_error(ResumeMkPrime(ck$cpFile, MatrixToPhyDat(mat), ck$tree),
               "differs from the data this chain was run on")
})

test_that("ResumeMkPrime accepts identical data, whatever the tip order", {
  ck <- .FingerprintCheckpoint()
  on.exit(unlink(ck$files), add = TRUE)
  mat <- .FingerprintMatrix()
  expect_no_error(
    ResumeMkPrime(ck$cpFile, MatrixToPhyDat(mat[6:1, ]), ck$tree))
})

test_that("A checkpoint without a fingerprint still resumes (#287)", {
  ck <- .FingerprintCheckpoint()
  on.exit(unlink(ck$files), add = TRUE)
  cp <- readRDS(ck$cpFile)
  expect_false(is.null(cp$mcmc$dataFingerprint))
  cp$mcmc$dataFingerprint <- NULL
  saveRDS(cp, ck$cpFile)
  pd <- MatrixToPhyDat(.FingerprintMatrix()[, c(2, 3, 1)])
  expect_no_error(ResumeMkPrime(ck$cpFile, pd, ck$tree))
})
