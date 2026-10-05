# #378: fixTopology = TRUE fixes every run to the supplied topology.
# #379: a resume aborts on a model that differs structurally from the
# checkpoint's, before the data are recoded or any chain state is rebuilt.

.ModelCheckTree <- function() {
  ape::read.tree(text = paste0("(((t1:0.1,t2:0.1):0.1,t3:0.1):0.1,",
                               "((t4:0.1,t5:0.1):0.1,(t6:0.1,t7:0.1):0.1):0.1);"))
}

# Three of the seven characters are parsimony-uninformative.
.ModelCheckData <- function() {
  mat <- cbind(c(0, 0, 0, 1, 1, 1, 1), c(0, 0, 1, 1, 2, 2, 2),
               c(0, 1, 0, 1, 0, 1, 1), c(0, 0, 0, 0, 0, 0, 1),
               c(1, 0, 0, 0, 0, 0, 0), c(0, 0, 1, 1, 1, 0, 0),
               c(0, 1, 1, 1, 1, 1, 1))
  rownames(mat) <- paste0("t", 1:7)
  MatrixToPhyDat(mat)
}

.ModelCheckMcmc <- function(cpFile, logFile) {
  MkPrimeMCMC(nIter = 200L, nRuns = 1L, checkEvery = 100L, autoTune = FALSE,
              maxWarmup = 20L, minWarmup = 0L, nCore = 1L, maxTime = 0.5,
              checkpointFile = cpFile, logFile = logFile)
}

.ModelCheckpoint <- function(model) {
  cpFile <- tempfile(fileext = ".ckp")
  logFile <- tempfile(fileext = ".log")
  setTimeLimit(elapsed = 120, transient = TRUE)
  on.exit(setTimeLimit(), add = TRUE)
  set.seed(379)
  suppressWarnings(suppressMessages(
    RunMkPrime(.ModelCheckData(), .ModelCheckTree(), model = model,
               mcmc = .ModelCheckMcmc(cpFile, logFile), verbosity = 0L)))
  list(cpFile = cpFile, logFile = logFile,
       files = c(cpFile, logFile, sub("\\.[^.]+$", "_1.log", logFile),
                 sub("\\.[^.]+$", "_trees.nwk", logFile)))
}

test_that("fixTopology = TRUE holds every run on the supplied tree (#378)", {
  tree <- .ModelCheckTree()
  setTimeLimit(elapsed = 120, transient = TRUE)
  on.exit(setTimeLimit(), add = TRUE)
  set.seed(378)
  result <- suppressWarnings(suppressMessages(
    RunMkPrime(.ModelCheckData(), tree, fixTopology = TRUE, verbosity = 0L,
               mcmc = MkPrimeMCMC(nIter = 200L, nRuns = 2L, autoTune = FALSE,
                                  maxWarmup = 20L, minWarmup = 0L,
                                  nCore = 1L, maxTime = 30))))
  expect_gt(length(result$trees), 0L)
  rf <- vapply(result$trees, function(tr) {
    TreeDist::RobinsonFoulds(tr, tree)
  }, numeric(1))
  expect_true(all(rf == 0))
})

test_that(".PerturbStart(topology = FALSE) keeps the topology (#378)", {
  tree <- .ModelCheckTree()
  set.seed(378)
  perturbed <- .PerturbStart(tree, topology = FALSE)
  expect_equal(TreeDist::RobinsonFoulds(perturbed, tree), 0)
  expect_false(isTRUE(all.equal(perturbed$edge.length, tree$edge.length)))
})

test_that("Resume aborts on informative -> variable coding (#379)", {
  ck <- .ModelCheckpoint(MkPrimeModel(coding = "informative"))
  on.exit(unlink(ck$files), add = TRUE)
  expect_error(
    ResumeMkPrime(ck$cpFile, .ModelCheckData(), model = MkPrimeModel(),
                  mcmc = list(nIter = 400L), verbosity = 0L),
    "coding.*informative.*checkpoint"
  )
})

test_that("Resume aborts on none -> informative coding (#379)", {
  ck <- .ModelCheckpoint(MkPrimeModel(coding = "none"))
  on.exit(unlink(ck$files), add = TRUE)
  expect_error(
    ResumeMkPrime(ck$cpFile, .ModelCheckData(),
                  model = MkPrimeModel(coding = "informative"),
                  mcmc = list(nIter = 400L), verbosity = 0L),
    "coding.*informative.*none.*checkpoint"
  )
})

test_that("RunMkPrime auto-resume checks the model too (#379)", {
  ck <- .ModelCheckpoint(MkPrimeModel(coding = "informative"))
  on.exit(unlink(ck$files), add = TRUE)
  expect_error(
    RunMkPrime(.ModelCheckData(), .ModelCheckTree(),
               model = MkPrimeModel(kPrimePrior = "geometric",
                                    coding = "informative"),
               mcmc = .ModelCheckMcmc(ck$cpFile, ck$logFile), verbosity = 0L),
    "kPrimePrior"
  )
})

test_that("A checkpoint that predates a model field still resumes (#379)", {
  ck <- .ModelCheckpoint(MkPrimeModel(coding = "informative"))
  on.exit(unlink(ck$files), add = TRUE)
  cp <- readRDS(ck$cpFile)
  cp$model$coding <- NULL
  saveRDS(cp, ck$cpFile)
  expect_no_error(suppressWarnings(suppressMessages(
    ResumeMkPrime(ck$cpFile, .ModelCheckData(),
                  model = MkPrimeModel(coding = "informative"),
                  mcmc = list(nIter = 400L), verbosity = 0L))))
})

test_that("An integer treeLengthShape keeps the checkpoint's rate", {
  model <- MkPrimeModel(coding = "informative", treeLengthRate = 0.37)
  ck <- .ModelCheckpoint(model)
  on.exit(unlink(ck$files), add = TRUE)
  res <- suppressWarnings(suppressMessages(
    ResumeMkPrime(ck$cpFile, .ModelCheckData(),
                  model = MkPrimeModel(coding = "informative",
                                       treeLengthShape = 2L),
                  mcmc = list(nIter = 400L), verbosity = 0L)))
  expect_equal(res$model$treeLengthRate, 0.37)
})

test_that("Resume aborts on changed neomorphic characters", {
  ck <- .ModelCheckpoint(MkPrimeModel())
  on.exit(unlink(ck$files), add = TRUE)
  expect_error(
    ResumeMkPrime(ck$cpFile, .ModelCheckData(), neomorphic = 1L,
                  mcmc = list(nIter = 400L), verbosity = 0L),
    "character types"
  )
})
