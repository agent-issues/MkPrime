# Checkpoints written before rate_neo, slice_rate_neo and joint_tl_rn were
# gated on .RateNeoLive() hold them for neomorphic-only data (#361).

.NeoOnlyCheckpoint <- function(td) {
  set.seed(1)
  nTip <- 8L
  nChar <- 10L
  mat <- matrix(sample(0:1, nTip * nChar, TRUE), nTip,
                dimnames = list(paste0("t", seq_len(nTip)), NULL))
  for (j in seq_len(nChar)) {
    if (length(unique(mat[, j])) < 2L) mat[1, j] <- 1L - mat[1, j]
  }
  pd <- MatrixToPhyDat(mat)
  ckp <- file.path(td, "neo.ckp")
  mcmc <- MkPrimeMCMC(nRuns = 1L, nChains = 1L, nIter = 100L, thin = 5L,
                      maxWarmup = 50L, minWarmup = 50L, autoTune = FALSE,
                      checkEvery = 50L, nCore = 1L, maxTime = 60,
                      checkpointFile = ckp, logFile = file.path(td, "neo.log"))
  setTimeLimit(elapsed = 120)
  on.exit(setTimeLimit(elapsed = Inf), add = TRUE)
  suppressWarnings(suppressMessages(RunMkPrime(
    pd, neomorphic = seq_len(nChar), mcmc = mcmc, verbosity = 0)))
  list(pd = pd, nChar = nChar, ckp = ckp)
}

.InjectMoves <- function(ckp, moves) {
  ck <- readRDS(ckp)
  ck$moveWeights[moves] <- 0.01
  ck$runs <- lapply(ck$runs, function(r) {
    for (counter in c("chain_accept", "chain_propose", "chain_time_ns",
                      "chain_slice_exp")) {
      if (!is.null(r[[counter]])) {
        r[[counter]] <- lapply(r[[counter]], function(x) {
          x[moves] <- 0L
          x
        })
      }
    }
    r
  })
  saveRDS(ck, ckp)
}

test_that("a neomorphic-only checkpoint holding rate_neo moves resumes (#361)", {
  skip_on_cran()
  skip_under_memcheck()
  td <- tempfile("mkp_inert_neo_")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  fx <- .NeoOnlyCheckpoint(td)
  .InjectMoves(fx$ckp, c("rate_neo", "slice_rate_neo", "joint_tl_rn"))

  setTimeLimit(elapsed = 120)
  on.exit(setTimeLimit(elapsed = Inf), add = TRUE)
  expect_message(
    suppressWarnings(res <- ResumeMkPrime(
      fx$ckp, fx$pd, neomorphic = seq_len(fx$nChar),
      mcmc = list(nIter = 150L), verbosity = 0)),
    "predates")
  expect_false(any(c("rate_neo", "slice_rate_neo", "joint_tl_rn") %in% names(readRDS(fx$ckp)$moveWeights)))
})

test_that("a checkpoint with any other extra move still aborts (#361)", {
  skip_on_cran()
  skip_under_memcheck()
  td <- tempfile("mkp_inert_neo_")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  fx <- .NeoOnlyCheckpoint(td)
  .InjectMoves(fx$ckp, c("rate_neo", "bogus_move"))

  setTimeLimit(elapsed = 120)
  on.exit(setTimeLimit(elapsed = Inf), add = TRUE)
  expect_error(
    suppressWarnings(suppressMessages(ResumeMkPrime(
      fx$ckp, fx$pd, neomorphic = seq_len(fx$nChar),
      mcmc = list(nIter = 150L), verbosity = 0))),
    "Resumed move set differs")
})
