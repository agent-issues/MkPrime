# #343: chain state numbers tips by data row, so a resume on the same taxa in
# another row order must pair them by label, not by position.

.TipOrderMatrix <- function() {
  set.seed(343)
  mat <- matrix(sample(0:1, 12 * 40, replace = TRUE), 12, 40,
                dimnames = list(paste0("t", 1:12), NULL))
  mat[1, ] <- 0
  mat[2, ] <- 1
  mat
}

.TipOrderCheckpoint <- function(mat) {
  cpFile <- tempfile(fileext = ".rds")
  logFile <- tempfile(fileext = ".log")
  set.seed(1343)
  tree <- ape::rtree(nrow(mat), tip.label = rownames(mat))
  setTimeLimit(elapsed = 120, transient = TRUE)
  on.exit(setTimeLimit(), add = TRUE)
  allow_warning(
    RunMkPrime(MatrixToPhyDat(mat), tree,
               mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 400L, thin = 5L,
                                  maxWarmup = 100L, minWarmup = 100L,
                                  autoTune = FALSE, checkEvery = 200L,
                                  checkpointFile = cpFile, logFile = logFile)),
    "stabilisation"
  )
  list(cpFile = cpFile,
       files = c(cpFile, logFile, sub("\\.[^.]+$", "_1.log", logFile),
                 sub("\\.[^.]+$", "_trees.nwk", logFile)))
}

# The chain's own log-likelihood, recomputed with its tips labelled in the
# order the chain was started on.
.ChainLogLik <- function(chain, mkd, model) {
  tree <- structure(
    list(edge = chain$edge,
         edge.length = chain$tree_length * chain$rel_br_lengths,
         Nnode = nrow(chain$edge) - nrow(mkd$matrix) + 1L,
         tip.label = rownames(mkd$matrix)),
    class = "phylo"
  )
  list(
    tree = tree,
    logLik = MkpLogLikelihood(
      tree, mkd, kPrime = chain$kPrime, rate_loss = chain$rate_loss,
      rate_log_sd = chain$rate_log_sd, nCat = model$nCat,
      coding = model$coding, rate_neo = chain$rate_neo %||% 1,
      relabel = model$relabel
    )
  )
}

test_that("ResumeMkPrime accepts identical data in another tip order, pairing tips by label (#343)", {
  mat <- .TipOrderMatrix()
  ck <- .TipOrderCheckpoint(mat)
  on.exit(unlink(ck$files), add = TRUE)
  before <- readRDS(ck$cpFile)
  mkd <- MkPrimeData(MatrixToPhyDat(mat))
  expect_equal(
    .ChainLogLik(before$runs[[1]]$chains[[1]], mkd, before$model)$logLik,
    before$runs[[1]]$chains[[1]]$log_lik, tolerance = 1e-8
  )

  setTimeLimit(elapsed = 120, transient = TRUE)
  on.exit(setTimeLimit(), add = TRUE)
  result <- allow_warning(
    ResumeMkPrime(ck$cpFile, MatrixToPhyDat(mat[12:1, ]),
                  mcmc = list(nIter = before$iter + 200L)),
    "stabilisation"
  )
  after <- readRDS(ck$cpFile)
  expect_equal(after$iter, before$iter + 200L)
  chain <- after$runs[[1]]$chains[[1]]
  labelled <- .ChainLogLik(chain, mkd, after$model)
  expect_equal(labelled$logLik, chain$log_lik, tolerance = 1e-8)

  lastTree <- result$trees[[length(result$trees)]]
  expect_equal(TreeDist::RobinsonFoulds(lastTree, labelled$tree), 0)
})

test_that("ResumeMkPrime refuses data on different taxa (#343)", {
  mat <- .TipOrderMatrix()
  ck <- .TipOrderCheckpoint(mat)
  on.exit(unlink(ck$files), add = TRUE)
  # Sorts where "t9" did, so the data fingerprint alone would not notice.
  renamed <- mat
  rownames(renamed)[9] <- "t9x"
  expect_error(ResumeMkPrime(ck$cpFile, MatrixToPhyDat(renamed)),
               "different taxa")
})

test_that("A checkpoint without tip labels pairs tips by position (#343)", {
  mat <- .TipOrderMatrix()
  ck <- .TipOrderCheckpoint(mat)
  on.exit(unlink(ck$files), add = TRUE)
  cp <- readRDS(ck$cpFile)
  expect_identical(cp$mcmc$tipLabels, rownames(mat))
  cp$mcmc$tipLabels <- NULL
  saveRDS(cp, ck$cpFile)
  local_mkp_verbosity()
  setTimeLimit(elapsed = 120, transient = TRUE)
  on.exit(setTimeLimit(), add = TRUE)
  messages <- capture_messages(allow_warning(
    ResumeMkPrime(ck$cpFile, MatrixToPhyDat(mat),
                  mcmc = list(nIter = cp$iter + 200L)),
    "stabilisation"
  ))
  expect_match(messages, "predates stored tip labels", all = FALSE)
  expect_identical(readRDS(ck$cpFile)$mcmc$tipLabels, rownames(mat))
})

# #346: pruning does not rescale, so a near-random many-state character on
# enough tips underflows the starting state to -Inf.
test_that("RunMkPrime aborts when the starting log-likelihood is -Inf (#346)", {
  set.seed(346)
  nTip <- 300
  mat <- matrix(sample(0:3, nTip, replace = TRUE), nTip, 1,
                dimnames = list(paste0("t", seq_len(nTip)), NULL))
  tree <- ape::rtree(nTip, tip.label = rownames(mat))
  tree$edge.length <- tree$edge.length * 20 / sum(tree$edge.length)
  err <- expect_error(
    RunMkPrime(MatrixToPhyDat(mat), tree, model = MkPrimeModel(nCat = 1L),
               mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 20L, minWarmup = 10L,
                                  maxWarmup = 10L)),
    "Starting log-likelihood is -Inf"
  )
  expect_match(conditionMessage(err), "Underflowing character: 1")
})

test_that("ResumeMkPrime aborts on a checkpoint whose chain is at -Inf (#346)", {
  mat <- .TipOrderMatrix()
  ck <- .TipOrderCheckpoint(mat)
  on.exit(unlink(ck$files), add = TRUE)
  cp <- readRDS(ck$cpFile)
  cp$runs[[1]]$chains[[1]]$log_lik <- -Inf
  saveRDS(cp, ck$cpFile)
  expect_error(ResumeMkPrime(ck$cpFile, MatrixToPhyDat(mat)),
               "Starting log-likelihood is -Inf")
})
