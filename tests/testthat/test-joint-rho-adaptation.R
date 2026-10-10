# Joint-move correlations are estimated at default warmup lengths (#301).

test_that(".AccumulateRhoSnapshot records rate_neo", {
  paramNames <- c("tree_length", "rate_loss", "rate_log_sd", "rate_neo")
  buf <- NULL
  set.seed(3011)
  for (i in 1:60) {
    state <- list(treeLength = rlnorm(1), rateLogSd = 0.8, rateLoss = 0.3,
                  rateNeo = rlnorm(1))
    buf <- MkPrime:::.AccumulateRhoSnapshot(buf, state, hasNeo = TRUE,
                                            paramNames)
  }
  expect_gt(min(buf[, "rate_neo"]), 0)
  expect_equal(unname(buf[60, "rate_neo"]), state$rateNeo)
  expect_true(is.finite(cor(log(buf[, "tree_length"]),
                            log(buf[, "rate_neo"]))))
  expect_equal(MkPrime:::.EstimateJointRhos(buf, hasNeo = TRUE)$rho_tl_rn,
               cor(log(buf[, "tree_length"]), log(buf[, "rate_neo"])))
})

test_that(".EstimateJointRhos keeps rhos it cannot re-estimate", {
  prior <- list(rho_tl_rls = 0.5, rho_tl_rl = -0.3, rho_tl_rn = 0.2)
  short <- cbind(tree_length = rlnorm(10), rate_log_sd = rlnorm(10))
  expect_equal(MkPrime:::.EstimateJointRhos(short, TRUE, prior = prior), prior)
  expect_equal(MkPrime:::.EstimateJointRhos(NULL, TRUE, prior = prior), prior)

  # Only the pair the samples can estimate moves.
  set.seed(3012)
  tl <- rlnorm(100)
  long <- cbind(tree_length = tl,
                rate_log_sd = exp(0.7 * log(tl) + rnorm(100, 0, 0.2)))
  rhos <- MkPrime:::.EstimateJointRhos(long, TRUE, prior = prior)
  expect_gt(rhos$rho_tl_rls, 0.5)
  expect_equal(rhos[c("rho_tl_rl", "rho_tl_rn")],
               prior[c("rho_tl_rl", "rho_tl_rn")])
})

test_that("a default-length run adapts every joint-move correlation", {
  skip_under_memcheck()
  skip_on_cran()
  set.seed(301)
  tree <- ape::rtree(20)
  tree$edge.length <- tree$edge.length / 5
  Simulate <- function(k) {
    as.character(ape::rTraitDisc(tree, k = k, states = as.character(seq_len(k) - 1L)))
  }
  Variable <- function(m) m[, apply(m, 2, function(x) length(unique(x)) > 1)]
  trans <- Variable(replicate(50, Simulate(sample(2:3, 1))))
  neo <- Variable(replicate(10, Simulate(2)))
  mat <- cbind(trans, neo)
  rownames(mat) <- tree$tip.label

  seen <- new.env()
  buildRhos <- MkPrime:::.BuildJointRhoMatrix
  local_mocked_bindings(
    .BuildJointRhoMatrix = function(chainRhos, moves, nChains) {
      seen$rhos <- chainRhos[[1]]
      buildRhos(chainRhos, moves, nChains)
    }
  )
  setTimeLimit(elapsed = 240, transient = TRUE)
  on.exit(setTimeLimit(elapsed = Inf, transient = TRUE), add = TRUE)
  result <- allow_warning(
    RunMkPrime(MatrixToPhyDat(mat), neomorphic = ncol(trans) + seq_len(ncol(neo)),
               mcmc = MkPrimeMCMC(nIter = 30000L, nRuns = 1L, nCore = 1L,
                                  maxTime = 180)),
    "without stabilisation|clamped"
  )
  expect_gt(nrow(result$samples), 0L)
  # Warmup ends in far fewer than the 50 batches its snapshots need, so the
  # estimate comes from Tuning; before #301 all three stayed at 0.
  expect_true(all(unlist(seen$rhos) != 0))
})

.RhoData <- function() {
  set.seed(302)
  tree <- ape::rtree(10)
  mat <- replicate(25, as.character(
    ape::rTraitDisc(tree, k = 3, rate = 2, states = c("0", "1", "2"))
  ))
  rownames(mat) <- tree$tip.label
  list(tree = tree, pd = MatrixToPhyDat(mat))
}

test_that("without tuning, rho is estimated from the first Sample rows (#405)", {
  skip_under_memcheck()
  skip_on_cran()
  data <- .RhoData()
  seen <- new.env()
  buildRhos <- MkPrime:::.BuildJointRhoMatrix
  local_mocked_bindings(
    .BuildJointRhoMatrix = function(chainRhos, moves, nChains) {
      seen$rhos <- chainRhos[[1]]
      buildRhos(chainRhos, moves, nChains)
    }
  )
  setTimeLimit(elapsed = 240, transient = TRUE)
  on.exit(setTimeLimit(elapsed = Inf, transient = TRUE), add = TRUE)
  set.seed(4051)
  # Default warmup, nIter / 2, holds at most 20 of the 50 snapshots it needs.
  result <- allow_warning(
    RunMkPrime(data$pd, data$tree, mcmc = MkPrimeMCMC(
      nIter = 20000L, nRuns = 1L, nCore = 1L, autoTune = FALSE
    )),
    "stabilis"
  )
  expect_gt(nrow(result$samples), 0L)
  expect_true(seen$rhos$rho_tl_rls != 0)
})

test_that("a rho pending in Sample survives a checkpoint (#405)", {
  skip_on_cran()
  data <- .RhoData()
  ckpFile <- tempfile(fileext = ".ckp")
  cancelFile <- tempfile()
  on.exit(unlink(c(ckpFile, cancelFile)), add = TRUE)
  setTimeLimit(elapsed = 240, transient = TRUE)
  on.exit(setTimeLimit(elapsed = Inf, transient = TRUE), add = TRUE)
  Job <- function(...) {
    MkPrimeMCMC(nIter = 4000L, minWarmup = 2000L, maxWarmup = 2000L,
                nRuns = 1L, nCore = 1L, thin = 20L, autoTune = FALSE,
                checkpointFile = ckpFile, cancelFile = cancelFile, ...)
  }
  # 500 iterations at thin 20 save 25 rows: half of what rho needs.
  Stopper <- function(info) if (info$iter >= 2500L) file.create(cancel)
  environment(Stopper) <- list2env(list(cancel = cancelFile),
                                   parent = baseenv())
  set.seed(4052)
  allow_warning(
    RunMkPrime(data$pd, data$tree,
               mcmc = Job(progressFn = Stopper, plotEvery = 500L)),
    "stabilis"
  )
  cp <- readRDS(ckpFile)
  expect_identical(cp$runs[[1]]$phase, "Sample")
  expect_true(cp$runs[[1]]$rhoPending)
  expect_equal(nrow(cp$runs[[1]]$rhoSampleBuf), 25L)
  expect_equal(cp$runs[[1]]$chain_rhos[[1]]$rho_tl_rls, 0)

  # The next 25 rows complete the estimate.
  cp$mcmc$progressFn <- NULL
  saveRDS(cp, ckpFile)
  file.create(cancelFile)
  allow_warning(
    ResumeMkPrime(ckpFile, data$pd, mcmc = list(plotEvery = 500L,
                                                progressFn = function(i) NULL)),
    "stabilis"
  )
  after <- readRDS(ckpFile)$runs[[1]]
  expect_false(after$rhoPending)
  expect_null(after$rhoSampleBuf)
  expect_true(after$chain_rhos[[1]]$rho_tl_rls != 0)
})

test_that("reported acceptance covers the Sample phase alone (#405)", {
  skip_on_cran()
  skip_under_memcheck()
  data <- .RhoData()
  ckpFile <- tempfile(fileext = ".ckp")
  on.exit(unlink(ckpFile), add = TRUE)
  setTimeLimit(elapsed = 240, transient = TRUE)
  on.exit(setTimeLimit(elapsed = Inf, transient = TRUE), add = TRUE)
  for (autoTune in c(FALSE, TRUE)) {
    set.seed(4053)
    allow_warning(
      RunMkPrime(data$pd, data$tree, overwrite = TRUE, mcmc = MkPrimeMCMC(
        nIter = 6000L, minWarmup = 2000L, maxWarmup = 2000L, nRuns = 1L,
        nCore = 1L, thin = 10L, autoTune = autoTune, checkpointFile = ckpFile
      )),
      "stabilis"
    )
    r <- readRDS(ckpFile)$runs[[1]]
    expect_identical(r$phase, "Sample")
    tuningIter <- if (autoTune) r$tuningIterUsed else 0L
    expect_equal(sum(r$chain_propose[[1]]), 6000L - 2000L - tuningIter)
  }
})
