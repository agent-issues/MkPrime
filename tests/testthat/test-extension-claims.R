# Extension ABI v3 (#441): integer parameters, core scalars in the tree view,
# and claims. toy_claim_extension_cpp() claims the neomorphic characters and
# scores them as MkPrime does (MkN at rate_loss, neo_scale, ACRV rates), so
# MkPrime's likelihood plus the extension's must equal MkPrime's alone. Its
# integer toy_switch has prior P(1) = 0.7 and a "flip" move; toy_mu ~
# Gamma(2, 4) enters only the prior.

ToyClaim <- function() MkPrimeExtension(function() toy_claim_extension_cpp())

ClaimFixture <- function(seed = 4422L) {
  set.seed(seed)
  ntax <- 8L
  mat <- matrix(sample.int(3, ntax * 14, replace = TRUE) - 1L,
                nrow = ntax, ncol = 14,
                dimnames = list(paste0("t", seq_len(ntax)), NULL))
  mat[, 1:6] <- pmin(mat[, 1:6], 1L)
  mat <- matrix(as.character(mat), nrow(mat), dimnames = dimnames(mat))
  mat[2, 3] <- "?"
  pd <- TreeTools::MatrixToPhyDat(mat)
  tree <- TreeTools::Preorder(
    ape::rtree(ntax, tip.label = rownames(mat),
               br = function(n) runif(n, 0.05, 0.3))
  )
  list(pd = pd, mkd = MkPrimeData(pd, neomorphic = 1:6), tree = tree)
}

# Data and state pointers as RunMkPrime() builds them.
ClaimPointers <- function(fx, model) {
  model <- MkPrime:::.FinalizeModel(model, fx$tree, fx$mkd)
  state <- MkPrime:::.InitState(fx$tree, fx$mkd, model)
  dataPtr <- MkPrime:::.InitMcmcData(fx$mkd, model)
  statePtr <- MkPrime:::.InitMcmcChain(state)
  fill_partition_cache(dataPtr, statePtr)
  allocate_cl_workspace(dataPtr, statePtr)
  list(dataPtr = dataPtr, statePtr = statePtr, state = state)
}

test_that("a v3 extension describes its integers and claims", {
  info <- extension_info_cpp(toy_claim_extension_cpp())
  expect_identical(info$abi_version, 3L)
  expect_true(info$claims_neomorphic)
  expect_true(info$log_lik)
  expect_identical(info$ints, "toy_switch")
  expect_identical(info$int_init, 0L)
  expect_identical(info$int_levels, 2L)
  ext <- ToyClaim()
  expect_identical(ext$ints, c(toy_switch = 0L))
  model <- MkPrimeModel(coding = "none", extensions = ext)
  expect_identical(MkPrime:::.ExtensionColumns(model),
                   c("toy_mu", "toy_switch"))
  expect_identical(MkPrime:::.ExtensionInitInts(model), 0L)
  expect_equal(MkPrime:::.ExtensionLogPrior(0.5, model, 1L),
               stats::dgamma(0.5, 2, 4, log = TRUE) + log(0.7),
               tolerance = 1e-12)
  expect_error(MkPrimeModel(extensions = list(ext, ext)), "defined twice")
  expect_error(extension_log_prior_cpp(toy_claim_extension_cpp(), 0.5),
               "1 integer parameters")
})

test_that("claimed characters are scored by the extension alone", {
  fx <- ClaimFixture()
  plain <- ClaimPointers(fx, MkPrimeModel(coding = "none"))
  claim <- ClaimPointers(fx, MkPrimeModel(coding = "none",
                                          extensions = ToyClaim()))
  expect_error(attach_extensions_cpp(claim$dataPtr,
                                     list(toy_claim_extension_cpp())),
               "already attached")
  s0 <- get_mcmc_state(claim$statePtr)
  plainLL <- get_state_log_lik(plain$statePtr)
  # The core no longer scores the neomorphic characters...
  expect_false(isTRUE(all.equal(s0$logLik, plainLL)))
  # ...but core plus extension is MkPrime's likelihood of all characters.
  expect_equal(s0$logLik + s0$extLogLik, plainLL, tolerance = 1e-9)
  expect_identical(s0$extInts, 0L)

  Move <- function(type, idx = 0L) {
    do_move_cpp(claim$dataPtr, claim$statePtr, type, idx, 0.6, 1.0, 1L, 1.0)
  }
  Check <- function(label) {
    s <- get_mcmc_state(claim$statePtr)
    full <- eval_full_loglik_cpp(plain$dataPtr, claim$statePtr)
    expect_equal(s$logLik + s$extLogLik, full, tolerance = 1e-9, label = label)
    expect_equal(s$logLik, eval_full_loglik_cpp(claim$dataPtr, claim$statePtr),
                 tolerance = 1e-9, label = label)
    expect_equal(s$logPrior, eval_log_prior_cpp(claim$dataPtr, claim$statePtr),
                 tolerance = 1e-10, label = label)
    s
  }
  set.seed(4423L)
  # Moves on the scalars the extension reads: tree_length, rate_loss,
  # rate_log_sd, rate_neo; then branch lengths and topology.
  for (i in 1:15) {
    for (type in c(0L, 1L, 2L, 3L, 4L, 5L)) {
      Move(type)
      Check(paste("move", type))
    }
  }
  s <- get_mcmc_state(claim$statePtr)
  expect_gt(s$rateLogSd, 0)

  # The integer move changes the prior only.
  before <- get_mcmc_state(claim$statePtr)
  accepted <- Move(41L, 0L)
  after <- Check("flip")
  if (accepted) {
    expect_identical(after$extInts, 1L - before$extInts)
    expect_equal(after$logPrior - before$logPrior,
                 log(c(0.3, 0.7)[after$extInts + 1L]) -
                   log(c(0.3, 0.7)[before$extInts + 1L]), tolerance = 1e-12)
  } else {
    expect_identical(after$extInts, before$extInts)
  }
  expect_identical(after$logLik, before$logLik)
  expect_identical(after$extLogLik, before$extLogLik)

  # Slice moves read only the core likelihood, so are refused.
  expect_false(do_move_cpp(claim$dataPtr, claim$statePtr, 19L,
                           MkPrime:::.ScalarParamIndex("rate_loss"), 0.6, 1.0,
                           1L, 1.0))
})

test_that("RunMkPrime samples and checkpoints integer parameters", {
  skip_on_cran()
  fx <- ClaimFixture(4424L)
  dir <- withr::local_tempdir()
  ckp <- file.path(dir, "run.ckp")
  model <- MkPrimeModel(coding = "none", extensions = list(ToyClaim()))
  set.seed(4425L)
  res <- allow_warning(RunMkPrime(
    fx$pd, fx$tree, neomorphic = 1:6, model = model,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 6000L, thin = 5L,
                       minWarmup = 1000L, maxWarmup = 1000L, autoTune = FALSE,
                       checkpointFile = ckp, maxTime = 120)
  ), "maxWarmup|stabilis")
  s <- res$samples
  expect_true(all(c("toy_mu", "toy_switch", "rate_loss", "rate_neo") %in%
                    colnames(s)))
  expect_true(all(s[, "toy_switch"] %in% 0:1))
  # Prior P(toy_switch = 1) = 0.7; nothing else involves it.
  expect_lt(abs(mean(s[, "toy_switch"]) - 0.7), 0.1)
  expect_true(all(is.finite(s[, "log_likelihood"])))
  expect_gt(length(unique(s[, "rate_loss"])), 10L)
  expect_false(any(c("slice_rate_loss", "slice_tree_length") %in%
                     names(res$acceptance)))

  ch <- readRDS(ckp)$runs[[1]]$chains[[1]]
  expect_length(ch$ext_ints, 1L)
  resumed <- allow_warning(ResumeMkPrime(
    ckp, fx$pd, fx$tree, neomorphic = 1:6,
    mcmc = list(progressFn = function(info) NULL, nIter = 7000L)
  ), "maxWarmup|stabilis")
  expect_gt(nrow(resumed$samples), nrow(s))
  expect_true(all(resumed$samples[, "toy_switch"] %in% 0:1))
})
