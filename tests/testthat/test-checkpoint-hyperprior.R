# Checkpoint round-trip for the pooled-σ hyperprior fields.
#
# Commit c6310ef ("feat(partition): pooled half-normal hyperprior on per-class
# σ_c") added two new chain-state fields — `hyper_tau` (population scale τ)
# and `class_rate_log_sd_z` (per-class non-centred z_c). `.InitRun` in
# R/RunMkPrime.R (~L669) writes them into the *initial* per-chain list, and
# `.RunMkPrimeSingleRun` (~L769) reads them on reconstruction. But the
# C++-state serialiser inside `.SaveCheckpoint` (~L2616) and the end-of-run
# serialiser inside `.RunMkPrimeSingleRun` (~L1625) hard-code the field list
# and DO NOT pull `hyperTau` / `classZ` out of `get_mcmc_state()`. External
# review flagged that no existing test exercises the on-disk serialise →
# load → restore path, so a missing-field regression here would slip through.
#
# This test runs a short hyperprior-active chain to a checkpoint and asserts:
#   (1) the on-disk checkpoint carries `hyper_tau` and `class_rate_log_sd_z`
#       on the (single) chain in the (single) run;
#   (2) both values are finite and positive (HN(1) support);
#   (3) `class_rate_log_sd_z` has length K = 2 matching the partition;
#   (4) ResumeMkPrime() reads that checkpoint and the resumed chain advances.
#
# K = 2 is achieved with partition = c(1,1,1,2,2,2,2) (3 chars in class 1,
# 4 in class 2). With `unlink = "shape"` and the default
# `priorOnClassRateLogSd = "hyperprior_pooled"`, `.InitStatePartitioned`
# activates the hyperprior branch (R/partition-api.R:278-284).

test_that("checkpoint round-trip preserves hyper_tau and class_rate_log_sd_z", {
  skip_on_cran()

  # 6 tips, 7 binary chars, all variable (so MkPrimeData keeps every column).
  set.seed(20260528L)
  nTip  <- 6L
  nChar <- 7L
  mat   <- matrix(sample(0:1, nTip * nChar, replace = TRUE),
                   nrow = nTip, ncol = nChar,
                   dimnames = list(paste0("t", seq_len(nTip)), NULL))
  for (j in seq_len(nChar)) {
    if (length(unique(mat[, j])) < 2L) mat[1L, j] <- 1L - mat[1L, j]
  }
  pd   <- MatrixToPhyDat(mat)
  tree <- Preorder(NJTree(pd, edgeLengths = TRUE))
  if (is.null(tree$edge.length) || any(tree$edge.length <= 0)) {
    tree$edge.length <- rep(0.1, nrow(tree$edge))
  }

  td <- tempfile("mkp_ckp_hyper_")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  ckp <- file.path(td, "run.ckp")
  log <- file.path(td, "run.log")

  partition <- c(1L, 1L, 1L, 2L, 2L, 2L, 2L)
  stopifnot(length(partition) == nChar)
  K <- length(unique(partition))
  stopifnot(K == 2L)

  # Streaming v3 checkpoint requires logFile != NULL; checkEvery = 50 inside
  # nIter = 200 with warmup = 100 means a checkpoint fires both during warmup
  # and during sample phase, so the hyperprior fields must survive both
  # serialiser paths.
  set.seed(8675309L)
  result1 <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 1L, nIter = 200L, thin = 5L,
                        maxWarmup = 100L, minWarmup = 100L, autoTune = FALSE,
                        checkEvery = 50L,
                        checkpointFile = ckp, logFile = log),
    partition = partition,
    unlink    = "shape"
  )

  expect_true(file.exists(ckp))

  cp <- readRDS(ckp)
  # Streaming checkpoints are v3 (R/RunMkPrime.R:2650).
  expect_equal(cp$version, 3L)
  expect_length(cp$runs, 1L)
  ch <- cp$runs[[1L]]$chains[[1L]]

  # (1) Fields must be present on the on-disk chain list. These two
  # expectations are the regression guard — they fail today because
  # .SaveCheckpoint's get_mcmc_state() block (R/RunMkPrime.R:2616) does
  # not copy s$hyperTau or s$classZ into the chain list.
  expect_false(is.null(ch$hyper_tau),
               label = "checkpoint chain list has hyper_tau field")
  expect_false(is.null(ch$class_rate_log_sd_z),
               label = "checkpoint chain list has class_rate_log_sd_z field")

  # (2) Half-normal support: finite and strictly positive.
  expect_true(is.finite(ch$hyper_tau))
  expect_gt(ch$hyper_tau, 0)
  expect_true(all(is.finite(ch$class_rate_log_sd_z)))
  expect_true(all(ch$class_rate_log_sd_z > 0))

  # (3) z_c length matches K = 2.
  expect_length(ch$class_rate_log_sd_z, K)

  # Capture pre-resume values for the post-resume comparison.
  tauBefore <- ch$hyper_tau
  zBefore   <- ch$class_rate_log_sd_z

  # (4) Resume must actually advance the chain. ResumeMkPrime does not
  # currently accept a `partition` argument (R/RunMkPrime.R:2802), so we
  # supply only the saved model; the partitioned state must come from the
  # checkpoint payload. If `init_mcmc_state` silently drops the hyperprior
  # because the saved chain list is missing the fields, the resumed run
  # advances but on the legacy single-σ path — that's still detected by
  # (1) above. If the fields are present, the resumed run is the same
  # model as the original.
  resumed <- suppressWarnings(ResumeMkPrime(ckp, pd, tree))
  expect_s3_class(resumed, "MkPosterior")
  expect_gt(resumed$nSamples, 0L)
  expect_true(file.exists(log) || file.exists(file.path(td, "run_1.log")))

  # Cross-check: after the resume completes the master checkpoint is
  # rewritten. The pooled-σ fields must still be present (and still
  # finite, positive, length K) — i.e. the resume path itself round-trips
  # them too.
  cp2 <- readRDS(ckp)
  ch2 <- cp2$runs[[1L]]$chains[[1L]]
  expect_false(is.null(ch2$hyper_tau),
               label = "post-resume checkpoint still carries hyper_tau")
  expect_false(is.null(ch2$class_rate_log_sd_z),
               label = "post-resume checkpoint still carries class_rate_log_sd_z")
  if (!is.null(ch2$hyper_tau) && !is.null(ch2$class_rate_log_sd_z)) {
    expect_true(is.finite(ch2$hyper_tau) && ch2$hyper_tau > 0)
    expect_length(ch2$class_rate_log_sd_z, K)
    expect_true(all(is.finite(ch2$class_rate_log_sd_z) &
                      ch2$class_rate_log_sd_z > 0))
  }

  # Sanity: the *initial* run's values lay in the same support; if the
  # serialiser dropped them, the pre-resume tauBefore/zBefore would be NULL
  # and these comparisons would error — guard with is.null() so the failure
  # mode stays attributable to (1) rather than this auxiliary check.
  if (!is.null(tauBefore)) expect_true(is.finite(tauBefore) && tauBefore > 0)
  if (!is.null(zBefore))   expect_length(zBefore, K)
})
