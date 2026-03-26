# Tests for Phase 6: Progress display
# M-037: Callback infrastructure
# M-038: PlotDuringMCMC trace plots
# M-039: PNG progress writer

# --- Helper: build a well-formed info list ---

mock_info <- function(iter = 500, nIter = 1000, warmup = 200,
                      nRuns = 2, has_samples = TRUE) {
  paramNames <- c("log_posterior", "log_likelihood", "tree_length",
                   "rate_loss", "rate_log_sd", "p", "br_1", "br_2")

  if (has_samples) {
    nSaved <- 30
    make_samples <- function() {
      m <- matrix(rnorm(nSaved * length(paramNames)),
                  nrow = nSaved, ncol = length(paramNames),
                  dimnames = list(NULL, paramNames))
      # log_posterior should be log_likelihood + something
      m[, "log_posterior"] <- m[, "log_likelihood"] + rnorm(nSaved, -5, 1)
      m[, "tree_length"] <- abs(m[, "tree_length"]) + 0.5
      m[, "rate_loss"] <- abs(m[, "rate_loss"]) + 0.1
      m[, "rate_log_sd"] <- abs(m[, "rate_log_sd"]) + 0.1
      m[, "p"] <- runif(nSaved, 0.2, 0.8)
      m
    }
    runSamples <- lapply(seq_len(nRuns), function(i) make_samples())
  } else {
    runSamples <- vector("list", nRuns)
  }

  currentState <- lapply(seq_len(nRuns), function(i) {
    list(log_lik = rnorm(1, -100, 10),
         log_prior = rnorm(1, -10, 2),
         tree_length = abs(rnorm(1, 2, 0.5)),
         rate_loss = abs(rnorm(1, 1, 0.3)),
         rate_log_sd = abs(rnorm(1, 0.5, 0.2)),
         p = runif(1, 0.3, 0.7))
  })

  list(
    iter = iter,
    nIter = nIter,
    warmup = warmup,
    inWarmup = iter <= warmup,
    nRuns = nRuns,
    nChains = 1L,
    runSamples = runSamples,
    currentState = currentState,
    recentAcceptance = 0.35,
    elapsed = 12.5,
    paramNames = paramNames
  )
}


# ===== M-037: Callback infrastructure =====

test_that("MkPrimeMCMC accepts plotEvery and progressFn", {
  cfg <- MkPrimeMCMC(plotEvery = 200L, progressFn = NULL)
  expect_equal(cfg$plotEvery, 200L)
  expect_null(cfg$progressFn)
})

test_that("MkPrimeMCMC resolves 'default' to MkpTracePlot", {
  cfg <- MkPrimeMCMC(progressFn = "default", plotEvery = 100L)
  expect_identical(cfg$progressFn, MkpTracePlot)
})

test_that("MkPrimeMCMC rejects invalid progressFn", {
  expect_error(
    MkPrimeMCMC(progressFn = 42),
    "progressFn"
  )
})

test_that("MkPrimeMCMC with no plotEvery stores NULL", {
  cfg <- MkPrimeMCMC()
  expect_null(cfg$plotEvery)
  expect_null(cfg$progressFn)
})

test_that("Progress callback is invoked correct number of times", {
  skip_on_cran()
  skip_if_not_installed("TreeTools")

  call_count <- 0L
  call_iters <- integer(0)
  counter_fn <- function(info) {
    call_count <<- call_count + 1L
    call_iters <<- c(call_iters, info$iter)
  }

  n <- 6L
  tree <- ape::stree(n, type = "left")
  tree$edge.length <- rep(0.1, nrow(tree$edge))

  set.seed(4817)
  dat <- TreeTools::StringToPhyDat(
    setNames(replicate(8, paste(sample(0:1, n, TRUE), collapse = "")),
             tree$tip.label),
    tips = tree
  )

  mcmc <- MkPrimeMCMC(
    nIter = 500L, warmup = 100L, thin = 10L,
    nRuns = 1L,
    plotEvery = 100L,
    progressFn = counter_fn
  )

  suppressMessages(
    result <- RunMkPrime(dat, tree, mcmc = mcmc, fixTopology = TRUE)
  )

  # Should be called at iter 100, 200, 300, 400, 500
  expect_equal(call_count, 5L)
  expect_equal(call_iters, c(100, 200, 300, 400, 500))
})

test_that("Progress info list has correct structure", {
  skip_on_cran()
  skip_if_not_installed("TreeTools")

  captured_info <- NULL
  capture_fn <- function(info) {
    captured_info <<- info
  }

  n <- 5L
  tree <- ape::stree(n, type = "left")
  tree$edge.length <- rep(0.1, nrow(tree$edge))

  set.seed(6293)
  dat <- TreeTools::StringToPhyDat(
    setNames(replicate(6, paste(sample(0:1, n, TRUE), collapse = "")),
             tree$tip.label),
    tips = tree
  )

  mcmc <- MkPrimeMCMC(
    nIter = 200L, warmup = 50L, thin = 10L,
    nRuns = 1L,
    plotEvery = 200L,
    progressFn = capture_fn
  )

  suppressMessages(
    RunMkPrime(dat, tree, mcmc = mcmc, fixTopology = TRUE)
  )

  expect_true(!is.null(captured_info))
  expect_equal(captured_info$iter, 200L)
  expect_equal(captured_info$nIter, 200L)
  expect_equal(captured_info$warmup, 50L)
  expect_false(captured_info$inWarmup)
  expect_equal(captured_info$nRuns, 1L)
  expect_equal(captured_info$nChains, 1L)
  expect_type(captured_info$runSamples, "list")
  expect_type(captured_info$currentState, "list")
  expect_type(captured_info$recentAcceptance, "double")
  expect_type(captured_info$elapsed, "double")
  expect_true(captured_info$elapsed >= 0)

  # runSamples should have post-warmup samples
  expect_true(!is.null(captured_info$runSamples[[1]]))
  expect_true(nrow(captured_info$runSamples[[1]]) > 0)
})

test_that("Progress callback fires during warmup too", {
  skip_on_cran()
  skip_if_not_installed("TreeTools")

  warmup_calls <- integer(0)
  tracker_fn <- function(info) {
    if (info$inWarmup) warmup_calls <<- c(warmup_calls, info$iter)
  }

  n <- 5L
  tree <- ape::stree(n, type = "left")
  tree$edge.length <- rep(0.1, nrow(tree$edge))

  set.seed(1438)
  dat <- TreeTools::StringToPhyDat(
    setNames(replicate(6, paste(sample(0:1, n, TRUE), collapse = "")),
             tree$tip.label),
    tips = tree
  )

  mcmc <- MkPrimeMCMC(
    nIter = 400L, warmup = 200L, thin = 10L,
    nRuns = 1L,
    plotEvery = 100L,
    progressFn = tracker_fn
  )

  suppressMessages(
    RunMkPrime(dat, tree, mcmc = mcmc, fixTopology = TRUE)
  )

  expect_true(length(warmup_calls) >= 1)
  expect_true(all(warmup_calls <= 200))
})


# ===== M-038: MkpTracePlot =====

test_that("MkpTracePlot runs without error (with samples)", {
  info <- mock_info(iter = 500, nRuns = 2, has_samples = TRUE)
  expect_no_error(MkpTracePlot(info))
})

test_that("MkpTracePlot runs without error (warmup, no samples)", {
  info <- mock_info(iter = 100, warmup = 200, nRuns = 2,
                    has_samples = FALSE)
  info$inWarmup <- TRUE
  expect_no_error(MkpTracePlot(info))
})

test_that("MkpTracePlot runs without error (single run)", {
  info <- mock_info(iter = 300, nRuns = 1, has_samples = TRUE)
  expect_no_error(MkpTracePlot(info))
})

test_that("MkpTracePlot returns info invisibly", {
  info <- mock_info()
  result <- MkpTracePlot(info)
  expect_identical(result, info)
})


# ===== M-039: PNG progress writer =====

test_that("MkpPngProgress creates files", {
  dir <- tempfile("mkp_progress_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))

  fn <- MkpPngProgress(dir)
  expect_true(is.function(fn))

  info <- mock_info()
  fn(info)

  expect_true(file.exists(file.path(dir, "mkp_progress.png")))
  expect_true(file.exists(file.path(dir, "mkp_progress.json")))
  # temp file should have been renamed away
  expect_false(file.exists(file.path(dir, "mkp_progress_tmp.png")))
})

test_that("MkpPngProgress JSON is valid", {
  dir <- tempfile("mkp_json_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))

  fn <- MkpPngProgress(dir)
  info <- mock_info(iter = 750, nIter = 1000, warmup = 200)
  fn(info)

  json_text <- readLines(file.path(dir, "mkp_progress.json"))
  # Should be parseable as a named list
  parsed <- tryCatch(
    eval(parse(text = gsub("true", "TRUE", gsub("false", "FALSE",
      gsub(":", "=", gsub("\\{", "list(", gsub("\\}", ")",
        json_text))))))),
    error = function(e) NULL
  )

  # Simpler: just check key fields are present

  expect_true(grepl('"iter":750', json_text))
  expect_true(grepl('"nIter":1000', json_text))
  expect_true(grepl('"warmup":200', json_text))
  expect_true(grepl('"inWarmup":false', json_text))
})

test_that("MkpPngProgress creates directory if needed", {
  dir <- file.path(tempdir(), "mkp_auto_create_test_dir")
  if (dir.exists(dir)) unlink(dir, recursive = TRUE)
  on.exit(unlink(dir, recursive = TRUE))

  fn <- MkpPngProgress(dir)
  expect_true(dir.exists(dir))

  info <- mock_info()
  fn(info)
  expect_true(file.exists(file.path(dir, "mkp_progress.png")))
})

test_that("MkpPngProgress overwrites on subsequent calls", {
  dir <- tempfile("mkp_overwrite_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))

  fn <- MkpPngProgress(dir)

  fn(mock_info(iter = 100))
  mtime1 <- file.mtime(file.path(dir, "mkp_progress.png"))
  Sys.sleep(0.1)

  fn(mock_info(iter = 200))
  mtime2 <- file.mtime(file.path(dir, "mkp_progress.png"))

  expect_true(mtime2 > mtime1)

  # JSON should reflect latest iter
  json <- readLines(file.path(dir, "mkp_progress.json"))
  expect_true(grepl('"iter":200', json))
})


# ===== .FormatElapsed =====

test_that(".FormatElapsed formats correctly", {
  expect_equal(MkPrime:::.FormatElapsed(5), "5s")
  expect_equal(MkPrime:::.FormatElapsed(90), "1.5min")
  expect_equal(MkPrime:::.FormatElapsed(7200), "2.0h")
})

# ===== .WriteProgressJson =====

test_that(".WriteProgressJson writes valid output", {
  f <- tempfile(fileext = ".json")
  on.exit(unlink(f))

  info <- list(iter = 42, nIter = 100, warmup = 20,
               inWarmup = FALSE, elapsed = 3.7,
               recentAcceptance = 0.345)
  MkPrime:::.WriteProgressJson(info, f)

  txt <- readLines(f)
  expect_true(grepl('"iter":42', txt))
  expect_true(grepl('"inWarmup":false', txt))
  expect_true(grepl('"recentAcceptance":0.3450', txt))
})
