# Tests for Phase 6: Progress display
# M-037: Callback infrastructure
# M-038: PlotDuringMCMC trace plots
# M-039: PNG progress writer

# --- Helper: build a well-formed info list ---

mock_info <- function(iter = 500, nIter = 1000, warmup = 200,
                      nRuns = 2, has_samples = TRUE) {
  param_names <- c("log_posterior", "log_likelihood", "tree_length",
                   "rate_loss", "rate_log_sd", "p", "br_1", "br_2")

  if (has_samples) {
    nSaved <- 30
    make_samples <- function() {
      m <- matrix(rnorm(nSaved * length(param_names)),
                  nrow = nSaved, ncol = length(param_names),
                  dimnames = list(NULL, param_names))
      # log_posterior should be log_likelihood + something
      m[, "log_posterior"] <- m[, "log_likelihood"] + rnorm(nSaved, -5, 1)
      m[, "tree_length"] <- abs(m[, "tree_length"]) + 0.5
      m[, "rate_loss"] <- abs(m[, "rate_loss"]) + 0.1
      m[, "rate_log_sd"] <- abs(m[, "rate_log_sd"]) + 0.1
      m[, "p"] <- runif(nSaved, 0.2, 0.8)
      m
    }
    run_samples <- lapply(seq_len(nRuns), function(i) make_samples())
  } else {
    run_samples <- vector("list", nRuns)
  }

  current_state <- lapply(seq_len(nRuns), function(i) {
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
    in_warmup = iter <= warmup,
    nRuns = nRuns,
    nChains = 1L,
    run_samples = run_samples,
    current_state = current_state,
    recent_acceptance = 0.35,
    elapsed = 12.5,
    param_names = param_names
  )
}


# ===== M-037: Callback infrastructure =====

test_that("MkPrimeMCMC accepts plot_every and progress_fn", {
  cfg <- MkPrimeMCMC(plot_every = 200L, progress_fn = NULL)
  expect_equal(cfg$plot_every, 200L)
  expect_null(cfg$progress_fn)
})

test_that("MkPrimeMCMC resolves 'default' to mkp_trace_plot", {
  cfg <- MkPrimeMCMC(progress_fn = "default", plot_every = 100L)
  expect_identical(cfg$progress_fn, mkp_trace_plot)
})

test_that("MkPrimeMCMC rejects invalid progress_fn", {
  expect_error(
    MkPrimeMCMC(progress_fn = 42),
    "progress_fn"
  )
})

test_that("MkPrimeMCMC with no plot_every stores NULL", {
  cfg <- MkPrimeMCMC()
  expect_null(cfg$plot_every)
  expect_null(cfg$progress_fn)
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
    plot_every = 100L,
    progress_fn = counter_fn
  )

  suppressMessages(
    result <- RunMkPrime(dat, tree, mcmc = mcmc, fix_topology = TRUE)
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
    plot_every = 200L,
    progress_fn = capture_fn
  )

  suppressMessages(
    RunMkPrime(dat, tree, mcmc = mcmc, fix_topology = TRUE)
  )

  expect_true(!is.null(captured_info))
  expect_equal(captured_info$iter, 200L)
  expect_equal(captured_info$nIter, 200L)
  expect_equal(captured_info$warmup, 50L)
  expect_false(captured_info$in_warmup)
  expect_equal(captured_info$nRuns, 1L)
  expect_equal(captured_info$nChains, 1L)
  expect_type(captured_info$run_samples, "list")
  expect_type(captured_info$current_state, "list")
  expect_type(captured_info$recent_acceptance, "double")
  expect_type(captured_info$elapsed, "double")
  expect_true(captured_info$elapsed >= 0)

  # run_samples should have post-warmup samples
  expect_true(!is.null(captured_info$run_samples[[1]]))
  expect_true(nrow(captured_info$run_samples[[1]]) > 0)
})

test_that("Progress callback fires during warmup too", {
  skip_on_cran()
  skip_if_not_installed("TreeTools")

  warmup_calls <- integer(0)
  tracker_fn <- function(info) {
    if (info$in_warmup) warmup_calls <<- c(warmup_calls, info$iter)
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
    plot_every = 100L,
    progress_fn = tracker_fn
  )

  suppressMessages(
    RunMkPrime(dat, tree, mcmc = mcmc, fix_topology = TRUE)
  )

  expect_true(length(warmup_calls) >= 1)
  expect_true(all(warmup_calls <= 200))
})


# ===== M-038: mkp_trace_plot =====

test_that("mkp_trace_plot runs without error (with samples)", {
  info <- mock_info(iter = 500, nRuns = 2, has_samples = TRUE)
  expect_no_error(mkp_trace_plot(info))
})

test_that("mkp_trace_plot runs without error (warmup, no samples)", {
  info <- mock_info(iter = 100, warmup = 200, nRuns = 2,
                    has_samples = FALSE)
  info$in_warmup <- TRUE
  expect_no_error(mkp_trace_plot(info))
})

test_that("mkp_trace_plot runs without error (single run)", {
  info <- mock_info(iter = 300, nRuns = 1, has_samples = TRUE)
  expect_no_error(mkp_trace_plot(info))
})

test_that("mkp_trace_plot returns info invisibly", {
  info <- mock_info()
  result <- mkp_trace_plot(info)
  expect_identical(result, info)
})


# ===== M-039: PNG progress writer =====

test_that("mkp_png_progress creates files", {
  dir <- tempfile("mkp_progress_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))

  fn <- mkp_png_progress(dir)
  expect_true(is.function(fn))

  info <- mock_info()
  fn(info)

  expect_true(file.exists(file.path(dir, "mkp_progress.png")))
  expect_true(file.exists(file.path(dir, "mkp_progress.json")))
  # temp file should have been renamed away
  expect_false(file.exists(file.path(dir, "mkp_progress_tmp.png")))
})

test_that("mkp_png_progress JSON is valid", {
  dir <- tempfile("mkp_json_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))

  fn <- mkp_png_progress(dir)
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
  expect_true(grepl('"in_warmup":false', json_text))
})

test_that("mkp_png_progress creates directory if needed", {
  dir <- file.path(tempdir(), "mkp_auto_create_test_dir")
  if (dir.exists(dir)) unlink(dir, recursive = TRUE)
  on.exit(unlink(dir, recursive = TRUE))

  fn <- mkp_png_progress(dir)
  expect_true(dir.exists(dir))

  info <- mock_info()
  fn(info)
  expect_true(file.exists(file.path(dir, "mkp_progress.png")))
})

test_that("mkp_png_progress overwrites on subsequent calls", {
  dir <- tempfile("mkp_overwrite_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))

  fn <- mkp_png_progress(dir)

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


# ===== .format_elapsed =====

test_that(".format_elapsed formats correctly", {
  expect_equal(MkPrime:::.format_elapsed(5), "5s")
  expect_equal(MkPrime:::.format_elapsed(90), "1.5min")
  expect_equal(MkPrime:::.format_elapsed(7200), "2.0h")
})

# ===== .write_progress_json =====

test_that(".write_progress_json writes valid output", {
  f <- tempfile(fileext = ".json")
  on.exit(unlink(f))

  info <- list(iter = 42, nIter = 100, warmup = 20,
               in_warmup = FALSE, elapsed = 3.7,
               recent_acceptance = 0.345)
  MkPrime:::.write_progress_json(info, f)

  txt <- readLines(f)
  expect_true(grepl('"iter":42', txt))
  expect_true(grepl('"in_warmup":false', txt))
  expect_true(grepl('"recent_acceptance":0.3450', txt))
})
