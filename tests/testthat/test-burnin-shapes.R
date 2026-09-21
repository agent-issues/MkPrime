# Shape regressions for the burnin API (#90): a streaming result leaves
# `samples` empty, and runs end at different iterations. The objects are
# built directly rather than sampled, so both shapes are reachable without
# the timing pressure a real run needs to produce them.

# A posterior with the given per-run sample counts, and no MCMC behind it.
.SyntheticPosterior <- function(perRunN, logFile = NULL, seed = 90L) {
  set.seed(seed)
  paramNames <- c("log_posterior", "log_likelihood", "tree_length",
                  "rate_log_sd")
  Draws <- function(n) {
    matrix(rnorm(n * length(paramNames)), nrow = n,
           ncol = length(paramNames), dimnames = list(NULL, paramNames))
  }
  empty <- Draws(0L)

  post <- MkPosterior(
    samples = if (is.null(logFile) && length(perRunN) == 1L) {
      Draws(perRunN)
    } else {
      empty
    },
    trees = list(), acceptance = numeric(0),
    model = NULL, data = NULL, mcmc = list(thin = 1L),
    warmup = 0L, tuning = list()
  )
  post$nSamples <- sum(perRunN)
  post$logFile <- logFile
  if (length(perRunN) > 1L) {
    post$nRuns <- length(perRunN)
    post$per_run <- lapply(perRunN, function(n) {
      list(samples = Draws(n), trees = list(), acceptance = numeric(0))
    })
  }
  post
}

# A Tracer-compatible log holding `n` samples, as a streaming run writes.
.WriteTestLog <- function(n, seed = 91L) {
  set.seed(seed)
  path <- withr::local_tempfile(fileext = ".log", .local_envir = parent.frame())
  utils::write.table(
    data.frame(Sample = seq_len(n), log_posterior = rnorm(n),
               log_likelihood = rnorm(n), tree_length = runif(n),
               rate_log_sd = runif(n)),
    path, sep = "\t", row.names = FALSE, quote = FALSE
  )
  path
}


test_that("SetBurnin loads samples a streaming run left on disk (#90)", {
  post <- .SyntheticPosterior(0L, logFile = .WriteTestLog(50L))
  expect_equal(nrow(post$samples), 0L)

  # The fraction resolves against the log's 50 rows, not the empty matrix.
  expect_equal(SetBurnin(post, 0.2)$burnin, 10L)
  expect_equal(SetBurnin(post, 5)$burnin, 5L)
})


test_that("AutoBurnin leaves a usable object when it cannot choose (#90)", {
  post <- .SyntheticPosterior(0L, logFile = .WriteTestLog(5L))

  # Every candidate is filtered out by the "leave at least 10 samples" rule,
  # so there is nothing to select; the object must stay usable regardless.
  expect_warning(out <- AutoBurnin(post), "Too few samples")
  expect_null(out$burnin)
  expect_equal(nrow(.PostBurninData(out)$samples), 5L)

  # An existing burnin is carried through rather than overwritten.
  expect_warning(kept <- AutoBurnin(SetBurnin(post, 2L)), "Too few samples")
  expect_equal(kept$burnin, 2L)
})


test_that("burnin is validated against the shortest run (#90)", {
  post <- .SyntheticPosterior(c(30L, 12L))

  # 20 is valid for run 1 and impossible for run 2.
  expect_error(SetBurnin(post, 20), "must be less than")
  expect_equal(SetBurnin(post, 0.5)$burnin, 6L)

  kept <- .PostBurninData(SetBurnin(post, 8))$per_run
  expect_equal(nrow(kept[[1]]$samples), 22L)
  expect_equal(nrow(kept[[2]]$samples), 4L)
  # Retained draws keep their original order.
  expect_equal(unname(kept[[2]]$samples[1, ]),
               unname(post$per_run[[2]]$samples[9, ]))
})


test_that(".KeepAfter is empty, never descending, once burnin exhausts a run", {
  expect_equal(.KeepAfter(0L, 5L), 1:5)
  expect_equal(.KeepAfter(3L, 5L), 4:5)
  expect_equal(.KeepAfter(5L, 5L), integer(0))
  expect_equal(.KeepAfter(9L, 5L), integer(0))
})


test_that("burnin resolves per-run logs for a streaming multi-run result", {
  post <- .SyntheticPosterior(c(0L, 0L),
                              logFile = c(.WriteTestLog(40L, seed = 1L),
                                          .WriteTestLog(16L, seed = 2L)))
  post$per_run <- lapply(post$per_run, function(r) {
    r$samples <- NULL
    r
  })

  # The shortest run holds 16, so half of it is 8 -- not half of run 1's 40.
  expect_equal(SetBurnin(post, 0.5)$burnin, 8L)
  expect_error(SetBurnin(post, 20), "must be less than")

  kept <- .PostBurninData(SetBurnin(post, 8))$per_run
  expect_equal(nrow(kept[[1]]$samples), 32L)
  expect_equal(nrow(kept[[2]]$samples), 8L)
})


test_that("SetBurnin declines a posterior with no samples at all (#90)", {
  # The PAR-006 zero-run result carries neither samples nor a log, so the
  # sample count is honestly 0 and there is nothing to discard.
  post <- .SyntheticPosterior(0L)

  expect_warning(out <- SetBurnin(post, 0.2), "No samples")
  expect_null(out$burnin)
})


test_that("AutoBurnin declines when no candidate is assessable (#90)", {
  # ESS is NA for a constant parameter, so every candidate ranks NA and
  # `which.max` returns integer(0) -- the second route to a 0-row `best`.
  post <- .SyntheticPosterior(100L)
  post$samples[] <- 1

  expect_warning(out <- AutoBurnin(post), "usable diagnostic")
  expect_null(out$burnin)
  expect_equal(nrow(.PostBurninData(out)$samples), 100L)
})


test_that(".PostBurninData says when a burnin empties a run", {
  post <- .SyntheticPosterior(c(30L, 12L))
  # Only reachable by assigning `$burnin` directly; SetBurnin rejects it.
  post$burnin <- 20L

  expect_warning(pb <- .PostBurninData(post), "discards every sample of run 2")
  expect_equal(nrow(pb$per_run[[2]]$samples), 0L)
})
