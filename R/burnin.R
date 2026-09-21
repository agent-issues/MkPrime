# Burnin selection and filtering for MkPosterior objects
#
# M-066: The posterior object stores a per-run sample count `burnin`.
# All accessors (summary, print, plot, ConvergenceDiagnostics) respect
# this value. AutoBurnin() searches for the optimal burnin that
# maximizes ESS while maintaining low R-hat.

#' Set burnin for a posterior object
#'
#' @param posterior `MkPosterior` object whose burnin to set.
#' @param burnin Numeric specifying the number of per-run samples to discard;
#'   a value < 1 is instead a fraction of the shortest run's sample count.
#'
#' @return A new `MkPosterior` with the updated burnin.
#' @export
SetBurnin <- function(posterior, burnin) {
  if (!inherits(posterior, "MkPosterior")) {
    cli::cli_abort("{.arg posterior} must be an {.cls MkPosterior} object.")
  }

  nPerRun <- .MinPerRunSamples(.ResolveSamples(posterior))

  if (burnin < 1) {
    burnin <- floor(burnin * nPerRun)
  }
  burnin <- as.integer(burnin)

  if (nPerRun == 0L) {
    cli::cli_warn(c(
      "No samples to apply a burnin to.",
      "i" = "Burnin left at {posterior$burnin %||% 0L}."
    ))
    # Return: unchanged
    return(posterior)
  }

  if (burnin < 0L) burnin <- 0L
  if (burnin >= nPerRun) {
    cli::cli_abort(
      "Burnin ({burnin}) must be less than the number of per-run samples ({nPerRun})."
    )
  }

  posterior$burnin <- burnin
  posterior
}


#' Automatically select optimal burnin
#'
#' Searches over candidate burnin fractions (0% to 50% of samples) and
#' selects the smallest burnin where R-hat is acceptable
#' (max < `rhatThreshold`), or the burnin that minimizes max(R-hat) if
#' convergence is not achieved.
#'
#' For single-run posteriors (no R-hat available), selects the burnin
#' that maximizes min(ESS) across key parameters.
#'
#' @inheritParams SetBurnin
#' @param rhatThreshold Maximum acceptable R-hat. Default 1.05.
#' @param fractions Candidate burnin fractions to evaluate. Default
#'   `seq(0, 0.5, by = 0.05)`.
#'
#' @return A new `MkPosterior` with the selected burnin set.
#' @export
AutoBurnin <- function(posterior,
                       rhatThreshold = 1.05,
                       fractions = seq(0, 0.5, by = 0.05)) {
  if (!inherits(posterior, "MkPosterior")) {
    cli::cli_abort("{.arg posterior} must be an {.cls MkPosterior} object.")
  }

  nRuns <- posterior$nRuns %||% 1L
  hasRhat <- nRuns >= 2L && !is.null(posterior$per_run)

  resolved <- .ResolveSamples(posterior)
  nPerRun <- .MinPerRunSamples(resolved)

  results <- data.frame(
    fraction = fractions,
    burnin = as.integer(floor(fractions * nPerRun)),
    minEss = NA_real_,
    maxRhat = NA_real_,
    stringsAsFactors = FALSE
  )

  # Ensure burnin values leave at least 10 samples
  results <- results[results$burnin < (nPerRun - 10L), , drop = FALSE]

  # Bail out rather than let the selections below yield a 0-row `best`, and
  # so `burnin <- integer(0)`: that value passes `%||%`, which guards only
  # NULL, and breaks every later accessor far from its cause.
  if (nrow(results) == 0L) {
    cli::cli_warn(c(
      "Too few samples ({nPerRun} per run) to select a burnin.",
      "i" = "Burnin left at {posterior$burnin %||% 0L}."
    ))
    # Return: unchanged
    return(posterior)
  }

  for (i in seq_len(nrow(results))) {
    bi <- results$burnin[i]
    pb <- .PostBurninData(resolved, bi)

    keyCols <- .KeyParamCols(pb$samples)
    if (length(keyCols) == 0L) next

    ess <- .ComputeEss(pb$samples[, keyCols, drop = FALSE])
    results$minEss[i] <- .MinOrNA(ess)

    if (hasRhat && length(pb$per_run) >= 2L) {
      rhat <- .ComputeRhat(pb$per_run, keyCols)
      if (!is.null(rhat) && length(rhat) > 0L) {
        results$maxRhat[i] <- .MaxOrNA(rhat)
      }
    }
  }

  if (hasRhat) {
    # Strategy: smallest burnin where max(R-hat) <= threshold
    converged <- results[!is.na(results$maxRhat) &
                         results$maxRhat <= rhatThreshold, , drop = FALSE]
    missedTarget <- nrow(converged) == 0L
    best <- if (missedTarget) {
      # No burnin achieves target R-hat; pick the one with lowest max(R-hat)
      results[which.min(results$maxRhat), ]
    } else {
      # Among converged, pick smallest burnin (preserves most samples / ESS)
      converged[which.min(converged$burnin), ]
    }
  } else {
    # Single run: maximize min(ESS)
    missedTarget <- FALSE
    best <- results[which.max(results$minEss), ]
  }

  # `which.min`/`which.max` return integer(0) when the column they rank is
  # all-NA -- every parameter constant over every candidate window.
  if (nrow(best) != 1L) {
    cli::cli_warn(c(
      "No burnin candidate yields a usable diagnostic.",
      "i" = "Every monitored parameter is constant over the retained samples.",
      "i" = "Burnin left at {posterior$burnin %||% 0L}."
    ))
    # Return: unchanged
    return(posterior)
  }

  if (missedTarget) {
    cli::cli_warn(c(
      "No burnin fraction achieves max(Rhat) <= {rhatThreshold}.",
      "i" = "Selected burnin = {best$burnin} \
             (max Rhat = {round(best$maxRhat, 3)}).",
      "i" = "Consider running the chain longer."
    ))
  }

  .Inform(c(
    "v" = "Auto burnin: {best$burnin} samples ({round(best$fraction * 100)}% of {nPerRun})",
    "i" = "min(ESS) = {round(best$minEss, 1)}{if (hasRhat) paste0(', max(Rhat) = ', round(best$maxRhat, 3)) else ''}"
  ))

  posterior$burnin <- best$burnin
  posterior$auto_burnin_results <- results
  posterior
}


# Internal: load samples a streaming run left on disk
#
# `.BuildResult` installs an empty `samples` matrix whenever a run streamed to
# `logFile` -- the normal state of a streaming result, not an error state.
# Every entry point that needs a sample count or the samples themselves must
# resolve them first.
#
# @param posterior MkPosterior object
# @return `posterior` with `$samples` (and per-run `$samples`) populated
.ResolveSamples <- function(posterior) {
  if (nrow(posterior$samples) == 0L && !is.null(posterior$logFile)) {
    posterior$samples <- ReadMkLog(posterior$logFile)
  }

  nRuns <- posterior$nRuns %||% 1L
  if (nRuns > 1L && !is.null(posterior$per_run) &&
      !is.null(posterior$logFile)) {
    for (i in seq_along(posterior$per_run)) {
      if (is.null(posterior$per_run[[i]]$samples) &&
          i <= length(posterior$logFile) &&
          file.exists(posterior$logFile[i])) {
        posterior$per_run[[i]]$samples <- ReadMkLog(posterior$logFile[i])
      }
    }
  }

  # Return:
  posterior
}


# Internal: how many samples each run holds
#
# Runs routinely differ in length: each stops on its own ESS criterion, and
# `maxTime` cuts them at different iterations.
#
# @param posterior MkPosterior object, already passed through `.ResolveSamples`
# @return Integer vector, one count per run
.PerRunSampleCounts <- function(posterior) {
  nRuns <- posterior$nRuns %||% 1L
  if (nRuns > 1L && !is.null(posterior$per_run) &&
      length(posterior$per_run) > 0L) {
    # Return:
    vapply(posterior$per_run, function(r) {
      as.integer(nrow(r$samples) %||% r$saved_idx %||%
                   ((posterior$nSamples %||% 0L) %/% nRuns))
    }, integer(1L))
  } else {
    # Return:
    nrow(posterior$samples)
  }
}


# Internal: samples available to the shortest run
#
# A single burnin is applied to every run, so it must be validated against
# the shortest.
#
# @inheritParams .PerRunSampleCounts
# @return Integer count
.MinPerRunSamples <- function(posterior) min(.PerRunSampleCounts(posterior))


# Internal: indices retained after discarding the first `bi` of `n`
#
# `seq(bi + 1L, n)` counts *downwards* when `bi >= n`, indexing out of bounds
# or -- worse -- silently reversing the draws.
.KeepAfter <- function(bi, n) seq_len(max(0L, n - bi)) + bi


# Internal: get post-burnin samples and trees
#
# @param posterior MkPosterior object
# @param burnin Override burnin (used by AutoBurnin grid search)
# @return List with $samples, $trees, $per_run (filtered)
.PostBurninData <- function(posterior, burnin = NULL) {
  posterior <- .ResolveSamples(posterior)
  nRuns <- posterior$nRuns %||% 1L

  bi <- burnin %||% (posterior$burnin %||% 0L)

  # Differential tree thinning: compute tree-side burnin
  thin     <- posterior$mcmc$thin %||% 1L
  treeThin <- posterior$treeThin %||% thin
  treeEvery <- max(1L, as.integer(treeThin / thin))
  treeBi <- as.integer(floor(bi / treeEvery))

  if (bi == 0L) {
    return(list(
      samples = posterior$samples,
      trees = posterior$trees,
      per_run = posterior$per_run
    ))
  }

  if (nRuns > 1L && !is.null(posterior$per_run)) {
    emptied <- which(.PerRunSampleCounts(posterior) <= bi)
    if (length(emptied) > 0L) {
      cli::cli_warn(
        "Burnin ({bi}) discards every sample of \
         {cli::qty(length(emptied))} run{?s} {.val {emptied}}."
      )
    }

    filteredRuns <- lapply(posterior$per_run, function(r) {
      keep <- .KeepAfter(bi, nrow(r$samples))
      treeKeep <- .KeepAfter(treeBi, length(r$trees))
      list(
        samples = r$samples[keep, , drop = FALSE],
        trees = r$trees[treeKeep],
        acceptance = r$acceptance
      )
    })

    allSamples <- do.call(rbind, lapply(filteredRuns, `[[`, "samples"))
    allTrees <- do.call(c, lapply(filteredRuns, `[[`, "trees"))

    list(
      samples = allSamples,
      trees = allTrees,
      per_run = filteredRuns
    )
  } else {
    keep <- .KeepAfter(bi, nrow(posterior$samples))
    treeKeep <- .KeepAfter(treeBi, length(posterior$trees))
    list(
      samples = posterior$samples[keep, , drop = FALSE],
      trees = posterior$trees[treeKeep],
      per_run = NULL
    )
  }
}


# Helper for print.MkPosterior: format sample count with burnin note
.PostBurninSampleCount <- function(posterior) {
  bi <- posterior$burnin %||% 0L
  pb <- .PostBurninData(posterior)
  n <- nrow(pb$samples)
  if (bi > 0L) {
    sprintf("%d (%d burnin discarded per run)", n, bi)
  } else {
    as.character(n)
  }
}
