# Burnin selection and filtering for MkPosterior objects
#
# M-066: The posterior object stores a per-run sample count `burnin`.
# All accessors (summary, print, plot, ConvergenceDiagnostics) respect
# this value. AutoBurnin() searches for the optimal burnin that
# maximizes ESS while maintaining low R-hat.

#' Set burnin for a posterior object
#'
#' @param posterior An `MkPosterior` object.
#' @param burnin Number of per-run samples to discard as burnin.
#'   If < 1, interpreted as a fraction of per-run sample count.
#'   If >= 1, interpreted as an absolute count.
#'
#' @return A new `MkPosterior` with the updated burnin.
#' @export
SetBurnin <- function(posterior, burnin) {
  if (!inherits(posterior, "MkPosterior")) {
    cli::cli_abort("{.arg posterior} must be an {.cls MkPosterior} object.")
  }

  nRuns <- posterior$nRuns %||% 1L
  if (nRuns > 1L && !is.null(posterior$per_run)) {
    nPerRun <- nrow(posterior$per_run[[1]]$samples) %||%
               posterior$per_run[[1]]$saved_idx %||%
               (posterior$nSamples %/% nRuns)
  } else {
    nPerRun <- nrow(posterior$samples)
  }

  if (burnin < 1) {
    burnin <- floor(burnin * nPerRun)
  }
  burnin <- as.integer(burnin)

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
#' @param posterior An `MkPosterior` object (typically with `nRuns >= 2`).
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

  if (hasRhat) {
    nPerRun <- nrow(posterior$per_run[[1]]$samples) %||%
               posterior$per_run[[1]]$saved_idx %||%
               (posterior$nSamples %/% nRuns)
  } else {
    nPerRun <- nrow(posterior$samples)
  }

  results <- data.frame(
    fraction = fractions,
    burnin = as.integer(floor(fractions * nPerRun)),
    minEss = NA_real_,
    maxRhat = NA_real_,
    stringsAsFactors = FALSE
  )

  # Ensure burnin values leave at least 10 samples
  results <- results[results$burnin < (nPerRun - 10L), , drop = FALSE]

  for (i in seq_len(nrow(results))) {
    bi <- results$burnin[i]
    pb <- .PostBurninData(posterior, bi)

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
    if (nrow(converged) > 0L) {
      # Among converged, pick smallest burnin (preserves most samples / ESS)
      best <- converged[which.min(converged$burnin), ]
    } else {
      # No burnin achieves target R-hat; pick the one with lowest max(R-hat)
      best <- results[which.min(results$maxRhat), ]
      cli::cli_warn(c(
        "No burnin fraction achieves max(Rhat) <= {rhatThreshold}.",
        "i" = "Selected burnin = {best$burnin} (max Rhat = {round(best$maxRhat, 3)}).",
        "i" = "Consider running the chain longer."
      ))
    }
  } else {
    # Single run: maximize min(ESS)
    best <- results[which.max(results$minEss), ]
  }

  cli::cli_inform(c(
    "v" = "Auto burnin: {best$burnin} samples ({round(best$fraction * 100)}% of {nPerRun})",
    "i" = "min(ESS) = {round(best$minEss, 1)}{if (hasRhat) paste0(', max(Rhat) = ', round(best$maxRhat, 3)) else ''}"
  ))

  posterior$burnin <- best$burnin
  posterior$auto_burnin_results <- results
  posterior
}


# Internal: get post-burnin samples and trees
#
# @param posterior MkPosterior object
# @param burnin Override burnin (used by AutoBurnin grid search)
# @return List with $samples, $trees, $per_run (filtered)
.PostBurninData <- function(posterior, burnin = NULL) {
  # Auto-load samples from log file when they haven't been loaded yet
  if (nrow(posterior$samples) == 0L && !is.null(posterior$logFile)) {
    posterior$samples <- ReadMkLog(posterior$logFile)
  }

  # Auto-load per-run samples from individual log files when NULL (streaming)
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
    filteredRuns <- lapply(posterior$per_run, function(r) {
      nSamp <- nrow(r$samples)
      keep <- seq(bi + 1L, nSamp)
      nTree <- length(r$trees)
      treeKeep <- if (treeBi < nTree) seq(treeBi + 1L, nTree) else integer(0)
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
    nSamp <- nrow(posterior$samples)
    keep <- seq(bi + 1L, nSamp)
    nTree <- length(posterior$trees)
    treeKeep <- if (treeBi < nTree) seq(treeBi + 1L, nTree) else integer(0)
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
