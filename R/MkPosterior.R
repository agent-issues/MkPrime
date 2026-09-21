# MkPosterior result object

#' @importFrom rlang `%||%`
NULL

#' Create an MkPosterior result object
#'
#' @param samples Matrix of posterior samples (nSaved x nParams).
#' @param trees List of `phylo` objects (sampled trees).
#' @param acceptance Named numeric vector of acceptance rates per move.
#' @param model The `MkPrimeModel` used.
#' @param data The `MkPrimeData` used.
#' @param mcmc The `MkPrimeMCMC` configuration.
#' @param warmup Number of warmup iterations.
#' @param tuning Final tuning parameters after adaptation.
#'
#' @return An S3 object of class `MkPosterior`.
#' @keywords internal
MkPosterior <- function(samples, trees, acceptance, model, data, mcmc,
                        warmup, tuning, warmup_trace = NULL) {
  structure(
    list(
      samples = samples,
      trees = trees,
      acceptance = acceptance,
      model = model,
      data = data,
      mcmc = mcmc,
      warmup = warmup,
      tuning = tuning,
      warmup_trace = warmup_trace
    ),
    class = "MkPosterior"
  )
}


#' @export
print.MkPosterior <- function(x, ...) {
  cli::cli_h1("MkPrime Posterior")

  nRuns   <- x$nRuns %||% 1L
  nChains <- x$mcmc$nChains %||% 1L
  thin    <- x$mcmc$thin %||% 1L

  treeThin <- x$treeThin %||% thin
  if (!is.null(x$mcmc)) {
    info <- c(
      "Iterations: {x$mcmc$nIter} ({x$warmup} warmup, thinned by {thin})"
    )
    if (treeThin > thin) {
      info <- c(info,
        "Tree samples: {length(x$trees)} (thinned by {treeThin})")
    }
  } else {
    info <- character(0)
  }

  if (nRuns > 1L) {
    perRunN <- .PerRunSampleCounts(x)
    info <- c(info, if (length(unique(perRunN)) == 1L) {
      "Runs: {nRuns} ({perRunN[[1]]} samples each)"
    } else {
      "Runs: {nRuns} ({paste(perRunN, collapse = ', ')} samples)"
    })
  }

  # PAR-009: surface run shrinkage when some runs were dropped
  drops <- x$dropped_runs
  if (!is.null(drops) && nrow(drops) > 0L) {
    reqN    <- x$requested_nRuns %||% (nRuns + nrow(drops))
    nDrop   <- nrow(drops)
    # Build per-drop description: "run 5 (unlaunched)", "run 6 (killed at 1.2s)"
    dropDesc <- vapply(seq_len(nDrop), function(i) {
      r <- drops[i, ]
      detail <- switch(r$reason,
        unlaunched = "unlaunched",
        killed     = paste0("killed at ", round(r$wait_s, 1), "s"),
        errored    = "errored",
        r$reason
      )
      paste0("run ", r$run, " (", detail, ")")
    }, character(1L))
    cli::cli_alert_warning(
      "{nDrop} of {reqN} run{?s} did not return: {paste(dropDesc, collapse = ', ')}."
    )
  }
  if (nChains > 1L) {
    info <- c(info, "Chains per run: {nChains} (cold + {nChains - 1L} heated)")
  }

  if (!is.null(x$logFile)) {
    info <- c(info,
      "Total samples: {x$nSamples} (streamed to disk)",
      "Parameters: {ncol(x$samples)}"
    )
  } else {
    info <- c(info,
      "Total samples: {(.PostBurninSampleCount(x))}",
      "Parameters: {ncol(x$samples)}"
    )
  }

  if (!is.null(x$data)) {
    info <- c(info,
      "Characters: {x$data$nChar} ({sum(x$data$type == 'transformational')} transformational, {sum(x$data$type == 'neomorphic')} neomorphic, {sum(x$data$type == 'known')} known)"
    )
  }

  if (!is.null(x$stop_reason)) {
    reason <- x$stop_reason
    if (!is.null(x$actual_iter)) {
      info <- c(info, "Stopped: {reason} (iter {x$actual_iter})")
    } else {
      info <- c(info, "Stopped: {reason}")
    }
  }

  cli::cli_ul(info)

  if (!is.null(x$logFile) && nrow(x$samples) == 0L) {
    cli::cli_alert_info(c(
      "Streaming mode: samples are on disk, not in memory.",
      "i" = "Load with: {.code result$samples <- ReadMkLog(result$logFile)}"
    ))
  }

  if (nRuns < 1L) {
    cli::cli_alert_warning(
      "No runs completed: this posterior carries no samples or diagnostics."
    )
  }

  if (length(x$acceptance) > 0L) {
    cli::cli_h2("Acceptance rates (cold chain)")
    for (nm in names(x$acceptance)) {
      cli::cli_li("{nm}: {format(round(x$acceptance[nm], 3), nsmall = 3)}")
    }
  }

  if (!is.null(x$swap_rates)) {
    cli::cli_h2("Swap rates")
    for (i in seq_along(x$swap_rates)) {
      cli::cli_li("chains {i}/{i+1}: {format(round(x$swap_rates[i], 3), nsmall = 3)}")
    }
  }

  if (nRuns >= 2L) {
    # trees = FALSE: topology ESS is expensive; call ConvergenceDiagnostics()
    # explicitly post-run if tree ESS is needed.
    diag <- tryCatch(ConvergenceDiagnostics(x, trees = FALSE),
                     error = function(e) NULL)
    if (!is.null(diag)) {
      print(diag)
    }
  }

  invisible(x)
}


#' @export
summary.MkPosterior <- function(object, ...) {
  pb <- .PostBurninData(object)
  s <- pb$samples
  if (nrow(s) == 0L) {
    cli::cli_abort(c(
      "No posterior samples available to summarise.",
      "i" = if (!is.null(object$logFile))
        "The log file{?s} {.file {object$logFile}} may be empty (header only)."
      else
        "The run may have been interrupted before sampling began."
    ))
  }
  # Use scalar params only (not individual kPrime_ or branch lengths)
  keyCols <- .PlotParamCols(s)
  key <- s[, keyCols, drop = FALSE]

  out <- data.frame(
    parameter = colnames(key),
    mean = colMeans(key),
    median = apply(key, 2, median),
    q025 = apply(key, 2, quantile, 0.025),
    q975 = apply(key, 2, quantile, 0.975),
    row.names = NULL
  )

  out$ESS <- .EssMatrix(key)[colnames(key)]

  nRuns <- object$nRuns %||% 1L
  if (nRuns >= 2L && !is.null(object$per_run)) {
    # trees = FALSE: topology ESS is expensive; omit from summary().
    diag <- tryCatch(ConvergenceDiagnostics(object, trees = FALSE),
                     error = function(e) NULL)
    if (!is.null(diag) && !is.null(diag$rhat)) {
      out$Rhat <- diag$rhat[out$parameter]
    }
  }

  out
}


#' @export
plot.MkPosterior <- function(x, ...) {
  pb <- .PostBurninData(x)
  s <- pb$samples
  if (nrow(s) == 0L) {
    cli::cli_abort(c(
      "No posterior samples available to plot.",
      "i" = if (!is.null(x$logFile))
        "The log file{?s} {.file {x$logFile}} may be empty (header only)."
      else
        "The run may have been interrupted before sampling began."
    ))
  }
  keyCols <- .PlotParamCols(s)
  nPanels <- length(keyCols)
  nCol <- min(3, nPanels)
  nRow <- ceiling(nPanels / nCol)

  nRuns <- x$nRuns %||% 1L

  oldpar <- par(mfrow = c(nRow, nCol), mar = c(3, 3, 2, 1))
  on.exit(par(oldpar))

  hasPerRunSamples <- nRuns > 1L && !is.null(pb$per_run) &&
    !is.null(pb$per_run[[1]]$samples)

  if (hasPerRunSamples) {
    colors <- grDevices::hcl.colors(nRuns, palette = "Set 2")

    for (colIdx in keyCols) {
      colName <- colnames(s)[colIdx]
      vals <- s[, colIdx]
      ylim <- range(vals, na.rm = TRUE)
      useLog <- colName %in% .LogScaleParams && all(vals > 0, na.rm = TRUE)
      logArg <- if (useLog) "y" else ""

      first <- TRUE
      for (run in seq_len(nRuns)) {
        if (nrow(pb$per_run[[run]]$samples) == 0L) next
        runData <- pb$per_run[[run]]$samples[, colIdx]
        iters <- seq_along(runData)
        if (first) {
          plot(iters, runData, type = "l", main = colName,
               xlab = "", ylab = "", col = colors[run], ylim = ylim,
               log = logArg)
          first <- FALSE
        } else {
          lines(iters, runData, col = colors[run])
        }
      }
    }
  } else {
    iters <- seq_len(nrow(s))
    for (colIdx in keyCols) {
      colName <- colnames(s)[colIdx]
      vals <- s[, colIdx]
      useLog <- colName %in% .LogScaleParams && all(vals > 0, na.rm = TRUE)
      logArg <- if (useLog) "y" else ""
      plot(iters, vals, type = "l", main = colName,
           xlab = "", ylab = "", col = "steelblue", log = logArg)
    }
  }
}
