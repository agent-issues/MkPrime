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
                        warmup, tuning) {
  structure(
    list(
      samples = samples,
      trees = trees,
      acceptance = acceptance,
      model = model,
      data = data,
      mcmc = mcmc,
      warmup = warmup,
      tuning = tuning
    ),
    class = "MkPosterior"
  )
}


#' @export
print.MkPosterior <- function(x, ...) {
  cli::cli_h1("MkPrime Posterior")

  nRuns <- x$nRuns %||% 1L
  nChains <- x$mcmc$nChains

  info <- c(
    "Iterations: {x$mcmc$nIter} ({x$warmup} warmup, thinned by {x$mcmc$thin})"
  )

  if (nRuns > 1L) {
    info <- c(info, "Runs: {nRuns} ({nrow(x$per_run[[1]]$samples)} samples each)")
  }
  if (nChains > 1L) {
    info <- c(info, "Chains per run: {nChains} (cold + {nChains - 1L} heated)")
  }

  if (!is.null(x$logFile)) {
    info <- c(info,
      "Total samples: {x$nSamples} (streamed to disk)",
      "Parameters: {ncol(x$samples)}",
      "Characters: {x$data$nChar} ({sum(x$data$type == 'transformational')} transformational, {sum(x$data$type == 'neomorphic')} neomorphic, {sum(x$data$type == 'known')} known)"
    )
  } else {
    info <- c(info,
      "Total samples: {(.PostBurninSampleCount(x))}",
      "Parameters: {ncol(x$samples)}",
      "Characters: {x$data$nChar} ({sum(x$data$type == 'transformational')} transformational, {sum(x$data$type == 'neomorphic')} neomorphic, {sum(x$data$type == 'known')} known)"
    )
  }

  if (!is.null(x$stop_reason)) {
    info <- c(info, "Stopped: {x$stop_reason} (iter {x$actual_iter})")
  }

  cli::cli_ul(info)

  if (!is.null(x$logFile)) {
    cli::cli_alert_info(c(
      "Streaming mode: samples are on disk, not in memory.",
      "i" = "Load with: {.code result$samples <- ReadMkLog(result$logFile)}"
    ))
  }

  cli::cli_h2("Acceptance rates (cold chain)")
  for (nm in names(x$acceptance)) {
    cli::cli_li("{nm}: {format(round(x$acceptance[nm], 3), nsmall = 3)}")
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

  if (requireNamespace("coda", quietly = TRUE)) {
    out$ESS <- apply(key, 2, function(col) {
      s <- sd(col, na.rm = TRUE); if (is.na(s) || s == 0) return(NA_real_)
      coda::effectiveSize(coda::mcmc(col))
    })
  }

  nRuns <- object$nRuns %||% 1L
  if (nRuns >= 2L && !is.null(object$per_run)) {
    # trees = FALSE: topology ESS is expensive; omit from summary().
    diag <- tryCatch(ConvergenceDiagnostics(object, trees = FALSE),
                     error = function(e) NULL)
    if (!is.null(diag) && !is.null(diag$psrf)) {
      out$PSRF <- diag$psrf[out$parameter]
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

  if (nRuns > 1L && !is.null(pb$per_run)) {
    colors <- grDevices::hcl.colors(nRuns, palette = "Set 2")

    for (colIdx in keyCols) {
      colName <- colnames(s)[colIdx]
      vals <- s[, colIdx]
      ylim <- range(vals, na.rm = TRUE)
      useLog <- colName %in% .LogScaleParams && all(vals > 0, na.rm = TRUE)
      logArg <- if (useLog) "y" else ""

      first <- TRUE
      for (run in seq_len(nRuns)) {
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
