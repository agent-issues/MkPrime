# Live progress plots for MkPrime MCMC
#
# Phase 6: PlotDuringMCMC trace plots and PNG output.

#' Draw live MCMC trace plots
#'
#' A progress callback that draws multi-panel base-R trace plots showing
#' log-posterior and key parameter values across iterations. Intended to be
#' passed as `progressFn` to [MkPrimeMCMC()], either directly or via the
#' convenience string `"default"`.
#'
#' @param info A named list produced by the MCMC loop, containing:
#'   \describe{
#'     \item{iter}{Current iteration number.}
#'     \item{nIter}{Total iterations.}
#'     \item{warmup}{Number of warmup iterations.}
#'     \item{inWarmup}{Logical; `TRUE` if still in warmup.}
#'     \item{nRuns}{Number of independent runs.}
#'     \item{nChains}{Number of chains per run.}
#'     \item{runSamples}{List of per-run sample matrices (may contain
#'       `NULL` if no post-warmup samples yet).}
#'     \item{currentState}{List of per-run cold chain state lists.}
#'     \item{recentAcceptance}{Rolling acceptance rate.}
#'     \item{elapsed}{Elapsed wall-clock seconds.}
#'     \item{paramNames}{Column names of the sample matrix.}
#'   }
#'
#' @details
#' During warmup (before post-warmup samples exist), a single-panel
#' log-posterior trace is drawn from `currentState` values accumulated
#' by the callback itself. After warmup, multi-panel traces show
#' log-posterior plus key model parameters from the saved samples.
#'
#' When `nRuns > 1`, each run is drawn in a distinct color using
#' [grDevices::hcl.colors()] with the `"Set 2"` palette.
#'
#' @return Called for its side effect (drawing a plot). Returns `info`
#'   invisibly.
#' @export
MkpTracePlot <- function(info) {
  nRuns <- info$nRuns

  # Determine which parameters to plot from saved samples
  hasSamples <- !is.null(info$runSamples[[1]])
  if (hasSamples) {
    cnames <- colnames(info$runSamples[[1]])
    keyParams <- intersect(
      c("log_posterior", "tree_length", "rate_loss", "rate_log_sd", "p"),
      cnames
    )
  } else {
    keyParams <- "log_posterior"
  }

  nPanels <- length(keyParams)
  nCol <- min(3L, nPanels)
  nRow <- ceiling(nPanels / nCol)

  colors <- if (nRuns > 1L) {
    grDevices::hcl.colors(nRuns, palette = "Set 2")
  } else {
    "steelblue"
  }

  oldpar <- par(mfrow = c(nRow, nCol),
                mar = c(3, 3.5, 2.5, 1),
                mgp = c(2, 0.6, 0))
  on.exit(par(oldpar))

  elapsedStr <- .FormatElapsed(info$elapsed)
  accStr <- format(round(info$recentAcceptance, 3), nsmall = 3)

  for (param in keyParams) {
    if (hasSamples) {
      .PlotTracePanel(param, info$runSamples, nRuns, colors,
                      info$warmup, info$nIter)
    } else {
      # During warmup, show current-state log-posterior
      .PlotWarmupPanel(info$currentState, nRuns, colors,
                       info$iter, info$nIter)
    }
  }

  # Overall title
  titleText <- sprintf("Iter %d / %d  |  acc: %s  |  %s",
                       info$iter, info$nIter, accStr, elapsedStr)
  mtext(titleText, outer = TRUE, line = -1.2, cex = 0.9)

  invisible(info)
}


#' Plot a single trace panel from saved samples
#' @keywords internal
.PlotTracePanel <- function(param, runSamples, nRuns, colors,
                             warmup, nIter) {
  # Collect data from all runs for y-axis range
  allVals <- unlist(lapply(runSamples, function(s) {
    if (!is.null(s) && param %in% colnames(s)) s[, param]
  }))

  if (length(allVals) == 0L || all(is.na(allVals))) {
    plot.new()
    title(main = param)
    return(invisible(NULL))
  }

  ylim <- range(allVals, na.rm = TRUE)
  if (diff(ylim) == 0) ylim <- ylim + c(-1, 1)

  first <- TRUE
  for (run in seq_len(nRuns)) {
    s <- runSamples[[run]]
    if (is.null(s) || !param %in% colnames(s)) next
    vals <- s[, param]
    n <- length(vals)
    if (n == 0) next
    iters <- seq_len(n)
    col <- if (nRuns > 1L) colors[run] else colors[1]

    if (first) {
      plot(iters, vals, type = "l", col = col, ylim = ylim,
           main = param, xlab = "sample", ylab = "",
           cex.main = 0.95, las = 1)
      first <- FALSE
    } else {
      lines(iters, vals, col = col)
    }
  }
}


#' Plot log-posterior during warmup (no saved samples yet)
#' @keywords internal
.PlotWarmupPanel <- function(currentState, nRuns, colors, iter, nIter) {
  logp <- vapply(currentState, function(s) {
    s$log_lik + s$log_prior
  }, numeric(1))

  plot(seq_len(nRuns), logp, pch = 19, col = colors,
       main = "log_posterior (current)", xlab = "run", ylab = "",
       xlim = c(0.5, nRuns + 0.5), las = 1, cex.main = 0.95)
}


#' Create a PNG-writing progress callback
#'
#' Returns a progress callback function suitable for `progressFn` in
#' [MkPrimeMCMC()]. Each invocation writes a trace-plot PNG and a JSON
#' status file to the specified directory. Files are written atomically
#' (write to temp, then rename) so that external consumers (e.g., a
#' Shiny app) never read a half-written file.
#'
#' @param dir Directory to write progress files to. Created if it does
#'   not exist.
#' @param width,height PNG dimensions in pixels.
#'
#' @return A function with signature `function(info)` suitable for use
#'   as `progressFn`.
#'
#' @details
#' Two files are written on each call:
#' \describe{
#'   \item{`mkp_progress.png`}{Trace plots (same as [MkpTracePlot()]).}
#'   \item{`mkp_progress.json`}{JSON object with fields `iter`, `nIter`,
#'     `warmup`, `inWarmup`, `elapsed`, and `recentAcceptance`.}
#' }
#' @export
MkpPngProgress <- function(dir, width = 800, height = 600) {
  force(dir)
  force(width)
  force(height)

  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)

  function(info) {
    tmpPng <- file.path(dir, "mkp_progress_tmp.png")
    finalPng <- file.path(dir, "mkp_progress.png")

    grDevices::png(tmpPng, width = width, height = height)
    on.exit(grDevices::dev.off(), add = TRUE)
    MkpTracePlot(info)
    grDevices::dev.off()
    on.exit(NULL)  # dev.off() already called
    file.rename(tmpPng, finalPng)

    .WriteProgressJson(info, file.path(dir, "mkp_progress.json"))
  }
}


#' Write a minimal JSON status file (no jsonlite dependency)
#' @keywords internal
.WriteProgressJson <- function(info, file) {
  json <- sprintf(
    paste0('{"iter":%d,"nIter":%d,"warmup":%d,"inWarmup":%s,',
           '"elapsed":%.1f,"recentAcceptance":%.4f}'),
    info$iter, info$nIter, info$warmup,
    if (info$inWarmup) "true" else "false",
    info$elapsed, info$recentAcceptance
  )
  writeLines(json, file)
}


#' Format elapsed seconds as human-readable string
#' @keywords internal
.FormatElapsed <- function(seconds) {
  if (seconds < 60) {
    sprintf("%.0fs", seconds)
  } else if (seconds < 3600) {
    sprintf("%.1fmin", seconds / 60)
  } else {
    sprintf("%.1fh", seconds / 3600)
  }
}
