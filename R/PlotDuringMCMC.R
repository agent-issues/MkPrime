# Live progress plots for MkPrime MCMC
#
# Phase 6: PlotDuringMCMC trace plots and PNG output.

# --- Stateful ESS history for the ESS-over-time panel ---
.tracePlotEnv <- new.env(parent = emptyenv())
.tracePlotEnv$essHistory <- NULL   # list of named numeric vectors
.tracePlotEnv$iterHistory <- NULL  # integer vector (iter at each snapshot)
.tracePlotEnv$lastIter <- 0L

#' Reset the ESS history used by [MkpTracePlot()]
#' @keywords internal
.ResetEssHistory <- function() {
  .tracePlotEnv$essHistory  <- list()
  .tracePlotEnv$iterHistory <- integer(0)
  .tracePlotEnv$lastIter    <- 0L
}


#' Per-parameter color palette
#'
#' Returns a named character vector of colors for key parameters.
#' @keywords internal
.ParamColors <- function(params) {
  # Fixed palette: visually distinct, legible on white background
  palette <- c(
    log_posterior = "#2166AC",   # blue
    tree_length   = "#B2182B",  # red
    rate_loss     = "#762A83",  # purple
    rate_log_sd   = "#1B7837",  # green
    rate_neo      = "#E08214",  # orange
    p             = "#D6604D",  # salmon
    beta_scale    = "#4393C3"   # light blue
  )
  matched <- palette[intersect(params, names(palette))]
  # Fall back to hcl.colors for any unknown params
  unknown <- setdiff(params, names(palette))
  if (length(unknown)) {
    extra <- grDevices::hcl.colors(length(unknown), "Dark 3")
    names(extra) <- unknown
    matched <- c(matched, extra)
  }
  matched[params]
}


#' Draw live MCMC trace plots
#'
#' A progress callback that draws multi-panel base-R trace plots showing
#' log-posterior and key parameter values across iterations, plus an
#' ESS-over-time panel showing effective sample size trajectories.
#' Intended to be passed as `progressFn` to [MkPrimeMCMC()], either
#' directly or via the convenience string `"default"`.
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
#' log-posterior plus key model parameters from the saved samples, with
#' each parameter in a distinct color. A final panel shows ESS over time
#' for all tracked parameters, using matching colors.
#'
#' When `nRuns > 1`, trace panels use per-run colors (from
#' [grDevices::hcl.colors()] with the `"Set 2"` palette) while the
#' ESS panel uses per-parameter colors with a legend.
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
      c("log_posterior", "tree_length", "rate_loss", "rate_log_sd",
        "rate_neo", "p"),
      cnames
    )
  } else {
    keyParams <- "log_posterior"
  }

  paramColors <- .ParamColors(keyParams)
  hasCoda <- requireNamespace("coda", quietly = TRUE)
  showEss <- hasSamples && hasCoda

  # Reset ESS history if new run detected (iter went backwards or in warmup)
  if (info$inWarmup || info$iter <= .tracePlotEnv$lastIter) {
    .ResetEssHistory()
  }

  # Compute and store ESS snapshot
  if (showEss) {
    combined <- do.call(rbind, info$runSamples)
    essSnap <- vapply(keyParams, function(p) {
      col <- combined[, p]
      s <- sd(col, na.rm = TRUE)
      if (is.na(s) || s == 0) return(NA_real_)
      as.numeric(coda::effectiveSize(coda::mcmc(col)))
    }, numeric(1))
    .tracePlotEnv$essHistory  <- c(.tracePlotEnv$essHistory, list(essSnap))
    .tracePlotEnv$iterHistory <- c(.tracePlotEnv$iterHistory, info$iter)
    .tracePlotEnv$lastIter    <- info$iter
  }

  # Panel layout: trace panels + ESS panel (if applicable)
  nPanels <- length(keyParams) + if (showEss) 1L else 0L
  nCol <- min(3L, nPanels)
  nRow <- ceiling(nPanels / nCol)

  runColors <- if (nRuns > 1L) {
    grDevices::hcl.colors(nRuns, palette = "Set 2")
  } else {
    NULL  # single-run: use per-parameter color
  }

  oldpar <- par(mfrow = c(nRow, nCol),
                mar = c(3, 3.5, 2.5, 1),
                mgp = c(2, 0.6, 0),
                oma = c(0, 0, 1.5, 0))
  on.exit(par(oldpar))

  elapsedStr <- .FormatElapsed(info$elapsed)
  accStr <- format(round(info$recentAcceptance, 3), nsmall = 3)

  for (param in keyParams) {
    if (hasSamples) {
      traceCol <- if (is.null(runColors)) paramColors[param] else runColors
      .PlotTracePanel(param, info$runSamples, nRuns, traceCol,
                      info$warmup, info$nIter)
    } else {
      wCol <- if (is.null(runColors)) paramColors[1] else runColors
      .PlotWarmupPanel(info$currentState, nRuns, wCol,
                       info$iter, info$nIter)
    }
  }

  # ESS-over-time panel
  if (showEss && length(.tracePlotEnv$essHistory) > 0L) {
    .PlotEssPanel(.tracePlotEnv$iterHistory,
                  .tracePlotEnv$essHistory,
                  paramColors)
  }

  # Overall title
  iterStr <- if (is.finite(info$nIter)) {
    sprintf("Iter %d / %d", info$iter, info$nIter)
  } else {
    sprintf("Iter %d", info$iter)
  }
  titleText <- sprintf("%s  |  acc: %s  |  %s", iterStr, accStr, elapsedStr)
  mtext(titleText, outer = TRUE, line = 0, cex = 0.9)

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


#' Plot ESS-over-time panel
#'
#' Draws ESS trajectories for all tracked parameters, with a compact
#' legend mapping colors to parameter names.
#' @param iters Integer vector of iteration numbers (x-axis).
#' @param essHistory List of named numeric vectors (one per snapshot).
#' @param paramColors Named character vector of colors per parameter.
#' @keywords internal
.PlotEssPanel <- function(iters, essHistory, paramColors) {
  params <- names(essHistory[[1]])
  nSnap <- length(essHistory)

  # Build matrix: rows = snapshots, cols = params
  essMat <- do.call(rbind, essHistory)

  # y-axis range (ignore NA)
  allEss <- as.numeric(essMat)
  finiteEss <- allEss[is.finite(allEss)]
  if (length(finiteEss) == 0L) {
    plot.new()
    title(main = "ESS", cex.main = 0.95)
    return(invisible(NULL))
  }
  ylim <- c(0, max(finiteEss) * 1.15)

  plot(iters, essMat[, 1], type = "n", ylim = ylim,
       main = "ESS", xlab = "iter", ylab = "",
       cex.main = 0.95, las = 1)

  for (p in params) {
    vals <- essMat[, p]
    ok <- is.finite(vals)
    if (any(ok)) {
      lines(iters[ok], vals[ok], col = paramColors[p], lwd = 1.5)
    }
  }

  # Compact legend
  legend("topleft", legend = params, col = paramColors[params],
         lwd = 1.5, cex = 0.65, bg = "white", seg.len = 1.2)
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
  nIterJson <- if (is.finite(info$nIter)) info$nIter else "null"
  json <- sprintf(
    paste0('{"iter":%d,"nIter":%s,"warmup":%d,"inWarmup":%s,',
           '"elapsed":%.1f,"recentAcceptance":%.4f}'),
    info$iter, nIterJson, info$warmup,
    if (info$inWarmup) "true" else "false",
    info$elapsed, info$recentAcceptance
  )
  writeLines(json, file)
}


#' Watch a MkPrime log file and draw live trace plots
#'
#' Polls one or more Tracer-compatible log files written by [RunMkPrime()]
#' (when `logFile` is set in [MkPrimeMCMC()]) and redraws trace plots in the
#' active graphics device at each interval. Intended to be run interactively
#' while a `RunMkPrime()` call executes in a background job or another R
#' session.
#'
#' Press **Ctrl+C** (or Escape in RStudio) to stop watching. The last-read
#' sample data is returned invisibly.
#'
#' @param logFiles Character vector of log file paths. For a single run,
#'   supply the value of `mcmc$logFile` directly. For multiple runs, supply
#'   the resolved per-run paths (e.g. `"run_1.log"`, `"run_2.log"`), or use
#'   [MkLogPaths()] to expand a base name.
#' @param interval Seconds between poll–redraw cycles. Default `2`.
#' @param params Character vector of parameter names to plot. Default
#'   `NULL` auto-selects key scalar parameters (`log_posterior`,
#'   `tree_length`, `rate_loss`, `rate_log_sd`, `p`, plus up to three
#'   `kPrime_*` columns). Branch lengths and `log_likelihood` are excluded.
#' @param maxSamples Maximum rows to read per file (tail of file). Set lower
#'   for faster redraws on large logs. Default `2000`.
#' @param warmup Integer. If supplied, draws a vertical dashed line at this
#'   sample index to mark the warmup boundary. Note that log files only
#'   contain post-warmup samples by default, so this is rarely needed.
#'   Default `NULL`.
#'
#' @return The last-read sample data, invisibly. For a single run, a numeric
#'   matrix (as returned by [ReadMkLog()]). For multiple runs, a named list
#'   of such matrices.
#'
#' @details
#' The function waits up to 30 seconds for `logFiles[1]` to appear on disk
#' (useful when called immediately after launching a background MCMC job).
#' If the file still does not exist after 30 seconds, the function aborts.
#'
#' @seealso [ReadMkLog()] to read a log file after the run has finished.
#'   [MkpTracePlot()] for the in-memory callback-based alternative.
#'
#' @export
MkpWatchLog <- function(logFiles,
                         interval   = 2,
                         params     = NULL,
                         maxSamples = 2000L,
                         warmup     = NULL) {
  .MkpWatchLogImpl(logFiles, interval, params, maxSamples, warmup,
                   maxUpdates = Inf)
}


# Internal implementation with maxUpdates for testability.
# @keywords internal
.MkpWatchLogImpl <- function(logFiles, interval, params, maxSamples, warmup,
                              maxUpdates = Inf) {
  logFiles <- as.character(logFiles)
  if (length(logFiles) == 0L)
    cli::cli_abort("{.arg logFiles} must be a non-empty character vector.")

  # Wait up to 30 s for the first file to appear
  waited <- 0L
  while (!file.exists(logFiles[1L])) {
    if (waited == 0L)
      cli::cli_alert_info("Waiting for {.file {logFiles[1L]}} to appear\u2026")
    Sys.sleep(1)
    waited <- waited + 1L
    if (waited >= 30L)
      cli::cli_abort("Log file {.file {logFiles[1L]}} did not appear after 30 s.")
  }

  lastData <- NULL
  nUpdates <- 0L

  tryCatch(
    repeat {
      nUpdates <- nUpdates + 1L

      # Read each file (skip silently if it doesn't exist yet)
      matList <- lapply(logFiles, function(f) {
        if (!file.exists(f)) return(NULL)
        tryCatch({
          m <- ReadMkLog(f)
          n <- nrow(m)
          if (n > maxSamples)
            m <- m[(n - maxSamples + 1L):n, , drop = FALSE]
          m
        }, error = function(e) NULL)  # file may be mid-write
      })
      matList <- Filter(Negate(is.null), matList)

      if (length(matList) > 0L) {
        lastData <- if (length(matList) == 1L) matList[[1L]] else matList
        .WatchLogPlot(matList, params, warmup)
      }

      if (nUpdates >= maxUpdates) break
      Sys.sleep(interval)
    },
    interrupt = function(e) {
      nSamples <- if (!is.null(lastData)) {
        if (is.list(lastData)) sum(vapply(lastData, nrow, integer(1L)))
        else nrow(lastData)
      } else 0L
      cli::cli_alert_info(
        "Stopped watching. {nSamples} sample{?s} read."
      )
    }
  )

  invisible(lastData)
}


#' Draw trace panels from log file data
#'
#' Internal helper called by [MkpWatchLog()]. Takes a list of sample matrices
#' (one per run, as returned by [ReadMkLog()]) and draws multi-panel trace
#' plots.
#'
#' @param matList Named list of numeric matrices (one per run).
#' @param params Character vector of parameters to plot, or `NULL` for
#'   auto-selection.
#' @param warmup Integer sample index at which to draw a warmup boundary
#'   line, or `NULL` to skip.
#' @keywords internal
.WatchLogPlot <- function(matList, params = NULL, warmup = NULL) {
  nRuns <- length(matList)
  # Use column names from the first non-null matrix
  cnames <- colnames(matList[[1L]])

  # Auto-select parameters if not specified
  if (is.null(params)) {
    scalars <- intersect(
      c("log_posterior", "tree_length", "rate_loss", "rate_log_sd", "p"),
      cnames
    )
    kpCols <- grep("^kPrime_", cnames, value = TRUE)
    if (length(kpCols) > 3L) kpCols <- kpCols[seq_len(3L)]
    params <- c(scalars, kpCols)
  }
  params <- intersect(params, cnames)
  if (length(params) == 0L) return(invisible(NULL))

  nPanels <- length(params)
  nCol    <- min(3L, nPanels)
  nRow    <- ceiling(nPanels / nCol)
  colors  <- if (nRuns > 1L) {
    grDevices::hcl.colors(nRuns, palette = "Set 2")
  } else {
    "steelblue"
  }

  oldpar <- par(mfrow = c(nRow, nCol),
                mar   = c(3, 3.5, 2.5, 1),
                mgp   = c(2, 0.6, 0),
                oma   = c(0, 0, 2, 0))
  on.exit(par(oldpar))

  # Compute total samples and last iteration for the title
  totalSamples <- sum(vapply(matList, nrow, integer(1L)))
  lastIter     <- max(vapply(matList, function(m) {
    rn <- suppressWarnings(as.integer(rownames(m)))
    if (all(is.na(rn))) 0L else max(rn, na.rm = TRUE)
  }, integer(1L)))

  for (param in params) {
    # Collect all values for y-axis range
    allVals <- unlist(lapply(matList, function(m) {
      if (param %in% colnames(m)) m[, param] else NULL
    }))
    if (length(allVals) == 0L || all(is.na(allVals))) {
      plot.new(); title(main = param, cex.main = 0.95); next
    }
    ylim <- range(allVals, na.rm = TRUE)
    if (diff(ylim) == 0) ylim <- ylim + c(-1, 1)

    first <- TRUE
    for (run in seq_len(nRuns)) {
      m <- matList[[run]]
      if (!param %in% colnames(m)) next
      vals <- m[, param]
      iters <- suppressWarnings(as.numeric(rownames(m)))
      if (all(is.na(iters))) iters <- seq_along(vals)

      col <- if (nRuns > 1L) colors[run] else colors[1L]
      if (first) {
        plot(iters, vals, type = "l", col = col, ylim = ylim,
             main = param, xlab = "iter", ylab = "",
             cex.main = 0.95, las = 1)
        if (!is.null(warmup) && warmup > min(iters))
          abline(v = warmup, lty = 2, col = "grey60")
        first <- FALSE
      } else {
        lines(iters, vals, col = col)
      }
    }
  }

  # Overall title
  mtext(
    sprintf("iter %s  |  %d sample%s  |  %s",
            format(lastIter, big.mark = ","),
            totalSamples,
            if (totalSamples == 1L) "" else "s",
            format(Sys.time(), "%H:%M:%S")),
    outer = TRUE, line = 0.5, cex = 0.9
  )
  invisible(NULL)
}


#' Expand a base log file name to per-run paths
#'
#' Given a base log file name and a number of runs, returns the resolved
#' per-run file paths that [RunMkPrime()] would create (e.g. `"run.log"` →
#' `c("run_1.log", "run_2.log")`). Useful when constructing the `logFiles`
#' argument to [MkpWatchLog()].
#'
#' @param logFile Character. Base log file name as passed to [MkPrimeMCMC()].
#' @param nRuns Integer. Number of runs.
#' @return Character vector of length `nRuns`.
#' @export
MkLogPaths <- function(logFile, nRuns) {
  MkPrime:::.LogFilePaths(logFile, as.integer(nRuns))
}


#' Construct the cancel-signal file path for a job directory
#'
#' Returns the conventional cancel-file path used by [RunMkPrime()] when
#' `cancelFile` is set. Passing a directory rather than a full path keeps
#' the cancel file alongside the log files in one place.
#'
#' @param jobDir Character. Directory that contains the log files for the run.
#' @return Character(1). Full path to the cancel-signal file
#'   (`<jobDir>/mkp_cancel.signal`).
#' @seealso [MkLogPaths()], [MkPrimeMCMC()]
#' @export
MkCancelPath <- function(jobDir) {
  file.path(jobDir, "mkp_cancel.signal")
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
