# Live progress plots for MkPrime MCMC
#
# Phase 6: PlotDuringMCMC trace plots and PNG output.

#' Draw live MCMC trace plots
#'
#' A progress callback that draws multi-panel base-R trace plots showing
#' log-posterior and key parameter values across iterations. Intended to be
#' passed as `progress_fn` to [MkPrimeMCMC()], either directly or via the
#' convenience string `"default"`.
#'
#' @param info A named list produced by the MCMC loop, containing:
#'   \describe{
#'     \item{iter}{Current iteration number.}
#'     \item{nIter}{Total iterations.}
#'     \item{warmup}{Number of warmup iterations.}
#'     \item{in_warmup}{Logical; `TRUE` if still in warmup.}
#'     \item{nRuns}{Number of independent runs.}
#'     \item{nChains}{Number of chains per run.}
#'     \item{run_samples}{List of per-run sample matrices (may contain
#'       `NULL` if no post-warmup samples yet).}
#'     \item{current_state}{List of per-run cold chain state lists.}
#'     \item{recent_acceptance}{Rolling acceptance rate.}
#'     \item{elapsed}{Elapsed wall-clock seconds.}
#'     \item{param_names}{Column names of the sample matrix.}
#'   }
#'
#' @details
#' During warmup (before post-warmup samples exist), a single-panel
#' log-posterior trace is drawn from `current_state` values accumulated
#' by the callback itself. After warmup, multi-panel traces show
#' log-posterior plus key model parameters from the saved samples.
#'
#' When `nRuns > 1`, each run is drawn in a distinct color using
#' [grDevices::hcl.colors()] with the `"Set 2"` palette.
#'
#' @return Called for its side effect (drawing a plot). Returns `info`
#'   invisibly.
#' @export
mkp_trace_plot <- function(info) {
  nRuns <- info$nRuns

  # Determine which parameters to plot from saved samples
  has_samples <- !is.null(info$run_samples[[1]])
  if (has_samples) {
    cnames <- colnames(info$run_samples[[1]])
    key_params <- intersect(
      c("log_posterior", "tree_length", "rate_loss", "rate_log_sd", "p"),
      cnames
    )
  } else {
    key_params <- "log_posterior"
  }

  nPanels <- length(key_params)
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

  elapsed_str <- .format_elapsed(info$elapsed)
  acc_str <- format(round(info$recent_acceptance, 3), nsmall = 3)

  for (param in key_params) {
    if (has_samples) {
      .plot_trace_panel(param, info$run_samples, nRuns, colors,
                        info$warmup, info$nIter)
    } else {
      # During warmup, show current-state log-posterior
      .plot_warmup_panel(info$current_state, nRuns, colors,
                         info$iter, info$nIter)
    }
  }

  # Overall title
  title_text <- sprintf("Iter %d / %d  |  acc: %s  |  %s",
                        info$iter, info$nIter, acc_str, elapsed_str)
  mtext(title_text, outer = TRUE, line = -1.2, cex = 0.9)

  invisible(info)
}


#' Plot a single trace panel from saved samples
#' @keywords internal
.plot_trace_panel <- function(param, run_samples, nRuns, colors,
                              warmup, nIter) {
  # Collect data from all runs for y-axis range
  all_vals <- unlist(lapply(run_samples, function(s) {
    if (!is.null(s) && param %in% colnames(s)) s[, param]
  }))

  if (length(all_vals) == 0L || all(is.na(all_vals))) {
    plot.new()
    title(main = param)
    return(invisible(NULL))
  }

  ylim <- range(all_vals, na.rm = TRUE)
  if (diff(ylim) == 0) ylim <- ylim + c(-1, 1)

  first <- TRUE
  for (run in seq_len(nRuns)) {
    s <- run_samples[[run]]
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
.plot_warmup_panel <- function(current_state, nRuns, colors, iter, nIter) {
  logp <- vapply(current_state, function(s) {
    s$log_lik + s$log_prior
  }, numeric(1))

  plot(seq_len(nRuns), logp, pch = 19, col = colors,
       main = "log_posterior (current)", xlab = "run", ylab = "",
       xlim = c(0.5, nRuns + 0.5), las = 1, cex.main = 0.95)
}


#' Create a PNG-writing progress callback
#'
#' Returns a progress callback function suitable for `progress_fn` in
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
#'   as `progress_fn`.
#'
#' @details
#' Two files are written on each call:
#' \describe{
#'   \item{`mkp_progress.png`}{Trace plots (same as [mkp_trace_plot()]).}
#'   \item{`mkp_progress.json`}{JSON object with fields `iter`, `nIter`,
#'     `warmup`, `in_warmup`, `elapsed`, and `recent_acceptance`.}
#' }
#' @export
mkp_png_progress <- function(dir, width = 800, height = 600) {
  force(dir)
  force(width)
  force(height)

  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)

  function(info) {
    tmp_png <- file.path(dir, "mkp_progress_tmp.png")
    final_png <- file.path(dir, "mkp_progress.png")

    grDevices::png(tmp_png, width = width, height = height)
    on.exit(grDevices::dev.off(), add = TRUE)
    mkp_trace_plot(info)
    grDevices::dev.off()
    on.exit(NULL)  # dev.off() already called
    file.rename(tmp_png, final_png)

    .write_progress_json(info, file.path(dir, "mkp_progress.json"))
  }
}


#' Write a minimal JSON status file (no jsonlite dependency)
#' @keywords internal
.write_progress_json <- function(info, file) {
  json <- sprintf(
    paste0('{"iter":%d,"nIter":%d,"warmup":%d,"in_warmup":%s,',
           '"elapsed":%.1f,"recent_acceptance":%.4f}'),
    info$iter, info$nIter, info$warmup,
    if (info$in_warmup) "true" else "false",
    info$elapsed, info$recent_acceptance
  )
  writeLines(json, file)
}


#' Format elapsed seconds as human-readable string
#' @keywords internal
.format_elapsed <- function(seconds) {
  if (seconds < 60) {
    sprintf("%.0fs", seconds)
  } else if (seconds < 3600) {
    sprintf("%.1fmin", seconds / 60)
  } else {
    sprintf("%.1fh", seconds / 3600)
  }
}
