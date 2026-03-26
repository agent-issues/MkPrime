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

  info <- c(info,
    "Total samples: {nrow(x$samples)}",
    "Parameters: {ncol(x$samples)}",
    "Characters: {x$data$nChar} ({sum(x$data$type == 'transformational')} transformational, {sum(x$data$type == 'neomorphic')} neomorphic, {sum(x$data$type == 'known')} known)"
  )

  if (!is.null(x$stop_reason)) {
    info <- c(info, "Stopped: {x$stop_reason} (iter {x$actual_iter})")
  }

  cli::cli_ul(info)

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
    diag <- tryCatch(convergence_diagnostics(x), error = function(e) NULL)
    if (!is.null(diag)) {
      cli::cli_h2("Convergence")
      cli::cli_li("Min ESS: {format(round(diag$min_ess, 1), nsmall = 1)}")
      if (!is.null(diag$psrf) && !is.na(diag$max_psrf)) {
        cli::cli_li("Max PSRF: {format(round(diag$max_psrf, 3), nsmall = 3)}")
        if (is.finite(diag$max_psrf) && diag$max_psrf > 1.05) {
          cli::cli_alert_warning(
            "PSRF > 1.05 suggests chains may not have converged. Consider running longer."
          )
        }
      }
      if (is.finite(diag$min_ess) && diag$min_ess < 200) {
        cli::cli_alert_warning(
          "ESS < 200 for some parameters. Consider running longer."
        )
      }
    }
  }

  invisible(x)
}


#' @export
summary.MkPosterior <- function(object, ...) {
  s <- object$samples
  key_cols <- .key_param_cols(s)
  key <- s[, key_cols, drop = FALSE]

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
      if (sd(col, na.rm = TRUE) == 0) return(NA_real_)
      coda::effectiveSize(coda::mcmc(col))
    })
  }

  nRuns <- object$nRuns %||% 1L
  if (nRuns >= 2L && !is.null(object$per_run)) {
    diag <- tryCatch(convergence_diagnostics(object), error = function(e) NULL)
    if (!is.null(diag) && !is.null(diag$psrf)) {
      out$PSRF <- diag$psrf[out$parameter]
    }
  }

  out
}


#' @export
plot.MkPosterior <- function(x, ...) {
  s <- x$samples
  key_cols <- .key_param_cols(s)
  nPanels <- length(key_cols)
  nCol <- min(3, nPanels)
  nRow <- ceiling(nPanels / nCol)

  nRuns <- x$nRuns %||% 1L

  oldpar <- par(mfrow = c(nRow, nCol), mar = c(3, 3, 2, 1))
  on.exit(par(oldpar))

  if (nRuns > 1L && !is.null(x$per_run)) {
    colors <- grDevices::hcl.colors(nRuns, palette = "Set 2")

    for (col_idx in key_cols) {
      col_name <- colnames(s)[col_idx]
      ylim <- range(s[, col_idx], na.rm = TRUE)

      first <- TRUE
      for (run in seq_len(nRuns)) {
        run_data <- x$per_run[[run]]$samples[, col_idx]
        iters <- seq_along(run_data)
        if (first) {
          plot(iters, run_data, type = "l", main = col_name,
               xlab = "", ylab = "", col = colors[run], ylim = ylim)
          first <- FALSE
        } else {
          lines(iters, run_data, col = colors[run])
        }
      }
    }
  } else {
    iters <- seq_len(nrow(s))
    for (col_idx in key_cols) {
      plot(iters, s[, col_idx], type = "l", main = colnames(s)[col_idx],
           xlab = "", ylab = "", col = "steelblue")
    }
  }
}
