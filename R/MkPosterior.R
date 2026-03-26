# MkPosterior result object

#' Create an MkPosterior result object
#'
#' @param samples Matrix of posterior samples (nSaved x nParams).
#' @param trees List of `phylo` objects (sampled trees with varying
#'   branch lengths, fixed topology).
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
  cli::cli_ul(c(
    "Iterations: {x$mcmc$nIter} ({x$warmup} warmup, thinned by {x$mcmc$thin})",
    "Samples: {nrow(x$samples)}",
    "Parameters: {ncol(x$samples)}",
    "Characters: {x$data$nChar} ({sum(x$data$type == 'transformational')} transformational, {sum(x$data$type == 'neomorphic')} neomorphic, {sum(x$data$type == 'known')} known)"
  ))

  cli::cli_h2("Acceptance rates")
  for (nm in names(x$acceptance)) {
    cli::cli_li("{nm}: {format(round(x$acceptance[nm], 3), nsmall = 3)}")
  }

  invisible(x)
}


#' @export
summary.MkPosterior <- function(object, ...) {
  s <- object$samples
  # Exclude branch length columns for summary
  key_cols <- grep("^(log_|tree_|rate_|p$|kPrime_)", colnames(s))
  key <- s[, key_cols, drop = FALSE]

  out <- data.frame(
    parameter = colnames(key),
    mean = colMeans(key),
    median = apply(key, 2, median),
    q025 = apply(key, 2, quantile, 0.025),
    q975 = apply(key, 2, quantile, 0.975),
    row.names = NULL
  )

  # ESS if coda is available
  if (requireNamespace("coda", quietly = TRUE)) {
    out$ESS <- apply(key, 2, function(col) {
      coda::effectiveSize(coda::mcmc(col))
    })
  }

  out
}


#' @export
plot.MkPosterior <- function(x, ...) {
  s <- x$samples
  key_cols <- grep("^(log_posterior|tree_|rate_|p$|kPrime_)", colnames(s))
  nPanels <- length(key_cols)
  nCol <- min(3, nPanels)
  nRow <- ceiling(nPanels / nCol)

  oldpar <- par(mfrow = c(nRow, nCol), mar = c(3, 3, 2, 1))
  on.exit(par(oldpar))

  iters <- seq_len(nrow(s))
  for (col in key_cols) {
    plot(iters, s[, col], type = "l", main = colnames(s)[col],
         xlab = "", ylab = "", col = "steelblue")
  }
}
