# Convergence monitoring for MkPrime MCMC
#
# Phase 5: ESS, R-hat (rank-normalized), convergence diagnostics.
# Replaced coda-based PSRF with native rank-normalized R-hat
# (Vehtari et al. 2021) and native ESS (Geyer 1992).

#' Compute convergence diagnostics for an MkPosterior
#'
#' Calculates effective sample size (ESS) per parameter and, when
#' `nRuns >= 2`, the rank-normalized split-chain R-hat
#' (Vehtari et al. 2021) across runs.  R-hat supersedes the classical
#' PSRF (Gelman-Rubin): it is robust to skewed and multimodal
#' posteriors, and also checks tail convergence.
#'
#' Optionally computes topology ESS using the Robinson-Foulds distance
#' via **TreeDist**.
#'
#' @param posterior An `MkPosterior` object.
#' @param trees Logical. If `TRUE` and **TreeDist** is installed, compute
#'   topology ESS (median pseudo-ESS) from the sampled trees.
#'   Defaults to `FALSE`.
#'   Each run is subsampled to at most 1,000 trees.
#' @param frechetESS Logical. If `TRUE`, also compute the Fréchet
#'   correlation ESS (Magee et al. 2021).  This requires the full
#'   n x n pairwise distance matrix rather than the partial cross-distance
#'   matrix used for median pseudo-ESS, so is substantially slower.
#'   Implies `trees = TRUE`.
#' @return An object of class `MkpDiagnostics`, a list with components:
#'   - `ess`: Named numeric vector of ESS per parameter (including kPrime).
#'   - `minEss`: Scalar minimum ESS across **scalar** parameters.
#'     Individual kPrime values are discrete nuisance parameters that are
#'     marginalized over, so they are excluded from the summary minimum
#'     to avoid blocking convergence (see M-098).
#'   - `rhat`: Named numeric vector of R-hat values per parameter
#'     (only if `nRuns >= 2`; `NULL` otherwise).
#'   - `maxRhat`: Scalar maximum R-hat across scalar parameters (or `NA`
#'     if single run). kPrime excluded for the same reason as `minEss`.
#'   - `treeEss`: Named numeric vector with `frechetCorrelationESS` and
#'     `medianPseudoESS` (or `NULL` if not computed).
#'     `frechetCorrelationESS` is `NA` unless `frechetESS = TRUE`.
#' @export
ConvergenceDiagnostics <- function(posterior, trees = FALSE,
                                   frechetESS = FALSE) {
  if (!inherits(posterior, "MkPosterior")) {
    cli::cli_abort("{.arg posterior} must be an {.cls MkPosterior} object.")
  }

  pb <- .PostBurninData(posterior)
  keyCols <- .KeyParamCols(pb$samples)
  nRuns <- posterior$nRuns %||% 1L

  # --- ESS (combined samples) ---
  ess <- .ComputeEss(pb$samples[, keyCols, drop = FALSE])

  # kPrime are discrete nuisance parameters — exclude from summary min/max
  # (M-098). Individual kPrime ESS/R-hat remain in the output for display.
  isConvParam <- !grepl("^kPrime_", names(ess)) & names(ess) != "log_likelihood"

  # --- R-hat across runs ---
  rhat <- NULL
  maxRhat <- NA_real_

  if (nRuns >= 2L && !is.null(posterior$per_run)) {
    rhat <- .ComputeRhat(pb$per_run, keyCols)
    maxRhat <- max(rhat[isConvParam[names(rhat) %in% names(ess)]],
                   na.rm = TRUE)
  }

  # --- Tree ESS ---
  if (isTRUE(frechetESS)) trees <- TRUE
  treeEss <- .ComputeTreeEss(pb, trees, frechet = isTRUE(frechetESS))

  structure(
    list(
      ess = ess,
      minEss = min(ess[isConvParam], na.rm = TRUE),
      rhat = rhat,
      maxRhat = maxRhat,
      treeEss = treeEss,
      nRuns = nRuns,
      nSamples = nrow(pb$samples),
      burnin = posterior$burnin %||% 0L
    ),
    class = "MkpDiagnostics"
  )
}


#' Print convergence diagnostics
#'
#' Displays a colour-coded table of ESS and R-hat values for each parameter.
#' ESS < 100 is shown in red, 100--199 in yellow, and \eqn{\ge}{>=} 200 in
#' plain text. Individual k' parameters are summarised as a single
#' min--median--max row.
#'
#' @param x An `MkpDiagnostics` object from [ConvergenceDiagnostics()].
#' @param ... Ignored.
#' @return `x` invisibly.
#' @export
print.MkpDiagnostics <- function(x, ...) {
  hasRhat <- !is.null(x$rhat)

  nms <- names(x$ess)
  # log_likelihood is redundant with log_posterior in display
  scalarNms <- nms[!grepl("^kPrime_", nms) & nms != "log_likelihood"]
  kPrimeNms <- nms[grepl("^kPrime_", nms)]

  cli::cli_rule(
    left = sprintf(
      "Convergence diagnostics  (%d run%s, %d samples)",
      x$nRuns, if (x$nRuns == 1L) "" else "s", x$nSamples
    )
  )

  # Column header
  if (hasRhat) {
    cat(sprintf("  %-20s  %6s  %7s\n", "Parameter", "ESS", "Rhat"))
  } else {
    cat(sprintf("  %-20s  %6s\n", "Parameter", "ESS"))
  }
  cat(sprintf("  %s\n", strrep("-", if (hasRhat) 38L else 28L)))

  # Scalar parameter rows
  for (nm in scalarNms) {
    essStr <- .FmtEss(x$ess[nm])
    if (hasRhat && nm %in% names(x$rhat)) {
      cat(sprintf("  %-20s  %s  %s\n", nm, essStr, .FmtRhat(x$rhat[nm])))
    } else {
      cat(sprintf("  %-20s  %s\n", nm, essStr))
    }
  }

  # kPrime summary row (min / median / max) — compact numbers, no fixed width
  if (length(kPrimeNms) > 0L) {
    kpEss <- x$ess[kPrimeNms]
    kpMin <- min(kpEss, na.rm = TRUE)
    kpMed <- median(kpEss, na.rm = TRUE)
    kpMax <- max(kpEss, na.rm = TRUE)
    essRange <- paste0(
      .ColorEss(kpMin), " / ",
      .ColorEss(kpMed), " / ",
      .ColorEss(kpMax),
      "  (min/med/max)"
    )
    label <- sprintf("kPrime (%d)", length(kPrimeNms))
    if (hasRhat && any(kPrimeNms %in% names(x$rhat))) {
      kpRhat <- x$rhat[kPrimeNms[kPrimeNms %in% names(x$rhat)]]
      rhatRange <- sprintf("%.3f\u2013%.3f",
        min(kpRhat, na.rm = TRUE), max(kpRhat, na.rm = TRUE))
      cat(sprintf("  %-20s  %s  Rhat %s\n", label, essRange, rhatRange))
    } else {
      cat(sprintf("  %-20s  %s\n", label, essRange))
    }
  }

  # Topology ESS row (always shown; NA when not computed)
  cat(sprintf("  %s\n", strrep("-", if (hasRhat) 38L else 28L)))
  if (!is.null(x$treeEss)) {
    frech <- x$treeEss[["frechetCorrelationESS"]]
    mdps  <- x$treeEss[["medianPseudoESS"]]
    cat(sprintf("  %-20s  %s\n", "topology (Fréchet)",  .FmtEss(frech)))
    cat(sprintf("  %-20s  %s\n", "topology (med.pseudo)", .FmtEss(mdps)))
  } else {
    cat(sprintf("  %-20s  %s\n", "topology (Fréchet)",   formatC("NA", width = 6)))
    cat(sprintf("  %-20s  %s\n", "topology (med.pseudo)", formatC("NA", width = 6)))
  }

  # Legend
  cat("\n")
  cat("  ESS: ",
      cli::col_green("\u2265 200"),
      "  ",
      cli::col_yellow("100\u2013199"),
      "  ",
      cli::col_red("< 100"),
      "\n", sep = "")
  if (hasRhat) {
    cat("  Rhat: \u2264 1.01 OK  ",
        cli::col_yellow("1.01\u20131.05"),
        "  ",
        cli::col_red("> 1.05"),
        "\n", sep = "")
  } else {
    cli::cli_alert_info(
      "Single run: R-hat not available. Use {.code nRuns >= 2} for convergence checking."
    )
  }

  invisible(x)
}


# Color an ESS value without fixed-width padding (used in range summaries)
.ColorEss <- function(ess) {
  if (is.na(ess) || !is.finite(ess)) return("NA")
  s <- as.character(round(ess))
  if (ess < 100) cli::col_red(s)
  else if (ess < 200) cli::col_yellow(s)
  else s
}


# Format ESS value with color coding
.FmtEss <- function(ess) {
  if (is.na(ess) || !is.finite(ess)) return(formatC("NA", width = 6))
  s <- formatC(round(ess), width = 6, format = "d")
  if (ess < 100) cli::col_red(s)
  else if (ess < 200) cli::col_yellow(s)
  else s
}


# Format R-hat value with color coding
# Thresholds: <= 1.01 OK, 1.01-1.05 caution, > 1.05 poor
# (Vehtari et al. 2021 recommend < 1.01 for reliable inference)
.FmtRhat <- function(rhat) {
  if (is.na(rhat) || !is.finite(rhat)) return(formatC("NA", width = 7))
  s <- formatC(rhat, width = 7, digits = 3, format = "f")
  if (rhat > 1.05) cli::col_red(s)
  else if (rhat > 1.01) cli::col_yellow(s)
  else s
}


#' Compute ESS for a matrix of samples
#' @keywords internal
.ComputeEss <- function(samples) {
  .EssMatrix(samples)
}


#' Compute R-hat across independent runs
#'
#' Computes rank-normalized split-chain R-hat (Vehtari et al. 2021)
#' for each parameter.  Each run provides one chain; the split-chain
#' procedure further splits each chain in half.
#'
#' @param perRun List of per-run objects, each with `$samples` matrix.
#' @param keyCols Integer vector of column indices.
#' @return Named numeric vector of R-hat values.
#' @keywords internal
.ComputeRhat <- function(perRun, keyCols) {
  nRuns <- length(perRun)
  paramNames <- colnames(perRun[[1]]$samples[, keyCols, drop = FALSE])
  nParams <- length(paramNames)

  rhat <- vapply(seq_len(nParams), function(j) {
    # Build nIter x nRuns matrix for this parameter
    chainMat <- do.call(cbind, lapply(perRun, function(r) {
      r$samples[, keyCols[j]]
    }))
    .Rhat(chainMat)
  }, numeric(1))

  names(rhat) <- paramNames
  rhat
}


#' Identify key parameter columns (exclude branch lengths)
#' @keywords internal
.KeyParamCols <- function(samples) {
  grep("^(log_|tree_|rate_|p$|kPrime_)", colnames(samples))
}


#' Identify scalar parameter columns for plotting (exclude kPrime_ and br_)
#' @keywords internal
.PlotParamCols <- function(samples) {
  grep("^(log_posterior$|tree_|rate_|p$)", colnames(samples))
}


#' Print a per-parameter ESS/R-hat progress table during MCMC
#'
#' Called at each convergence-check interval to show a table of per-parameter
#' ESS and R-hat (when available). Mirrors the format of
#' [print.MkpDiagnostics()].
#'
#' On dynamic terminals, consecutive tables overwrite each other using ANSI
#' cursor-up codes. Pass `prevLines` (the return value of the previous call)
#' to enable overwriting.
#'
#' @param diagCheck Return value of `.CheckConvergence()`.
#' @param nRuns Number of independent runs.
#' @param iter Current iteration number.
#' @param nSamples Total saved samples across all runs.
#' @param prevLines Number of lines printed by the previous call (0 on first
#'   call). Used to overwrite the old table on dynamic terminals.
#' @return Number of lines printed (invisibly), for passing as `prevLines`
#'   to the next call.
#' @keywords internal
.PrintProgressTable <- function(diagCheck, nRuns, iter, nSamples,
                                prevLines = 0L) {
  ess  <- diagCheck$ess
  rhat <- diagCheck$rhat
  hasRhat <- !is.null(rhat)

  nms       <- names(ess)
  scalarNms <- nms[!grepl("^kPrime_", nms) & nms != "log_likelihood"]
  kPrimeNms <- nms[grepl("^kPrime_", nms)]

  # Build all output as a character vector (one element per line)
  out <- character()
  out <- c(out, cli::rule(left = sprintf(
    "Progress diagnostics  (iter %d | %d run%s | %d samples)",
    iter, nRuns, if (nRuns == 1L) "" else "s", nSamples
  )))

  if (hasRhat) {
    out <- c(out, sprintf("  %-20s  %6s  %7s", "Parameter", "ESS", "Rhat"))
  } else {
    out <- c(out, sprintf("  %-20s  %6s", "Parameter", "ESS"))
  }
  out <- c(out, sprintf("  %s", strrep("-", if (hasRhat) 38L else 28L)))

  for (nm in scalarNms) {
    essStr <- .FmtEss(ess[[nm]])
    if (hasRhat && nm %in% names(rhat)) {
      out <- c(out, sprintf("  %-20s  %s  %s", nm, essStr,
                            .FmtRhat(rhat[[nm]])))
    } else {
      out <- c(out, sprintf("  %-20s  %s", nm, essStr))
    }
  }

  # kPrime summary row (min / median / max)
  if (length(kPrimeNms) > 0L) {
    kpEss    <- ess[kPrimeNms]
    kpFinite <- kpEss[is.finite(kpEss)]
    if (length(kpFinite) > 0L) {
      kpMin <- min(kpFinite)
      kpMed <- stats::median(kpFinite)
      kpMax <- max(kpFinite)
    } else {
      kpMin <- kpMed <- kpMax <- NA_real_
    }
    essRange <- paste0(
      .ColorEss(kpMin), " / ", .ColorEss(kpMed), " / ", .ColorEss(kpMax),
      "  (min/med/max)"
    )
    label <- sprintf("kPrime (%d)", length(kPrimeNms))
    if (hasRhat && any(kPrimeNms %in% names(rhat))) {
      kpRhat    <- rhat[kPrimeNms[kPrimeNms %in% names(rhat)]]
      kpRhatFin <- kpRhat[is.finite(kpRhat)]
      if (length(kpRhatFin) > 0L) {
        rhatRange <- sprintf("%.3f\u2013%.3f", min(kpRhatFin), max(kpRhatFin))
        out <- c(out, sprintf("  %-20s  %s  Rhat %s", label, essRange,
                              rhatRange))
      } else {
        out <- c(out, sprintf("  %-20s  %s", label, essRange))
      }
    } else {
      out <- c(out, sprintf("  %-20s  %s", label, essRange))
    }
  }

  out <- c(out, "", paste0(
    "  ESS: ",
    cli::col_green("\u2265 200"), "  ",
    cli::col_yellow("100\u2013199"), "  ",
    cli::col_red("< 100")
  ))
  if (hasRhat) {
    out <- c(out, paste0(
      "  Rhat: \u2264 1.01 OK  ",
      cli::col_yellow("1.01\u20131.05"), "  ",
      cli::col_red("> 1.05")
    ))
  }

  nLines <- length(out)

  # On dynamic terminals, erase the previous table before printing the new one
  if (prevLines > 0L && cli::is_dynamic_tty()) {
    cat(sprintf("\x1b[%dA\x1b[0J", prevLines))
  }
  cat(paste(out, collapse = "\n"), "\n", sep = "")
  invisible(nLines)
}


#' Build the full ticker page sequence for the progress bar
#'
#' Returns a character vector of page content strings (without the
#' `logP:` prefix, which the caller prepends).  Summary pages
#' (minESS/R-hat) are interleaved with detail pages that show 2
#' parameters each with their full names and colour-coded ESS.
#'
#' Colour matches the table conventions: red < 100, yellow 100--199,
#' plain >= 200.  Non-finite ESS values (e.g. `rate_loss` when no
#' neomorphic characters are present) are silently dropped.
#'
#' @param diagCheck Return value of [.CheckConvergence()].
#' @return Character vector of page strings.
#' @keywords internal
.BuildTickerPages <- function(diagCheck) {
  # Single page: minESS (+ R-hat when multi-run).
  # Per-parameter detail removed (M-139): live trace + ESS panel
  # in the plot callback provides richer information.
  .TickerSummaryStr(diagCheck)
}


#' Build headline summary string for progress ticker
#'
#' Produces `"minESS: 42"` (single run) or
#' `"minESS: 42 | Rhat: 1.003"` (multi-run) with colour matching
#' the table conventions.  Used on the dominant ticker page (M-097).
#'
#' @param diagCheck Return value of `.CheckConvergence()`.
#' @return Single string.
#' @keywords internal
.TickerSummaryStr <- function(diagCheck) {
  minEss <- diagCheck$minEss
  if (!is.finite(minEss)) {
    essStr <- "?"
  } else {
    rval <- round(minEss)
    sval <- as.character(rval)
    essStr <- if (rval < 100) cli::col_red(sval)
              else if (rval < 200) cli::col_yellow(sval)
              else cli::col_green(sval)
  }
  s <- paste0("minESS: ", essStr)

  maxRhat <- diagCheck$maxRhat
  if (!is.na(maxRhat) && is.finite(maxRhat)) {
    rhatFmt <- sprintf("%.3f", maxRhat)
    rhatStr <- if (maxRhat > 1.05) cli::col_red(rhatFmt)
               else if (maxRhat > 1.01) cli::col_yellow(rhatFmt)
               else cli::col_green(rhatFmt)
    s <- paste0(s, " \u2502 Rhat: ", rhatStr)
  }
  s
}


# Compute topology ESS from sampled trees using internal .TreeESS + TreeDist RF.
# Returns named numeric vector (frechetCorrelationESS, medianPseudoESS)
# as the minimum across runs (conservative), or NULL if skipped or failed.
#
# When frechet = FALSE (default), only the median pseudo-ESS is computed
# using a cross-distance matrix (maxRows x n), which is much cheaper than
# the full n x n pairwise matrix required for Fréchet ESS.
.ComputeTreeEss <- function(pb, trees, frechet = FALSE) {
  if (isFALSE(trees)) return(NULL)
  if (!requireNamespace("TreeDist", quietly = TRUE)) return(NULL)

  # Build per-run tree lists (post-burnin)
  perRunTrees <- if (!is.null(pb$per_run)) {
    lapply(pb$per_run, `[[`, "trees")
  } else {
    list(pb$trees)
  }
  perRunTrees <- Filter(function(x) length(x) >= 5L, perRunTrees)
  if (length(perRunTrees) == 0L) return(NULL)

  # Subsample each run to at most 1000 trees
  maxPerRun <- 1000L
  totalTrees <- sum(vapply(perRunTrees, length, integer(1L)))
  if (any(vapply(perRunTrees, length, integer(1L)) > maxPerRun)) {
    if (interactive()) {
      cli::cli_alert_info(
        "Tree ESS: subsampling to {maxPerRun} trees/run \\
         ({totalTrees} total available)."
      )
    }
    perRunTrees <- lapply(perRunTrees, function(tr) {
      n <- length(tr)
      if (n > maxPerRun) tr[round(seq(1, n, length.out = maxPerRun))] else tr
    })
  }

  if (interactive()) cli::cli_progress_message("Computing tree ESS\u2026")

  tryCatch({
    chainRows <- lapply(perRunTrees, function(chain) {
      .TreeESS(chain, dist_fn = TreeDist::RobinsonFoulds,
               frechet = frechet)
    })
    essMat <- do.call(rbind, chainRows)
    # Minimum across runs — conservative multi-chain estimate.
    apply(essMat, 2, min, na.rm = TRUE)
  }, error = function(e) {
    cli::cli_warn("Tree ESS computation failed: {conditionMessage(e)}")
    NULL
  })
}
