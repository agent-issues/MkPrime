# Convergence monitoring for MkPrime MCMC
#
# Phase 5: ESS, PSRF (Gelman-Rubin), convergence diagnostics.

#' Compute convergence diagnostics for an MkPosterior
#'
#' Calculates effective sample size (ESS) per parameter and, when
#' `nRuns >= 2`, the potential scale reduction factor (PSRF, Gelman-Rubin
#' diagnostic) across runs. Optionally computes topology ESS using the
#' Robinson-Foulds distance via **TreeDist**.
#'
#' @param posterior An `MkPosterior` object.
#' @param trees Logical. If `TRUE` and **TreeDist** is installed, compute
#'   topology ESS (Fréchet correlation ESS and median pseudo-ESS) from the
#'   sampled trees. Defaults to `FALSE` because Robinson-Foulds distance
#'   computation is O(n^2) and can be slow for large posteriors — reserve
#'   for deliberate post-run calls once parameter ESS has been satisfied.
#'   Each run is subsampled to at most 1,000 trees.
#' @return An object of class `MkpDiagnostics`, a list with components:
#'   - `ess`: Named numeric vector of ESS per parameter.
#'   - `minEss`: Scalar minimum ESS across parameters.
#'   - `psrf`: Named numeric vector of PSRF point estimates per parameter
#'     (only if `nRuns >= 2`; `NULL` otherwise).
#'   - `maxPsrf`: Scalar maximum PSRF (or `NA` if single run).
#'   - `treeEss`: Named numeric vector with `frechetCorrelationESS` and
#'     `medianPseudoESS` (or `NULL` if not computed).
#' @export
ConvergenceDiagnostics <- function(posterior, trees = FALSE) {
  if (!inherits(posterior, "MkPosterior")) {
    cli::cli_abort("{.arg posterior} must be an {.cls MkPosterior} object.")
  }

  pb <- .PostBurninData(posterior)
  keyCols <- .KeyParamCols(pb$samples)
  nRuns <- posterior$nRuns %||% 1L

  # --- ESS (combined samples) ---
  ess <- .ComputeEss(pb$samples[, keyCols, drop = FALSE])

  # --- PSRF across runs ---
  psrf <- NULL
  maxPsrf <- NA_real_

  if (nRuns >= 2L && !is.null(posterior$per_run)) {
    psrf <- .ComputePsrf(pb$per_run, keyCols)
    maxPsrf <- max(psrf, na.rm = TRUE)
  }

  # --- Tree ESS ---
  treeEss <- .ComputeTreeEss(pb, trees)

  structure(
    list(
      ess = ess,
      minEss = min(ess),
      psrf = psrf,
      maxPsrf = maxPsrf,
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
#' Displays a colour-coded table of ESS and PSRF values for each parameter.
#' ESS < 100 is shown in red, 100--199 in yellow, and ≥ 200 in plain text.
#' Individual k' parameters are summarised as a single min–median–max row.
#'
#' @param x An `MkpDiagnostics` object from [ConvergenceDiagnostics()].
#' @param ... Ignored.
#' @return `x` invisibly.
#' @export
print.MkpDiagnostics <- function(x, ...) {
  hasPsrf <- !is.null(x$psrf)

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
  if (hasPsrf) {
    cat(sprintf("  %-20s  %6s  %7s\n", "Parameter", "ESS", "PSRF"))
  } else {
    cat(sprintf("  %-20s  %6s\n", "Parameter", "ESS"))
  }
  cat(sprintf("  %s\n", strrep("-", if (hasPsrf) 38L else 28L)))

  # Scalar parameter rows
  for (nm in scalarNms) {
    essStr <- .FmtEss(x$ess[nm])
    if (hasPsrf && nm %in% names(x$psrf)) {
      cat(sprintf("  %-20s  %s  %s\n", nm, essStr, .FmtPsrf(x$psrf[nm])))
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
    if (hasPsrf && any(kPrimeNms %in% names(x$psrf))) {
      kpPsrf <- x$psrf[kPrimeNms[kPrimeNms %in% names(x$psrf)]]
      psrfRange <- sprintf("%.3f\u2013%.3f",
        min(kpPsrf, na.rm = TRUE), max(kpPsrf, na.rm = TRUE))
      cat(sprintf("  %-20s  %s  PSRF %s\n", label, essRange, psrfRange))
    } else {
      cat(sprintf("  %-20s  %s\n", label, essRange))
    }
  }

  # Topology ESS row (always shown; NA when not computed)
  cat(sprintf("  %s\n", strrep("-", if (hasPsrf) 38L else 28L)))
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
  if (hasPsrf) {
    cat("  PSRF: \u2264 1.05 OK  ",
        cli::col_yellow("1.05\u20131.1"),
        "  ",
        cli::col_red("> 1.1"),
        "\n", sep = "")
  } else {
    cli::cli_alert_info(
      "Single run: PSRF not available. Use {.code nRuns >= 2} for convergence checking."
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


# Format PSRF value with color coding
.FmtPsrf <- function(psrf) {
  if (is.na(psrf) || !is.finite(psrf)) return(formatC("NA", width = 7))
  s <- formatC(psrf, width = 7, digits = 3, format = "f")
  if (psrf > 1.1) cli::col_red(s)
  else if (psrf > 1.05) cli::col_yellow(s)
  else s
}


#' Compute ESS for a matrix of samples
#' @keywords internal
.ComputeEss <- function(samples) {
  if (!requireNamespace("coda", quietly = TRUE)) {
    cli::cli_warn("Package {.pkg coda} needed for ESS; returning NA.")
    return(rep(NA_real_, ncol(samples)))
  }

  ess <- apply(samples, 2, function(col) {
    s <- sd(col, na.rm = TRUE)
    if (is.na(s) || s == 0) return(NA_real_)
    coda::effectiveSize(coda::mcmc(col))
  })
  names(ess) <- colnames(samples)
  ess
}


#' Compute PSRF (Gelman-Rubin) across independent runs
#'
#' Uses `coda::gelman.diag()` on the cold chain samples from each run.
#' @keywords internal
.ComputePsrf <- function(perRun, keyCols) {
  if (!requireNamespace("coda", quietly = TRUE)) {
    cli::cli_warn("Package {.pkg coda} needed for PSRF; returning NA.")
    return(NULL)
  }

  # Build mcmc.list: one mcmc object per run
  chainList <- lapply(perRun, function(r) {
    coda::mcmc(r$samples[, keyCols, drop = FALSE])
  })
  mcmcList <- coda::mcmc.list(chainList)

  # gelman.diag returns a list with $psrf (matrix: point est + upper CI)
  gd <- tryCatch(
    coda::gelman.diag(mcmcList, multivariate = FALSE),
    error = function(e) NULL
  )

  if (is.null(gd)) return(NULL)

  # Extract point estimates
  psrf <- gd$psrf[, 1]
  names(psrf) <- colnames(perRun[[1]]$samples[, keyCols, drop = FALSE])
  psrf
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


#' Print a per-parameter ESS/PSRF progress table during MCMC
#'
#' Called at each convergence-check interval to show a table of per-parameter
#' ESS and PSRF (when available). Mirrors the format of
#' [print.MkpDiagnostics()].
#'
#' @param diagCheck Return value of `.CheckConvergence()`.
#' @param nRuns Number of independent runs.
#' @param iter Current iteration number.
#' @param nSamples Total saved samples across all runs.
#' @keywords internal
.PrintProgressTable <- function(diagCheck, nRuns, iter, nSamples) {
  ess  <- diagCheck$ess
  psrf <- diagCheck$psrf
  hasPsrf <- !is.null(psrf)

  nms       <- names(ess)
  scalarNms <- nms[!grepl("^kPrime_", nms) & nms != "log_likelihood"]
  kPrimeNms <- nms[grepl("^kPrime_", nms)]

  cli::cli_rule(
    left = sprintf(
      "Progress diagnostics  (iter %d | %d run%s | %d samples)",
      iter, nRuns, if (nRuns == 1L) "" else "s", nSamples
    )
  )

  if (hasPsrf) {
    cat(sprintf("  %-20s  %6s  %7s\n", "Parameter", "ESS", "PSRF"))
  } else {
    cat(sprintf("  %-20s  %6s\n", "Parameter", "ESS"))
  }
  cat(sprintf("  %s\n", strrep("-", if (hasPsrf) 38L else 28L)))

  for (nm in scalarNms) {
    essStr <- .FmtEss(ess[[nm]])
    if (hasPsrf && nm %in% names(psrf)) {
      cat(sprintf("  %-20s  %s  %s\n", nm, essStr, .FmtPsrf(psrf[[nm]])))
    } else {
      cat(sprintf("  %-20s  %s\n", nm, essStr))
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
    if (hasPsrf && any(kPrimeNms %in% names(psrf))) {
      kpPsrf    <- psrf[kPrimeNms[kPrimeNms %in% names(psrf)]]
      kpPsrfFin <- kpPsrf[is.finite(kpPsrf)]
      if (length(kpPsrfFin) > 0L) {
        psrfRange <- sprintf("%.3f\u2013%.3f", min(kpPsrfFin), max(kpPsrfFin))
        cat(sprintf("  %-20s  %s  PSRF %s\n", label, essRange, psrfRange))
      } else {
        cat(sprintf("  %-20s  %s\n", label, essRange))
      }
    } else {
      cat(sprintf("  %-20s  %s\n", label, essRange))
    }
  }

  cat("\n  ESS: ",
      cli::col_green("\u2265 200"), "  ",
      cli::col_yellow("100\u2013199"), "  ",
      cli::col_red("< 100"),
      "\n", sep = "")
  if (hasPsrf) {
    cat("  PSRF: \u2264 1.05 OK  ",
        cli::col_yellow("1.05\u20131.1"), "  ",
        cli::col_red("> 1.1"),
        "\n", sep = "")
  }
  invisible(NULL)
}


# Compute topology ESS from sampled trees using internal .TreeESS + TreeDist RF.
# Returns named numeric vector (frechetCorrelationESS, medianPseudoESS)
# as the minimum across runs (conservative), or NULL if skipped or failed.
#
# NOTE: tree ESS via RF distances is expensive (O(n^2) distance matrix) and
# is intentionally excluded from print()/summary() and from the checkEvery
# polling callback.  Call ConvergenceDiagnostics(posterior, trees = TRUE)
# explicitly after a run has finished — or after parameter ESS has been
# satisfied — to obtain topology ESS.
.ComputeTreeEss <- function(pb, trees) {
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
      .TreeESS(chain, dist_fn = TreeDist::RobinsonFoulds)
    })
    essMat <- do.call(rbind, chainRows)
    # Minimum across runs — conservative multi-chain estimate.
    apply(essMat, 2, min, na.rm = TRUE)
  }, error = function(e) {
    cli::cli_warn("Tree ESS computation failed: {conditionMessage(e)}")
    NULL
  })
}
