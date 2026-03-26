# Convergence monitoring for MkPrime MCMC
#
# Phase 5: ESS, PSRF (Gelman-Rubin), convergence diagnostics.

#' Compute convergence diagnostics for an MkPosterior
#'
#' Calculates effective sample size (ESS) per parameter and, when
#' `nRuns >= 2`, the potential scale reduction factor (PSRF, Gelman-Rubin
#' diagnostic) across runs.
#'
#' @param posterior An `MkPosterior` object.
#' @return A list with components:
#'   - `ess`: Named numeric vector of ESS per parameter.
#'   - `min_ess`: Scalar minimum ESS across parameters.
#'   - `psrf`: Named numeric vector of PSRF point estimates per parameter
#'     (only if `nRuns >= 2`; `NULL` otherwise).
#'   - `max_psrf`: Scalar maximum PSRF (or `NA` if single run).
#' @export
convergence_diagnostics <- function(posterior) {
  if (!inherits(posterior, "MkPosterior")) {
    cli::cli_abort("{.arg posterior} must be an {.cls MkPosterior} object.")
  }

  key_cols <- .key_param_cols(posterior$samples)
  nRuns <- posterior$nRuns %||% 1L

  # --- ESS (combined samples) ---
  ess <- .compute_ess(posterior$samples[, key_cols, drop = FALSE])

  # --- PSRF across runs ---
  psrf <- NULL
  max_psrf <- NA_real_

  if (nRuns >= 2L && !is.null(posterior$per_run)) {
    psrf <- .compute_psrf(posterior$per_run, key_cols)
    max_psrf <- max(psrf, na.rm = TRUE)
  }

  list(
    ess = ess,
    min_ess = min(ess),
    psrf = psrf,
    max_psrf = max_psrf
  )
}


#' Compute ESS for a matrix of samples
#' @keywords internal
.compute_ess <- function(samples) {
  if (!requireNamespace("coda", quietly = TRUE)) {
    cli::cli_warn("Package {.pkg coda} needed for ESS; returning NA.")
    return(rep(NA_real_, ncol(samples)))
  }

  ess <- apply(samples, 2, function(col) {
    if (all(is.na(col)) || sd(col, na.rm = TRUE) == 0) return(NA_real_)
    coda::effectiveSize(coda::mcmc(col))
  })
  names(ess) <- colnames(samples)
  ess
}


#' Compute PSRF (Gelman-Rubin) across independent runs
#'
#' Uses `coda::gelman.diag()` on the cold chain samples from each run.
#' @keywords internal
.compute_psrf <- function(per_run, key_cols) {
  if (!requireNamespace("coda", quietly = TRUE)) {
    cli::cli_warn("Package {.pkg coda} needed for PSRF; returning NA.")
    return(NULL)
  }

  # Build mcmc.list: one mcmc object per run
  chain_list <- lapply(per_run, function(r) {
    coda::mcmc(r$samples[, key_cols, drop = FALSE])
  })
  mcmc_list <- coda::mcmc.list(chain_list)

  # gelman.diag returns a list with $psrf (matrix: point est + upper CI)
  gd <- tryCatch(
    coda::gelman.diag(mcmc_list, multivariate = FALSE),
    error = function(e) NULL
  )

  if (is.null(gd)) return(NULL)

  # Extract point estimates
  psrf <- gd$psrf[, 1]
  names(psrf) <- colnames(per_run[[1]]$samples[, key_cols, drop = FALSE])
  psrf
}


#' Identify key parameter columns (exclude branch lengths)
#' @keywords internal
.key_param_cols <- function(samples) {
  grep("^(log_|tree_|rate_|p$|kPrime_)", colnames(samples))
}
