# Model specification for MkPrime MCMC
#
# Defines priors and model options. Used by RunMkPrime() and log_prior().

#' Specify an MkPrime model
#'
#' @param coding Ascertainment bias correction: `"variable"` (default) or
#'   `"none"`.
#' @param nCat Number of ACRV rate categories (default 6).
#' @param relabel Apply Mk' relabelling correction for transformational
#'   characters? Default `TRUE`.
#' @param tree_length_shape,tree_length_rate Shape and rate for the Gamma prior
#'   on tree length. Defaults: shape = 2, rate = 2 / `exp_steps`.
#' @param exp_steps Expected number of character state changes. Used to set
#'   the tree length prior scale. Default `NULL` (computed from parsimony score
#'   of the starting tree).
#' @param rate_loss_meanlog,rate_loss_sdlog Parameters for the LogNormal prior
#'   on `rate_loss` (neomorphic asymmetry). Defaults: meanlog = 0, sdlog = 2.
#' @param rate_log_sd_shape,rate_log_sd_rate Shape and rate for the Gamma prior
#'   on `rate_log_sd` (ACRV dispersion). Defaults: shape = 1, rate = 1.
#' @param kprime_hyper_a,kprime_hyper_b Parameters for the Beta prior on the
#'   geometric hyperprior parameter `p`. Defaults: a = 1, b = 1 (uniform).
#'
#' @return An S3 object of class `MkPrimeModel`.
#' @export
MkPrimeModel <- function(
    coding = "variable",
    nCat = 6L,
    relabel = TRUE,
    tree_length_shape = 2,
    tree_length_rate = NULL,
    exp_steps = NULL,
    rate_loss_meanlog = 0,
    rate_loss_sdlog = 2,
    rate_log_sd_shape = 1,
    rate_log_sd_rate = 1,
    kprime_hyper_a = 1,
    kprime_hyper_b = 1
) {
  coding <- match.arg(coding, c("variable", "none"))

  # Derive tree_length_rate from exp_steps if not provided

  if (is.null(tree_length_rate) && !is.null(exp_steps)) {
    tree_length_rate <- 2 / exp_steps
  }

  structure(
    list(
      coding = coding,
      nCat = as.integer(nCat),
      relabel = relabel,
      tree_length_shape = tree_length_shape,
      tree_length_rate = tree_length_rate,
      exp_steps = exp_steps,
      rate_loss_meanlog = rate_loss_meanlog,
      rate_loss_sdlog = rate_loss_sdlog,
      rate_log_sd_shape = rate_log_sd_shape,
      rate_log_sd_rate = rate_log_sd_rate,
      kprime_hyper_a = kprime_hyper_a,
      kprime_hyper_b = kprime_hyper_b
    ),
    class = "MkPrimeModel"
  )
}


#' Finalize model with data-derived defaults
#'
#' Sets `exp_steps` and `tree_length_rate` if not user-specified.
#' Called internally by [RunMkPrime()] before MCMC starts.
#'
#' @param model An `MkPrimeModel` object.
#' @param tree A `phylo` object (starting tree).
#' @param mkd An `MkPrimeData` object.
#' @return Updated `MkPrimeModel` with all defaults resolved.
#' @keywords internal
.finalize_model <- function(model, tree, mkd) {
  if (is.null(model$exp_steps)) {
    # Default: nChar * mean_kObs / 4 gives a rough expected changes estimate.
    # A better default would use parsimony score, but that requires the
    # original phyDat. Users can set exp_steps explicitly for more control.
    mean_kObs <- mean(mkd$kObs)
    model$exp_steps <- max(1, mkd$nChar * (mean_kObs - 1) / 2)
  }
  if (is.null(model$tree_length_rate)) {
    model$tree_length_rate <- 2 / model$exp_steps
  }
  model
}


#' Compute total log-prior density
#'
#' @param state A list with current parameter values:
#'   `tree_length`, `rel_br_lengths`, `rate_loss`, `rate_log_sd`,
#'   `kPrime` (integer vector), `p` (hyperprior).
#' @param model An `MkPrimeModel` object (finalized).
#' @param mkd An `MkPrimeData` object (for kObs and character types).
#' @return Scalar log-prior density.
#' @keywords internal
log_prior <- function(state, model, mkd) {
  # Boundary checks — return -Inf for out-of-support values

  if (state$tree_length <= 0) return(-Inf)
  if (state$rate_log_sd < 0) return(-Inf)
  if (any(state$rel_br_lengths <= 0)) return(-Inf)

  has_neo <- any(mkd$type == "neomorphic")
  if (has_neo && state$rate_loss <= 0) return(-Inf)

  trans_idx <- which(mkd$type == "transformational")
  if (length(trans_idx)) {
    if (state$p <= 0 || state$p >= 1) return(-Inf)
    if (any(state$kPrime[trans_idx] < mkd$kObs[trans_idx])) return(-Inf)
  }

  lp <- 0.0

  # Tree length: Gamma prior
  lp <- lp + dgamma(state$tree_length,
                     shape = model$tree_length_shape,
                     rate = model$tree_length_rate,
                     log = TRUE)

  # Relative branch lengths: Dirichlet(1, ..., 1) = uniform on simplex
  # log-density is constant: log((n-1)!). Doesn't affect MH ratios but
  # included for correct log-posterior reporting.
  n_edges <- length(state$rel_br_lengths)
  lp <- lp + lfactorial(n_edges - 1L)

  # rate_loss: LogNormal prior (neomorphic characters only)
  if (has_neo) {
    lp <- lp + dlnorm(state$rate_loss,
                      meanlog = model$rate_loss_meanlog,
                      sdlog = model$rate_loss_sdlog,
                      log = TRUE)
  }

  # rate_log_sd: Gamma prior (rate_log_sd = 0 is a boundary; dgamma(0) = 0
  # for shape >= 1, but log(0) = -Inf. Treat 0 specially as a valid point.)
  if (state$rate_log_sd > 0) {
    lp <- lp + dgamma(state$rate_log_sd,
                       shape = model$rate_log_sd_shape,
                       rate = model$rate_log_sd_rate,
                       log = TRUE)
  }
  # When rate_log_sd == 0 and shape == 1, the density is finite (rate);
  # when shape > 1, density is 0. Handle both:
  if (state$rate_log_sd == 0 && model$rate_log_sd_shape > 1) return(-Inf)

  # k'_i: Geometric(p) shifted by kObs_i
  # P(k'_i = kObs_i + u) = p * (1-p)^u, u = 0, 1, 2, ...
  if (length(trans_idx)) {
    u <- state$kPrime[trans_idx] - mkd$kObs[trans_idx]
    lp <- lp + length(trans_idx) * log(state$p) + sum(u) * log1p(-state$p)

    # p: Beta hyperprior
    lp <- lp + dbeta(state$p,
                     shape1 = model$kprime_hyper_a,
                     shape2 = model$kprime_hyper_b,
                     log = TRUE)
  }

  lp
}


#' @export
print.MkPrimeModel <- function(x, ...) {
  cli::cli_h1("MkPrime Model")
  cli::cli_ul(c(
    "Coding: {x$coding}",
    "ACRV categories: {x$nCat}",
    "Relabelling correction: {x$relabel}",
    "Tree length prior: Gamma({x$tree_length_shape}, {x$tree_length_rate %||% 'auto'})",
    "rate_loss prior: LogNormal({x$rate_loss_meanlog}, {x$rate_loss_sdlog})",
    "rate_log_sd prior: Gamma({x$rate_log_sd_shape}, {x$rate_log_sd_rate})",
    "k' hyperprior p: Beta({x$kprime_hyper_a}, {x$kprime_hyper_b})"
  ))
  invisible(x)
}
