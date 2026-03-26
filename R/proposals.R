# MCMC proposal functions for MkPrime
#
# Each proposal returns list(value = ..., log_hastings = ...).
# The Hastings ratio is log q(current | proposed) - log q(proposed | current).

#' Scale proposal for positive scalars
#'
#' Proposes x' = x * exp(tuning * (U - 0.5)) where U ~ Uniform(0,1).
#' Hastings ratio: log(x'/x) = tuning * (U - 0.5).
#'
#' @param x Current value (positive scalar).
#' @param tuning Scale parameter controlling proposal width.
#' @return `list(value, log_hastings)`.
#' @keywords internal
propose_scale <- function(x, tuning = 1.0) {
  u <- runif(1)
  m <- exp(tuning * (u - 0.5))
  list(value = x * m, log_hastings = log(m))
}


#' BetaSimplex proposal for simplex vectors
#'
#' Picks element `index` and a random other element, proposes redistributing
#' mass between them using a Beta draw. Other elements are unchanged.
#'
#' @param x Current simplex vector (sums to 1, all positive).
#' @param index Index of the element to perturb. If NULL, picks randomly.
#' @param tuning Concentration parameter (higher = more conservative).
#'   Effective concentration is `tuning * 2`.
#' @return `list(value, log_hastings)`.
#' @keywords internal
propose_beta_simplex <- function(x, index = NULL, tuning = 10.0) {
  n <- length(x)
  if (n < 2L) {
    return(list(value = x, log_hastings = 0))
  }

  if (is.null(index)) {
    index <- sample.int(n, 1L)
  }
  # Pick a random other element
  other <- sample.int(n - 1L, 1L)
  if (other >= index) other <- other + 1L

  # Current values of the two elements
  old_a <- x[index]
  old_b <- x[other]
  total <- old_a + old_b

  if (total <= 0) {
    return(list(value = x, log_hastings = 0))
  }

  # Current fraction: f = old_a / total
  old_f <- old_a / total

  # Propose new fraction from Beta centered on old_f
  alpha <- old_f * tuning + 1
  beta_param <- (1 - old_f) * tuning + 1
  new_f <- rbeta(1, alpha, beta_param)

  # Compute new values
  new_a <- new_f * total
  new_b <- (1 - new_f) * total
  x_new <- x
  x_new[index] <- new_a
  x_new[other] <- new_b

  # Hastings ratio: q(old|new) / q(new|old)
  # Forward: Beta(old_f * tuning + 1, (1-old_f) * tuning + 1) at new_f
  # Reverse: Beta(new_f * tuning + 1, (1-new_f) * tuning + 1) at old_f
  log_fwd <- dbeta(new_f, alpha, beta_param, log = TRUE)
  rev_alpha <- new_f * tuning + 1
  rev_beta <- (1 - new_f) * tuning + 1
  log_rev <- dbeta(old_f, rev_alpha, rev_beta, log = TRUE)

  list(value = x_new, log_hastings = log_rev - log_fwd)
}


#' BoundedIntegerWalk proposal
#'
#' Proposes x' = x + delta where delta ~ Uniform(-window, ..., window).
#' Rejects (returns current value with log_hastings = -Inf) if x' < lower.
#'
#' @param x Current integer value.
#' @param lower Lower bound (inclusive).
#' @param window Half-width of the proposal window.
#' @return `list(value, log_hastings)`.
#' @keywords internal
propose_bounded_int_walk <- function(x, lower, window = 1L) {
  delta <- sample(-window:window, 1L)
  x_new <- x + delta

  if (x_new < lower) {
    return(list(value = x, log_hastings = -Inf))
  }

  # Symmetric proposal: log_hastings = 0
  # (Both x→x' and x'→x have same number of valid proposals)
  list(value = x_new, log_hastings = 0)
}
