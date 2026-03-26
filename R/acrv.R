# Among-Character Rate Variation (ACRV)
#
# Discretized lognormal distribution following Wagner (2012).
# Each rate category gets a rate multiplier; the rates are normalized
# so their mean is 1.0 (rates are relative).

#' Compute discretized lognormal rate categories
#'
#' @param rate_log_sd Standard deviation on the log scale (sigma parameter
#'   of the lognormal). When 0, returns a single category with rate 1.
#' @param nCat Number of discrete rate categories (default 6).
#'
#' @return Numeric vector of `nCat` rate multipliers, normalized to mean 1.
#' @keywords internal
discrete_lognormal_rates <- function(rate_log_sd, nCat = 6L) {
  if (rate_log_sd <= 0) {
    return(rep(1.0, nCat))
  }

  # Lognormal with mean 1: mu = -sigma^2/2
  mu <- -rate_log_sd^2 / 2

  # Category boundaries at quantiles
  boundaries <- qnorm(seq(0, 1, length.out = nCat + 1L))

  # Mean of truncated normal in each interval, then exponentiate
  # For a normal(mu, sigma), the mean in [a, b] is:
  #   mu + sigma * (phi(a') - phi(b')) / (Phi(b') - Phi(a'))
  # where a' = (a - mu)/sigma, b' = (b - mu)/sigma
  # But we're working on the log scale, so the rate for category i is
  # the mean of the lognormal within quantile interval [q_i, q_{i+1}].

  # Simpler approach: use the midpoint quantile of each category
  midpoints <- (seq_len(nCat) - 0.5) / nCat
  log_rates <- qnorm(midpoints, mean = mu, sd = rate_log_sd)
  rates <- exp(log_rates)

  # Normalize to mean 1
  rates / mean(rates)
}
