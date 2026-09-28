# Among-Character Rate Variation (ACRV)
#
# Discretized lognormal distribution following Wagner (2012).
# Each rate category gets a rate multiplier; the rates are normalized
# so their mean is 1.0 (rates are relative).

#' Compute discretized lognormal rate categories
#'
#' @param rateLogSd Standard deviation on the log scale (sigma parameter
#'   of the lognormal). When 0, returns `nCat` categories, all with rate 1.
#' @param nCat Number of discrete rate categories (default 6).
#'
#' @return Numeric vector of `nCat` rate multipliers, normalized to mean 1.
#' @keywords internal
DiscreteLognormalRates <- function(rateLogSd, nCat = 6L) {
  if (rateLogSd <= 0) {
    return(rep(1.0, nCat))
  }

  # Lognormal with mean 1: mu = -sigma^2/2
  mu <- -rateLogSd^2 / 2

  # Each category's rate is the lognormal quantile at its interval's
  # midpoint probability -- the intended approximation.
  midpoints <- (seq_len(nCat) - 0.5) / nCat
  logRates <- qnorm(midpoints, mean = mu, sd = rateLogSd)
  rates <- exp(logRates)

  # Normalize to mean 1
  rates / mean(rates)
}
