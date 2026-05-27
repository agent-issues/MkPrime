#!/usr/bin/env Rscript
# Sanity check: at posterior sigma values, how much does the ACRV fix shift the
# effective mean rate?
#
# Discrete-mean error of three implementations:
#  BUG    : dnLognormal(0, sigma)            -- pre-fix RB
#  MU_FIX : dnLognormal(-sigma^2/2, sigma)   -- mu fix only, no renorm
#  HL2014 : MU_FIX rates / mean(MU_FIX rates) -- Harrison & Larsson 2014 eq. 2

nCat <- 6
midpoints <- (seq_len(nCat) - 0.5) / nCat

for (sigma in c(0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0)) {
  rates_bug <- exp(qnorm(midpoints, mean = 0, sd = sigma))
  rates_mu  <- exp(qnorm(midpoints, mean = -sigma^2 / 2, sd = sigma))
  rates_hl  <- rates_mu / mean(rates_mu)
  cat(sprintf("sigma=%.2f  BUG=%.4f  MU_FIX=%.4f  HL2014=%.4f\n",
              sigma, mean(rates_bug), mean(rates_mu), mean(rates_hl)))
}
