# Step 5: ESS estimator sanity check.
#
# Test whether .Ess() is robust on slow-mixing integer-drift traces by
# generating AR(1) traces with known autocorrelation rho, thinning by
# 1/10/100, and recomputing ESS.
#
# Expectation: true_ess = n * (1 - rho) / (1 + rho).
# Thinning by t: rho_t = rho^t, true_ess_t = (n/t) * (1-rho^t) / (1+rho^t).
# Stable across thinning means the estimator is sound for that rho.

suppressPackageStartupMessages({
  devtools::load_all(quiet = TRUE)
})

set.seed(2026)

.AR1 <- function(n, rho, sd = 1) {
  x <- numeric(n)
  x[1] <- rnorm(1, sd = sd / sqrt(1 - rho^2))
  noise <- rnorm(n - 1, sd = sd)
  for (i in 2:n) x[i] <- rho * x[i-1] + noise[i-1]
  x
}

.IntegerDrift <- function(n, rho, levels = 2:5, level_prob = c(0.7, 0.2, 0.07, 0.03)) {
  # Latent AR(1) -> discretise via CDF inversion.
  z <- .AR1(n, rho)
  u <- pnorm(z)
  cdf <- cumsum(level_prob) / sum(level_prob)
  k <- sapply(u, function(x) levels[which(x <= cdf)[1]])
  k
}

n <- 30000
rhos <- c(0.5, 0.9, 0.99, 0.999)

cat("=== Continuous AR(1) ===\n")
cat(sprintf("%-7s %-9s %-11s %-11s %-11s %-11s\n",
            "rho", "true_ess", "ess_thin=1", "ess_thin=10", "ess_thin=100", "verdict"))
for (rho in rhos) {
  x <- .AR1(n, rho)
  true_ess <- n * (1 - rho) / (1 + rho)
  e1   <- .Ess(x)
  e10  <- .Ess(x[seq(1, n, by = 10)])
  e100 <- .Ess(x[seq(1, n, by = 100)])
  verdict <- if (abs(e1 - e10) / e1 < 0.2 && abs(e1 - e100) / e1 < 0.4) "STABLE"
             else if (e1 < e10 && e1 < e100) "PESSIMISTIC"
             else if (e1 > e10 && e1 > e100) "OPTIMISTIC"
             else "MIXED"
  cat(sprintf("%-7.3f %-9.0f %-11.0f %-11.0f %-11.0f %-11s\n",
              rho, true_ess, e1, e10, e100, verdict))
}

cat("\n=== Integer-drift trace (mimics k') ===\n")
cat(sprintf("%-7s %-9s %-11s %-11s %-11s %-11s\n",
            "rho", "—", "ess_thin=1", "ess_thin=10", "ess_thin=100", "verdict"))
for (rho in rhos) {
  k <- .IntegerDrift(n, rho)
  e1   <- .Ess(k)
  e10  <- .Ess(k[seq(1, n, by = 10)])
  e100 <- .Ess(k[seq(1, n, by = 100)])
  verdict <- if (abs(e1 - e10) / e1 < 0.2 && abs(e1 - e100) / e1 < 0.4) "STABLE"
             else if (e1 < e10 && e1 < e100) "PESSIMISTIC"
             else if (e1 > e10 && e1 > e100) "OPTIMISTIC"
             else "MIXED"
  cat(sprintf("%-7.3f %-9s %-11.0f %-11.0f %-11.0f %-11s\n",
              rho, "-", e1, e10, e100, verdict))
}

cat("\n=== Bimodal trace (mimics topology switches) ===\n")
cat(sprintf("%-7s %-11s %-11s %-11s %-11s\n",
            "rho", "ess_thin=1", "ess_thin=10", "ess_thin=100", "verdict"))
for (rho in c(0.99, 0.999, 0.9995)) {
  z <- .AR1(n, rho)
  x <- ifelse(z > 0, 1, 0) + rnorm(n, sd = 0.01)  # noisy bimodal
  e1   <- .Ess(x)
  e10  <- .Ess(x[seq(1, n, by = 10)])
  e100 <- .Ess(x[seq(1, n, by = 100)])
  verdict <- if (abs(e1 - e10) / e1 < 0.2 && abs(e1 - e100) / e1 < 0.4) "STABLE"
             else if (e1 < e10 && e1 < e100) "PESSIMISTIC"
             else if (e1 > e10 && e1 > e100) "OPTIMISTIC"
             else "MIXED"
  cat(sprintf("%-7.4f %-11.0f %-11.0f %-11.0f %-11s\n",
              rho, e1, e10, e100, verdict))
}
