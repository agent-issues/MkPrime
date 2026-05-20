## Verify detailed balance for case 30 logit-MH on p with EG-like prior.
## Self-contained: implements the C++ acceptance ratio and checks the
## stationary distribution matches the analytic target.

set.seed(7)
options(scipen = 99)

## Target ∝ π(p) with EG prior at fixed k'=3, body=(0.6, 0.3, 0.1), tail=0.
## P(k'=3) = 0.6 p (1-p) + 0.3 p
log_target <- function(p) {
  if (p <= 0 || p >= 1) return(-Inf)
  log(0.6 * p * (1 - p) + 0.3 * p) + dbeta(p, 1, 1, log = TRUE)
}

## Bactrian-like symmetric perturbation on logit-scale.
sigma <- 1.0  # step
bactrian <- function() {
  m <- 0.95
  sdv <- sqrt(1 - m * m)
  z <- rnorm(1, 0, sdv)
  raw <- if (runif(1) < 0.5) m + z else -m + z
  raw / sqrt(12)
}

step <- function(p) {
  lp <- log(p / (1 - p))
  lp_new <- lp + sigma * bactrian()
  pn <- plogis(lp_new)
  logH <- log(pn) + log1p(-pn) - log(p) - log1p(-p)
  logA <- log_target(pn) - log_target(p) + logH
  if (log(runif(1)) < logA) pn else p
}

nIter <- 2e5
samples <- numeric(nIter)
p <- 0.5
for (i in seq_len(nIter)) {
  p <- step(p)
  samples[i] <- p
}

## Analytic target via numerical integration.
xx <- seq(1e-4, 1 - 1e-4, length.out = 10000)
dens_un <- exp(vapply(xx, log_target, numeric(1)))
Z <- sum(dens_un) * (xx[2] - xx[1])
target_pdf <- dens_un / Z

## Compare CDF: KS-like distance + mean
emp_mean <- mean(samples)
ana_mean <- sum(xx * target_pdf) * (xx[2] - xx[1])
cat(sprintf("empirical mean = %.4f\n", emp_mean))
cat(sprintf("analytic  mean = %.4f\n", ana_mean))
cat(sprintf("diff = %.4f\n", emp_mean - ana_mean))

## KS distance between CDFs.
emp_cdf <- ecdf(samples)
ana_cdf <- cumsum(target_pdf) * (xx[2] - xx[1])
ks <- max(abs(emp_cdf(xx) - ana_cdf))
cat(sprintf("KS distance = %.4f\n", ks))

stopifnot(abs(emp_mean - ana_mean) < 0.01)
stopifnot(ks < 0.02)
cat("DETAILED BALANCE OK\n")
