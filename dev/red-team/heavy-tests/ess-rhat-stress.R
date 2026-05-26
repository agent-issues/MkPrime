# dev/red-team/heavy-tests/ess-rhat-stress.R
#
# Lane D3 — ESS / R-hat validity stress test.
#
# Validates MkPrime's `R/ess.R` (Geyer 1992 + Vehtari et al. 2021 truncated
# estimator) and `R/Convergence.R` (rank-normalised split-chain R-hat) against
# adversarial processes whose true integrated autocorrelation time (IAT) or
# convergence diagnostic value is analytically known.
#
# Usage:
#   Rscript dev/red-team/heavy-tests/ess-rhat-stress.R [--quick]
#
#   --quick    Small-N mode (<= 90s on a laptop). One AR(1) rho, small N,
#              fast bimodal & sticky scenarios. Used to confirm execution.
#   (no flag)  Full grid: AR(1) over {0.1,0.5,0.9,0.99,0.999} x {1e3,1e4,1e5},
#              AR(2) cases, bimodal, near-non-stationary, sticky, constant
#              chain, unequal-length CONV-002 regression.
#
# Output:   dev/red-team/heavy-tests/ess-rhat-stress-results/
#   - summary.rds         results table (data.frame)
#   - results.csv         human-readable CSV of all scenarios
#   - verdict.txt         per-row PASS/FAIL/WARN + overall verdict
#
# Reproducibility: set.seed() set explicitly per scenario; full grid replicates
# AR(1) cases 20x to estimate measurement variance.
#
# References:
#   Geyer CJ (1992). Statistical Science 7, 473-483.
#   Vehtari A et al. (2021). Bayesian Analysis 16, 667-718.

suppressWarnings(suppressMessages(pkgload::load_all(".", quiet = TRUE)))

args <- commandArgs(trailingOnly = TRUE)
QUICK <- "--quick" %in% args

OUTDIR <- "dev/red-team/heavy-tests/ess-rhat-stress-results"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# Pull internal functions (intentionally not exported by MkPrime).
.Ess           <- MkPrime:::.Ess
.EssVector     <- MkPrime:::.EssVector
.EssMatrix     <- MkPrime:::.EssMatrix
.EssMultiChain <- MkPrime:::.EssMultiChain
.Rhat          <- MkPrime:::.Rhat
.ComputeRhat   <- MkPrime:::.ComputeRhat

# ---- Theoretical truths --------------------------------------------------

# AR(1):  x_t = rho * x_{t-1} + eps_t,  eps_t ~ N(0,1)
#   stationary var  = 1 / (1 - rho^2)
#   IAT             = (1 + rho) / (1 - rho)
#   true ESS        = N / IAT = N * (1 - rho) / (1 + rho)
SimAR1 <- function(n, rho, sd = 1) {
  x <- numeric(n)
  e <- rnorm(n, sd = sd)
  x[1] <- e[1] / sqrt(1 - rho^2)
  for (i in 2:n) x[i] <- rho * x[i - 1] + e[i]
  x
}

IatAR1 <- function(rho) (1 + rho) / (1 - rho)

# AR(2):  x_t = phi1 * x_{t-1} + phi2 * x_{t-2} + eps_t
# Closed-form IAT = (1 + sum_{k=1}^infty rho_k) where rho_k = ACF at lag k.
# We compute it numerically by summing many lags of the analytical ACF.
SimAR2 <- function(n, phi1, phi2) {
  as.numeric(arima.sim(model = list(ar = c(phi1, phi2)), n = n))
}

IatAR2 <- function(phi1, phi2, lagMax = 2000L) {
  # ACF from Yule-Walker / recursion. acf-rho-k satisfies:
  #   rho_k = phi1 rho_{k-1} + phi2 rho_{k-2}
  # rho_0 = 1; rho_1 = phi1 / (1 - phi2).
  rho <- numeric(lagMax + 1L)
  rho[1] <- 1
  if (lagMax >= 1L) rho[2] <- phi1 / (1 - phi2)
  if (lagMax >= 2L) {
    for (k in 3:(lagMax + 1L)) rho[k] <- phi1 * rho[k - 1L] + phi2 * rho[k - 2L]
  }
  1 + 2 * sum(rho[-1L])
}


# ---- Result accumulator --------------------------------------------------

results <- data.frame(
  scenario   = character(),
  parameter  = character(),
  n          = integer(),
  measured   = numeric(),
  truth      = numeric(),
  tolerance  = numeric(),
  rel_err    = numeric(),
  status     = character(),
  note       = character(),
  stringsAsFactors = FALSE
)

AddRow <- function(scenario, parameter, n, measured, truth, tolerance,
                   status, note = "") {
  rel_err <- if (is.finite(truth) && truth != 0 && is.finite(measured)) {
    abs(measured - truth) / abs(truth)
  } else NA_real_
  results[nrow(results) + 1L, ] <<- list(
    scenario, parameter, as.integer(n), measured, truth, tolerance,
    rel_err, status, note
  )
}

# Decide PASS/FAIL using relative-error tolerance.
Classify <- function(measured, truth, tol) {
  if (!is.finite(measured) || !is.finite(truth)) return("FAIL")
  if (truth == 0) return(if (abs(measured) <= tol) "PASS" else "FAIL")
  if (abs(measured - truth) / abs(truth) <= tol) "PASS" else "FAIL"
}


# =========================================================================
# Scenario 1: AR(1) — true ESS = N (1-rho) / (1+rho)
# =========================================================================
cat("[1/7] AR(1) ESS validity\n")

if (QUICK) {
  ar1_grid <- expand.grid(rho = c(0.5, 0.9), n = 1e4, KEEP.OUT.ATTRS = FALSE)
  nReps <- 3L
} else {
  ar1_grid <- expand.grid(rho = c(0.1, 0.5, 0.9, 0.99, 0.999),
                          n = c(1e3, 1e4, 1e5),
                          KEEP.OUT.ATTRS = FALSE)
  nReps <- 20L
}

# Per-rho tolerance. Geyer's truncated estimator has a known finite-N bias
# at high rho (truncation cuts the slow-decaying tail). Tolerances widen
# accordingly, justified in companion .md.
RhoTol <- function(rho) {
  if      (rho <= 0.5)  0.10
  else if (rho <= 0.9)  0.10
  else if (rho <= 0.99) 0.30
  else                  0.60   # rho = 0.999: estimator strongly biased low
}

for (i in seq_len(nrow(ar1_grid))) {
  rho <- ar1_grid$rho[i]
  n   <- as.integer(ar1_grid$n[i])
  set.seed(1000L + i)
  ess_meas <- replicate(nReps, .Ess(SimAR1(n, rho)))
  measured <- mean(ess_meas)
  truth    <- n * (1 - rho) / (1 + rho)
  tol      <- RhoTol(rho)
  status   <- Classify(measured, truth, tol)
  AddRow("AR(1)", sprintf("rho=%g", rho), n, measured, truth, tol, status,
         sprintf("mean over %d reps", nReps))
}


# =========================================================================
# Scenario 2: AR(2) — analytical IAT via lag-recursion
# =========================================================================
cat("[2/7] AR(2) ESS validity\n")

if (QUICK) {
  ar2_cases <- list(
    list(phi1 = 0.5,  phi2 = -0.2, n = 1e4, label = "mild")
  )
} else {
  ar2_cases <- list(
    list(phi1 = 0.5,  phi2 = -0.2, n = 1e4, label = "mild"),
    list(phi1 = 0.8,  phi2 = -0.3, n = 1e4, label = "oscillatory"),
    list(phi1 = 1.2,  phi2 = -0.4, n = 1e4, label = "persistent")
  )
}

for (cs in ar2_cases) {
  set.seed(2000L)
  reps <- if (QUICK) 3L else 10L
  ess_meas <- replicate(reps, .Ess(SimAR2(cs$n, cs$phi1, cs$phi2)))
  measured <- mean(ess_meas)
  iat      <- IatAR2(cs$phi1, cs$phi2)
  truth    <- cs$n / iat
  # AR(2) tolerance: 15% (less sticky than AR(1)@0.9 typically)
  tol <- 0.15
  status <- Classify(measured, truth, tol)
  AddRow("AR(2)", sprintf("%s phi=(%.2f,%.2f)", cs$label, cs$phi1, cs$phi2),
         cs$n, measured, truth, tol, status,
         sprintf("IAT=%.3f", iat))
}


# =========================================================================
# Scenario 3: Bimodal target — chains stuck in different modes
#   Two chains, each pinned to a different mode (no proposal crosses 0).
#   Each chain is Normal(+/-mu, 1) -> per-chain ESS reasonable.
#   Cross-chain R-hat should be huge (>> 1.1).
# =========================================================================
cat("[3/7] Bimodal: per-mode ESS + cross-chain R-hat\n")

# Average ESS over replicates to reduce Geyer-estimator variance at small N.
# Larger N also makes per-mode-ESS recovery much tighter.
mu <- 5.0
n  <- if (QUICK) 2000L else 5000L
nRepsBi <- if (QUICK) 5L else 10L

set.seed(3000L)
ess1_reps <- replicate(nRepsBi, .Ess(rnorm(n, mean = +mu, sd = 1)))
set.seed(3001L)
ess2_reps <- replicate(nRepsBi, .Ess(rnorm(n, mean = -mu, sd = 1)))
ess1 <- mean(ess1_reps)
ess2 <- mean(ess2_reps)
# Each per-mode chain is iid Normal -> true ESS = N (no autocorrelation).
# Tolerance widened to 20% because Geyer estimator variance at N=2000 is ~10%.
truth_essMode <- n
tolMode <- 0.20
AddRow("Bimodal", "per-mode-ESS chain1", n, ess1, truth_essMode, tolMode,
       Classify(ess1, truth_essMode, tolMode),
       sprintf("mean over %d reps", nRepsBi))
AddRow("Bimodal", "per-mode-ESS chain2", n, ess2, truth_essMode, tolMode,
       Classify(ess2, truth_essMode, tolMode),
       sprintf("mean over %d reps", nRepsBi))

# R-hat must diagnose non-mixing. Bimodal split between chains -> R-hat huge.
# Single draws here -- Rhat is so large that variance is irrelevant.
set.seed(3100L)
chain1 <- rnorm(n, mean = +mu, sd = 1)
chain2 <- rnorm(n, mean = -mu, sd = 1)
mat <- cbind(chain1, chain2)
rhatBimodal <- .Rhat(mat)
status <- if (is.finite(rhatBimodal) && rhatBimodal >= 1.1) "PASS" else "FAIL"
AddRow("Bimodal", "cross-chain Rhat", n, rhatBimodal, NA_real_, 1.1, status,
       "expect Rhat >= 1.1 (separated modes)")


# =========================================================================
# Scenario 4: Near-non-stationary random walk (acceptance ~0%)
#   Random walk Metropolis on N(0,1) with absurdly large proposal sd ->
#   almost every proposal rejected -> chain is nearly constant.
#   Expected: ESS << N; R-hat across two such chains should be large
#   (chains stuck at different positions).
# =========================================================================
cat("[4/7] Near-non-stationary chain\n")

RWMetropolis <- function(n, propSd, target_logpdf, seed) {
  set.seed(seed)
  x <- numeric(n)
  x[1] <- 0
  for (i in 2:n) {
    y <- x[i - 1L] + rnorm(1, sd = propSd)
    if (log(runif(1)) < target_logpdf(y) - target_logpdf(x[i - 1L])) {
      x[i] <- y
    } else {
      x[i] <- x[i - 1L]
    }
  }
  x
}

n <- if (QUICK) 1000L else 5000L
# propSd=100 vs target sd=1 -> acceptance ~ 1% (chain mostly stuck).
chainA <- RWMetropolis(n, 100, function(z) dnorm(z, log = TRUE), seed = 4001L)
chainB <- RWMetropolis(n, 100, function(z) dnorm(z, log = TRUE), seed = 4002L)
essA <- .Ess(chainA)
# Pass criterion: ESS must be <much less than> N. We use ESS / N < 0.1.
status <- if (is.finite(essA) && essA / n < 0.1) "PASS" else "FAIL"
AddRow("NearNonStationary", "single-chain ESS", n, essA, n * 0.05, 1.0, status,
       sprintf("ESS/N = %.3f; expect << 0.1", essA / n))

# R-hat across two such chains: each is stuck at a different random walk
# trajectory; they should *not* falsely pass the Rhat <= 1.01 criterion.
matNS <- cbind(chainA, chainB)
rhatNS <- .Rhat(matNS)
status <- if (is.finite(rhatNS) && rhatNS > 1.01) "PASS" else "FAIL"
AddRow("NearNonStationary", "cross-chain Rhat", n, rhatNS, NA_real_, 1.01,
       status, "expect Rhat > 1.01 (non-mixed)")


# =========================================================================
# Scenario 5: Sticky chain (AR(1) at rho=0.99 / 0.999)
#   Same family as Scenario 1 but pulled out as an explicit "sticky" case
#   matching the brief. The brief mentions Pólya urn; there is no closed-
#   form IAT for the Pólya urn, so we use very-high-rho AR(1) as the
#   adversarial sticky proxy with the known IAT.
# =========================================================================
cat("[5/7] Sticky chain (rho=0.99)\n")

set.seed(5000L)
n <- if (QUICK) 2e4L else 1e5L
rho <- 0.99
ess_meas <- mean(replicate(if (QUICK) 2L else 5L,
                           .Ess(SimAR1(n, rho))))
truth <- n * (1 - rho) / (1 + rho)
tol   <- 0.30
status <- Classify(ess_meas, truth, tol)
AddRow("Sticky AR(1)", sprintf("rho=%g", rho), n, ess_meas, truth, tol, status,
       sprintf("ESS << N (true=%.1f, N=%d)", truth, n))


# =========================================================================
# Scenario 6: Constant chain — regression test for CONV-001 family.
#   A constant chain has undefined R-hat (between-chain var = 0, within = 0).
#   The implementation must return NA / Inf, NOT a finite ratio that would
#   pass a "Rhat <= 1.01" criterion.
# =========================================================================
cat("[6/7] Constant chain regression (CONV-001 family)\n")

# Single constant chain
constVec <- rep(3.14, 1000L)
ess_c <- .EssVector(constVec)
status <- if (is.na(ess_c)) "PASS" else "FAIL"
AddRow("Constant", ".EssVector single chain", 1000L, ess_c, NA_real_,
       NA_real_, status, "expect NA on constant input")

# Multi-chain constant (all chains identical constant)
matC <- cbind(constVec, constVec)
rhat_c <- .Rhat(matC)
status <- if (is.na(rhat_c) || !is.finite(rhat_c) || rhat_c > 1.05) "PASS" else "FAIL"
AddRow("Constant", ".Rhat two constant chains", 1000L, rhat_c, NA_real_,
       NA_real_, status, "expect NA/Inf, NOT < 1.01")

# Multi-chain with two different constants (each chain const, different value)
matD <- cbind(rep(1, 1000L), rep(2, 1000L))
rhat_d <- .Rhat(matD)
# within-chain var = 0 -> .RhatClassical returns NA. Acceptable answer is
# NA (current behaviour) — anything > 1.01 also OK, but NOT <= 1.01.
status <- if (is.na(rhat_d) || !is.finite(rhat_d) || rhat_d > 1.01) "PASS" else "FAIL"
AddRow("Constant", ".Rhat two diff-constant chains", 1000L, rhat_d, NA_real_,
       NA_real_, status, "different constants per chain -> must not pass")

# `.ComputeRhat` on a per-run constant-chain regression:
# Build a posterior-shaped per_run list and check ConvergenceDiagnostics
# behaves (no falsely-low Rhat).
fakePerRun_const <- list(
  list(samples = matrix(rep(0, 1000L), ncol = 1L,
                        dimnames = list(NULL, "rate_x"))),
  list(samples = matrix(rep(0, 1000L), ncol = 1L,
                        dimnames = list(NULL, "rate_x")))
)
rhatVec_const <- .ComputeRhat(fakePerRun_const, keyCols = 1L)
status <- if (all(is.na(rhatVec_const) | !is.finite(rhatVec_const) |
                  rhatVec_const > 1.01)) "PASS" else "FAIL"
AddRow("Constant", ".ComputeRhat const per_run", 1000L,
       as.numeric(rhatVec_const[1]), NA_real_, NA_real_, status,
       "regression: must not yield Rhat <= 1.01 on constants")


# =========================================================================
# Scenario 7: Unequal-length per-run chains — regression test for CONV-002.
#   Per the brief, .ComputeRhat does:
#     do.call(cbind, lapply(perRun, function(r) r$samples[, keyCols[j]]))
#   without length-equalisation. If chains differ in length, cbind() will
#   either recycle silently (R 3.x style) or emit "number of rows of result
#   is not a multiple of vector length" warning. Either is a regression.
# =========================================================================
cat("[7/7] Unequal-length per_run regression (CONV-002)\n")

set.seed(7000L)
nA <- 1000L
nB <- 600L
perRunUneq <- list(
  list(samples = matrix(rnorm(nA), ncol = 1L,
                        dimnames = list(NULL, "rate_x"))),
  list(samples = matrix(rnorm(nB), ncol = 1L,
                        dimnames = list(NULL, "rate_x")))
)

warningEmitted <- FALSE
errorEmitted   <- FALSE
warnMsg        <- ""
result <- tryCatch(
  withCallingHandlers(
    .ComputeRhat(perRunUneq, keyCols = 1L),
    warning = function(w) {
      warningEmitted <<- TRUE
      warnMsg <<- conditionMessage(w)
      invokeRestart("muffleWarning")
    }
  ),
  error = function(e) { errorEmitted <<- TRUE; NA }
)

# Expected (post-fix): no warning, no error, finite Rhat from harmonised
# (tail-equalised) chains. Current (buggy) state: warning emitted OR error.
status <- if (!warningEmitted && !errorEmitted && is.finite(result[1])) "PASS" else "FAIL"
AddRow("CONV-002", "unequal-length per_run", min(nA, nB),
       if (is.numeric(result) && length(result) >= 1L) as.numeric(result[1]) else NA_real_,
       NA_real_, NA_real_, status,
       sprintf("warn=%s err=%s msg=%s",
               warningEmitted, errorEmitted,
               substr(warnMsg, 1, 60)))


# =========================================================================
# Sanity probe (always run): does .Ess inflate above N when the safety cap
# fires?  White noise, small N. If measured > 1.5*N, raise WARN (potential
# ESS-CAP-001 finding) but don't FAIL outright (cap is documented).
# =========================================================================
cat("[probe] ESS safety-cap on white noise (small N)\n")

set.seed(9000L)
n <- 500L
ess_probe <- mean(replicate(if (QUICK) 5L else 30L,
                            .Ess(rnorm(n))))
status <- if (is.finite(ess_probe) && ess_probe <= 1.5 * n) "PASS" else "WARN"
AddRow("ESS cap probe", "white noise small N", n, ess_probe, n, 0.5,
       status, "watch for ESS > 1.5*N from safety cap")


# ---- Write outputs -------------------------------------------------------

cat("\nWriting results to ", OUTDIR, "\n", sep = "")

saveRDS(results, file = file.path(OUTDIR, "summary.rds"))
write.csv(results, file = file.path(OUTDIR, "results.csv"), row.names = FALSE)

# verdict.txt: stable, grep-friendly.
overall <- if (any(results$status == "FAIL")) {
  "FAIL"
} else if (any(results$status == "WARN")) {
  "WARN"
} else {
  "PASS"
}

verdictLines <- c(
  sprintf("ESS/R-hat stress test — %s mode", if (QUICK) "QUICK" else "FULL"),
  sprintf("run date: %s", Sys.time()),
  sprintf("nScenarios: %d", nrow(results)),
  "",
  sprintf("%-22s  %-32s  %8s  %12s  %12s  %9s  %6s",
          "scenario", "parameter", "n", "measured", "truth",
          "tolerance", "status"),
  strrep("-", 110)
)
for (i in seq_len(nrow(results))) {
  verdictLines <- c(verdictLines,
    sprintf("%-22s  %-32s  %8d  %12.4g  %12.4g  %9.3g  %6s",
            results$scenario[i],
            substr(results$parameter[i], 1, 32),
            results$n[i],
            results$measured[i],
            results$truth[i],
            results$tolerance[i],
            results$status[i]))
}
verdictLines <- c(verdictLines,
  strrep("-", 110),
  sprintf("OVERALL: %s", overall))

writeLines(verdictLines, file.path(OUTDIR, "verdict.txt"))

cat(paste(verdictLines, collapse = "\n"), "\n", sep = "")
cat("\nOverall verdict:", overall, "\n")

invisible(NULL)
