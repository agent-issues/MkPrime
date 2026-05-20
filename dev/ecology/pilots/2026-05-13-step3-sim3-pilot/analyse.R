# analyse.R -------------------------------------------------------------------
# Post-processing for Step 3 pilot: reads the log produced by run-pilot.R
# and computes the five diagnostic metrics vs pre-fix baseline.
#
# Pre-fix headline numbers (from 200k-iter chain, n=194 post-burnin):
#   tree_length : median=22, IQR=[11.06, 124.78], max=7496
#   gamma_e (theta=0.5 proxy): median=1.130, IQR=[1.006, 2.174]
#   tl/gamma_e  : median=14.49, IQR=[10.43, 41.21], max=4009.94
#   cor(log tl, log gamma_e) = +0.629
#
# Truth: tree_length ≈ 9 (phi=4, nEco=60 ecomorphic, Goldilocks)

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
  library("coda")
})

PILOT_DIR <- "dev/pilots/2026-05-13-step3-sim3-pilot"
logFile   <- file.path(PILOT_DIR, "sim3-pilot.log")

if (!file.exists(logFile)) stop("Log file not found: ", logFile)

samp_all <- ReadMkLog(logFile)
cat("Total rows in log:", nrow(samp_all), "\n")
cat("Columns:", paste(colnames(samp_all), collapse = ", "), "\n\n")

# Discard first 25% as burnin (mirrors pre-fix convention)
n <- nrow(samp_all)
discard_n <- ceiling(n / 4)
samp <- samp_all[seq.int(discard_n + 1L, n), , drop = FALSE]
cat("Retained (post 25% discard):", nrow(samp), "rows\n\n")

# Helper: quantile shortcut
qs <- function(x, p) stats::quantile(x, p, names = FALSE)

# ---- Metric 1: tree_length trace summary ------------------------------------
tl <- samp[, "tree_length"]
cat("=== 1. tree_length ===\n")
cat(sprintf("  median  = %.3f   (pre-fix: 22.0,  truth ≈ 9)\n",   median(tl)))
cat(sprintf("  IQR     = [%.3f, %.3f]  (pre-fix: [11.06, 124.78])\n",
            qs(tl, 0.25), qs(tl, 0.75)))
cat(sprintf("  max     = %.3f   (pre-fix: 7496)\n\n",              max(tl)))

# ---- Metric 2: Posterior covariance -----------------------------------------
phi_col  <- if ("phi" %in% colnames(samp)) samp[, "phi"] else NULL
pi0_col  <- if ("pi0" %in% colnames(samp)) samp[, "pi0"] else NULL
th_col   <- if ("theta_1" %in% colnames(samp)) samp[, "theta_1"] else NULL

cat("=== 2. Posterior covariance (log/logit-transformed) ===\n")
if (!is.null(phi_col) && !is.null(pi0_col) && !is.null(th_col)) {
  lp   <- log(phi_col)
  lpi0 <- qlogis(pmin(pmax(pi0_col, 1e-8), 1 - 1e-8))
  lth  <- qlogis(pmin(pmax(th_col, 1e-8), 1 - 1e-8))
  ltl  <- log(tl)
  covMat <- cov(cbind(log_phi = lp, logit_pi0 = lpi0,
                      logit_theta1 = lth, log_tree_length = ltl))
  print(round(covMat, 4))
  cat("\n")
} else {
  cat("  WARNING: phi/pi0/theta_1 columns not all present.\n")
  cat("  Available cols:", paste(colnames(samp), collapse = ", "), "\n\n")
}

# ---- Metric 3: ESS ----------------------------------------------------------
cat("=== 3. ESS ===\n")
ess_cols <- c("phi", "pi0", "theta_1", "tree_length")
for (nm in ess_cols) {
  if (nm %in% colnames(samp)) {
    ess_val <- coda::effectiveSize(samp[, nm])
    cat(sprintf("  %-14s ESS = %.1f  (n=%d)\n", nm, ess_val, nrow(samp)))
  } else {
    cat(sprintf("  %-14s MISSING\n", nm))
  }
}
cat("\n")

# ---- Metric 4: gamma_e per sample using logged theta ------------------------
cat("=== 4. gamma_e (using actual logged theta) ===\n")
phi_truth <- 4  # from simulation config
if (!is.null(phi_col) && !is.null(pi0_col) && !is.null(th_col)) {
  # gamma_e = pi0 + (1 - pi0) * [theta * phi + (1 - theta) / phi]
  gE_actual <- pi0_col + (1 - pi0_col) * (th_col * phi_col + (1 - th_col) / phi_col)
  cat(sprintf("  median  = %.3f   (pre-fix proxy: 1.130)\n",  median(gE_actual)))
  cat(sprintf("  IQR     = [%.3f, %.3f]  (pre-fix: [1.006, 2.174])\n",
              qs(gE_actual, 0.25), qs(gE_actual, 0.75)))
  cat(sprintf("  max     = %.3f\n\n", max(gE_actual)))

  # tl / gamma_e
  tl_over_gE <- tl / gE_actual
  cat("  tl / gamma_e:\n")
  cat(sprintf("    median = %.3f  (pre-fix: 14.49, truth ≈ 9)\n", median(tl_over_gE)))
  cat(sprintf("    IQR    = [%.3f, %.3f]  (pre-fix: [10.43, 41.21])\n",
              qs(tl_over_gE, 0.25), qs(tl_over_gE, 0.75)))
  cat(sprintf("    max    = %.3f  (pre-fix: 4010)\n\n", max(tl_over_gE)))

  # Also the theta=0.5 proxy for apples-to-apples comparison
  theta_proxy <- 0.5
  gE_proxy <- pi0_col + (1 - pi0_col) * (theta_proxy * phi_col + (1 - theta_proxy) / phi_col)
  cat(sprintf("  gamma_e (theta=0.5 proxy, for pre-fix comparison):\n"))
  cat(sprintf("    median = %.3f\n", median(gE_proxy)))
} else {
  cat("  Cannot compute: ecology columns missing.\n\n")
}
cat("\n")

# ---- Metric 5: cor(log tl, log gamma_e) ------------------------------------
cat("=== 5. cor(log tl, log gamma_e) ===\n")
if (!is.null(phi_col) && !is.null(pi0_col) && !is.null(th_col)) {
  r_actual <- cor(log(tl), log(gE_actual))
  r_proxy  <- cor(log(tl), log(gE_proxy))
  cat(sprintf("  cor(log tl, log gamma_e [actual theta])  = %.4f  (pre-fix: +0.629)\n", r_actual))
  cat(sprintf("  cor(log tl, log gamma_e [theta=0.5 proxy]) = %.4f\n\n", r_proxy))
} else {
  cat("  Cannot compute: ecology columns missing.\n\n")
}

# ---- Decision rule ----------------------------------------------------------
cat("=== Decision rule ===\n")
if (!is.null(phi_col) && !is.null(pi0_col) && !is.null(th_col)) {
  r <- cor(log(tl), log(gE_actual))
  tl_max <- max(tl)
  tl_med <- median(tl)
  cat(sprintf("  tree_length max=%.1f, median=%.3f\n", tl_max, tl_med))
  cat(sprintf("  cor(log tl, log gamma_e)=%.4f\n", r))
  if (tl_max < 500 && abs(r) < 0.4) {
    cat("  --> BRANCH: Tail GONE, cov benign. Ship simple Bactrians; proceed to multi-rep.\n")
  } else if (tl_max < 5000 && abs(r) >= 0.4) {
    cat("  --> BRANCH: Tail ATTENUATED but gamma-tl cov still > 0.4. Recommend Step 4 (v2-no-gamma).\n")
  } else if (tl_max >= 5000) {
    cat("  --> BRANCH: Tail UNCHANGED. Fundamental problem; reconsider model.\n")
  } else {
    cat("  --> BRANCH: INTERMEDIATE. max=", tl_max, " cor=", round(r, 4), ". Assess manually.\n")
  }
}
cat("\n")

# ---- phi / pi0 / theta summaries -------------------------------------------
cat("=== Hyperparameter summaries ===\n")
for (nm in c("phi", "pi0", "theta_1", "log_likelihood")) {
  if (nm %in% colnames(samp)) {
    x <- samp[, nm]
    cat(sprintf("  %-14s median=%.4f  IQR=[%.4f, %.4f]  CI95=[%.4f, %.4f]\n",
                nm, median(x), qs(x, 0.25), qs(x, 0.75), qs(x, 0.025), qs(x, 0.975)))
  }
}
cat("\n")

cat("Analysis complete.\n")
