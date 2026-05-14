# truth-init-diag.R — analyse truth-init aware chain.
# Goal: report the trajectory of phi, pi0, theta_1, tree_length, log_lik,
# log_post across the chain so we can see whether truth holds or drifts.
logFile <- "inst/simulations/ecology/truth-init-results/truth-init.log"
log <- read.table(logFile, sep = "\t", header = TRUE, comment.char = "#",
                  check.names = FALSE)
cat("Rows:", nrow(log), "  Columns:", ncol(log), "\n")

keepCols <- intersect(c("Sample", "log_posterior", "log_likelihood",
                        "tree_length", "rate_loss", "rate_log_sd",
                        "rate_neo", "p", "phi", "pi0", "theta_1",
                        "swap_cold"),
                      colnames(log))
subDf <- log[, keepCols]
cat("\n=== First 5 samples ===\n"); print(round(head(subDf, 5), 3))
cat("\n=== Last 5 samples ===\n");  print(round(tail(subDf, 5), 3))

# Trajectory milestones: first sample, sample at 25%, 50%, 75%, last
n <- nrow(subDf)
if (n >= 5) {
  idx <- unique(c(1, floor(n * c(0.10, 0.25, 0.5, 0.75)), n))
  cat("\n=== Milestone snapshots ===\n")
  ms <- subDf[idx, , drop = FALSE]
  rownames(ms) <- paste0(round(100 * (idx - 1) / max(1, n - 1)), "%")
  print(round(ms, 3))
}

cat("\n=== Median of last 50% ===\n")
keep <- subDf[seq.int(ceiling(n / 2), n), , drop = FALSE]
print(round(apply(keep[, -1], 2, median), 3))

cat("\n=== Truth: phi=4, pi0=0.75, theta_1=1, TL=", 13.5,
    ", rate_loss=1, rate_neo=1 ===\n", sep = "")
