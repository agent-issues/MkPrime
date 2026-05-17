.libPaths(c("/nobackup/pjjg18/mkp-sim3-multirep-v3/lib", .libPaths()))
library(MkPrime)
res <- readRDS("/nobackup/pjjg18/mkp-sim3-phi2/results/aware-result.rds")
res <- RelabelEcology(res)

# Post-relabel samples (in-memory object, already relabeled)
res$samples <- ReadMkLog(res$logFile)
s <- as.data.frame(res$samples)

cat("=== POST-relabel (from RDS after RelabelEcology) ===\n")
# phi in the RDS samples — check if already relabeled
if ("phi" %in% names(s)) {
  cat("phi (from log): min=", round(min(s$phi),3),
      " med=", round(median(s$phi),3),
      " max=", round(max(s$phi),3), "\n")
}

# Also check the trees object for any phi stored there
# The relabeled samples are in res$samples after ReadMkLog+RelabelEcology
# RelabelEcology operates on res$samples matrix
res2 <- readRDS("/nobackup/pjjg18/mkp-sim3-phi2/results/aware-result.rds")
res2 <- RelabelEcology(res2)
# res2$samples is the in-memory matrix (empty in streaming mode)
# but res2 has been relabeled — check what's available
cat("\nnrow(res2$samples):", nrow(as.data.frame(res2$samples)), "\n")
cat("res2$nSamples:", res2$nSamples, "\n")

# RelabelEcology on streaming result — does it relabel the log file?
# Check run_phi2.R: res <- RelabelEcology(res) then saveRDS
# The log file is NOT modified by RelabelEcology (it relabels in-memory samples)
# So: read log -> raw phi; read RDS samples -> also raw phi (streaming)
# RelabelEcology only flips the in-memory matrix which is empty in streaming mode

# The correct approach: read log, then apply relabel logic manually
s_raw <- as.data.frame(ReadMkLog(res$logFile))
if ("phi" %in% names(s_raw)) {
  phi_raw <- s_raw$phi
  phi_relabeled <- ifelse(phi_raw < 1, 1/phi_raw, phi_raw)
  cat("\n=== Manually relabeled phi (from log, flip phi<1 -> 1/phi) ===\n")
  cat("phi_raw:       med=", round(median(phi_raw),3),
      " mean=", round(mean(phi_raw),3), "\n")
  cat("phi_relabeled: med=", round(median(phi_relabeled),3),
      " mean=", round(mean(phi_relabeled),3),
      " min=", round(min(phi_relabeled),3),
      " max=", round(max(phi_relabeled),3), "\n")
  cat("Fraction phi_raw < 1:", round(mean(phi_raw < 1),3), "\n")
  cat("Truth phi=2\n")
}
