# Step 1 + Step 3: decompose log_posterior into log_lik + log_prior on
# Hamilton 17140607/17141180 tail samples, compute swap rate, and check
# per-character k' trajectories.

LOCAL <- "inst/hamilton/mkp-study-17140607-17141180/mkp_eg_sample"
dirs  <- list.dirs(LOCAL, recursive = FALSE)
cat(sprintf("Analysing %d subdirs.\n\n", length(dirs)))

summary_rows <- list()
for (d in dirs) {
  tag <- basename(d)
  f1 <- file.path(d, "mkp_eg_run_1.log")
  if (!file.exists(f1) || file.info(f1)$size < 1000) next

  # Read TSV; first row is column headers, # lines were stripped on remote.
  trace <- read.table(f1, header = TRUE, sep = "\t", comment.char = "#",
                       check.names = FALSE)
  if (!nrow(trace)) next

  lp <- trace$log_posterior
  ll <- trace$log_likelihood
  lpr <- lp - ll

  # k' columns
  kp_cols <- grep("^kPrime_", colnames(trace), value = TRUE)
  kp_mat <- as.matrix(trace[, kp_cols, drop = FALSE])

  # swap_cold (per-sample swap count): rate = mean / (attempts per sample)
  # We don't know attempts per sample directly — report mean count as
  # rough swap activity.
  swap_mean <- mean(trace$swap_cold)

  # topology hash autocorrelation: number of distinct hashes seen
  topo_uniq <- length(unique(trace$topo_hash))

  # Largest k' seen across all chars; max range per char
  max_k <- max(kp_mat)
  per_char_range <- apply(kp_mat, 2, function(x) diff(range(x)))
  q90_range <- quantile(per_char_range, 0.9)
  q50_range <- quantile(per_char_range, 0.5)

  summary_rows[[tag]] <- data.frame(
    tag      = tag,
    n        = nrow(trace),
    logL_sd  = round(sd(ll), 2),
    logPr_sd = round(sd(lpr), 2),
    logP_sd  = round(sd(lp), 2),
    cor_L_Pr = round(cor(ll, lpr), 3),
    swap_avg = round(swap_mean, 3),
    n_topo   = topo_uniq,
    p_min    = round(min(trace$p), 3),
    p_max    = round(max(trace$p), 3),
    max_k    = max_k,
    k_range_med = q50_range,
    k_range_q90 = q90_range
  )
}

all_sum <- do.call(rbind, summary_rows)
rownames(all_sum) <- NULL
cat("=== Per-task summary (run_1.log tail = last 50k samples) ===\n")
print(all_sum, row.names = FALSE)

# Aggregate sense-check
cat(sprintf(
  "\nAcross %d tasks, logPrior sd ranges %.1f-%.1f (median %.1f)\n",
  nrow(all_sum), min(all_sum$logPr_sd), max(all_sum$logPr_sd), median(all_sum$logPr_sd)
))
cat(sprintf("Across %d tasks, logLik   sd ranges %.1f-%.1f (median %.1f)\n",
  nrow(all_sum), min(all_sum$logL_sd), max(all_sum$logL_sd), median(all_sum$logL_sd)
))
cat(sprintf("Ratio logPr_sd / logL_sd: median %.2f\n",
            median(all_sum$logPr_sd / all_sum$logL_sd)))
cat(sprintf("cor(logL, logPrior): range %.2f to %.2f, median %.2f\n",
            min(all_sum$cor_L_Pr), max(all_sum$cor_L_Pr),
            median(all_sum$cor_L_Pr)))
cat(sprintf("swap_cold mean per sample: range %.3f to %.3f, median %.3f\n",
            min(all_sum$swap_avg), max(all_sum$swap_avg),
            median(all_sum$swap_avg)))

saveRDS(all_sum, "data-raw/step1_hamilton_summary.rds")

# Trajectory plot for one task — kPrime drift for the worst-range char.
d <- dirs[1]
f1 <- file.path(d, "mkp_eg_run_1.log")
trace <- read.table(f1, header = TRUE, sep = "\t", comment.char = "#",
                     check.names = FALSE)
kp_cols <- grep("^kPrime_", colnames(trace), value = TRUE)
kp_mat <- as.matrix(trace[, kp_cols, drop = FALSE])
per_char_range <- apply(kp_mat, 2, function(x) diff(range(x)))
worst_chars <- order(per_char_range, decreasing = TRUE)[1:5]
cat(sprintf("\nIn %s, chars with widest k' range across last 50k samples:\n",
            basename(d)))
for (i in worst_chars) {
  cat(sprintf("  %-12s  range = %d (min=%d, max=%d), median=%.0f\n",
              kp_cols[i], per_char_range[i],
              min(kp_mat[, i]), max(kp_mat[, i]),
              median(kp_mat[, i])))
}
