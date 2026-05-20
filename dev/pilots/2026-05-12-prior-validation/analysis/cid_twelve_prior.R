# 12-arm CID comparison: the fixed-k ramp + Mk' prior-shape sweep + oracle.
#
# Arms grouped by family:
#   mk family (k = function of kObs):        mk, mk_kp1, mk_kp2, mk_k9, mk_k15, mk_k24, mk_k40
#   Mk' family (per-character learned k'):   mkp_eg, mkp_geo, mkp_highk, mkp_logs
#   Oracle (k = k_true per character):       mk_ktrue
#
# Hypotheses tested in this analysis:
#   H1: mk_ktrue sits below mk_k40 — oracle is the ceiling
#   H2: mkp_highk (strong high-k prior) matches mk_k40 — paradox is prior-driven
#   H3: mkp_logs (heavy 1/k tail) sits between mkp_geo and mk_k40
#
# Output: cid_twelve_prior.rds, .png, and a text report.

SUMMARY_DIR <- "dev/pilots/2026-05-12-prior-validation/summary"
OUT_DIR     <- "dev/pilots/2026-05-12-prior-validation/analysis"

ARMS <- c("mk", "mk_kp1", "mk_kp2", "mk_k9", "mk_k15", "mk_k24", "mk_k40",
          "mkp_eg", "mkp_geo", "mkp_highk", "mkp_logs", "mk_ktrue")

# Colour palette: greys/blues for the mk ramp, distinct for Mk', gold for oracle
COLS <- c(mk        = "#999999",
          mk_kp1    = "#cccccc",
          mk_kp2    = "#aaaaaa",
          mk_k9     = "#888888",
          mk_k15    = "#666666",
          mk_k24    = "#444444",
          mk_k40    = "#222222",
          mkp_eg    = "#3cb44b",
          mkp_geo   = "#911eb4",
          mkp_highk = "#e6194B",
          mkp_logs  = "#f58231",
          mk_ktrue  = "#daa520")

# Implicit k under each arm (for the ramp plot)
RAMP_K <- c(mk = 2.0, mk_kp1 = 3.0, mk_kp2 = 4.0,
            mk_k9 = 9, mk_k15 = 15, mk_k24 = 24, mk_k40 = 40,
            mkp_eg = 3.0, mkp_geo = 3.0, mkp_highk = 23, mkp_logs = 15,
            mk_ktrue = 3.78)  # approximate mean k_true across simulation

rows <- list()
for (arm in ARMS) {
  pattern <- paste0("^", arm, "_t[0-9]+_r[0-9]+\\.rds$")
  files <- list.files(SUMMARY_DIR, pattern, full.names = TRUE)
  cat(sprintf("%-12s %3d files\n", arm, length(files)))
  for (f in files) {
    x <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(x) || !length(x$cid)) next
    rows[[length(rows) + 1L]] <- data.frame(
      arm        = arm,
      tag        = x$tag,
      tree_idx   = x$tree_idx,
      rep_idx    = x$rep_idx,
      n_trees    = length(x$cid),
      mean_cid   = mean(x$cid, na.rm = TRUE),
      median_cid = median(x$cid, na.rm = TRUE),
      sd_cid     = sd(x$cid, na.rm = TRUE)
    )
  }
}
df <- do.call(rbind, rows)
df$arm <- factor(df$arm, levels = ARMS)
cat(sprintf("\nTotal rows: %d  (tasks per arm: %s)\n",
            nrow(df),
            paste(table(df$arm), collapse = "/")))

saveRDS(df, file.path(OUT_DIR, "cid_twelve_prior.rds"))

# ---- Marginal task-mean CID by arm ------------------------------------------
cat("\n=== Marginal task-mean CID by arm (all available tasks) ===\n")
agg <- aggregate(mean_cid ~ arm, data = df,
                 FUN = function(v) c(n      = length(v),
                                     mean   = mean(v),
                                     median = median(v),
                                     sd     = sd(v),
                                     se     = sd(v) / sqrt(length(v))))
print(agg)

# ---- Paired comparisons of interest -----------------------------------------
.paired <- function(arm_a, arm_b) {
  common <- Reduce(intersect, list(df$tag[df$arm == arm_a],
                                   df$tag[df$arm == arm_b]))
  if (length(common) == 0) return(NULL)
  paired <- df[df$arm %in% c(arm_a, arm_b) & df$tag %in% common, ]
  wide <- reshape(paired[, c("tag", "arm", "mean_cid")],
                  idvar = "tag", timevar = "arm", direction = "wide")
  names(wide) <- sub("^mean_cid\\.", "", names(wide))
  diff <- wide[[arm_b]] - wide[[arm_a]]
  cat(sprintf("\n=== %s vs %s (paired, n=%d) ===\n",
              arm_a, arm_b, length(common)))
  cat(sprintf("  mean(%s - %s)   = %+.5f\n", arm_b, arm_a,
              mean(diff, na.rm = TRUE)))
  cat(sprintf("  median(%s - %s) = %+.5f\n", arm_b, arm_a,
              median(diff, na.rm = TRUE)))
  cat(sprintf("  %s wins (CID lower): %d / %d\n", arm_b,
              sum(diff < 0, na.rm = TRUE), sum(!is.na(diff))))
  p <- tryCatch(binom.test(sum(diff < 0, na.rm = TRUE),
                           sum(!is.na(diff)))$p.value, error = function(e) NA)
  cat(sprintf("  binomial sign-test p = %.3g\n", p))
}

# H1: oracle test
.paired("mk_k40", "mk_ktrue")
# H2: high-k prior test
.paired("mkp_eg", "mkp_highk")
.paired("mk_k40", "mkp_highk")
# H3: heavy-tail prior test
.paired("mkp_geo", "mkp_logs")
.paired("mk_k40", "mkp_logs")
# Also compare oracle to all Mk' priors
.paired("mkp_highk", "mk_ktrue")

# ---- Headline plot ----------------------------------------------------------
png(file.path(OUT_DIR, "cid_twelve_prior.png"),
    width = 2400, height = 800, res = 110)
op <- par(mfrow = c(1, 2), mar = c(5, 4.5, 3, 1))

# Box plot, ordered by mean CID
arm_order <- names(sort(tapply(df$mean_cid, df$arm, mean, na.rm = TRUE),
                        decreasing = TRUE))
boxplot(mean_cid ~ factor(arm, levels = arm_order), data = df,
        col = COLS[arm_order], border = "grey30",
        ylab = "Task-mean CID to truth",
        xlab = "",
        main = "All 12 arms (ordered by mean CID)",
        las = 2, cex.axis = 0.85)

# "Effective k" scatter
arm_means <- tapply(df$mean_cid, df$arm, mean, na.rm = TRUE)
plot(RAMP_K[ARMS], arm_means[ARMS],
     pch = 19, cex = 1.5, col = COLS[ARMS],
     xlab = "Implicit / mean k", ylab = "Mean CID to truth",
     main = "Family-by-family CID vs k",
     xlim = c(0, 42), las = 1)
text(RAMP_K[ARMS], arm_means[ARMS],
     labels = ARMS, pos = 3, cex = 0.7, col = COLS[ARMS])
legend("topright", legend = c("mk family", "Mk' family", "oracle"),
       col = c("#444444", "#3cb44b", "#daa520"),
       pch = 19, bty = "n", cex = 0.9)

par(op); dev.off()
cat("\nSaved: cid_twelve_prior.png and cid_twelve_prior.rds\n")
