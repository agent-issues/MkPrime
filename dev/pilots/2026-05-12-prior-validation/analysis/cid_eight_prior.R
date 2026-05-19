# 8-way CID comparison, adding mk_k24 to the 7-way arms.
# Tests whether the flexibility ramp continues past k=15 or begins to plateau.

SUMMARY_DIR <- "dev/pilots/2026-05-12-prior-validation/summary"
OUT_DIR     <- "dev/pilots/2026-05-12-prior-validation/analysis"

ARMS <- c("mk", "mk_kp1", "mk_kp2", "mk_k9", "mk_k15", "mk_k24", "mkp_eg", "mkp_geo")
COLS <- c(mk      = "#7f7f7f",
          mk_kp1  = "#e6194B",
          mk_kp2  = "#f58231",
          mk_k9   = "#4363d8",
          mk_k15  = "#1f3a8f",
          mk_k24  = "#00aacc",
          mkp_eg  = "#3cb44b",
          mkp_geo = "#911eb4")

rows <- list()
for (arm in ARMS) {
  pattern <- paste0("^", arm, "_t[0-9]+_r[0-9]+\\.rds$")
  files <- list.files(SUMMARY_DIR, pattern, full.names = TRUE)
  cat(sprintf("%-10s %d files\n", arm, length(files)))
  for (f in files) {
    x <- readRDS(f)
    if (!length(x$cid)) next
    rows[[length(rows) + 1L]] <- data.frame(
      arm       = arm,
      tag       = x$tag,
      tree_idx  = x$tree_idx,
      rep_idx   = x$rep_idx,
      n_trees   = length(x$cid),
      mean_cid  = mean(x$cid, na.rm = TRUE),
      median_cid = median(x$cid, na.rm = TRUE),
      sd_cid    = sd(x$cid, na.rm = TRUE)
    )
  }
}
df <- do.call(rbind, rows)
df$arm <- factor(df$arm, levels = ARMS)
cat(sprintf("\nTotal rows: %d  (tasks per arm: %s)\n",
            nrow(df),
            paste(table(df$arm), collapse = "/")))

saveRDS(df, file.path(OUT_DIR, "cid_eight_prior.rds"))

cat("\n=== Marginal task-mean CID by arm (all available tasks) ===\n")
agg <- aggregate(mean_cid ~ arm, data = df,
                 FUN = function(v) c(n      = length(v),
                                     mean   = mean(v),
                                     median = median(v),
                                     sd     = sd(v),
                                     se     = sd(v) / sqrt(length(v))))
print(agg)

# ---- Paired mk_k15 vs mk_k24 ------------------------------------------------
fixedK_arms   <- c("mk_k15", "mk_k24")
fixedK_common <- Reduce(intersect,
                        split(df$tag[df$arm %in% fixedK_arms],
                              droplevels(df$arm[df$arm %in% fixedK_arms])))
cat(sprintf("\nTasks common to mk_k15 + mk_k24: %d\n", length(fixedK_common)))

paired_fk <- df[df$arm %in% fixedK_arms & df$tag %in% fixedK_common, ]
wide_fk   <- reshape(paired_fk[, c("tag", "arm", "mean_cid")],
                     idvar = "tag", timevar = "arm", direction = "wide")
names(wide_fk) <- sub("^mean_cid\\.", "", names(wide_fk))

cat("\n=== mk_k15 vs mk_k24 (paired) ===\n")
diff_fk <- wide_fk$mk_k24 - wide_fk$mk_k15
cat(sprintf("  mean(k24 - k15)  = %+.5f\n", mean(diff_fk, na.rm = TRUE)))
cat(sprintf("  median(k24 - k15) = %+.5f\n", median(diff_fk, na.rm = TRUE)))
cat(sprintf("  k24 wins (k24 < k15): %d / %d\n",
            sum(diff_fk < 0, na.rm = TRUE), sum(!is.na(diff_fk))))
cat(sprintf("  k24 loses:            %d\n", sum(diff_fk > 0, na.rm = TRUE)))
cat(sprintf("  ties:                 %d\n", sum(diff_fk == 0, na.rm = TRUE)))
sg <- tryCatch(binom.test(sum(diff_fk < 0, na.rm = TRUE),
                          sum(!is.na(diff_fk)))$p.value,
               error = function(e) NA)
cat(sprintf("  binomial sign-test p = %.3g\n", sg))

# ---- Monotonicity ramp summary ----------------------------------------------
ramp_arms <- c("mk_kp1", "mk_kp2", "mk_k9", "mk_k15", "mk_k24")
ramp_common <- Reduce(intersect,
                      split(df$tag[df$arm %in% ramp_arms],
                            droplevels(df$arm[df$arm %in% ramp_arms])))
cat(sprintf("\n=== Ramp arms (n=%d common tasks) ===\n", length(ramp_common)))
ramp_df <- df[df$arm %in% ramp_arms & df$tag %in% ramp_common, ]
ramp_means <- tapply(ramp_df$mean_cid, ramp_df$arm, mean, na.rm = TRUE)
for (a in ramp_arms) cat(sprintf("  %-10s  %.4f\n", a, ramp_means[a]))
cat(sprintf("\n  Deltas:\n"))
cat(sprintf("  kp1->kp2 = %+.4f\n", ramp_means["mk_kp2"] - ramp_means["mk_kp1"]))
cat(sprintf("  kp2->k9  = %+.4f\n", ramp_means["mk_k9"]  - ramp_means["mk_kp2"]))
cat(sprintf("  k9->k15  = %+.4f\n", ramp_means["mk_k15"] - ramp_means["mk_k9"]))
cat(sprintf("  k15->k24 = %+.4f\n", ramp_means["mk_k24"] - ramp_means["mk_k15"]))

# ---- All-arm paired (mkp_geo-limited common subset) -------------------------
tasks_all <- Reduce(intersect, split(df$tag, df$arm))
cat(sprintf("\n=== All-arm paired (n=%d tasks, mkp_geo-limited) ===\n",
            length(tasks_all)))
paired_all <- df[df$tag %in% tasks_all, ]
wide_all   <- reshape(paired_all[, c("tag", "arm", "mean_cid")],
                      idvar = "tag", timevar = "arm", direction = "wide")
names(wide_all) <- sub("^mean_cid\\.", "", names(wide_all))
common_summary <- data.frame(
  arm    = ARMS,
  mean   = sapply(ARMS, function(a) mean(wide_all[[a]], na.rm = TRUE)),
  median = sapply(ARMS, function(a) median(wide_all[[a]], na.rm = TRUE))
)
print(common_summary, row.names = FALSE)

# ---- Headline plot ----------------------------------------------------------
png(file.path(OUT_DIR, "cid_eight_prior.png"),
    width = 2000, height = 700, res = 110)
op <- par(mfrow = c(1, 3), mar = c(4, 4.5, 3, 1))

boxplot(mean_cid ~ arm, data = df,
        col = COLS[ARMS], border = "grey30",
        ylab = "Task-mean CID to truth",
        main = "All available tasks per arm",
        las = 2, cex.axis = 0.85)
mtext(side = 1, line = 3.5,
      text = sprintf("n_tasks: %s",
                     paste(paste0(ARMS, "=", table(df$arm)),
                           collapse = ", ")),
      cex = 0.6, col = "grey40")

# k15 vs k24 paired-diff histogram
hist(diff_fk, breaks = 30,
     col = "#00aacc", border = "grey30",
     xlab = "CID(mk_k24) - CID(mk_k15)",
     main = sprintf("k24 vs k15 (n=%d common tasks)", length(fixedK_common)))
abline(v = 0, lty = 2, col = "grey40")
abline(v = mean(diff_fk, na.rm = TRUE), col = "red", lwd = 2)

# Ramp plot
ramp_k <- c(mk_kp1 = 1, mk_kp2 = 2, mk_k9 = 9, mk_k15 = 15, mk_k24 = 24)
plot(ramp_k[ramp_arms], ramp_means[ramp_arms],
     type = "b", pch = 19, col = COLS[ramp_arms],
     xlab = "k (fixed upper bound)", ylab = "Mean CID to truth",
     main = sprintf("Flexibility ramp (n=%d tasks)", length(ramp_common)),
     xlim = c(0, 26), las = 1)
text(ramp_k[ramp_arms], ramp_means[ramp_arms],
     labels = sprintf("%.4f", ramp_means[ramp_arms]),
     pos = 3, cex = 0.7)

par(op); dev.off()
cat("\nSaved: cid_eight_prior.png and cid_eight_prior.rds\n")
