# 7-way CID comparison, adding mk_k15 to the 6-way arms.
# Tests whether the flexibility ramp continues past k=9 or plateaus.

SUMMARY_DIR <- "dev/pilots/2026-05-12-prior-validation/summary"
OUT_DIR     <- "dev/pilots/2026-05-12-prior-validation/analysis"

ARMS <- c("mk", "mk_kp1", "mk_kp2", "mk_k9", "mk_k15", "mkp_eg", "mkp_geo")
COLS <- c(mk      = "#7f7f7f",
          mk_kp1  = "#e6194B",
          mk_kp2  = "#f58231",
          mk_k9   = "#4363d8",
          mk_k15  = "#1f3a8f",
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

saveRDS(df, file.path(OUT_DIR, "cid_seven_prior.rds"))

cat("\n=== Marginal task-mean CID by arm (all available tasks) ===\n")
print(aggregate(mean_cid ~ arm, data = df,
                FUN = function(v) c(n = length(v),
                                    mean = mean(v),
                                    median = median(v),
                                    sd = sd(v),
                                    se = sd(v) / sqrt(length(v)))))

# ---- Paired apples-to-apples ---------------------------------------------
# 1. All-arms common subset (constrained by mkp_geo's 26 tasks)
tasks_all <- Reduce(intersect, split(df$tag, df$arm))
cat(sprintf("\nTasks common to all %d arms: %d\n",
            length(ARMS), length(tasks_all)))

# 2. mk_k9 vs mk_k15 head-to-head (the new comparison)
fixedK_arms <- c("mk_k9", "mk_k15")
fixedK_common <- Reduce(intersect, split(df$tag[df$arm %in% fixedK_arms],
                                         droplevels(df$arm[df$arm %in% fixedK_arms])))
cat(sprintf("Tasks common to mk_k9 + mk_k15: %d\n", length(fixedK_common)))

paired_fk <- df[df$arm %in% fixedK_arms & df$tag %in% fixedK_common, ]
wide_fk <- reshape(paired_fk[, c("tag", "arm", "mean_cid")],
                   idvar = "tag", timevar = "arm", direction = "wide")
names(wide_fk) <- sub("^mean_cid\\.", "", names(wide_fk))

cat("\n=== mk_k9 vs mk_k15 (paired) ===\n")
diff_fk <- wide_fk$mk_k15 - wide_fk$mk_k9
cat(sprintf("  mean(k15 - k9)  = %+.5f\n", mean(diff_fk, na.rm = TRUE)))
cat(sprintf("  median(k15 - k9) = %+.5f\n", median(diff_fk, na.rm = TRUE)))
cat(sprintf("  k15 wins (k15 < k9): %d / %d\n",
            sum(diff_fk < 0, na.rm = TRUE), sum(!is.na(diff_fk))))
cat(sprintf("  k15 loses:           %d\n", sum(diff_fk > 0, na.rm = TRUE)))
cat(sprintf("  ties:                %d\n", sum(diff_fk == 0, na.rm = TRUE)))
sg <- tryCatch(binom.test(sum(diff_fk < 0, na.rm = TRUE),
                          sum(!is.na(diff_fk)))$p.value,
               error = function(e) NA)
cat(sprintf("  binomial sign-test p = %.3g\n", sg))

# ---- Per-arm summary on the all-arms common subset -----------------------
paired_all <- df[df$tag %in% tasks_all, ]
wide_all <- reshape(paired_all[, c("tag", "arm", "mean_cid")],
                    idvar = "tag", timevar = "arm", direction = "wide")
names(wide_all) <- sub("^mean_cid\\.", "", names(wide_all))
cat(sprintf("\n=== All-arm paired (n=%d tasks, mkp_geo-limited) ===\n",
            length(tasks_all)))
common_summary <- data.frame(
  arm    = ARMS,
  mean   = sapply(ARMS, function(a) mean(wide_all[[a]], na.rm = TRUE)),
  median = sapply(ARMS, function(a) median(wide_all[[a]], na.rm = TRUE))
)
print(common_summary, row.names = FALSE)

# ---- Headline plot -------------------------------------------------------
png(file.path(OUT_DIR, "cid_seven_prior.png"),
    width = 1600, height = 600, res = 110)
op <- par(mfrow = c(1, 2), mar = c(4, 4.5, 3, 1))

boxplot(mean_cid ~ arm, data = df,
        col = COLS[ARMS], border = "grey30",
        ylab = "Task-mean CID to truth",
        main = "All available tasks per arm",
        las = 1)
mtext(side = 1, line = 2.5,
      text = sprintf("n_tasks: %s",
                     paste(paste0(ARMS, "=", table(df$arm)),
                           collapse = ", ")),
      cex = 0.7, col = "grey40")

# mk_k9 vs mk_k15 paired-diff histogram
hist(diff_fk, breaks = 30,
     col = "#1f3a8f", border = "grey30",
     xlab = "CID(mk_k15) - CID(mk_k9)",
     main = sprintf("k15 vs k9 (n=%d common tasks)", length(fixedK_common)))
abline(v = 0, lty = 2, col = "grey40")
abline(v = mean(diff_fk, na.rm = TRUE), col = "red", lwd = 2)

par(op); dev.off()
cat("\nSaved: cid_seven_prior.png and cid_seven_prior.rds\n")
