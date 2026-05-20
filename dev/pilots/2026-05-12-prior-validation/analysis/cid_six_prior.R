# 6-way CID comparison across prior arms.
# For each task (tree x rep), each arm's CID-to-truth is averaged across
# its thinned posterior (1000 trees per RDS). We compare arms on:
#   1. Marginal task-mean CID distribution (all available tasks per arm)
#   2. Paired apples-to-apples: tasks that all arms ran (rep_01 only,
#      since mkp_geo was only re-run for rep_01 with thin=100)
#
# Arms:
#   mk       — kObs plug-in (no hidden states)
#   mk_kp1   — kObs+1 plug-in
#   mk_kp2   — kObs+2 plug-in
#   mk_k9    — k=9 across all variable chars
#   mkp_eg   — Mk' with empirical_geometric prior on k'
#   mkp_geo  — Mk' with plain geometric prior on k' (rep_01 only)

SUMMARY_DIR <- "dev/pilots/2026-05-12-prior-validation/summary"
OUT_DIR     <- "dev/pilots/2026-05-12-prior-validation/analysis"

ARMS <- c("mk", "mk_kp1", "mk_kp2", "mk_k9", "mkp_eg", "mkp_geo")
COLS <- c(mk      = "#7f7f7f",
          mk_kp1  = "#e6194B",
          mk_kp2  = "#f58231",
          mk_k9   = "#4363d8",
          mkp_eg  = "#3cb44b",
          mkp_geo = "#911eb4")

# ---- Load all summaries -------------------------------------------------
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

saveRDS(df, file.path(OUT_DIR, "cid_six_prior.rds"))

# ---- 1. Marginal per-arm summary ----------------------------------------
cat("\n=== Marginal task-mean CID by arm (all available tasks) ===\n")
print(aggregate(mean_cid ~ arm, data = df,
                FUN = function(v) c(n = length(v),
                                    mean = mean(v),
                                    median = median(v),
                                    sd = sd(v),
                                    se = sd(v) / sqrt(length(v)))))

# ---- 2. Paired apples-to-apples (rep_01 only) ---------------------------
tasks_all <- Reduce(intersect, split(df$tag, df$arm))
cat(sprintf("\nTasks common to all %d arms: %d\n",
            length(ARMS), length(tasks_all)))

paired <- df[df$tag %in% tasks_all, ]
wide <- reshape(paired[, c("tag", "arm", "mean_cid")],
                idvar = "tag", timevar = "arm", direction = "wide")
names(wide) <- sub("^mean_cid\\.", "", names(wide))

cat("\n=== Per-task CID on common subset (head) ===\n"); print(head(wide))

cat("\n=== Paired summary on common subset ===\n")
common_summary <- data.frame(
  arm    = ARMS,
  n      = colSums(!is.na(wide[, ARMS])),
  mean   = sapply(ARMS, function(a) mean(wide[[a]], na.rm = TRUE)),
  median = sapply(ARMS, function(a) median(wide[[a]], na.rm = TRUE))
)
print(common_summary)

cat("\n=== Paired wins vs mk_kp1 (n tasks where arm beats mk_kp1) ===\n")
for (a in setdiff(ARMS, "mk_kp1")) {
  diff <- wide[[a]] - wide$mk_kp1
  cat(sprintf("  %-8s  mean_diff=%+.4f  wins=%2d/%d  loss=%2d  ties=%d\n",
              a, mean(diff, na.rm = TRUE),
              sum(diff < 0, na.rm = TRUE), sum(!is.na(diff)),
              sum(diff > 0, na.rm = TRUE),
              sum(diff == 0, na.rm = TRUE)))
}

# ---- Plots --------------------------------------------------------------
png(file.path(OUT_DIR, "cid_six_prior.png"),
    width = 1500, height = 600, res = 110)
op <- par(mfrow = c(1, 2), mar = c(4, 4.5, 3, 1))

# Left: marginal distribution per arm (all tasks)
boxplot(mean_cid ~ arm, data = df,
        col = COLS[ARMS], border = "grey30",
        ylab = "Task-mean CID to truth",
        main = "All available tasks per arm",
        las = 1)
mtext(side = 1, line = 2.5,
      text = sprintf("n_tasks: %s",
                     paste(paste0(ARMS, "=", table(df$arm)),
                           collapse = ", ")),
      cex = 0.75, col = "grey40")

# Right: paired diff vs mk_kp1 (common subset)
diff_long <- do.call(rbind, lapply(setdiff(ARMS, "mk_kp1"), function(a) {
  data.frame(arm = a, d = wide[[a]] - wide$mk_kp1)
}))
diff_long$arm <- factor(diff_long$arm,
                        levels = setdiff(ARMS, "mk_kp1"))
boxplot(d ~ arm, data = diff_long,
        col = COLS[setdiff(ARMS, "mk_kp1")], border = "grey30",
        ylab = "CID(arm) - CID(mk_kp1)",
        main = sprintf("Paired diff vs mk_kp1 (n=%d tasks)",
                       length(tasks_all)),
        las = 1)
abline(h = 0, lty = 2, col = "grey40")

par(op); dev.off()
cat("\nSaved: cid_six_prior.png and cid_six_prior.rds\n")
