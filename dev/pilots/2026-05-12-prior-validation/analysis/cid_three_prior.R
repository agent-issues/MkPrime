# CID three-prior comparison: mk floor vs geo vs EG
#
# Each summary RDS has a length-1000 cid vector (ClusteringInfoDistance from
# thinned posterior trees to the true generating tree). Summarize per task and
# compare distributions across priors.

SUMMARY_DIR <- "dev/pilots/2026-05-12-prior-validation/summary"
EG_CID      <- "C:/Users/pjjg18/GitHub/mkprime/report-data/mkp-eg-260/cid_compare_260.rds"
OUT_DIR     <- "dev/pilots/2026-05-12-prior-validation/analysis"

per_task_cid <- function(prefix, label) {
  files <- list.files(SUMMARY_DIR, paste0("^", prefix, "_.*\\.rds$"), full.names = TRUE)
  do.call(rbind, lapply(files, function(f) {
    x <- readRDS(f)
    cid <- x$cid
    if (!length(cid) || all(is.na(cid))) return(NULL)
    data.frame(prior = label, task = x$tag, tree_idx = x$tree_idx,
               rep_idx = x$rep_idx,
               cid_mean = mean(cid, na.rm = TRUE),
               cid_median = median(cid, na.rm = TRUE),
               cid_q05 = quantile(cid, 0.05, na.rm = TRUE),
               cid_q95 = quantile(cid, 0.95, na.rm = TRUE),
               n_trees = sum(!is.na(cid)))
  }))
}

mk  <- per_task_cid("mk",      "mk")
geo <- per_task_cid("mkp_geo", "geo")

# EG: reuse existing per-task means from cid_compare_260.rds
cc <- readRDS(EG_CID)
eg <- data.frame(prior = "EG",
                 task = cc$task,
                 tree_idx = cc$tree,
                 rep_idx = as.integer(sub(".*_r([0-9]+)$", "\\1", cc$task)),
                 cid_mean = cc$cid_eg,
                 cid_median = NA, cid_q05 = NA, cid_q95 = NA,
                 n_trees = cc$n_eg)
eg <- eg[!is.na(eg$cid_mean), ]

all_arms <- rbind(mk, geo, eg)

cat("=== Per-arm task counts ===\n")
print(table(all_arms$prior))

cat("\n=== Aggregate cid_mean by prior (all tasks) ===\n")
agg_all <- aggregate(cid_mean ~ prior, data = all_arms,
                     FUN = function(v) c(n_tasks = length(v),
                                         mean = mean(v),
                                         median = median(v),
                                         q25 = quantile(v, 0.25),
                                         q75 = quantile(v, 0.75),
                                         sd = sd(v)))
print(agg_all)

# Fair like-for-like: restrict to the 26 tasks present in all three arms (rep_01 only)
common_tasks <- Reduce(intersect, list(mk$task, geo$task, eg$task))
cat("\n=== Common tasks across all 3 arms: n =", length(common_tasks), "===\n")
print(common_tasks)

sub <- all_arms[all_arms$task %in% common_tasks, ]
cat("\n=== Aggregate cid_mean on common task set ===\n")
agg_common <- aggregate(cid_mean ~ prior, data = sub,
                        FUN = function(v) c(n = length(v),
                                            mean = mean(v),
                                            median = median(v),
                                            sd = sd(v)))
print(agg_common)

# Paired test: per common task, compare arms
cat("\n=== Paired Δcid by task (geo − EG, geo − mk, mk − EG) ===\n")
wide <- reshape(sub[, c("task", "prior", "cid_mean")],
                idvar = "task", timevar = "prior", direction = "wide")
wide$d_geo_EG <- wide$cid_mean.geo - wide$cid_mean.EG
wide$d_geo_mk <- wide$cid_mean.geo - wide$cid_mean.mk
wide$d_mk_EG  <- wide$cid_mean.mk  - wide$cid_mean.EG
for (cn in c("d_geo_EG", "d_geo_mk", "d_mk_EG")) {
  v <- wide[[cn]]
  cat(sprintf("  %-9s  mean=%+.4f  median=%+.4f  signed-rank p=%.3g\n",
              cn, mean(v, na.rm = TRUE), median(v, na.rm = TRUE),
              wilcox.test(v)$p.value))
}

# --- save ---
saveRDS(list(all_arms = all_arms, common = sub, paired_wide = wide),
        file.path(OUT_DIR, "cid_three_prior.rds"))
cat("\nSaved:", file.path(OUT_DIR, "cid_three_prior.rds"), "\n")

# --- plot: boxplots on common task set ---
png(file.path(OUT_DIR, "cid_three_prior.png"),
    width = 900, height = 500, res = 110)
op <- par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1))
boxplot(cid_mean ~ prior, data = sub,
        ylab = "Posterior mean CID to truth",
        main = paste0("Common tasks (n=", length(common_tasks), ")"),
        col = c("steelblue", "tomato", "darkgreen"))
mtext("EG vs geo vs mk floor", side = 3, line = 0.2, cex = 0.85)
# all tasks per arm
boxplot(cid_mean ~ prior, data = all_arms,
        ylab = "Posterior mean CID to truth",
        main = "All tasks per arm",
        col = c("steelblue", "tomato", "darkgreen"))
mtext(paste0("mk n=", sum(all_arms$prior=="mk"),
             ", geo n=", sum(all_arms$prior=="geo"),
             ", EG n=",  sum(all_arms$prior=="EG")),
      side = 3, line = 0.2, cex = 0.85)
par(op); dev.off()
cat("Saved:", file.path(OUT_DIR, "cid_three_prior.png"), "\n")
