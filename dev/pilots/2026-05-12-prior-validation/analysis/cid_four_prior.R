# CID four-arm comparison: mk floor (M(kObs)) vs geo vs EG vs mk_kp1 (M(kObs+1))
#
# Extends cid_three_prior.R with the mk_kp1 reference arm.

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

# `mk` prefix would collide with `mk_kp1`; load mk explicitly via
# anchored regex on the filename.
mk_files <- list.files(SUMMARY_DIR, "^mk_t[0-9]+_r[0-9]+\\.rds$",
                       full.names = TRUE)
mk <- do.call(rbind, lapply(mk_files, function(f) {
  x <- readRDS(f); cid <- x$cid
  if (!length(cid) || all(is.na(cid))) return(NULL)
  data.frame(prior = "mk", task = x$tag, tree_idx = x$tree_idx,
             rep_idx = x$rep_idx,
             cid_mean = mean(cid, na.rm = TRUE),
             cid_median = median(cid, na.rm = TRUE),
             cid_q05 = quantile(cid, 0.05, na.rm = TRUE),
             cid_q95 = quantile(cid, 0.95, na.rm = TRUE),
             n_trees = sum(!is.na(cid)))
}))

geo     <- per_task_cid("mkp_geo", "geo")
mk_kp1  <- per_task_cid("mk_kp1",  "mk_kp1")

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

all_arms <- rbind(mk, geo, eg, mk_kp1)

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

# Fair like-for-like: restrict to tasks present in all four arms
common_tasks <- Reduce(intersect, list(mk$task, geo$task, eg$task, mk_kp1$task))
cat("\n=== Common tasks across all 4 arms: n =", length(common_tasks), "===\n")
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
cat("\n=== Paired Δcid by task ===\n")
wide <- reshape(sub[, c("task", "prior", "cid_mean")],
                idvar = "task", timevar = "prior", direction = "wide")
contrasts <- list(
  d_mk_kp1_EG  = c("mk_kp1", "EG"),
  d_mk_kp1_geo = c("mk_kp1", "geo"),
  d_mk_kp1_mk  = c("mk_kp1", "mk"),
  d_geo_EG     = c("geo",    "EG"),
  d_geo_mk     = c("geo",    "mk"),
  d_mk_EG      = c("mk",     "EG")
)
for (nm in names(contrasts)) {
  a <- contrasts[[nm]][1]; b <- contrasts[[nm]][2]
  v <- wide[[paste0("cid_mean.", a)]] - wide[[paste0("cid_mean.", b)]]
  cat(sprintf("  %-13s  mean=%+.4f  median=%+.4f  signed-rank p=%.3g\n",
              nm, mean(v, na.rm = TRUE), median(v, na.rm = TRUE),
              wilcox.test(v)$p.value))
}

# --- save ---
saveRDS(list(all_arms = all_arms, common = sub, paired_wide = wide),
        file.path(OUT_DIR, "cid_four_prior.rds"))
cat("\nSaved:", file.path(OUT_DIR, "cid_four_prior.rds"), "\n")

# --- plot: boxplots on common task set ---
arm_levels <- c("mk", "mk_kp1", "geo", "EG")
arm_cols   <- c(mk = "tomato", mk_kp1 = "orange",
                geo = "steelblue", EG = "darkgreen")
sub$prior      <- factor(sub$prior,      levels = arm_levels)
all_arms$prior <- factor(all_arms$prior, levels = arm_levels)

png(file.path(OUT_DIR, "cid_four_prior.png"),
    width = 1000, height = 500, res = 110)
op <- par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1))
boxplot(cid_mean ~ prior, data = sub,
        ylab = "Posterior mean CID to truth",
        main = paste0("Common tasks (n=", length(common_tasks), ")"),
        col = arm_cols[arm_levels])
mtext("mk floor vs mk_kp1 vs geo vs EG", side = 3, line = 0.2, cex = 0.85)
boxplot(cid_mean ~ prior, data = all_arms,
        ylab = "Posterior mean CID to truth",
        main = "All tasks per arm",
        col = arm_cols[arm_levels])
mtext(paste0("mk n=", sum(all_arms$prior=="mk"),
             ", mk_kp1 n=", sum(all_arms$prior=="mk_kp1"),
             ", geo n=", sum(all_arms$prior=="geo"),
             ", EG n=",  sum(all_arms$prior=="EG")),
      side = 3, line = 0.2, cex = 0.85)
par(op); dev.off()
cat("Saved:", file.path(OUT_DIR, "cid_four_prior.png"), "\n")
