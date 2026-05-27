## Red-team: stratify CID by tree shape (J1)
## Goal: does mk_k40's advantage over mk depend on tree balance?

suppressPackageStartupMessages({
  library(ggplot2)
})

out_dir <- "C:/Users/pjjg18/GitHub/mkp/dev/pilots/2026-05-12-prior-validation/analysis"

cid  <- readRDS(file.path(out_dir, "cid_nine_prior.rds"))
meta <- read.csv("C:/Users/pjjg18/GitHub/mkprime/tree-inference/trees_meta.csv",
                 stringsAsFactors = FALSE)
meta$tree_idx <- as.integer(sub("tree_", "", meta$tree_id))

cid <- merge(cid, meta[, c("tree_idx", "j1_index", "tree_length")],
             by = "tree_idx", all.x = TRUE)

cat("J1 distribution across 26 trees:\n")
print(summary(meta$j1_index))
cat("Quartiles:\n")
print(quantile(meta$j1_index, c(0, .25, .5, .75, 1)))

## ---- focus arms: drop mkp_geo (1 rep), keep ramp + mkp_eg ----
arms_keep <- c("mk", "mk_kp1", "mk_kp2", "mk_k9", "mk_k15",
               "mk_k24", "mk_k40", "mkp_eg")
cid <- droplevels(cid[cid$arm %in% arms_keep, ])
cid$arm <- factor(cid$arm, levels = arms_keep)

## ============================================================
## 1. CID vs J1, one curve per arm (mean across reps per tree)
## ============================================================
per_tree <- aggregate(mean_cid ~ arm + tree_idx + j1_index,
                      data = cid, FUN = mean)
per_tree_med <- aggregate(median_cid ~ arm + tree_idx + j1_index,
                          data = cid, FUN = median)

p1 <- ggplot(per_tree, aes(j1_index, mean_cid, colour = arm)) +
  geom_point(alpha = 0.7) +
  geom_smooth(se = FALSE, method = "loess", span = 1.2, linewidth = 0.7) +
  labs(title = "Mean CID-to-truth vs J1 (per-tree mean across 10 reps)",
       x = "J1 (balance index; 0 = pectinate, 1 = balanced)",
       y = "Mean CID across reps") +
  theme_minimal()
ggsave(file.path(out_dir, "redteam_stratify_1_cid_vs_j1.png"),
       p1, width = 9, height = 6, dpi = 130)

## ============================================================
## 2. Paired Δ(mk_k40 − mk) per tree×rep, vs J1
## ============================================================
mk    <- cid[cid$arm == "mk",   c("tree_idx", "rep_idx", "mean_cid", "median_cid", "j1_index")]
mk40  <- cid[cid$arm == "mk_k40", c("tree_idx", "rep_idx", "mean_cid", "median_cid")]
names(mk)[3:4]   <- c("cid_mk",   "med_mk")
names(mk40)[3:4] <- c("cid_mk40", "med_mk40")
pair <- merge(mk, mk40, by = c("tree_idx", "rep_idx"))
pair$delta <- pair$cid_mk40 - pair$cid_mk           # negative = mk_k40 better

cat("\nOverall paired Δ(mk_k40 - mk), mean CID, across all task-reps:\n")
print(summary(pair$delta))
cat("n pairs:", nrow(pair), "  pct where mk_k40 better:",
    round(mean(pair$delta < 0) * 100, 1), "%\n")

## per-tree mean diff
diff_tree <- aggregate(delta ~ tree_idx + j1_index, data = pair, FUN = mean)
diff_tree$n <- aggregate(delta ~ tree_idx, data = pair, FUN = length)$delta

p2 <- ggplot(diff_tree, aes(j1_index, delta)) +
  geom_hline(yintercept = 0, lty = 2, colour = "grey50") +
  geom_point(size = 2.5) +
  geom_smooth(method = "loess", se = TRUE, span = 1.2) +
  labs(title = "Per-tree mean Δ CID = mk_k40 − mk (negative = mk_k40 wins)",
       x = "J1 balance index",
       y = "mean(CID_mk_k40 − CID_mk) over 10 reps") +
  theme_minimal()
ggsave(file.path(out_dir, "redteam_stratify_2_delta_vs_j1.png"),
       p2, width = 9, height = 6, dpi = 130)

## also a Bland-style: paired CIDs
p2b <- ggplot(pair, aes(cid_mk, cid_mk40, colour = j1_index)) +
  geom_abline(slope = 1, intercept = 0, lty = 2, colour = "grey50") +
  geom_point(alpha = 0.7) +
  scale_colour_viridis_c() +
  labs(title = "Paired CID per task: mk vs mk_k40 (below y=x means k40 better)",
       x = "CID mk", y = "CID mk_k40", colour = "J1") +
  theme_minimal()
ggsave(file.path(out_dir, "redteam_stratify_4_scatter.png"),
       p2b, width = 8, height = 7, dpi = 130)

## ============================================================
## 3. Median-based ramp (robust)
## ============================================================
## per task: take median_cid already in data.  Then median across reps per
## (arm, tree), then median across trees per arm -- the 'all-median' ramp.
med_per_tree <- aggregate(median_cid ~ arm + tree_idx,
                          data = cid, FUN = median)
ramp_median <- aggregate(median_cid ~ arm, data = med_per_tree,
                         FUN = median)
ramp_mean   <- aggregate(mean_cid ~ arm,
                         data = aggregate(mean_cid ~ arm + tree_idx,
                                          data = cid, FUN = mean),
                         FUN = mean)
ramp <- merge(ramp_mean, ramp_median, by = "arm")
ramp <- ramp[match(arms_keep, ramp$arm), ]
cat("\nRamp comparison (mean-of-means vs median-of-medians):\n")
print(ramp, row.names = FALSE)

## ============================================================
## 4. J1-binned summary of (mk_k40 vs mk) gap
## ============================================================
diff_tree$j1_bin <- cut(diff_tree$j1_index,
                       breaks = quantile(meta$j1_index,
                                         c(0, 1/3, 2/3, 1)),
                       include.lowest = TRUE,
                       labels = c("pectinate (low J1)",
                                  "intermediate",
                                  "balanced (high J1)"))
cat("\nMean Δ(mk_k40 - mk) by J1 tercile (per-tree means):\n")
binagg <- aggregate(delta ~ j1_bin, data = diff_tree,
                    FUN = function(d) c(mean = mean(d),
                                         median = median(d),
                                         n = length(d),
                                         pct_neg = mean(d < 0)))
print(binagg)

## a paired t test per bin
cat("\nPaired t-test of Δ within each J1 tercile (per-tree means):\n")
for (b in levels(diff_tree$j1_bin)) {
  d <- diff_tree$delta[diff_tree$j1_bin == b]
  tt <- t.test(d)
  cat(sprintf("  %s: n=%d  mean=%.4f  t=%.2f  p=%.3g\n",
              b, length(d), mean(d), tt$statistic, tt$p.value))
}

## ============================================================
## 5. Per-arm CID-vs-J1 regression slope (does mk_k40 win more on
##    balanced trees because pectinate trees are uniformly harder?)
## ============================================================
cat("\nLinear slope of per-tree mean CID on J1 by arm:\n")
slopes <- do.call(rbind, lapply(split(per_tree, per_tree$arm), function(df) {
  fit <- lm(mean_cid ~ j1_index, data = df)
  cf <- summary(fit)$coefficients
  data.frame(arm = unique(df$arm),
             intercept = cf[1, 1],
             slope = cf[2, 1],
             slope_se = cf[2, 2],
             p = cf[2, 4])
}))
print(slopes, row.names = FALSE)

## tile plot of per-tree paired delta
p5 <- ggplot(diff_tree, aes(reorder(factor(tree_idx), j1_index),
                            delta, fill = j1_index)) +
  geom_col() +
  scale_fill_viridis_c() +
  geom_hline(yintercept = 0) +
  labs(title = "Per-tree paired gap (k40 - mk); trees ordered by J1",
       x = "tree_idx (sorted by J1)",
       y = "mean Δ CID across 10 reps") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, size = 7))
ggsave(file.path(out_dir, "redteam_stratify_3_pertree_bars.png"),
       p5, width = 10, height = 5, dpi = 130)

## save the per-tree diff for inspection
write.csv(diff_tree[order(diff_tree$j1_index), ],
          file.path(out_dir, "redteam_stratify_diff_per_tree.csv"),
          row.names = FALSE)

cat("\nDone. Plots & csv saved to: ", out_dir, "\n")
