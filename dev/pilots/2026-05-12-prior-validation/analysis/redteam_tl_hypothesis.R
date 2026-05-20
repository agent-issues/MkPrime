#!/usr/bin/env Rscript
# Red-team test: does posterior tree length (TL) inflate with k as JC(k) predicts,
# and does that inflation track the mk_k40 > Mk' CID advantage?
#
# Mechanism: under JC(k), at small T, off-diag transition prob per unit T ~ T/(k-1)*(k/k).
# More precisely, the substitution prob per branch over time T is
#     p_subst(T,k) = (k-1)/k * (1 - exp(-k*T/(k-1)))
# For matched observed substitution rate p between two k values,
#     T_k = -((k-1)/k) * log(1 - k*p/(k-1))
# so T_k / T_2 grows roughly linearly with k at small p (factor ~ k/2 at p->0).

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ape)
})

summary_dir <- "dev/pilots/2026-05-12-prior-validation/summary"
out_dir     <- "dev/pilots/2026-05-12-prior-validation/analysis"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Load all summaries -------------------------------------------------
files <- list.files(summary_dir, pattern = "\\.rds$", full.names = TRUE)
cat("Total summary files:", length(files), "\n")

extract_one <- function(f) {
  x <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(x)) return(NULL)
  sm <- x$scalar_means
  tl <- if ("tree_length" %in% names(sm)) unname(sm["tree_length"]) else NA_real_
  ll <- if ("log_likelihood" %in% names(sm)) unname(sm["log_likelihood"]) else NA_real_
  lp <- if ("log_posterior" %in% names(sm)) unname(sm["log_posterior"]) else NA_real_
  # CID mean
  cid_mean <- if (!is.null(x$cid)) mean(x$cid, na.rm = TRUE) else NA_real_
  # tree size (n tips) from first thinned tree (cheap)
  n_tip <- if (!is.null(x$thinned_trees) && length(x$thinned_trees) > 0) {
    length(x$thinned_trees[[1]]$tip.label)
  } else NA_integer_
  data.table(
    arm       = x$arm,
    tree_idx  = x$tree_idx,
    rep_idx   = x$rep_idx,
    tag       = x$tag,
    n_tip     = n_tip,
    tree_length = tl,
    mean_cid    = cid_mean,
    log_lik     = ll,
    log_post    = lp,
    n_samples   = x$n_samples
  )
}

t0 <- Sys.time()
rows <- lapply(files, extract_one)
dt   <- rbindlist(rows, fill = TRUE)
cat(sprintf("Loaded %d rows in %.1fs\n", nrow(dt),
            as.numeric(Sys.time() - t0, units = "secs")))

# Drop the freak few that have no TL
dt <- dt[!is.na(tree_length)]
cat("After dropping NA TL:", nrow(dt), "rows\n")

# Map arm -> k (fixed-k arms only; non-fixed arms get NA)
arm_k <- c(mk = 2L, mk_kp1 = NA, mk_kp2 = NA, mk_k9 = 9L, mk_k15 = 15L,
           mk_k24 = 24L, mk_k40 = 40L, mk_ktrue = NA, mk_tlshrink = NA,
           mkp_eg = NA, mkp_geo = NA, mkp_highk = NA, mkp_logs = NA)
dt[, k_fixed := arm_k[arm]]

# Order arms for plotting (full 12-arm cohort + tlshrink)
arm_order <- c("mk", "mk_kp1", "mk_kp2", "mk_k9", "mk_k15", "mk_k24", "mk_k40",
               "mk_ktrue", "mk_tlshrink",
               "mkp_eg", "mkp_geo", "mkp_highk", "mkp_logs")
dt[, arm := factor(arm, levels = arm_order)]

# ---- 2. Per-arm summary ----------------------------------------------------
arm_summ <- dt[, .(
  n_runs   = .N,
  mean_TL  = mean(tree_length),
  med_TL   = median(tree_length),
  sd_TL    = sd(tree_length),
  mean_CID = mean(mean_cid, na.rm = TRUE)
), keyby = arm]
print(arm_summ)

fwrite(arm_summ, file.path(out_dir, "redteam_tl_arm_summary.csv"))

# ---- 3. Plot 1: boxplot per arm --------------------------------------------
p1 <- ggplot(dt, aes(x = arm, y = tree_length, fill = arm)) +
  geom_boxplot(outlier.size = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.15, size = 0.4) +
  scale_y_log10() +
  labs(title = "Posterior mean tree length by arm",
       subtitle = "log10 axis; one point per (task,rep)",
       x = NULL, y = "Posterior mean tree length") +
  theme_bw() + theme(legend.position = "none",
                     axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(out_dir, "redteam_tl_1_boxplot.png"), p1,
       width = 7, height = 4.5, dpi = 150)

# ---- 4. Plot 2: TL ramp vs k for the fixed-k arms --------------------------
dt_fixed <- dt[!is.na(k_fixed)]
ramp <- dt_fixed[, .(mean_TL = mean(tree_length),
                     sd_TL   = sd(tree_length),
                     n       = .N), keyby = .(arm, k_fixed)]
print(ramp)

# theoretical: assume some observed substitution prob p_obs and compute T(k)
# We pick p_obs such that T(k=2) matches the empirical mean_TL[k=2].
# p_subst(T,k) = (k-1)/k * (1 - exp(-k T/(k-1)))
# Solve for p_obs at k=2 using T2 = mean_TL[k=2] (TL is summed over branches;
# we use the mean *branch* substitution prob proxy — TL/n_branches roughly
# scales with per-branch T, so we use TL as proxy and look at ratios).
T2 <- ramp[k_fixed == 2L, mean_TL]
# substitution prob implied at k=2 with branch length = T2 / n_internal_branches
# But TL totals over all branches; ratio across k is what matters. For a
# constant data pattern, the implied *per-branch* T scales with k via:
#   T(k) such that p_subst(T(k),k) = p_subst(T2,2)
# At small p, T(k)/T2 -> k/2 . Plot theoretical curve.
ks <- c(2, 9, 15, 24, 40)
# We compute matched-T for a representative per-branch T2/n_branches.
# But TL itself is sum of branches; ratios of TL across k = ratios of per-branch T
# (same tree structure on average). So plot theoretical TL = T2 * ratio(k).
# ratio_theo(k) via solving p_subst equality:
p_subst <- function(T, k) (k - 1) / k * (1 - exp(-k * T / (k - 1)))
solve_T <- function(p, k) {
  # invert: -((k-1)/k) * log(1 - k*p/(k-1))
  arg <- 1 - k * p / (k - 1)
  if (arg <= 0) return(NA_real_)
  -((k - 1) / k) * log(arg)
}
# Use mean per-branch T at k=2 to back out p_obs.
# For TL ratios, the per-branch p_subst is constant across k by assumption.
# A reasonable proxy: pick p_obs such that mean per-branch T at k=2 = T2 / (2 n_tip - 3)
# We approximate using median n_tip in dataset.
n_tip_med <- median(dt$n_tip, na.rm = TRUE)
n_br_med  <- 2 * n_tip_med - 3
T2_per_branch <- T2 / n_br_med
p_obs <- p_subst(T2_per_branch, 2)
cat(sprintf("Median n_tip=%g, n_branches=%g, TL(k=2)=%g => per-branch T=%g, implied p_obs=%g\n",
            n_tip_med, n_br_med, T2, T2_per_branch, p_obs))

theo_TL <- vapply(ks, function(k) solve_T(p_obs, k) * n_br_med, numeric(1))
theo <- data.table(k_fixed = ks, mean_TL_theo = theo_TL)
print(theo)

p2 <- ggplot(ramp, aes(x = k_fixed, y = mean_TL)) +
  geom_errorbar(aes(ymin = mean_TL - sd_TL/sqrt(n),
                    ymax = mean_TL + sd_TL/sqrt(n)),
                width = 0.6, colour = "grey50") +
  geom_point(size = 3, colour = "steelblue") +
  geom_line(colour = "steelblue") +
  geom_line(data = theo, aes(y = mean_TL_theo),
            colour = "firebrick", linetype = "dashed") +
  geom_point(data = theo, aes(y = mean_TL_theo), colour = "firebrick", shape = 4) +
  labs(title = "TL vs k: empirical (blue) vs JC(k) theory anchored at k=2 (red)",
       subtitle = sprintf("Theory: matched per-branch substitution prob p=%.3f (n_br=%.0f)",
                          p_obs, n_br_med),
       x = "k (fixed-k arm)", y = "Posterior mean tree length") +
  theme_bw()
ggsave(file.path(out_dir, "redteam_tl_2_ramp_vs_k.png"), p2,
       width = 6.5, height = 4.5, dpi = 150)

# ---- 5. Plot 3: per-task scatter TL(mk_k40) vs TL(mkp_eg) ------------------
wide <- dcast(dt, tag + tree_idx + rep_idx ~ arm,
              value.var = "tree_length")
# Some arms (mkp_geo) only have 26 runs; the join keeps NAs there.
p3 <- ggplot(wide[!is.na(mk_k40) & !is.na(mkp_eg)],
             aes(x = mkp_eg, y = mk_k40)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey50", linetype = "dashed") +
  geom_point(alpha = 0.45, size = 0.9) +
  scale_x_log10() + scale_y_log10() +
  labs(title = "Per-task posterior TL: mk_k40 vs mkp_eg",
       subtitle = "Above the y=x line means mk_k40 used longer branches than Mk' (eg)",
       x = "TL (mkp_eg)", y = "TL (mk_k40)") +
  theme_bw()
ggsave(file.path(out_dir, "redteam_tl_3_scatter_k40_vs_eg.png"), p3,
       width = 5.5, height = 5, dpi = 150)

# Test if TL_k40 > TL_eg systematically.
paired <- wide[!is.na(mk_k40) & !is.na(mkp_eg)]
cat(sprintf("\nPaired n=%d (mk_k40 vs mkp_eg):\n", nrow(paired)))
cat(sprintf("  median TL_k40 = %.3f, median TL_eg = %.3f, ratio = %.2fx\n",
            median(paired$mk_k40), median(paired$mkp_eg),
            median(paired$mk_k40) / median(paired$mkp_eg)))
wt <- wilcox.test(paired$mk_k40, paired$mkp_eg, paired = TRUE,
                  alternative = "greater")
cat(sprintf("  Wilcoxon paired (k40 > eg): V=%g, p=%.3g\n",
            wt$statistic, wt$p.value))

# ---- 6. Does TL inflation correlate with CID improvement? ------------------
# Per task: delta_CID = CID(mkp_eg) - CID(mk_k40) (positive = mk_k40 wins)
#          ratio_TL = TL(mk_k40) / TL(mkp_eg)
cid_wide <- dcast(dt, tag + tree_idx + rep_idx ~ arm, value.var = "mean_cid")
both <- merge(wide[, .(tag, tree_idx, rep_idx, TL_k40 = mk_k40, TL_eg = mkp_eg)],
              cid_wide[, .(tag, CID_k40 = mk_k40, CID_eg = mkp_eg)],
              by = "tag")
both <- both[!is.na(TL_k40) & !is.na(TL_eg) & !is.na(CID_k40) & !is.na(CID_eg)]
both[, delta_CID := CID_eg - CID_k40]   # positive = mk_k40 better
both[, log_TL_ratio := log(TL_k40 / TL_eg)]

ct <- cor.test(both$log_TL_ratio, both$delta_CID, method = "spearman")
cat(sprintf("\nSpearman cor (log TL_k40/TL_eg vs delta_CID): rho=%.3f, p=%.3g\n",
            ct$estimate, ct$p.value))

p4 <- ggplot(both, aes(x = log_TL_ratio, y = delta_CID)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_point(alpha = 0.5, size = 0.9) +
  geom_smooth(method = "lm", colour = "steelblue", se = TRUE) +
  labs(title = "Per-task: does TL inflation predict mk_k40's CID advantage?",
       subtitle = sprintf("Spearman rho=%.3f, p=%.2g (n=%d)",
                          ct$estimate, ct$p.value, nrow(both)),
       x = "log( TL_k40 / TL_eg )",
       y = "CID_eg - CID_k40 (positive = mk_k40 closer to truth)") +
  theme_bw()
ggsave(file.path(out_dir, "redteam_tl_4_dtl_vs_dcid.png"), p4,
       width = 6, height = 4.5, dpi = 150)

# ---- 7. Stratify by tree size ---------------------------------------------
# Split TL inflation by quartile of n_tip
if (!all(is.na(dt$n_tip))) {
  paired_nt <- merge(paired,
                     dt[arm == "mk_k40", .(tag, n_tip)],
                     by = "tag", all.x = TRUE)
  qb <- unique(quantile(paired_nt$n_tip, c(0, .25, .5, .75, 1), na.rm = TRUE))
  if (length(qb) >= 2) {
    paired_nt[, n_tip_q := cut(n_tip, breaks = qb, include.lowest = TRUE)]
  } else {
    paired_nt[, n_tip_q := factor("all")]
  }
  by_size <- paired_nt[, .(
    n = .N,
    med_TL_k40 = median(mk_k40),
    med_TL_eg  = median(mkp_eg),
    ratio      = median(mk_k40 / mkp_eg)
  ), keyby = n_tip_q]
  print(by_size)
  fwrite(by_size, file.path(out_dir, "redteam_tl_by_treesize.csv"))
}

# ---- 8. Save data ----------------------------------------------------------
saveRDS(list(dt = dt, ramp = ramp, theo = theo, paired = paired, both = both),
        file.path(out_dir, "redteam_tl_hypothesis.rds"))

cat("\nDone. Outputs in", out_dir, "\n")
