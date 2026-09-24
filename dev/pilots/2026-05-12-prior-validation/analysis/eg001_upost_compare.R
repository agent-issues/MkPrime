# EG-001 u_post comparison: geo (26 trees × rep_01) vs EG (260 tasks)
#
# Hypothesis: EG prior over-shrinks u via missing per-character truncation
# normaliser, producing a u_post ≈ 0 anchor. The geo prior, having no such
# truncation issue, should produce u_post that tracks u_true more closely.

SUMMARY_DIR <- "dev/pilots/2026-05-12-prior-validation/summary"
EG_POST     <- "C:/Users/pjjg18/GitHub/mkprime/report-data/mkp-eg-260/post_means.rds"
GT_ROOT     <- "C:/Users/pjjg18/GitHub/mkprime/tree-inference"
OUT_DIR     <- "dev/pilots/2026-05-12-prior-validation/analysis"

# Ground truth, reordered into the character order the sampler actually used.
#
# CORRECTED 2026-09-19 (#54). `run_one.R` builds its matrix from
# `sort(list.files(..., "^chr[0-9]+\\.nex$"))`, which is LEXICAL --
# chr1, chr10, chr11, ..., chr2 -- so `kPrime_i` is the i-th lexically sorted
# file. `ground_truth.csv` is in NUMERIC order. This function used to return
# the numeric-order rows and every caller below paired them positionally
# against lexical-order posteriors, so each character's posterior was compared
# against a different character's truth.
#
# The published rho of -0.11 was therefore a permutation null by construction:
# the analysis could not have detected tracking had there been any. Corrected,
# both arms track at about +0.31. Marginal statistics -- mean, median, the
# whole u_post distribution -- are exactly invariant, because the error permutes
# a multiset; only anything paired per character changes.
# `charIdx` is the summary's recorded character order; summaries without one
# came from runs that sorted characters lexically.
gt_for <- function(task, charIdx = NULL) {
  tree <- as.integer(sub("^t([0-9]+).*", "\\1", task))
  rep  <- as.integer(sub(".*_r([0-9]+)$", "\\1", task))
  tr   <- sprintf("tree_%02d", tree)
  rp   <- sprintf("rep_%02d",  rep)
  gt   <- read.csv(file.path(GT_ROOT, tr, rp, "ground_truth.csv"))

  if (is.null(charIdx)) {
    files <- list.files(file.path(GT_ROOT, tr, rp), "^chr[0-9]+\\.nex$")
    charIdx <- as.integer(sub("^chr([0-9]+)\\.nex$", "\\1", sort(files)))
  }
  stopifnot(setequal(charIdx, gt$char_idx), !anyDuplicated(charIdx))
  gt[match(charIdx, gt$char_idx), , drop = FALSE]
}

# k' >= kObs holds by construction, so u_post < 0 is impossible for a correctly
# aligned pair. Counting impossible values is a falsifier that needs no
# modelling assumption: it was 115 (EG) and 248 (geo) of 1300 char-tasks under
# the old pairing, and 0 under this one. Asserting it here means the pairing
# cannot silently regress.
check_alignment <- function(d) {
  bad <- sum(d$u_post < -1e-9)
  if (bad > 0L) {
    stop(sprintf(
      "%d impossible u_post < 0 values: characters are misaligned (see #54)",
      bad))
  }
  d
}

# --- geo arm: read all mkp_geo_*.rds ---
geo_files <- list.files(SUMMARY_DIR, "^mkp_geo_.*\\.rds$", full.names = TRUE)
cat("geo files:", length(geo_files), "\n")
geo_per <- do.call(rbind, lapply(geo_files, function(f) {
  x <- readRDS(f)
  gt <- gt_for(x$tag, x$char_idx)
  k_post <- unname(x$kp_means)
  if (length(k_post) != nrow(gt)) return(NULL)
  data.frame(
    prior  = "geo",
    task   = x$tag,
    char   = seq_len(nrow(gt)),
    k_obs  = gt$kObs,
    k_true = gt$k_true,
    u_true = gt$u_true,
    k_post = k_post,
    u_post = k_post - gt$kObs
  )
}))

# --- EG arm: existing post_means.rds, subset to rep_01 only for fair compare ---
eg_post <- readRDS(EG_POST)
eg_per_full <- do.call(rbind, lapply(seq_along(eg_post), function(i) {
  r    <- eg_post[[i]]
  task <- if (!is.null(r$task)) r$task else names(eg_post)[i]
  k_post <- r$k_post
  gt   <- gt_for(task, r$char_idx)
  if (length(k_post) != nrow(gt)) return(NULL)
  data.frame(
    prior  = "EG",
    task   = task,
    char   = seq_len(nrow(gt)),
    k_obs  = gt$kObs,
    k_true = gt$k_true,
    u_true = gt$u_true,
    k_post = unname(k_post),
    u_post = unname(k_post) - gt$kObs
  )
}))
eg_per_r01 <- eg_per_full[grepl("_r01$", eg_per_full$task), ]

geo_per    <- check_alignment(geo_per)
eg_per_r01 <- check_alignment(eg_per_r01)

# Combined
both <- rbind(geo_per, eg_per_r01)

# --- summary tables ---
cat("\n=== Marginal summary of u_post by prior (rep_01 only) ===\n")
summ <- aggregate(u_post ~ prior, data = both, FUN = function(v) {
  c(n = length(v), mean = mean(v), median = median(v),
    p_eq_0 = mean(abs(v) < 0.05), p_lt_05 = mean(v < 0.5),
    q25 = quantile(v, 0.25), q75 = quantile(v, 0.75),
    sd = sd(v))
})
print(summ)

cat("\n=== Spearman corr(u_post, u_true) by prior ===\n")
sp <- aggregate(cbind(u_post, u_true) ~ prior, data = both,
                FUN = identity, simplify = FALSE)
for (p in unique(both$prior)) {
  sub <- both[both$prior == p, ]
  rho <- cor(sub$u_post, sub$u_true, method = "spearman")
  pe  <- cor(sub$u_post, sub$u_true, method = "pearson")
  rmse <- sqrt(mean((sub$u_post - sub$u_true)^2))
  bias <- mean(sub$u_post - sub$u_true)
  cat(sprintf("  %-3s  n=%4d  rho=%+.3f  pearson=%+.3f  RMSE(u)=%.3f  bias=%+.3f\n",
              p, nrow(sub), rho, pe, rmse, bias))
}

cat("\n=== Per-task Spearman (mean ± sd over tasks) ===\n")
per_task <- do.call(rbind, lapply(split(both, list(both$prior, both$task), drop = TRUE),
                                  function(d) {
  data.frame(prior = d$prior[1], task = d$task[1],
             rho   = suppressWarnings(cor(d$u_post, d$u_true, method = "spearman")),
             n     = nrow(d))
}))
agg <- aggregate(rho ~ prior, data = per_task,
                 FUN = function(v) c(mean = mean(v, na.rm = TRUE),
                                     sd   = sd(v,   na.rm = TRUE),
                                     n_tasks = sum(!is.na(v))))
print(agg)

cat("\n=== Where is mass at u_post = 0  vs  u_true > 0  (false-floor)? ===\n")
# A character with u_true > 0 should NOT have u_post ≈ 0 under a fair prior.
both$u_true_pos <- both$u_true > 0
ff <- aggregate(u_post ~ prior + u_true_pos, data = both, FUN = function(v)
  c(n = length(v), p_post_near_0 = mean(v < 0.5), median = median(v)))
print(ff)

# --- save outputs ---
saveRDS(list(geo_per = geo_per, eg_per_r01 = eg_per_r01, both = both,
             per_task = per_task),
        file.path(OUT_DIR, "eg001_upost_compare.rds"))
cat("\nSaved:", file.path(OUT_DIR, "eg001_upost_compare.rds"), "\n")

# --- plot ---
png(file.path(OUT_DIR, "eg001_upost_vs_utrue.png"), width = 1100, height = 500, res = 110)
op <- par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1))
for (p in c("EG", "geo")) {
  sub <- both[both$prior == p, ]
  plot(jitter(sub$u_true, 0.5), sub$u_post,
       pch = 16, cex = 0.5, col = adjustcolor("steelblue", 0.35),
       xlab = "u_true (= k_true - kObs)",
       ylab = "u_post (= E[k_prime|data] - kObs)",
       main = sprintf("%s prior  (n=%d chars, 26 trees × rep_01)", p, nrow(sub)),
       ylim = range(both$u_post), xlim = range(both$u_true))
  abline(0, 1, col = "red", lty = 2)
  abline(h = 0, col = "grey50", lty = 3)
  rho <- cor(sub$u_post, sub$u_true, method = "spearman")
  mtext(sprintf("Spearman = %+.3f", rho), side = 3, line = 0.3, cex = 0.9)
}
par(op); dev.off()
cat("Saved:", file.path(OUT_DIR, "eg001_upost_vs_utrue.png"), "\n")
