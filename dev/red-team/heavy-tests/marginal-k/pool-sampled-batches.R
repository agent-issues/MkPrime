#!/usr/bin/env Rscript
# pool-sampled-batches.R — apply the PRE-REGISTERED 2-batch SBC criterion to the
# sampled_k geometric arm (Stage 2). Reads the two seed-batch sims.rds files and
# computes the operative verdict — NOT the per-run all-3-AD>0.4 strict gate,
# which rejects a perfect sampler ~78% of the time (0.6^3).
#
# Criterion (MARGINAL-K-CACHE-002-resume.md §POST-FIX PASS CRITERION; the same
# bar that marginal_k passed). DECLARE PASS iff ALL of:
#   (a) GLOBAL: pooled (N=400) AD > 0.05 on each of tree_length, rate_log_sd, p,
#       AND no single batch has any param with AD < 0.01.
#   (b) TARGETED low-p (the actual Stage-1 fix): among p_true < 0.10 sims pooled,
#       fraction p-rank <= 2 is <= 0.15 AND their mean normalised p-rank in
#       [0.30, 0.70].
#   (c) NO FREEZE REGRESSION: tree_length & rate_log_sd extreme(0|L) fraction
#       <= ~0.03 in each batch.
#
# Usage (after both aggregators (17333101 / 17333103) finish):
#   Rscript pool-sampled-batches.R [b1_dir] [b2_dir]
# Defaults to sbc-results-sampled-b1 / -b2 under this script's directory.

MKDIR <- "dev/red-team/heavy-tests/marginal-k"
args  <- commandArgs(trailingOnly = TRUE)
b1dir <- if (length(args) >= 1) args[[1]] else file.path(MKDIR, "sbc-results-sampled-b1")
b2dir <- if (length(args) >= 2) args[[2]] else file.path(MKDIR, "sbc-results-sampled-b2")

readSims <- function(d) {
  f <- file.path(d, "sims.rds")
  if (!file.exists(f)) stop("missing sims.rds in ", d, " (has the aggregator run?)")
  s <- readRDS(f)
  good <- Filter(function(x) !isTRUE(x$skipped) && !is.null(x$ranks), s)
  if (!length(good)) stop("no good sims in ", f)
  good
}
b1 <- readSims(b1dir)
b2 <- readSims(b2dir)
cat(sprintf("batch1 good sims: %d   batch2 good sims: %d   pooled: %d\n",
            length(b1), length(b2), length(b1) + length(b2)))

# Common L_used across BOTH batches (mirrors the driver's single-L normalisation)
L_used <- max(vapply(c(b1, b2), function(s) s$L, numeric(1)))
normRank <- function(r, L = L_used) (r + 0.5) / (L + 1)

adP <- function(x) {                       # AD vs Uniform(0,1), as in the driver
  x <- x[is.finite(x)]
  if (length(x) < 4L) return(NA_real_)
  if (requireNamespace("goftest", quietly = TRUE))
    suppressWarnings(goftest::ad.test(x, null = "punif")$p.value)
  else suppressWarnings(ks.test(x, "punif")$p.value)
}
pull <- function(sims, param) vapply(sims, function(s) s$ranks[[param]], numeric(1))

params <- c("tree_length", "rate_log_sd", "p")

# ---- (a) GLOBAL: pooled AD + per-batch AD ---------------------------------
pooled <- c(b1, b2)
pooledAD  <- sapply(params, function(p) adP(normRank(pull(pooled, p))))
batch1AD  <- sapply(params, function(p) adP(normRank(pull(b1, p))))
batch2AD  <- sapply(params, function(p) adP(normRank(pull(b2, p))))

cat("\n--- (a) GLOBAL pooled AD (>0.05) + per-batch AD (no <0.01) ---\n")
for (p in params)
  cat(sprintf("  %-12s pooled=%.4f  b1=%.4f  b2=%.4f\n",
              p, pooledAD[[p]], batch1AD[[p]], batch2AD[[p]]))
clause_a <- all(pooledAD > 0.05, na.rm = TRUE) &&
            all(c(batch1AD, batch2AD) >= 0.01, na.rm = TRUE)

# ---- (b) TARGETED low-p (the actual fix) ----------------------------------
p_true  <- vapply(pooled, function(s) s$p_true, numeric(1))
p_rank  <- pull(pooled, "p")
lowp    <- p_true < 0.10
n_lowp  <- sum(lowp)
frac_le2  <- if (n_lowp) mean(p_rank[lowp] <= 2L) else NA_real_
mean_norm <- if (n_lowp) mean(normRank(p_rank[lowp])) else NA_real_
cat(sprintf("\n--- (b) TARGETED low-p (p_true<0.10, n=%d) ---\n", n_lowp))
cat(sprintf("  frac p-rank<=2 = %.3f (<=0.15)   mean normRank = %.3f (in [0.30,0.70])\n",
            frac_le2, mean_norm))
clause_b <- isTRUE(n_lowp >= 4L) && isTRUE(frac_le2 <= 0.15) &&
            isTRUE(mean_norm >= 0.30 && mean_norm <= 0.70)
if (n_lowp < 4L)
  cat("  ! too few low-p sims to judge (b); inspect manually / draw more low-p sims\n")

# ---- (c) NO FREEZE REGRESSION ---------------------------------------------
extremeFrac <- function(sims, param)
  mean(vapply(sims, function(s) s$ranks[[param]] %in% c(0L, s$L), logical(1)))
cat("\n--- (c) freeze: tl/rls extreme(0|L) fraction (<=~0.03) ---\n")
fr <- list(
  b1_tl = extremeFrac(b1, "tree_length"), b1_rls = extremeFrac(b1, "rate_log_sd"),
  b2_tl = extremeFrac(b2, "tree_length"), b2_rls = extremeFrac(b2, "rate_log_sd"))
for (nm in names(fr)) cat(sprintf("  %-7s %.3f\n", nm, fr[[nm]]))
clause_c <- all(unlist(fr) <= 0.03 + 1e-9)

verdict <- if (isTRUE(clause_a && clause_b && clause_c)) "PASS" else "FAIL"
cat(sprintf("\n==== PRE-REGISTERED 2-BATCH VERDICT: %s ====\n", verdict))
cat(sprintf("  (a) global AD       : %s\n", if (clause_a) "PASS" else "FAIL"))
cat(sprintf("  (b) targeted low-p  : %s\n", if (isTRUE(clause_b)) "PASS" else "FAIL/NA"))
cat(sprintf("  (c) no freeze       : %s\n", if (clause_c) "PASS" else "FAIL"))

out <- file.path(MKDIR, "sampled-pooled-verdict.txt")
writeLines(c(
  sprintf("PRE-REGISTERED 2-BATCH VERDICT (sampled_k): %s", verdict),
  sprintf("pooled good sims: %d (b1=%d b2=%d), L_used=%d", length(pooled), length(b1), length(b2), L_used),
  sprintf("(a) pooled AD  tl=%.4f rls=%.4f p=%.4f  [>0.05]", pooledAD[["tree_length"]], pooledAD[["rate_log_sd"]], pooledAD[["p"]]),
  sprintf("    batch1 AD  tl=%.4f rls=%.4f p=%.4f  [no <0.01]", batch1AD[["tree_length"]], batch1AD[["rate_log_sd"]], batch1AD[["p"]]),
  sprintf("    batch2 AD  tl=%.4f rls=%.4f p=%.4f", batch2AD[["tree_length"]], batch2AD[["rate_log_sd"]], batch2AD[["p"]]),
  sprintf("(b) low-p n=%d frac<=2=%.3f mean=%.3f", n_lowp, frac_le2, mean_norm),
  sprintf("(c) freeze b1(tl/rls)=%.3f/%.3f b2(tl/rls)=%.3f/%.3f", fr$b1_tl, fr$b1_rls, fr$b2_tl, fr$b2_rls)
), out)
cat(sprintf("\nwrote %s\n", out))
quit(status = if (verdict == "PASS") 0L else 2L)
