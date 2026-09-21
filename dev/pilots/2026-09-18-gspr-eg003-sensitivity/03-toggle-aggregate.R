# Aggregate the gibbsSpr on/off grid.
SP <- Sys.getenv("GSPR_DIR", ".")
fs <- list.files(file.path(SP, "grid"), "\\.rds$", full.names = TRUE)
cat("cells:", length(fs), "\n")

d <- do.call(rbind, lapply(fs, function(f) {
  x <- readRDS(f)
  u <- x$u_post
  data.frame(tree = x$tree_idx, arm = x$arm, seed = x$seed,
             nsamp = x$n_samples, elapsed = x$elapsed,
             mean_u = mean(u), median_u = median(u),
             p_lt_05 = mean(u < 0.5),
             rho = suppressWarnings(cor(u, x$u_true, method = "spearman")),
             p_mean = x$p_mean, tree_length = x$tree_length,
             log_lik = x$log_lik,
             mcse_k = median(x$kp_se, na.rm = TRUE))
}))
d$arm <- factor(d$arm, levels = c("on", "off"))
cat("\ncells per (tree, arm):\n"); print(table(d$tree, d$arm))

METRICS <- c("mean_u", "median_u", "p_lt_05", "rho", "p_mean", "tree_length", "log_lik")

cat("\n=== within-arm seed-to-seed SD (Monte Carlo error at this chain length) ===\n")
for (m in METRICS) {
  s <- tapply(d[[m]], list(d$tree, d$arm), sd)
  cat(sprintf("  %-12s  median over cells = %.4f   max = %.4f\n",
              m, median(s, na.rm = TRUE), max(s, na.rm = TRUE)))
}

cat("\n=== paired arm contrast (on - off), per dataset then pooled ===\n")
res <- list()
for (m in METRICS) {
  on  <- tapply(d[[m]][d$arm == "on"],  d$tree[d$arm == "on"],  mean)
  off <- tapply(d[[m]][d$arm == "off"], d$tree[d$arm == "off"], mean)
  sd_on  <- tapply(d[[m]][d$arm == "on"],  d$tree[d$arm == "on"],  sd)
  sd_off <- tapply(d[[m]][d$arm == "off"], d$tree[d$arm == "off"], sd)
  k  <- intersect(names(on), names(off))
  df <- on[k] - off[k]
  nD <- length(k)
  # SE of the pooled paired difference: between-dataset spread of the
  # (already seed-averaged) differences.
  se <- sd(df) / sqrt(nD)
  tt <- t.test(df)
  # pure MC-noise floor: what SE would we expect if the arms were identical?
  se_mc <- sqrt(mean((sd_on[k]^2 + sd_off[k]^2) / 3)) / sqrt(nD)
  cat(sprintf("  %-12s  on=%+.4f  off=%+.4f  diff=%+.4f  SE=%.4f  95%%CI [%+.4f,%+.4f]  p=%.3f  (MC-only SE %.4f)\n",
              m, mean(on[k]), mean(off[k]), mean(df), se,
              tt$conf.int[1], tt$conf.int[2], tt$p.value, se_mc))
  res[[m]] <- list(on = on[k], off = off[k], diff = df, ci = tt$conf.int,
                   p = tt$p.value, se = se, se_mc = se_mc, nD = nD)
}

cat("\n=== per-dataset paired differences ===\n")
tab <- data.frame(tree = as.integer(names(res$rho$diff)),
                  d_mean_u = round(as.numeric(res$mean_u$diff), 4),
                  d_rho    = round(as.numeric(res$rho$diff), 4),
                  d_plt05  = round(as.numeric(res$p_lt_05$diff), 4),
                  d_p      = round(as.numeric(res$p_mean$diff), 4))
print(tab)

# --- equivalence verdict -----------------------------------------------------
# Reference scales: the character-ordering artefact moves rho by +0.427;
# the geo-vs-EG prior contrast moves mean u_post by 1.18.
EQ <- list(rho = 0.10, mean_u = 0.10, p_lt_05 = 0.05)
cat("\n=== equivalence verdict (TOST-style: is the 95% CI inside +/- margin?) ===\n")
verdict <- "PASS"
for (m in names(EQ)) {
  ci <- res[[m]]$ci; mar <- EQ[[m]]
  ok <- ci[1] > -mar && ci[2] < mar
  if (!ok) verdict <- "FAIL"
  cat(sprintf("  %-10s margin +/-%.2f  CI [%+.4f,%+.4f]  -> %s\n",
              m, mar, ci[1], ci[2], if (ok) "EQUIVALENT" else "NOT SHOWN EQUIVALENT"))
}
cat("\nVERDICT:", verdict, "\n")
saveRDS(list(cells = d, contrasts = res, verdict = verdict),
        file.path(SP, "toggle-summary.rds"))
