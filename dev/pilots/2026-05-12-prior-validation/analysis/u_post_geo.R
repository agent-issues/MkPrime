# Per-character u_post audit for mkp_geo (plain geometric prior on k').
# Twin of u_post_eg.R; tests whether geo's posterior also concentrates
# near kObs+1 (which would explain why mk_kp1 beats geo on CID).

SUMMARY_DIR <- "dev/pilots/2026-05-12-prior-validation/summary"
OUT_DIR     <- "dev/pilots/2026-05-12-prior-validation/analysis"

files <- list.files(SUMMARY_DIR, "^mkp_geo_t[0-9]+_r[0-9]+\\.rds$",
                    full.names = TRUE)
cat("Found", length(files), "mkp_geo summary files\n")

rows <- list()
for (f in files) {
  x <- readRDS(f)
  if (is.null(x$kObs) || !length(x$kp_means)) next
  kobs <- as.integer(x$kObs)
  kp   <- as.numeric(x$kp_means)
  if (length(kp) != length(kobs)) next
  rows[[length(rows) + 1L]] <- data.frame(
    tag      = x$tag,
    tree_idx = x$tree_idx,
    rep_idx  = x$rep_idx,
    char_idx = seq_along(kobs),
    kobs     = kobs,
    kp_mean  = kp,
    u_post   = kp - kobs,
    p_mean   = if (length(x$p_mean)) x$p_mean else NA_real_
  )
}
df <- do.call(rbind, rows)
cat("Loaded", length(unique(df$tag)), "tasks /", nrow(df), "char-rows\n\n")

cat("=== u_post pooled distribution ===\n"); print(summary(df$u_post))
cat("\nQuantiles:\n"); print(quantile(df$u_post, c(.01,.05,.25,.5,.75,.95,.99)))

cat("\n=== u_post by kObs (mean kp_mean should be > kObs) ===\n")
print(aggregate(cbind(kp_mean, u_post) ~ kobs, data = df,
                FUN = function(v) c(n = length(v),
                                    mean = mean(v),
                                    median = median(v),
                                    sd = sd(v))))

cat("\n=== Violations (kp < kObs)? ===\n")
cat("n char-rows with u_post < 0:", sum(df$u_post < 0), "/", nrow(df), "\n")

cat("\n=== Where does u_post fall? ===\n")
brk <- c(-Inf, 0.5, 1.5, 2.5, 3.5, Inf)
lab <- c("<0.5", "[0.5,1.5)", "[1.5,2.5)", "[2.5,3.5)", ">=3.5")
print(round(prop.table(table(cut(df$u_post, breaks=brk, labels=lab,
                                  right=FALSE))), 3))

cat("\n=== p_mean (geometric hyperparameter) across tasks ===\n")
print(summary(unique(df[, c("tag","p_mean")])$p_mean))

saveRDS(df, file.path(OUT_DIR, "u_post_geo.rds"))

png(file.path(OUT_DIR, "u_post_geo.png"),
    width = 1200, height = 500, res = 110)
op <- par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1))
hist(df$u_post, breaks = 50, col = "steelblue", border = "white",
     xlab = "u_post = E[k'|data] - kObs",
     main = paste0("mkp_geo per-char u_post (n=", nrow(df), ")"))
abline(v = 1, lty = 2, col = "tomato", lwd = 2)
abline(v = 2, lty = 2, col = "orange", lwd = 2)
legend("topright", legend = c("kp1 (k=kObs+1)", "kp2 (k=kObs+2)"),
       col = c("tomato","orange"), lty = 2, lwd = 2, bty = "n")
boxplot(u_post ~ kobs, data = df, xlab = "kObs", ylab = "u_post",
        main = "u_post by kObs", col = "lightblue")
abline(h = 1, lty = 2, col = "tomato", lwd = 2)
abline(h = 2, lty = 2, col = "orange", lwd = 2)
par(op); dev.off()
cat("\nSaved: u_post_geo.png\n")
