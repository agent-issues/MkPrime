suppressPackageStartupMessages(library(MkPrime))
res <- readRDS("/nobackup/pjjg18/mkp-sim3-v4cross-b/results/aware-result.rds")
S <- res$samples
phi <- S[, "phi"]; ll <- S[, "log_likelihood"]; th <- S[, "topo_hash"]
n <- length(phi)
cat("n=", n, "\n")
cat("phi by quarter:\n")
qrt <- cut(seq_len(n), 4, labels = FALSE)
for (q in 1:4) cat(sprintf("  Q%d: median=%.2f  IQR=[%.2f, %.2f]  ll=%.1f\n", q,
  median(phi[qrt == q]), quantile(phi[qrt == q], .25),
  quantile(phi[qrt == q], .75), mean(ll[qrt == q])))
cat("unique topo_hash:", length(unique(th)), "\n")
print(head(sort(table(th), decreasing = TRUE), 6))
cat("cor(phi, log_lik)=", round(cor(phi, ll), 3), "\n")
cat("min(phi) =", min(phi), " (post-relabel: phi<1 would mean no reflection)\n")

log <- tryCatch(MkPrime::ReadMkLog("/nobackup/pjjg18/mkp-sim3-v4cross-b/results/aware.log"),
                error = function(e) NULL)
if (!is.null(log)) {
  cat("raw-log cols:", head(colnames(log), 20), "\n")
  cat("raw-log nrows:", nrow(log), "\n")
  if ("phi" %in% colnames(log)) {
    rp <- log[, "phi"]
    cat(sprintf("raw phi: n=%d  min=%.3f  max=%.3f  median=%.3f  frac<1=%.3f\n",
                length(rp), min(rp), max(rp), median(rp), mean(rp < 1)))
  }
}

## Are there trees in the posterior where B IS monophyletic? Compare their phi.
trees <- res$trees; class(trees) <- "multiPhylo"
tn <- res$data$taxon_names
B <- tn[9:12]
mono_idx <- vapply(trees, function(tr) ape::is.monophyletic(tr, B), logical(1))
cat(sprintf("\nB-mono trees: %d / %d\n", sum(mono_idx), length(trees)))
cat(sprintf("phi | B-mono:    median=%.2f  IQR=[%.2f, %.2f]\n",
            median(phi[mono_idx]), quantile(phi[mono_idx], .25), quantile(phi[mono_idx], .75)))
cat(sprintf("phi | B-not-mono: median=%.2f  IQR=[%.2f, %.2f]\n",
            median(phi[!mono_idx]), quantile(phi[!mono_idx], .25), quantile(phi[!mono_idx], .75)))
cat(sprintf("log_lik | B-mono:    mean=%.2f\n", mean(ll[mono_idx])))
cat(sprintf("log_lik | B-not-mono: mean=%.2f\n", mean(ll[!mono_idx])))
## Mixing: do B-mono samples cluster temporally?
cat("B-mono indices (sample order):", which(mono_idx), "\n")
