## Diagnose aware-chain B-clade asymmetry on sim3-v4cross-b
## Usage: Rscript diagnose.R <result.rds>
suppressPackageStartupMessages({
  library(MkPrime)
  library(ape)
})

args <- commandArgs(trailingOnly = TRUE)
path <- if (length(args) >= 1) args[1] else
  "/nobackup/pjjg18/mkp-sim3-v4cross-b/results/aware-result.rds"

res <- readRDS(path)
cat("\n=== file:", path, "===\n")
cat("relabelled attr:", isTRUE(attr(res, "relabelled")), "\n")
cat("ecologyAware:", isTRUE(res$model$ecologyAware), "\n")
cat("magnitudeMode:", res$model$magnitudeMode, "\n")
cat("nSamples:", res$nSamples, " stop:", res$stop_reason, "\n")

## --- taxa / ecology layout ---
tn <- res$data$taxon_names
eco <- res$data$ecology
cat("\nTaxa:", tn, "\n")
if (length(eco) == length(tn)) print(data.frame(tip = tn, eco = eco))

## --- samples matrix ---
S <- res$samples
cat("\nSample cols (first 20):", head(colnames(S), 20), "\n")
phi_col <- grep("^phi", colnames(S), value = TRUE)
pi0_col <- grep("^pi0|^rho0", colnames(S), value = TRUE)
cat("phi cols:", phi_col, "  pi0/rho0 cols:", pi0_col, "\n")

if (length(phi_col)) {
  phi <- S[, phi_col[1]]
  cat(sprintf("\nphi: median=%.3f  mean=%.3f  q=[%.3f, %.3f, %.3f, %.3f, %.3f]  range=[%.3f, %.3f]\n",
              median(phi), mean(phi),
              quantile(phi, 0.025), quantile(phi, 0.25),
              quantile(phi, 0.5),   quantile(phi, 0.75), quantile(phi, 0.975),
              min(phi), max(phi)))
  cat("phi <1 fraction:", mean(phi < 1), "\n")
  cat("phi histogram bins (0,.25,.5,.75,1,1.5,2,3,5,Inf):\n")
  print(table(cut(phi, c(0, .25, .5, .75, 1, 1.5, 2, 3, 5, Inf))))
}
if (length(pi0_col)) {
  p0 <- S[, pi0_col[1]]
  cat(sprintf("\n%s: median=%.3f  q=[%.3f, %.3f, %.3f]\n",
              pi0_col[1], median(p0), quantile(p0,.025), quantile(p0,.5), quantile(p0,.975)))
}

## --- trees: clade monophyly & alternative B partners ---
trees <- res$trees
class(trees) <- "multiPhylo"
## NOTE: taxa are stored in order A,C,B,D (matches actual data layout).
## User-labels:  A=tip 1:4 (A4..A1), B=tip 9:12 (B4..B1), C=tip 5:8 (C4..C1), D=tip 13:16 (D4..D1)
A <- tn[1:4]; C <- tn[5:8]; B <- tn[9:12]; D <- tn[13:16]
cat("\nA:", A, " B:", B, " C:", C, " D:", D, "\n")

is_mono <- function(tr, set) ape::is.monophyletic(tr, set)
mono <- function(set) mean(vapply(trees, is_mono, logical(1), set = set))
cat(sprintf("\nMonophyly: A=%.3f  B=%.3f  C=%.3f  D=%.3f\n",
            mono(A), mono(B), mono(C), mono(D)))

sister <- function(tr, s1, s2) {
  all_t <- tr$tip.label
  out  <- setdiff(all_t, c(s1, s2))
  ape::is.monophyletic(tr, c(s1, s2)) &&
    ape::is.monophyletic(tr, s1) && ape::is.monophyletic(tr, s2)
}
cat(sprintf("(A,C) sister: %.3f   (B,D) sister: %.3f\n",
            mean(vapply(trees, sister, logical(1), s1 = A, s2 = C)),
            mean(vapply(trees, sister, logical(1), s1 = B, s2 = D))))

## Where does each B tip go? For trees where B is NOT mono, find the smallest
## clade containing each B tip and report composition.
notmono_idx <- which(!vapply(trees, is_mono, logical(1), set = B))
cat("\nTrees where B non-mono:", length(notmono_idx), "of", length(trees), "\n")

## Top alternative splits involving B tips
pp <- ape::prop.part(trees)
freq <- attr(pp, "number") / length(trees)
labs <- attr(pp, "labels")
splits <- lapply(pp, function(i) sort(labs[i]))
contains_B <- vapply(splits, function(s) any(s %in% B), logical(1))
not_pure_B <- vapply(splits, function(s) !setequal(s, B), logical(1))
small <- lengths(splits) <= 8 & lengths(splits) >= 2
keep <- contains_B & not_pure_B & small
ord <- order(-freq[keep])
cat("\nTop splits containing B tips (excluding pure {B1..B4}):\n")
for (k in head(ord, 15)) {
  idx <- which(keep)[k]
  cat(sprintf("  freq=%.3f  size=%d  {%s}\n",
              freq[idx], length(splits[[idx]]),
              paste(splits[[idx]], collapse = ",")))
}

## Per-B-tip sister composition: smallest non-trivial clade containing each B tip
cat("\nPer-B-tip nearest neighbours (mode across posterior):\n")
for (bt in B) {
  partners <- character(length(trees))
  for (i in seq_along(trees)) {
    tr <- trees[[i]]
    desc <- prop.part(list(tr))
    s <- lapply(desc, function(ix) sort(attr(desc, "labels")[ix]))
    s <- s[vapply(s, function(x) bt %in% x && length(x) >= 2 && length(x) <= 6, logical(1))]
    if (!length(s)) { partners[i] <- "NA"; next }
    smallest <- s[[which.min(lengths(s))]]
    partners[i] <- paste(setdiff(smallest, bt), collapse = ",")
  }
  tab <- sort(table(partners), decreasing = TRUE)
  cat(sprintf("  %s : ", bt))
  for (k in head(names(tab), 4)) cat(sprintf("[%s]=%d  ", k, tab[[k]]))
  cat("\n")
}

## --- z: per-character mean across posterior ---
zs <- res$z_samples
if (length(zs)) {
  Z <- do.call(cbind, lapply(zs, function(z) as.integer(z[, 1] != 0)))  # was eco-affected (any non-zero)
  zmean <- rowMeans(Z)
  ## Bin characters by which-clade-distinguishing they are.
  ## Heuristic: for each char, compute symmetric diff between A-pattern and B-pattern.
  M <- res$data$matrix  # 16 x nChar, rownames = taxa
  rn <- rownames(M)
  Aidx <- match(A, rn); Bidx <- match(B, rn)
  Cidx <- match(C, rn); Didx <- match(D, rn)
  A_var <- apply(M[Aidx, , drop = FALSE], 2, function(x) length(unique(x[!is.na(x)])) > 1)
  B_var <- apply(M[Bidx, , drop = FALSE], 2, function(x) length(unique(x[!is.na(x)])) > 1)
  ## Characters that VARY within A but not B (so "B-distinguishing-supportive" — keep B together)
  ## vs characters that vary within B but not A.
  only_A_var <- A_var & !B_var
  only_B_var <- B_var & !A_var
  both_var   <- A_var & B_var
  neither    <- !A_var & !B_var
  cat(sprintf("\nChar partitions (n=%d):\n  vary-within-A-only: %d  vary-within-B-only: %d  both: %d  neither: %d\n",
              length(A_var), sum(only_A_var), sum(only_B_var), sum(both_var), sum(neither)))
  cat(sprintf("Mean z(=eco-affected) by partition:\n"))
  cat(sprintf("  vary-A-only:  %.3f\n", mean(zmean[only_A_var])))
  cat(sprintf("  vary-B-only:  %.3f\n", mean(zmean[only_B_var])))
  cat(sprintf("  vary-both:    %.3f\n", mean(zmean[both_var])))
  cat(sprintf("  vary-neither: %.3f\n", mean(zmean[neither])))

  ## Per-tip "informativeness" inflation: for each tip, fraction of characters where
  ## that tip is the sole 1 (autapomorphy) and z==eco for it.
  cat("\nz mean overall:", mean(zmean), "\n")
}

## Save phi trace for plot
out_dir <- dirname(path)
saveRDS(list(phi = if (length(phi_col)) S[, phi_col[1]] else NULL,
             zmean = if (exists("zmean")) zmean else NULL,
             mono = c(A = mono(A), B = mono(B), C = mono(C), D = mono(D))),
        file.path(out_dir, "diagnose-summary.rds"))
cat("\nSaved diagnose-summary.rds\n")
