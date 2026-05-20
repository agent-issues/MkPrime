## audit-rerun.R --------------------------------------------------------------
## Re-examine the v4cross-b aware-chain "mode-trap" diagnosis against
## root-invariant (corrected) bipartition scoring.
##
## Usage:
##   Rscript audit-rerun.R <result.rds> <sim3-scoring.R> <out-dir>
##
## Outputs:
##   <out-dir>/audit-rerun.txt   - human-readable log of all findings
##   <out-dir>/audit-rerun.rds   - per-sample data frame
suppressPackageStartupMessages({
  library(MkPrime)
  library(ape)
  library(TreeTools)
  library(TreeDist)
})

args <- commandArgs(trailingOnly = TRUE)
path <- if (length(args) >= 1) args[1] else
  "/nobackup/pjjg18/mkp-sim3-v4cross-b/results/aware-result.rds"
scoringSrc <- if (length(args) >= 2) args[2] else
  "/nobackup/pjjg18/mkp-sim3-v4cross-b/diagnose/sim3-scoring.R"
outDir <- if (length(args) >= 3) args[3] else
  "/nobackup/pjjg18/mkp-sim3-v4cross-b/diagnose"

dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
source(scoringSrc, local = TRUE)

logFile <- file.path(outDir, "audit-rerun.txt")
sink(logFile, split = TRUE)
on.exit(sink(), add = TRUE)

cat("=== audit-rerun.R ===\n")
cat("date:", format(Sys.time()), "\n")
cat("file:", path, "\n")

res <- readRDS(path)
cat("relabelled attr:", isTRUE(attr(res, "relabelled")), "\n")
cat("ecologyAware:", isTRUE(res$model$ecologyAware), "\n")
cat("magnitudeMode:", res$model$magnitudeMode, "\n")
cat("nSamples:", res$nSamples, " stop:", res$stop_reason, "\n\n")

tn <- res$data$taxon_names
cat("taxa:", tn, "\n")
A <- tn[1:4]; C <- tn[5:8]; B <- tn[9:12]; D <- tn[13:16]
cat("A:", A, "\nB:", B, "\nC:", C, "\nD:", D, "\n\n")

trees <- res$trees
class(trees) <- "multiPhylo"
S <- res$samples
phi_col <- grep("^phi", colnames(S), value = TRUE)[1]
ll_col  <- grep("^log_lik|^ll$|log_likelihood", colnames(S), value = TRUE)[1]
phi <- if (!is.na(phi_col)) S[, phi_col] else rep(NA_real_, length(trees))
ll  <- if (!is.na(ll_col))  S[, ll_col]  else rep(NA_real_, length(trees))
cat("phi col:", phi_col, "   ll col:", ll_col, "\n")
cat("nSamples (trees):", length(trees), "  nrow(S):", nrow(S), "\n\n")

n <- length(trees)

## -----------------------------------------------------------------------------
## Per-sample diagnostics
## -----------------------------------------------------------------------------
cat("Computing per-sample bipart tests + root configurations ...\n")
legacy_B  <- HasBipartLegacy(trees, B)
correct_B <- HasBipartSplits(trees, B)

## Root config key per sample
root_key <- vapply(trees, function(tr) {
  tipLab <- tr$tip.label
  nT <- length(tipLab)
  rootNode <- nT + 1L
  rootChildren <- tr$edge[tr$edge[, 1L] == rootNode, 2L]
  if (length(rootChildren) < 2L) return(NA_character_)
  descTips <- function(node) {
    if (node <= nT) return(tipLab[node])
    kids <- tr$edge[tr$edge[, 1L] == node, 2L]
    unlist(lapply(kids, descTips))
  }
  side1 <- sort(descTips(rootChildren[1L]))
  side2 <- sort(descTips(rootChildren[2L]))
  key1 <- paste(side1, collapse = ",")
  key2 <- paste(side2, collapse = ",")
  if (key1 < key2) key1 else key2
}, character(1))

root_inside_B <- RootInsideSet(trees, B)

df <- data.frame(
  idx        = seq_len(n),
  phi        = phi,
  ll         = ll,
  legacy_B   = legacy_B,
  correct_B  = correct_B,
  root_inB   = root_inside_B,
  root_key   = root_key,
  stringsAsFactors = FALSE
)
saveRDS(df, file.path(outDir, "audit-rerun.rds"))

cat("\n=== Step 1: Reproduce 'contiguous block 1..83' check ===\n")
cat(sprintf("Total samples: %d\n", n))
cat(sprintf("LEGACY    B-mono count: %d  (%.3f)\n",   sum(legacy_B),  mean(legacy_B)))
cat(sprintf("CORRECTED B-mono count: %d  (%.3f)\n",   sum(correct_B), mean(correct_B)))

## Time-series block structure: where do True/False segments sit?
legacy_runs  <- rle(legacy_B)
correct_runs <- rle(correct_B)
cat("\nLegacy RLE (value, length):\n")
print(data.frame(value = legacy_runs$values, length = legacy_runs$lengths))

cat("\nCorrected RLE (value, length):\n")
print(data.frame(value = correct_runs$values, length = correct_runs$lengths))

## Show first/last 5 indices of each "block" for legacy
cat("\nLegacy B-mono indices (first 20):", head(which(legacy_B), 20), "\n")
cat("Legacy B-mono indices (last 5):",   tail(which(legacy_B), 5),  "\n")
cat("Legacy B-broken indices (first 5):", head(which(!legacy_B), 5), "\n")

cat("\nCorrected B-mono indices (first 20):", head(which(correct_B), 20), "\n")
cat("Corrected B-mono indices (last 5):",   tail(which(correct_B), 5),  "\n")
cat("Corrected B-broken indices (first 5):", head(which(!correct_B), 5), "\n")

## -----------------------------------------------------------------------------
## Where do legacy and corrected disagree? -> root inside B sets
## -----------------------------------------------------------------------------
cat("\n=== Discrepancy analysis: legacy vs corrected B-mono ===\n")
disagree <- legacy_B != correct_B
cat(sprintf("Samples where legacy != corrected: %d / %d (%.3f)\n",
            sum(disagree), n, mean(disagree)))
cat("In disagreement set: legacy=TRUE only:", sum(legacy_B  & !correct_B), "\n")
cat("In disagreement set: corrected=TRUE only:", sum(!legacy_B & correct_B), "\n")
cat("Root inside B (TRUE) by case:\n")
print(table(legacy = legacy_B, corrected = correct_B, root_in_B = root_inside_B,
            useNA = "ifany"))

## -----------------------------------------------------------------------------
## Phi / ll conditional on the TWO scorings
## -----------------------------------------------------------------------------
cat("\n=== Phi/ll conditional on monophyly (LEGACY) ===\n")
cat(sprintf("LEGACY  B-mono   : n=%4d  phi.med=%.3f  phi.IQR=[%.3f,%.3f]  ll.mean=%.2f\n",
            sum(legacy_B),
            median(phi[legacy_B]),  quantile(phi[legacy_B],  .25), quantile(phi[legacy_B],  .75),
            mean(ll[legacy_B])))
cat(sprintf("LEGACY  B-broken : n=%4d  phi.med=%.3f  phi.IQR=[%.3f,%.3f]  ll.mean=%.2f\n",
            sum(!legacy_B),
            median(phi[!legacy_B]), quantile(phi[!legacy_B], .25), quantile(phi[!legacy_B], .75),
            mean(ll[!legacy_B])))

cat("\n=== Phi/ll conditional on monophyly (CORRECTED) ===\n")
if (sum(correct_B) > 0L && sum(!correct_B) > 0L) {
  cat(sprintf("CORR    B-mono   : n=%4d  phi.med=%.3f  phi.IQR=[%.3f,%.3f]  ll.mean=%.2f\n",
              sum(correct_B),
              median(phi[correct_B]),  quantile(phi[correct_B],  .25), quantile(phi[correct_B],  .75),
              mean(ll[correct_B])))
  cat(sprintf("CORR    B-broken : n=%4d  phi.med=%.3f  phi.IQR=[%.3f,%.3f]  ll.mean=%.2f\n",
              sum(!correct_B),
              median(phi[!correct_B]), quantile(phi[!correct_B], .25), quantile(phi[!correct_B], .75),
              mean(ll[!correct_B])))
} else {
  cat("CORR: only one class present.\n")
  cls <- if (sum(correct_B) > 0) "all B-mono" else "all B-broken"
  cat("  ", cls, "\n")
  cat(sprintf("    n=%d  phi.med=%.3f  ll.mean=%.2f\n",
              n, median(phi), mean(ll)))
}

## -----------------------------------------------------------------------------
## Root config breakdown
## -----------------------------------------------------------------------------
cat("\n=== Root configurations (canonicalised smaller side) ===\n")
rootTab <- sort(table(root_key), decreasing = TRUE)
cat("Number of distinct root configs:", length(rootTab), "\n")
print(rootTab)

cat("\n=== Root config vs CORRECTED B-mono ===\n")
ct <- table(root_key, corrected_B = correct_B)
print(ct)

cat("\n=== Root config vs LEGACY B-mono ===\n")
ctL <- table(root_key, legacy_B = legacy_B)
print(ctL)

## Root-config transition: find the index of root-config change
cat("\n=== Root config over time (first 90 samples) ===\n")
rkey_short <- substr(root_key, 1, 40)
cat("Sample 1..10 root keys:\n")
print(rkey_short[1:10])
cat("Sample 80..95 root keys:\n")
print(rkey_short[80:min(95, n)])

## Where does the root config "switch"?
diff_pos <- which(root_key[-1] != root_key[-n]) + 1L  # idx where a switch begins
cat("\nNumber of root-config switches:", length(diff_pos), "\n")
cat("First 20 switch positions:", head(diff_pos, 20), "\n")

## -----------------------------------------------------------------------------
## Unique topology count: legacy vs root-invariant
## -----------------------------------------------------------------------------
cat("\n=== Unique tree counts ===\n")
topo_hash_col <- grep("^topo_hash", colnames(S), value = TRUE)
if (length(topo_hash_col)) {
  th <- S[, topo_hash_col[1]]
  cat("Unique topo_hash values (legacy chain field):", length(unique(th)), "\n")
}

## Root-invariant: enumerate distinct unrooted bipartition sets per tree
cat("Computing unique unrooted topologies (via Splits sets) ...\n")
sigKeys <- vapply(trees, function(tr) {
  spl <- TreeTools::as.Splits(tr, tipLabels = trees[[1]]$tip.label)
  m <- as.logical(spl)
  if (is.null(dim(m))) m <- matrix(m, nrow = 1L)
  # Canonicalise each row: prefer the half containing tip 1
  canonRows <- t(apply(m, 1L, function(r) if (r[1L]) r else !r))
  # Drop trivial splits (all TRUE, all FALSE, or single TRUE)
  keep <- apply(canonRows, 1L, function(r) {
    s <- sum(r)
    s > 1L && s < length(r) - 0L && !all(r) && !all(!r)
  })
  canonRows <- canonRows[keep, , drop = FALSE]
  # Sort rows for stable hash
  rowKeys <- apply(canonRows, 1L, function(r) paste(as.integer(r), collapse = ""))
  paste(sort(rowKeys), collapse = "|")
}, character(1))
cat("Unique unrooted topologies (Splits-based):", length(unique(sigKeys)), "\n")

## -----------------------------------------------------------------------------
## Cross-tab summary
## -----------------------------------------------------------------------------
cat("\n=== Cross-tab: root_in_B vs corrected B-mono ===\n")
print(table(root_in_B = root_inside_B, corrected_B = correct_B))

cat("\n=== Cross-tab: root_in_B vs legacy B-mono ===\n")
print(table(root_in_B = root_inside_B, legacy_B = legacy_B))

cat("\nDone.\n")
