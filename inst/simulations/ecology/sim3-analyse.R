# sim3-analyse.R --------------------------------------------------------------
# Post-hoc analysis of an Sim 3 (or Sim 3b) MCMC result file.
# Uses TreeDist::ClusteringInfoDist (CID) as the primary tree-distance
# measure: it is a generalised information-theoretic distance that
# degrades gracefully on near-misses (unlike RF, which is binary at the
# bipartition level).  Reports per-chain CID distributions, plus exact
# bipartition probabilities and a "closer to truth vs closer to wrong"
# tally.
suppressPackageStartupMessages({
  library("TreeTools")
  library("phangorn")
  library("TreeDist")
})

resultFile <- if (length(commandArgs(trailingOnly = TRUE)) > 0) {
  commandArgs(trailingOnly = TRUE)[1]
} else {
  "inst/simulations/ecology/sim3-mcmc-result.rds"
}
cat("Reading result file:", resultFile, "\n")
res <- readRDS(resultFile)

trueTree  <- res$tree
wrongTree <- Preorder(ape::read.tree(text =
  paste0("(((A1,A2),(A3,A4)),((B1,B2),(B3,B4)),",
         "(((C1,C2),(C3,C4)),((D1,D2),(D3,D4))));")))

discard <- function(x) x[seq.int(ceiling(length(x) / 4) + 1L, length(x))]

scoreTrees <- function(trees, label) {
  if (length(trees) == 0L) {
    cat(label, ": no trees\n"); return(invisible())
  }
  class(trees) <- "multiPhylo"
  # CID is normalized to [0, 1]: 0 = identical, 1 = maximally different.
  cidTrue  <- as.numeric(TreeDist::ClusteringInfoDist(
    trees, trueTree, normalize = TRUE))
  cidWrong <- as.numeric(TreeDist::ClusteringInfoDist(
    trees, wrongTree, normalize = TRUE))
  closer <- ifelse(cidTrue < cidWrong, "truth",
            ifelse(cidTrue > cidWrong, "wrong", "tied"))
  cat(sprintf("\n== %s (n = %d) ==\n", label, length(trees)))
  cat(sprintf("  CID to TRUE  : mean=%.3f  median=%.3f  min=%.3f  max=%.3f\n",
              mean(cidTrue), median(cidTrue), min(cidTrue), max(cidTrue)))
  cat(sprintf("  CID to WRONG : mean=%.3f  median=%.3f  min=%.3f  max=%.3f\n",
              mean(cidWrong), median(cidWrong), min(cidWrong), max(cidWrong)))
  cat("  Closer to:\n"); print(table(closer))
  cat(sprintf("  Trees at CID=0 to TRUE : %d  (%.1f%%)\n",
              sum(cidTrue == 0), 100 * mean(cidTrue == 0)))
  cat(sprintf("  Trees at CID=0 to WRONG: %d  (%.1f%%)\n",
              sum(cidWrong == 0), 100 * mean(cidWrong == 0)))
  cat(sprintf("  Mean CID gap (wrong - truth): %.4f\n",
              mean(cidWrong - cidTrue)))
}

scoreTrees(discard(res$resBlind$trees), "BLIND")
scoreTrees(discard(res$resAware$trees), "AWARE")

# Bipartition-level: how often does each reference split appear in the
# posterior?
trueSplitTips  <- c(paste0("A", 1:4), paste0("C", 1:4))
wrongSplitTips <- c(paste0("A", 1:4), paste0("B", 1:4))
otherSplits <- list(
  trueB = c(paste0("B", 1:4), paste0("D", 1:4)),
  wrongCD = c(paste0("C", 1:4), paste0("D", 1:4))
)
# Use root-invariant HasBipartSplits() — see sim3-scoring.R. The old
# prop.part-based hasBipart() in this file was root-dependent.
source("inst/simulations/ecology/sim3-scoring.R")
hasBipart <- HasBipartSplits

reportSplits <- function(trees, label) {
  if (length(trees) == 0L) return(invisible())
  cat(sprintf("\n== %s bipartition support (n=%d) ==\n", label, length(trees)))
  cat(sprintf("  P(A,C) [TRUE half]                : %.3f\n",
              mean(hasBipart(trees, trueSplitTips))))
  cat(sprintf("  P(B,D) [TRUE other half]          : %.3f\n",
              mean(hasBipart(trees, otherSplits$trueB))))
  cat(sprintf("  P(A,B) [WRONG, convergent halves] : %.3f\n",
              mean(hasBipart(trees, wrongSplitTips))))
  cat(sprintf("  P(C,D) [WRONG other half]         : %.3f\n",
              mean(hasBipart(trees, otherSplits$wrongCD))))
}
reportSplits(discard(res$resBlind$trees), "BLIND")
reportSplits(discard(res$resAware$trees), "AWARE")
