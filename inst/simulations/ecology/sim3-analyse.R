# sim3-analyse.R --------------------------------------------------------------
# Post-hoc analysis of sim3-mcmc-result.rds: RF and clustering-info
# distances from each posterior tree to the true and the convergent-
# misled reference topologies.  Reports per-chain distributions and
# "closer to truth vs closer to wrong" counts as a softer metric than
# exact-bipartition matching.
suppressPackageStartupMessages({
  library("TreeTools")
  library("phangorn")
})

res <- readRDS("inst/simulations/ecology/sim3-mcmc-result.rds")

trueTree  <- res$tree
wrongTree <- Preorder(ape::read.tree(text =
  paste0("(((A1,A2),(A3,A4)),((B1,B2),(B3,B4)),",
         "(((C1,C2),(C3,C4)),((D1,D2),(D3,D4))));")))

discard <- function(x) x[seq.int(ceiling(length(x) / 4) + 1L, length(x))]

scoreTrees <- function(trees, label) {
  if (length(trees) == 0L) {
    cat(label, ": no trees\n"); return(invisible())
  }
  rfTrue  <- vapply(trees, function(tr) phangorn::RF.dist(tr, trueTree),
                    numeric(1))
  rfWrong <- vapply(trees, function(tr) phangorn::RF.dist(tr, wrongTree),
                    numeric(1))
  closer <- ifelse(rfTrue < rfWrong, "truth",
            ifelse(rfTrue > rfWrong, "wrong", "tied"))
  cat(sprintf("\n== %s (n = %d) ==\n", label, length(trees)))
  cat(sprintf("  RF to TRUE  : mean=%.2f  median=%.1f  min=%d  max=%d\n",
              mean(rfTrue), median(rfTrue), min(rfTrue), max(rfTrue)))
  cat(sprintf("  RF to WRONG : mean=%.2f  median=%.1f  min=%d  max=%d\n",
              mean(rfWrong), median(rfWrong), min(rfWrong), max(rfWrong)))
  cat("  Closer to:\n"); print(table(closer))
  cat(sprintf("  Trees at RF=0 to TRUE : %d  (%.1f%%)\n",
              sum(rfTrue == 0), 100 * mean(rfTrue == 0)))
  cat(sprintf("  Trees at RF=0 to WRONG: %d  (%.1f%%)\n",
              sum(rfWrong == 0), 100 * mean(rfWrong == 0)))
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
hasBipart <- function(treeList, splitTips) {
  vapply(treeList, function(tr) {
    cl <- ape::prop.part(tr)
    tips <- attr(cl, "labels")
    splitSet <- which(tips %in% splitTips)
    any(vapply(cl, function(p) setequal(p, splitSet), logical(1)))
  }, logical(1))
}

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
