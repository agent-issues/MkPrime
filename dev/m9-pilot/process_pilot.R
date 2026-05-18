#!/usr/bin/env Rscript
# Process the t_kv pilot output: load trees, prune to wcTree taxa, compute CID.
# Confirms the pilot -> CID-vs-WCT pipeline works on partial samples.

suppressPackageStartupMessages({
  library(ape)
  library(TreeTools)
  library(TreeDist)
})

pid <- "07203"
treesDir <- file.path("dev", "m9-pilot", "syab07203")

wcTrees <- setNames(unclass(ape::read.tree(
  "../neotrans/inst/wct/wellCorroboratedTrees.nwk")),
  paste0("0720", 0:6))
wcTree <- wcTrees[[pid]]
outgroup <- c("Didelphis", "Macropus")

read_rb_trees <- function(f) {
  # RB .trees is TSV: Iteration Posterior Likelihood Prior phylogeny
  # (no NEXUS header). The phylogeny column carries [&index=N] node
  # comments that ape::read.tree tolerates if we leave them in.
  tab <- read.delim(f, stringsAsFactors = FALSE)
  trs <- lapply(tab$phylogeny, function(s) ape::read.tree(text = s))
  structure(trs, class = "multiPhylo")
}

t1 <- read_rb_trees(file.path(treesDir, "t_kv_run_1.trees"))
t2 <- read_rb_trees(file.path(treesDir, "t_kv_run_2.trees"))

cat(sprintf("Loaded %d trees (run_1) + %d trees (run_2)\n",
            length(t1), length(t2)))
cat(sprintf("First tree has %d tips\n", length(t1[[1]]$tip.label)))

commonTips <- intersect(t1[[1]]$tip.label, wcTree$tip.label)
cat(sprintf("Common tips with asher WCT: %d / wc=%d / sample=%d\n",
            length(commonTips), length(wcTree$tip.label),
            length(t1[[1]]$tip.label)))

# Prune both sides to common leaves, root MCMC trees on outgroup
wcPruned <- KeepTip(wcTree, commonTips) |> RootTree(outgroup)
prune_and_root <- function(trees) {
  lapply(trees, function(tr) {
    KeepTip(tr, commonTips) |> RootTree(outgroup)
  }) |> structure(class = "multiPhylo")
}
t1p <- prune_and_root(t1)
t2p <- prune_and_root(t2)

# Drop the first 25% of each run as informal burnin (RB's burnin
# directive applies before mcmc starts; we have only post-burnin samples
# anyway, but a residual burnin is conventional)
burn <- function(trs) {
  n <- length(trs)
  trs[seq(max(2, floor(n * 0.25)), n)]
}
t1b <- burn(t1p)
t2b <- burn(t2p)

cat(sprintf("Post-burnin: %d (run_1), %d (run_2)\n", length(t1b), length(t2b)))

cid1 <- ClusteringInfoDistance(t1b, wcPruned, normalize = TRUE)
cid2 <- ClusteringInfoDistance(t2b, wcPruned, normalize = TRUE)

cat(sprintf("CID vs asher (run_1): mean=%.4f sd=%.4f min=%.4f max=%.4f\n",
            mean(cid1), sd(cid1), min(cid1), max(cid1)))
cat(sprintf("CID vs asher (run_2): mean=%.4f sd=%.4f min=%.4f max=%.4f\n",
            mean(cid2), sd(cid2), min(cid2), max(cid2)))

# CID between runs (chain convergence sanity check)
pooled_cid <- ClusteringInfoDistance(
  c(sample(t1b, min(20, length(t1b))),
    sample(t2b, min(20, length(t2b)))),
  normalize = TRUE
)
cat(sprintf("\nInter-tree CID (run_1 vs run_2 pooled, 20+20 sample):\n"))
cat(sprintf("  Range: %.4f - %.4f, median %.4f\n",
            min(pooled_cid[lower.tri(pooled_cid)]),
            max(pooled_cid[lower.tri(pooled_cid)]),
            median(pooled_cid[lower.tri(pooled_cid)])))
