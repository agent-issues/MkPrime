#!/usr/bin/env Rscript
# Process the 07203 pilot output: load trees for each model, prune to wcTree
# taxa, compute CID-to-WCT. Confirms the pilot -> CID-vs-WCT pipeline.

suppressPackageStartupMessages({
  library(ape)
  library(TreeTools)
  library(TreeDist)
})

pid <- "07203"
models <- c("by_nt_9v", "t_9v", "t_kv")
treesDir <- file.path("dev", "m9-pilot", "syab07203")

wcTrees <- setNames(unclass(ape::read.tree(
  "../neotrans/inst/wct/wellCorroboratedTrees.nwk")),
  paste0("0720", 0:6))
wcTree <- wcTrees[[pid]]
outgroup <- c("Didelphis", "Macropus")

read_rb_trees <- function(f) {
  # RB .trees is TSV: Iteration Posterior Likelihood Prior phylogeny
  tab <- read.delim(f, stringsAsFactors = FALSE)
  trs <- lapply(tab$phylogeny, function(s) ape::read.tree(text = s))
  structure(trs, class = "multiPhylo")
}

burn <- function(trs) {
  n <- length(trs)
  trs[seq(max(2, floor(n * 0.25)), n)]
}

process_model <- function(model) {
  f1 <- file.path(treesDir, paste0(model, "_run_1.trees"))
  f2 <- file.path(treesDir, paste0(model, "_run_2.trees"))
  if (!file.exists(f1) || !file.exists(f2)) {
    cat(sprintf("\n[%s] missing tree files, skipping\n", model))
    return(NULL)
  }
  t1 <- read_rb_trees(f1)
  t2 <- read_rb_trees(f2)
  commonTips <- intersect(t1[[1]]$tip.label, wcTree$tip.label)
  wcPruned <- KeepTip(wcTree, commonTips) |> RootTree(outgroup)
  prune_and_root <- function(trees) {
    lapply(trees, function(tr) {
      KeepTip(tr, commonTips) |> RootTree(outgroup)
    }) |> structure(class = "multiPhylo")
  }
  t1b <- burn(prune_and_root(t1))
  t2b <- burn(prune_and_root(t2))
  cid1 <- ClusteringInfoDistance(t1b, wcPruned, normalize = TRUE)
  cid2 <- ClusteringInfoDistance(t2b, wcPruned, normalize = TRUE)
  cidAll <- c(cid1, cid2)
  cat(sprintf("\n[%s] n=%d+%d post-burnin trees, common tips=%d/%d\n",
              model, length(t1b), length(t2b),
              length(commonTips), length(wcTree$tip.label)))
  cat(sprintf("  CID vs asher: mean=%.4f sd=%.4f min=%.4f max=%.4f\n",
              mean(cidAll), sd(cidAll), min(cidAll), max(cidAll)))
  data.frame(
    matrix = pid, model = model,
    n = length(cidAll),
    cid_mean = mean(cidAll), cid_sd = sd(cidAll),
    cid_min = min(cidAll), cid_max = max(cidAll),
    stringsAsFactors = FALSE
  )
}

results <- do.call(rbind, lapply(models, process_model))
cat("\n=== Summary table ===\n")
print(results, row.names = FALSE)

saveRDS(results, file.path(treesDir, "cid_summary.rds"))
cat(sprintf("\nSaved: %s\n", file.path(treesDir, "cid_summary.rds")))
