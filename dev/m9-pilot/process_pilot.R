#!/usr/bin/env Rscript
# Process M9 pilot/long-run output for one matrix: load trees, prune to wcTree
# taxa, compute CID-to-WCT for each model.
#
# Usage: Rscript dev/m9-pilot/process_pilot.R <pid>
# where <pid> is e.g. "07200" (matches both syab07200/ and the WCT list).

suppressPackageStartupMessages({
  library(ape)
  library(TreeTools)
  library(TreeDist)
})

args <- commandArgs(trailingOnly = TRUE)
pid <- if (length(args) >= 1) args[[1]] else "07203"
models <- c("by_nt_9v", "t_9v", "t_kv")
treesDir <- file.path("dev", "m9-pilot", paste0("syab", pid))

# Outgroup mapping mirrors ../neotrans/R/AsherSmithSetup.R
outgroups <- list(
  "07200" = "Ichthyornis",
  "07201" = "aaCrocodylia",
  "07202" = c("Tupaia", "Dermoptera"),
  "07203" = c("Didelphis", "Macropus"),
  "07204" = c("Monodelphis", "Sarcophilus"),
  "07205" = c("Monodelphis", "Sarcophilus"),
  "07206" = "Ornithorhynchus"
)
outgroup <- outgroups[[pid]]
stopifnot(!is.null(outgroup))

# Trees are named by position, in `outgroups` order; the outgroup check below
# catches a reordered file.
wcTreesRaw <- unclass(ape::read.tree(
  "../neotrans/inst/wct/wellCorroboratedTrees.nwk"))
stopifnot(length(wcTreesRaw) == length(outgroups))
wcTrees <- setNames(wcTreesRaw, names(outgroups))
wcTree <- wcTrees[[pid]]
stopifnot("outgroup missing from the WCT; trees out of order?" =
            all(outgroup %in% wcTree[["tip.label"]]))

de_zz <- function(tr) {
  # Mirrors neotrans::DeZZ: strip the "zz" prefix only. The species suffix is
  # kept until ToGenusTree, after the tips have been selected.
  lab <- tr[["tip.label"]]
  zz <- startsWith(lab, "zz")
  lab[zz] <- substr(lab[zz], 3, nchar(lab[zz]))
  tr[["tip.label"]] <- lab
  tr
}

GenusOf <- function(lab) sub("([^_]+)_.*", "\\1", lab, perl = TRUE)

# Keep one tip (the alphabetically first label) per genus in `genera`, then
# relabel to genus so tips match the WCT's genus-only labels. Collapsing labels
# first would give KeepTip duplicated tips when congeners are present.
ToGenusTree <- function(tr, genera) {
  lab <- tr[["tip.label"]]
  g <- GenusOf(lab)
  ord <- order(lab)
  keep <- lab[ord][g[ord] %in% genera & !duplicated(g[ord])]
  tr <- KeepTip(tr, keep)
  tr[["tip.label"]] <- GenusOf(tr[["tip.label"]])
  tr
}

read_rb_trees <- function(f) {
  tab <- read.delim(f, stringsAsFactors = FALSE)
  trs <- lapply(tab$phylogeny, function(s) de_zz(ape::read.tree(text = s)))
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
  commonTips <- intersect(GenusOf(t1[[1]]$tip.label), wcTree$tip.label)
  wcPruned <- KeepTip(wcTree, commonTips) |> RootTree(outgroup)
  prune_and_root <- function(trees) {
    lapply(trees, function(tr) {
      ToGenusTree(tr, commonTips) |> RootTree(outgroup)
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
  cat(sprintf("  CID vs WCT: mean=%.4f sd=%.4f min=%.4f max=%.4f\n",
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
