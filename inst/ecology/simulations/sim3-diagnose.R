suppressPackageStartupMessages({library("TreeTools"); library("TreeDist")})
res <- readRDS("inst/ecology/simulations/sim3-mcmc-result.rds")
samp <- read.table("sim3-aware.log", header = TRUE, sep = "\t", comment.char = "#")
cat("== Sim3 aware chain diagnostics ==\n")
cat(sprintf("nSamples=%d  iters %d to %d\n", nrow(samp),
            samp$generation[1], tail(samp$generation, 1)))
cat(sprintf("tree_length: median=%.2f IQR=[%.2f,%.2f] min=%.2f max=%.2f (truth ~9.1)\n",
            median(samp$tree_length),
            stats::quantile(samp$tree_length, 0.25, names = FALSE),
            stats::quantile(samp$tree_length, 0.75, names = FALSE),
            min(samp$tree_length), max(samp$tree_length)))
cat(sprintf("logL: median=%.0f IQR=[%.0f,%.0f]\n",
            median(samp$log_likelihood),
            stats::quantile(samp$log_likelihood, 0.25, names = FALSE),
            stats::quantile(samp$log_likelihood, 0.75, names = FALSE)))
cat(sprintf("phi: median=%.2f IQR=[%.2f,%.2f]\n",
            median(samp$phi),
            stats::quantile(samp$phi, 0.25, names = FALSE),
            stats::quantile(samp$phi, 0.75, names = FALSE)))
cat(sprintf("pi0: median=%.3f IQR=[%.3f,%.3f]\n",
            median(samp$pi0),
            stats::quantile(samp$pi0, 0.25, names = FALSE),
            stats::quantile(samp$pi0, 0.75, names = FALSE)))
trees <- res$resAware$trees[
  seq.int(ceiling(length(res$resAware$trees) / 4) + 1L,
          length(res$resAware$trees))]
class(trees) <- "multiPhylo"
cat(sprintf("Distinct topologies (rough): %d of %d samples\n",
            length(unique(vapply(trees, function(t) digest::digest(
              t$edge), character(1)))), length(trees)))
# Per-bipartition frequency for each TRUE split, plus the wrong split.
trueT <- res$tree
trueSplits <- TreeTools::as.Splits(trueT)
postSplits <- TreeTools::as.Splits(trees, tipLabels = trueT$tip.label)
sf <- TreeTools::SplitFrequency(trueT, postSplits)
cat("\nPosterior frequency of each TRUE bipartition (",
    length(trueSplits), " internal splits in truth):\n")
print(round(sf / length(trees), 3))
cat("\nMean: ", round(mean(sf / length(trees)), 3), "\n")
# Same on blind for comparison.
treesB <- res$resBlind$trees[
  seq.int(ceiling(length(res$resBlind$trees) / 4) + 1L,
          length(res$resBlind$trees))]
class(treesB) <- "multiPhylo"
postSplitsB <- TreeTools::as.Splits(treesB, tipLabels = trueT$tip.label)
sfB <- TreeTools::SplitFrequency(trueT, postSplitsB)
cat("\nBLIND posterior frequency of each TRUE bipartition:\n")
print(round(sfB / length(treesB), 3))
cat("Mean: ", round(mean(sfB / length(treesB)), 3), "\n")
# Also tree-length on blind log.
sampB <- read.table("sim3-blind.log", header = TRUE, sep = "\t",
                    comment.char = "#")
cat(sprintf("\nBLIND tree_length: median=%.2f IQR=[%.2f,%.2f]\n",
            median(sampB$tree_length),
            stats::quantile(sampB$tree_length, 0.25, names = FALSE),
            stats::quantile(sampB$tree_length, 0.75, names = FALSE)))
cat(sprintf("BLIND logL: median=%.0f IQR=[%.0f,%.0f]\n",
            median(sampB$log_likelihood),
            stats::quantile(sampB$log_likelihood, 0.25, names = FALSE),
            stats::quantile(sampB$log_likelihood, 0.75, names = FALSE)))
