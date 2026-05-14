# Diagnose mixing under empirical_geometric prior on a small fast dataset.
# Reports per-move acceptance rates so we can confirm whether the `mh_p`
# multiplicative move is the bottleneck.

suppressPackageStartupMessages({
  devtools::load_all(quiet = TRUE)
  library(ape)
})

set.seed(2026)

nTip <- 12L
nCharTotal <- 120L
true_tree <- rtree(nTip, tip.label = paste0("t", seq_len(nTip)))
true_tree$edge.length <- runif(nrow(true_tree$edge), 0.05, 0.25)
true_tree <- TreeTools::Preorder(true_tree)

kTrueAll <- sample(c(2L, 3L, 4L, 5L), nCharTotal, replace = TRUE,
                    prob = c(0.40, 0.30, 0.20, 0.10))
sim_mat <- matrix(NA_integer_, nTip, nCharTotal,
                   dimnames = list(true_tree$tip.label, NULL))
for (ch in seq_len(nCharTotal)) {
  kTrue <- kTrueAll[ch]
  node_states <- integer(2 * nTip - 1)
  rootIdx <- nTip + 1L
  node_states[rootIdx] <- sample.int(kTrue, 1L) - 1L
  edges <- true_tree$edge
  el    <- true_tree$edge.length
  for (e in rev(seq_len(nrow(edges)))) {
    pa <- edges[e, 1]; ch2 <- edges[e, 2]; t <- el[e]
    pSame <- 1 / kTrue + (1 - 1 / kTrue) * exp(-kTrue * t)
    if (runif(1) < pSame) {
      node_states[ch2] <- node_states[pa]
    } else {
      node_states[ch2] <- sample(setdiff(seq.int(0, kTrue - 1L),
                                           node_states[pa]), 1L)
    }
  }
  sim_mat[, ch] <- node_states[seq_len(nTip)]
}
variable <- apply(sim_mat, 2, function(x) length(unique(x)) > 1)
sim_mat <- sim_mat[, variable, drop = FALSE]
pd <- TreeTools::MatrixToPhyDat(sim_mat)
mkd <- MkPrimeData(pd)
start_tree <- TreeTools::NJTree(pd, edgeLengths = TRUE)

cat(sprintf("Variable characters: %d\n", ncol(sim_mat)))
cat("kObs distribution: "); print(table(mkd$kObs))

cat("\n=== Running empirical_geometric arm (long enough to develop ESS) ===\n")
set.seed(7)
t0 <- Sys.time()
res <- RunMkPrime(
  mkd, start_tree,
  model = MkPrimeModel(coding = "variable",
                        kPrimePrior = "empirical_geometric",
                        expSteps = sum(true_tree$edge.length)),
  mcmc = MkPrimeMCMC(nIter = 30000L, thin = 30L,
                     maxWarmup = 8000L, minWarmup = 8000L,
                     autoTune = FALSE,
                     nRuns = 2L, nChains = 2L,
                     progressFn = function(...) invisible())
)
cat(sprintf("Wall time: %s\n", format(Sys.time() - t0, digits = 3)))
saveRDS(res, "data-raw/diagnose_eg_res.rds")

cat("\n--- Per-move acceptance ---\n")
print(res$acceptance)

cat("\n--- p posterior ---\n")
pcols <- grep("^p$|^p_", colnames(res$samples), value = TRUE)
if (length(pcols)) {
  p_samples <- res$samples[, pcols[1], drop = TRUE]
  cat(sprintf("  mean=%.4f  sd=%.4f  range=[%.4f, %.4f]\n",
              mean(p_samples), sd(p_samples),
              min(p_samples), max(p_samples)))
} else {
  cat("  p column not found in samples\n")
}

cat("\n--- Convergence ---\n")
cat(sprintf("  Stop reason: %s\n", res$stop_reason %||% "(unknown)"))
if (!is.null(res$ess_min)) {
  cat(sprintf("  Min ESS: %.1f\n", res$ess_min))
}
cat(sprintf("  Samples kept: %d\n", nrow(res$samples)))
