# Diagnostic: is Sigma u = 32 on truth=0 a mixing or structural issue?
suppressPackageStartupMessages({
  devtools::load_all(quiet = TRUE)
  library(ape)
})

set.seed(7)
nTip <- 7L
nChar <- 60L
tree <- rtree(nTip, tip.label = paste0("t", seq_len(nTip)))
tree$edge.length <- runif(nrow(tree$edge), 0.05, 0.15)
tree <- TreeTools::Preorder(tree)

kTrue <- 2L
buildBinary <- function() {
  node_states <- integer(2 * nTip - 1)
  rootIdx <- nTip + 1L
  node_states[rootIdx] <- sample.int(kTrue, 1L) - 1L
  edges <- tree$edge
  el    <- tree$edge.length
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
  node_states[seq_len(nTip)]
}
sim_mat <- matrix(NA_integer_, nTip, 0L,
                   dimnames = list(tree$tip.label, NULL))
while (ncol(sim_mat) < nChar) {
  col <- buildBinary()
  if (length(unique(col)) == kTrue) {
    sim_mat <- cbind(sim_mat, col)
  }
}
sim_mat <- sim_mat[, seq_len(nChar)]
colnames(sim_mat) <- NULL
pd <- TreeTools::MatrixToPhyDat(sim_mat)
mkd <- MkPrimeData(pd)

set.seed(202)
suppressWarnings(suppressMessages({
  res <- RunMkPrime(
    mkd, tree,
    model = MkPrimeModel(coding = "variable",
                         kPrimePrior = "empirical_geometric",
                         expSteps = sum(tree$edge.length)),
    mcmc = MkPrimeMCMC(nIter = 10000L, thin = 20L,
                       maxWarmup = 5000L, minWarmup = 5000L,
                       autoTune = FALSE,
                       nRuns = 1L, nChains = 1L)
  )
}))
kpCols <- grep("^kPrime_", colnames(res$samples), value = TRUE)
kPrimeIdx <- as.integer(sub("kPrime_", "", kpCols))
ord <- order(kPrimeIdx)

p_trace <- res$samples[, "p"]
sigmaU_trace <- rowSums(res$samples[, kpCols, drop = FALSE]) -
                sum(mkd$kObs[mkd$type == "transformational"])

cat("Total iterations sampled:", length(p_trace), "\n")
cat("p_mean (overall):       ", mean(p_trace), "\n")
cat("p_mean (last half):     ", mean(p_trace[(length(p_trace)/2):length(p_trace)]), "\n")
cat("p_mean (first quarter): ", mean(p_trace[1:(length(p_trace)/4)]), "\n")
cat("p_mean (last quarter):  ", mean(p_trace[(3*length(p_trace)/4):length(p_trace)]), "\n")
cat("\n")
cat("Sigma u (overall):      ", mean(sigmaU_trace), "\n")
cat("Sigma u (last half):    ", mean(sigmaU_trace[(length(sigmaU_trace)/2):length(sigmaU_trace)]), "\n")
cat("Sigma u (first quarter):", mean(sigmaU_trace[1:(length(sigmaU_trace)/4)]), "\n")
cat("Sigma u (last quarter): ", mean(sigmaU_trace[(3*length(sigmaU_trace)/4):length(sigmaU_trace)]), "\n")
cat("\n")
cat("mh_logit_p acceptance:  ", res$acceptance[["mh_logit_p"]], "\n")
cat("joint_p_kprime acc:     ", res$acceptance[["joint_p_kprime"]], "\n")
