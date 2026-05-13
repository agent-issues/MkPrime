# Quick sweep over (a, b) Beta hyperprior to calibrate Fix B.
suppressPackageStartupMessages({
  devtools::load_all(quiet = TRUE)
  library(ape)
})

set.seed(2026)
nTip <- 8L
nCharTotal <- 80L
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
kTrueAll <- kTrueAll[variable]

pd <- TreeTools::MatrixToPhyDat(sim_mat)
mkd <- MkPrimeData(pd)
unseen <- sum(kTrueAll - mkd$kObs)
cat("Truth Sigma u =", unseen, " nChar =", ncol(sim_mat), "\n")

start_tree <- TreeTools::NJTree(pd, edgeLengths = TRUE)

runOne <- function(a, b, seed = 101) {
  set.seed(seed)
  suppressWarnings(suppressMessages({
    res <- RunMkPrime(
      mkd, start_tree,
      model = MkPrimeModel(coding = "variable",
                           kPrimePrior = "empirical_geometric",
                           kprimeHyperA = a, kprimeHyperB = b,
                           expSteps = sum(true_tree$edge.length)),
      mcmc = MkPrimeMCMC(nIter = 3000L, thin = 20L,
                         maxWarmup = 1500L, minWarmup = 1500L,
                         autoTune = FALSE,
                         nRuns = 1L, nChains = 1L)
    )
  }))
  kpCols <- grep("^kPrime_", colnames(res$samples), value = TRUE)
  kPostMean <- colMeans(res$samples[, kpCols, drop = FALSE])
  kPrimeIdx <- as.integer(sub("kPrime_", "", kpCols))
  ord <- order(kPrimeIdx)
  uPostMean <- (kPostMean - mkd$kObs)[ord]
  pPostMean <- mean(res$samples[, "p"])
  list(SumU = sum(uPostMean), pMean = pPostMean,
       mh_logit_p_acc = res$acceptance[["mh_logit_p"]])
}

cands <- list(
  c(4, 1),
  c(7, 1),
  c(10, 1),
  c(15, 1),
  c(20, 1),
  c(8, 2),
  c(12, 2)
)
cat(sprintf("%6s %6s %10s %10s %12s\n", "a", "b", "SumU", "p_mean", "mh_logit_p"))
for (ab in cands) {
  r <- runOne(ab[1], ab[2])
  cat(sprintf("%6.1f %6.1f %10.2f %10.3f %12.3f\n",
              ab[1], ab[2], r$SumU, r$pMean, r$mh_logit_p_acc))
}
