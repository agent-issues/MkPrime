# Convergence dry-run for Fix B.  Re-uses the original smoke dataset
# (8 tips, ~65 variable chars, truth Sigma u = 50) but pushes warmup
# 5x longer to see if the chain plateaus.

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
  ns <- integer(2 * nTip - 1); ns[nTip + 1L] <- sample.int(kTrue, 1L) - 1L
  for (e in rev(seq_len(nrow(true_tree$edge)))) {
    pa <- true_tree$edge[e, 1]; ch2 <- true_tree$edge[e, 2]
    t  <- true_tree$edge.length[e]
    pSame <- 1 / kTrue + (1 - 1 / kTrue) * exp(-kTrue * t)
    ns[ch2] <- if (runif(1) < pSame) ns[pa] else
                sample(setdiff(0:(kTrue - 1L), ns[pa]), 1L)
  }
  sim_mat[, ch] <- ns[seq_len(nTip)]
}
variable <- apply(sim_mat, 2, function(x) length(unique(x)) > 1)
sim_mat <- sim_mat[, variable, drop = FALSE]
kTrueAll <- kTrueAll[variable]
pd <- TreeTools::MatrixToPhyDat(sim_mat)
mkd <- MkPrimeData(pd)
unseen <- sum(kTrueAll - mkd$kObs)
start_tree <- TreeTools::NJTree(pd, edgeLengths = TRUE)

cat("Truth Sigma u =", unseen, "; nChar =", ncol(sim_mat), "\n\n")

set.seed(101)
t0 <- Sys.time()
suppressWarnings(suppressMessages({
  res <- RunMkPrime(
    mkd, start_tree,
    model = MkPrimeModel(coding = "variable",
                         kPrimePrior = "empirical_geometric",
                         expSteps = sum(true_tree$edge.length)),
    mcmc = MkPrimeMCMC(nIter = 15000L, thin = 20L,
                       maxWarmup = 10000L, minWarmup = 10000L,
                       autoTune = FALSE,
                       nRuns = 1L, nChains = 1L)
  )
}))
cat("Wall time:", format(Sys.time() - t0, digits = 3), "\n\n")

kpCols <- grep("^kPrime_", colnames(res$samples), value = TRUE)
kPrimeIdx <- as.integer(sub("kPrime_", "", kpCols))
ord <- order(kPrimeIdx)
sumU <- rowSums(res$samples[, kpCols, drop = FALSE]) -
        sum(mkd$kObs[mkd$type == "transformational"])
pTrace <- res$samples[, "p"]
logPostCol <- intersect(c("log_post", "logPost", "logP"),
                       colnames(res$samples))
logPost <- if (length(logPostCol) > 0L) res$samples[, logPostCol[[1]]] else NULL

# Quartile-based stability check
n <- length(sumU)
qs <- list(q1 = 1:(n %/% 4),
           q2 = (n %/% 4 + 1):(n %/% 2),
           q3 = (n %/% 2 + 1):(3 * n %/% 4),
           q4 = (3 * n %/% 4 + 1):n)

cat("Sample columns:\n  ", paste(colnames(res$samples), collapse = ", "), "\n\n")
cat("Quartile means (n =", n, "samples):\n")
cat(sprintf("           Q1       Q2       Q3       Q4\n"))
if (!is.null(logPost)) {
  cat(sprintf("logPost:   %7.2f  %7.2f  %7.2f  %7.2f\n",
              mean(logPost[qs$q1]), mean(logPost[qs$q2]),
              mean(logPost[qs$q3]), mean(logPost[qs$q4])))
}
cat(sprintf("Sigma u:   %7.2f  %7.2f  %7.2f  %7.2f  (truth %d)\n",
            mean(sumU[qs$q1]), mean(sumU[qs$q2]),
            mean(sumU[qs$q3]), mean(sumU[qs$q4]), unseen))
cat(sprintf("p:         %7.3f  %7.3f  %7.3f  %7.3f\n",
            mean(pTrace[qs$q1]), mean(pTrace[qs$q2]),
            mean(pTrace[qs$q3]), mean(pTrace[qs$q4])))

cat("\nFinal acceptance:\n")
print(round(res$acceptance, 3))
