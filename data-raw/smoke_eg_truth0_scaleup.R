# Structural-floor scaling test: 12 tips, 200 binary chars, short
# branches.  Truth Sigma u = 0.  Pre-scale-up check that Fix B's
# residual floor (~0.5/char on the 60-char diagnostic) doesn't blow up
# on a representative-sized binary dataset.

suppressPackageStartupMessages({
  devtools::load_all(quiet = TRUE)
  library(ape)
})

set.seed(13)
nTip <- 12L
nChar <- 200L
kTrue <- 2L
tree <- rtree(nTip, tip.label = paste0("t", seq_len(nTip)))
tree$edge.length <- runif(nrow(tree$edge), 0.05, 0.15)
tree <- TreeTools::Preorder(tree)

buildBinary <- function() {
  ns <- integer(2 * nTip - 1); ns[nTip + 1L] <- sample.int(kTrue, 1L) - 1L
  for (e in rev(seq_len(nrow(tree$edge)))) {
    pa <- tree$edge[e, 1]; ch <- tree$edge[e, 2]; t <- tree$edge.length[e]
    pSame <- 1 / kTrue + (1 - 1 / kTrue) * exp(-kTrue * t)
    ns[ch] <- if (runif(1) < pSame) ns[pa] else
               sample(setdiff(0:(kTrue - 1L), ns[pa]), 1L)
  }
  ns[seq_len(nTip)]
}

sim <- matrix(NA_integer_, nTip, 0L,
              dimnames = list(tree$tip.label, NULL))
attempts <- 0L
while (ncol(sim) < nChar && attempts < 20L * nChar) {
  attempts <- attempts + 1L
  col <- buildBinary()
  if (length(unique(col)) == kTrue) sim <- cbind(sim, col)
}
sim <- sim[, seq_len(nChar)]
colnames(sim) <- NULL
pd <- TreeTools::MatrixToPhyDat(sim)
mkd <- MkPrimeData(pd)
cat("Dataset: ", nTip, " tips, ", nChar, " binary chars (truth Sigma u = 0)\n",
    sep = "")

set.seed(202)
t0 <- Sys.time()
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
cat("Wall time:", format(Sys.time() - t0, digits = 3), "\n\n")

kpCols <- grep("^kPrime_", colnames(res$samples), value = TRUE)
sumU <- rowSums(res$samples[, kpCols, drop = FALSE]) -
        sum(mkd$kObs[mkd$type == "transformational"])
pTrace <- res$samples[, "p"]

n <- length(sumU)
qs <- list(q1 = 1:(n %/% 4),
           q4 = (3 * n %/% 4 + 1):n)

cat("Quartile means (n =", n, "samples):\n")
cat(sprintf("Sigma u:   Q1 = %7.2f   Q4 = %7.2f   (truth 0)\n",
            mean(sumU[qs$q1]), mean(sumU[qs$q4])))
cat(sprintf("Sigma u / nChar at Q4: %.3f\n",
            mean(sumU[qs$q4]) / nChar))
cat(sprintf("p:         Q1 = %7.3f   Q4 = %7.3f\n",
            mean(pTrace[qs$q1]), mean(pTrace[qs$q4])))

cat("\nAcceptance:\n")
print(round(res$acceptance, 3))
