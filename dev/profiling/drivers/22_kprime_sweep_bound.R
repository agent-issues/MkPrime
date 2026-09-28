# #22 (MK-09): does the k' Gibbs sweep still evaluate candidates beyond the
# per-partition K cap?
#
# The bound #22 proposes landed with #66 (commit 1638678, partKoMax in
# compute_per_kprime_log_lik). This driver measures, on Sun2018 under the
# plain geometric prior (the only arm truncated at K), (a) candidates
# evaluated per character against nEff = K - kObs + 1, across the (K, p) grid,
# and (b) one sweep's wall time at the capped K against K = 256
# (kMaxKprimeCand, i.e. what the pre-#66 loop enumerated).
#
# Usage: Rscript dev/profiling/drivers/22_kprime_sweep_bound.R [lib]

args <- commandArgs(TRUE)
if (length(args)) .libPaths(c(args[1], .libPaths()))
suppressPackageStartupMessages({
  library(MkPrime)
  library(TreeTools)
})

nexFile <- system.file("datasets/Sun2018.nex", package = "TreeSearch")
if (!nzchar(nexFile)) {
  nexFile <- file.path(Sys.getenv("TREESEARCH_SRC"),
                       "inst/datasets/Sun2018.nex")
}
sun <- ReadAsPhyDat(nexFile)
mkd <- suppressWarnings(MkPrimeData(sun))
tree <- Preorder(RenumberTips(NJTree(sun, edgeLengths = TRUE), names(sun)))
tree$edge.length[tree$edge.length <= 0] <- 1e-8
kObs <- mkd$kObs[mkd$type == "transformational"]

Setup <- function(truncK, p) {
  model <- MkPrime:::.FinalizeModel(
    MkPrimeModel(kPrimePrior = "geometric", kprimeTruncK = truncK),
    tree, mkd)
  state <- MkPrime:::.InitState(tree, mkd, model)
  state$p <- p
  dataPtr <- MkPrime:::.InitMcmcData(mkd, model)
  statePtr <- MkPrime:::.InitMcmcChain(state)
  MkPrime:::fill_partition_cache(dataPtr, statePtr)
  MkPrime:::allocate_cl_workspace(dataPtr, statePtr)
  list(data = dataPtr, state = statePtr)
}

grid <- expand.grid(K = c(10L, 30L, 100L, 200L), p = c(0.01, 0.1, 0.5))
grid$evaluated <- grid$nEff <- grid$over <- NA_integer_
for (i in seq_len(nrow(grid))) {
  s <- Setup(grid$K[i], grid$p[i])
  cand <- MkPrime:::kprime_sweep_candidates(s$data, s$state, 1)
  nEff <- pmax(grid$K[i] - kObs + 1L, 0L)
  grid$evaluated[i] <- sum(cand)
  grid$nEff[i] <- sum(nEff)
  grid$over[i] <- sum(pmax(cand - nEff, 0L))
}
print(grid)
cat(sprintf("characters evaluated past their cap, whole grid: %d\n",
            sum(grid$over)))

SweepSec <- function(truncK, p) {
  s <- Setup(truncK, p)
  median(replicate(5, system.time(
    MkPrime:::kprime_sweep_candidates(s$data, s$state, 1))[[3]]))
}
for (p in c(0.01, 0.1)) {
  cat(sprintf("p = %.2f: sweep at K = 30: %.3f s; at K = 256: %.3f s\n",
              p, SweepSec(30L, p), SweepSec(256L, p)))
}
