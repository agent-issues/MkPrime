# #346: cost of conditional CL rescaling on the ordinary path.
#
# Rescaling fires only for characters whose CLs fall below 2^-256, which
# Sun2018 never reaches, so everything timed here is the overhead of the
# check alone. Times each kernel family on Sun2018 (54 taxa, 225 chars),
# then a short RunMkPrime() run (driver 01's settings) for iter/s.
#
# Usage: Rscript dev/profiling/drivers/23_cl_rescaling.R <lib> <label> <rep>
# Prints CSV rows: label,rep,metric,value. Run it against two libraries in
# alternation, in the same job, to compare builds. MKP_PKG names a renamed
# build (build-agent.sh).

setTimeLimit(elapsed = 900, transient = FALSE)
args <- commandArgs(TRUE)
lib <- args[1]
label <- if (length(args) > 1) args[2] else "lib"
rep <- if (length(args) > 2) args[3] else "1"
pkg <- Sys.getenv("MKP_PKG", "MkPrime")
suppressPackageStartupMessages({
  library(pkg, lib.loc = lib, character.only = TRUE)
  library("TreeTools")
})
ns <- asNamespace(pkg)

Out <- function(metric, value) {
  cat(sprintf("%s,%s,%s,%.6g\n", label, rep, metric, value))
}

# Milliseconds per call of `fn()`, repeated until at least `minTime` s.
MsPerCall <- function(fn, minTime = 3) {
  fn()
  n <- 0L
  t0 <- proc.time()[["elapsed"]]
  repeat {
    for (i in 1:10) fn()
    n <- n + 10L
    el <- proc.time()[["elapsed"]] - t0
    if (el >= minTime) break
  }
  1000 * el / n
}

set.seed(5813)
pd <- ReadAsPhyDat(system.file("datasets/Sun2018.nex", package = "TreeSearch"))
tree <- NJTree(pd, edgeLengths = TRUE)
tree$edge.length[tree$edge.length <= 0] <- 1e-8
tree <- Preorder(tree)
tokens <- PhyDatToMatrix(pd)
states <- suppressWarnings(matrix(as.integer(tokens), nrow(tokens)))
states[is.na(states)] <- -1L
p <- tree$edge[, 1]
ch <- tree$edge[, 2]
el <- tree$edge.length
rates <- ns$header_lognormal_rates(0.6, qnorm((1:4 - 0.5) / 4))
k <- max(states) + 2L

Out("flat_jc_ms", MsPerCall(function() {
  ns$test_flat_jc_uncollapsed(p, ch, el, states, k, FALSE, 1)
}))
Out("flat_jc_acrv_ms", MsPerCall(function() {
  ns$test_flat_jc_uncollapsed(p, ch, el, states, k, TRUE, rates)
}))
Out("flat_jc_acrv_collapsed_ms", MsPerCall(function() {
  ns$test_flat_jc_collapsed(p, ch, el, states, k + 3L, k, TRUE, rates)
}))
Out("persite_ms", MsPerCall(function() {
  ns$test_persite_uncollapsed(p, ch, el, states, k, rates)
}))
Out("pruner_acrv_ms", MsPerCall(function() {
  ns$pruning_jc_acrv(p, ch, el, states, k, rep(1 / k, k), rates)
}))

mkd <- MkPrimeData(pd)
model <- ns$.FinalizeModel(MkPrimeModel(nCat = 4L), tree, mkd)
NewChain <- function() {
  dataPtr <- ns$.InitMcmcData(mkd, model)
  statePtr <- ns$.InitMcmcChain(ns$.InitState(tree, mkd, model))
  ns$fill_partition_cache(dataPtr, statePtr)
  ns$allocate_cl_workspace(dataPtr, statePtr)
  list(dataPtr = dataPtr, statePtr = statePtr)
}
chain <- NewChain()
moves <- c(nni = 5L, beta_simplex = 4L)
for (i in seq_along(moves)) {
  Out(paste0("move_", names(moves)[i], "_ms"), MsPerCall(function() {
    ns$do_move_cpp(chain$dataPtr, chain$statePtr, moves[[i]], 0L, 0.5, 2, 1L,
                   1.0)
  }))
}
chain <- NewChain()
edge <- ns$get_mcmc_state(chain$statePtr)$edge
pruneRow <- which(edge[, 1] != NTip(tree) + 1L & edge[, 2] > NTip(tree))[3]
Out("gibbs_spr_enumerate_ms", MsPerCall(function() {
  ns$gibbs_spr_enumerate_cpp(chain$dataPtr, chain$statePtr, pruneRow)
}))

nIter <- 800L
t0 <- proc.time()[["elapsed"]]
res <- RunMkPrime(
  mkd, tree,
  model = MkPrimeModel(coding = "variable", kPrimePrior = "empirical_geometric"),
  mcmc = MkPrimeMCMC(nIter = nIter, thin = 50L, maxWarmup = 200L,
                     minWarmup = 200L, autoTune = FALSE, nRuns = 1L,
                     nChains = 1L, maxTime = 600,
                     progressFn = function(...) invisible())
)
Out("runmkprime_iter_per_s", nIter / (proc.time()[["elapsed"]] - t0))
