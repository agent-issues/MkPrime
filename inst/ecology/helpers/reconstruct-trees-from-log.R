# reconstruct-trees-from-log.R ----------------------------------------------
# Salvage helper for runs created BEFORE the brColStart fix (commit be1bdbb).
#
# Bug background: when ecologyAware=TRUE, the brColStart formula in
# RunMkPrime() / ResumeMkPrime() omitted the ecology column block
# (phi, pi0, theta_*). The MCMC log file (tree_length + br_*) is correct,
# but the stored phylo trees in res$trees have their first nEco edges set
# to kPrime-derived garbage; only the remaining edges are correct.
#
# This helper rebuilds res$trees with edge.length = tree_length * br_*
# read directly from the log. Use it on any res object produced before
# the fix. Trees produced AFTER the fix already satisfy this invariant
# (enforced by test STREAM-006).
#
# Usage:
#   source("inst/ecology/helpers/reconstruct-trees-from-log.R")
#   res$trees <- ReconstructTreesFromLog(res)
#
# Or, for a stored RDS wrapper:
#   r <- readRDS(path)
#   r$res$trees <- ReconstructTreesFromLog(r$res)
#   if (!is.null(r$relab)) r$relab$trees <- ReconstructTreesFromLog(r$relab)
#   saveRDS(r, path)

ReconstructTreesFromLog <- function(res, logFile = res$logFile) {
  stopifnot(file.exists(logFile))
  samp <- MkPrime::ReadMkLog(logFile)
  br_cols <- grep("^br_", colnames(samp))
  nEdge_log <- length(br_cols)
  nEdge_tree <- nrow(res$trees[[1]]$edge)
  if (nEdge_log != nEdge_tree) {
    stop(sprintf("br_* count (%d) does not match tree edge count (%d)",
                 nEdge_log, nEdge_tree))
  }
  trees <- res$trees
  n <- min(length(trees), nrow(samp))
  for (i in seq_len(n)) {
    tl <- samp[i, "tree_length"]
    br <- as.numeric(samp[i, br_cols])
    trees[[i]]$edge.length <- tl * br
  }
  trees
}
