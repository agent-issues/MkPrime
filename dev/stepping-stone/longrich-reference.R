# Long stepping-stone runs on Longrich2010 (20 tips, Mk with k = kObs) to pin
# its log marginal likelihood, against which shorter budgets are judged (#66).
# Each task runs one descending and one ascending ladder; with long enough
# stones the two converge.
#
# Usage, from the source root: Rscript dev/stepping-stone/longrich-reference.R
#   <seed> <nIter> [outDir]
suppressPackageStartupMessages({
  source("dev/red-team/heavy-tests/load-mkprime.R")
  LoadMkPrime(".")
})
args <- commandArgs(trailingOnly = TRUE)
seed <- as.integer(args[1])
nIter <- as.integer(args[2])
outDir <- if (length(args) > 2L) args[3] else "."

pd <- TreeSearch::inapplicable.phyData[["Longrich2010"]]
mkd <- suppressWarnings(MkPrimeData(pd))
trans <- which(mkd$type == "transformational")
knownStates <- setNames(mkd$kObs[trans], mkd$orig_col[trans])
mkd <- suppressWarnings(MkPrimeData(pd, knownStates = knownStates))

set.seed(seed)
elapsed <- system.time(
  ss <- mkp_stepping_stone(mkd, nStones = 20L, nIter = nIter,
                           warmup = nIter %/% 10L, nRuns = 2L,
                           verbose = FALSE)
)[["elapsed"]]
row <- data.frame(seed = seed, nIter = nIter,
                  descending = ss$run_log_marginal[1],
                  ascending = ss$run_log_marginal[2],
                  seconds = elapsed)
print(row)
write.csv(row, file.path(outDir, sprintf("longrich-%d-%d.csv", nIter, seed)),
          row.names = FALSE)
