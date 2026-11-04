#!/usr/bin/env Rscript
# Exact arm of the empirical-vs-geometric k' prior study (#179).
#
# For each simulated dataset, scores both priors' exact k' posterior on the
# true tree (score.R), with no MCMC: the effect of the prior alone, free of
# sampler noise and topology uncertainty. Seconds per dataset.
#
# Usage, from the repository root:
#   Rscript dev/red-team/heavy-tests/empirical-benefit/exact.R [nRep] [out]
#
# Writes one row per dataset x arm to `out` (default
# dev/red-team/heavy-tests/empirical-benefit/results/exact.rds).

suppressPackageStartupMessages({
  source("dev/red-team/heavy-tests/load-mkprime.R")
  LoadMkPrime(".")
})
source("dev/red-team/heavy-tests/empirical-benefit/forward.R")
source("dev/red-team/heavy-tests/empirical-benefit/score.R")

args <- commandArgs(trailingOnly = TRUE)
nRep <- if (length(args) >= 1L) as.integer(args[[1]]) else 20L
outFile <- if (length(args) >= 2L) args[[2]] else
  "dev/red-team/heavy-tests/empirical-benefit/results/exact.rds"
dir.create(dirname(outFile), recursive = TRUE, showWarnings = FALSE)

K <- 30L
EMP <- MkPrime::empiricalNObs
ARMS <- c("empirical_geometric", "geometric")
configs <- expand.grid(
  truth = c("empirical", "geometric"),
  pTrue = c(0.3, 0.6, 0.9),
  nTip = 30L,
  nChar = 100L,
  expSteps = c(3, 8),
  stringsAsFactors = FALSE
)

.OneDataset <- function(cfg, seed) {
  set.seed(seed)
  tree <- .SimTree(cfg$nTip, shape = 2, expSteps = cfg$expSteps)
  d <- .SimDataset(tree, cfg$nChar, cfg$truth, cfg$pTrue, EMP, K)
  logLik <- .PerCharLogLik(tree, d$states, d$kObs, K)
  logPV <- .LogPVariable(tree, K)
  rows <- lapply(ARMS, function(arm) {
    post <- .ExactKPosterior(logLik, logPV, arm, EMP, K)
    s <- .ScoreK(post$kPost, d$kTrue)
    data.frame(
      cfg, seed = seed, arm = arm, treeLength = sum(tree$edge.length),
      logScore = mean(s$logScore), cover = mean(s$cover),
      bias = mean(s$err), rmse = sqrt(mean(s$err^2)),
      pMean = sum(post$pGrid * post$pPost),
      meanKTrue = mean(d$kTrue), meanKObs = mean(d$kObs),
      fracKObs2 = mean(d$kObs == 2L), fracKObs3 = mean(d$kObs == 3L)
    )
  })
  # Return:
  do.call(rbind, rows)
}

results <- list()
for (ci in seq_len(nrow(configs))) {
  for (r in seq_len(nRep)) {
    seed <- 179000L + 1000L * ci + r
    results[[length(results) + 1L]] <- .OneDataset(configs[ci, ], seed)
  }
  cat(sprintf("config %d/%d done\n", ci, nrow(configs)))
  saveRDS(do.call(rbind, results), outFile)
}
