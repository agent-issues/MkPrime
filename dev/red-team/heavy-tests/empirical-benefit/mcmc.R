#!/usr/bin/env Rscript
# MCMC arm of the empirical-vs-geometric k' prior study (#179).
#
# One task = one simulated dataset, fitted with free topology under both k'
# priors. Scores k' (score.R: log score, 90% HPD coverage, bias, RMSE) from
# the draws, the trees by mean normalized CID and RF to the true tree, and
# also the exact known-tree posterior as the reference the draws should
# approach. The paired difference between arms on one dataset is the unit
# of analysis.
#
# Usage, from the repository root:
#   Rscript dev/red-team/heavy-tests/empirical-benefit/mcmc.R <task> [--pilot]
#
# <task> is 1-based and indexes CONFIGS x replicate (see .Task()). Writes
# results/mcmc-<task>.rds; a task whose file exists is skipped, so an array
# can be resubmitted after a walltime kill.

stopifnot(requireNamespace("TreeDist", quietly = TRUE))
suppressPackageStartupMessages({
  source("dev/red-team/heavy-tests/load-mkprime.R")
  LoadMkPrime(".")
})
source("dev/red-team/heavy-tests/empirical-benefit/forward.R")
source("dev/red-team/heavy-tests/empirical-benefit/score.R")

args <- commandArgs(trailingOnly = TRUE)
task <- as.integer(args[[1]])
pilot <- "--pilot" %in% args

K <- 30L
EMP <- MkPrime::empiricalNObs
ARMS <- c("empirical_geometric", "geometric")
nTip <- 30L
nChar <- 100L
expSteps <- 3
# Row 1 is the headline: its recorded kObs match the corpus (two thirds
# binary). Row 2 has the most unobserved states; rows 3-4 are the geometric
# controls, row 4 the mostly-binary truth on which the empirical floor loses.
CONFIGS <- data.frame(
  truth = c("empirical", "empirical", "geometric", "geometric"),
  pTrue = c(0.6, 0.3, 0.6, 0.9),
  stringsAsFactors = FALSE
)
nRep <- 50L
MCMC <- if (pilot) {
  list(nIter = 20000L, thin = 10L, warmup = 5000L, maxTime = 300)
} else {
  list(nIter = 100000L, thin = 20L, warmup = 20000L, maxTime = 3000)
}
outDir <- Sys.getenv(
  "MKP_EB_OUT", "dev/red-team/heavy-tests/empirical-benefit/results"
)
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
outFile <- file.path(outDir, sprintf("mcmc-%04d.rds", task))
if (file.exists(outFile)) {
  cat("Exists, skipping:", outFile, "\n")
  quit(status = 0L)
}

.Task <- function(task) {
  ci <- (task - 1L) %/% nRep + 1L
  if (ci > nrow(CONFIGS)) stop("Task ", task, " is beyond the design")
  # Return:
  list(cfg = CONFIGS[ci, ], rep = (task - 1L) %% nRep + 1L)
}

.Fit <- function(mkd, arm) {
  model <- MkPrimeModel(
    coding = "variable", nCat = 1L, kPrimePrior = arm,
    kprimeTruncK = K, empiricalNObs = if (arm == ARMS[[1]]) EMP,
    expSteps = expSteps
  )
  mcmc <- MkPrimeMCMC(
    nIter = MCMC$nIter, thin = MCMC$thin, minWarmup = MCMC$warmup %/% 2L,
    maxWarmup = MCMC$warmup, autoTune = !pilot, nRuns = 2L,
    maxTime = MCMC$maxTime
  )
  t0 <- Sys.time()
  res <- RunMkPrime(mkd, model = model, mcmc = mcmc, overwrite = TRUE,
                    verbosity = 0L)
  res$wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  # Return:
  res
}

tk <- .Task(task)
seed <- 179000L + task
set.seed(seed)
tree <- .SimTree(nTip, shape = 2, expSteps = expSteps)
d <- .SimDataset(tree, nChar, tk$cfg$truth, tk$cfg$pTrue, EMP, K)
mkd <- MkPrimeData(TreeTools::MatrixToPhyDat(d$states))
# Column order matters (#54): compare each character's tip partition.
Partition <- function(m) apply(m, 2, function(x) match(x, unique(x)))
stopifnot(
  identical(Partition(mkd$matrix[rownames(d$states), ]), Partition(d$states)),
  all(mkd$type == "transformational")
)
logLik <- .PerCharLogLik(tree, d$states, d$kObs, K)
logPV <- .LogPVariable(tree, K)

rows <- lapply(ARMS, function(arm) {
  res <- .Fit(mkd, arm)
  kCols <- paste0("kPrime_", seq_len(nChar))
  sampled <- .ScoreK(.SampledKPosterior(res$samples[, kCols], K), d$kTrue,
                     floor = 0.5 / nrow(res$samples))
  exact <- .ScoreK(.ExactKPosterior(logLik, logPV, arm, EMP, K)$kPost,
                   d$kTrue)
  trees <- .ScoreTrees(res$trees, tree)
  data.frame(
    tk$cfg, rep = tk$rep, task = task, seed = seed, arm = arm,
    nTip = nTip, nChar = nChar, treeLength = sum(tree$edge.length),
    logScore = mean(sampled$logScore), cover = mean(sampled$cover),
    bias = mean(sampled$err), rmse = sqrt(mean(sampled$err^2)),
    exactLogScore = mean(exact$logScore), exactBias = mean(exact$err),
    cid = trees[["cid"]], rf = trees[["rf"]],
    nDraw = nrow(res$samples), wall = res$wall,
    meanKTrue = mean(d$kTrue), meanKObs = mean(d$kObs)
  )
})
out <- do.call(rbind, rows)
print(out)
saveRDS(list(rows = out, kTrue = d$kTrue, kObs = d$kObs, tree = tree),
        outFile)
