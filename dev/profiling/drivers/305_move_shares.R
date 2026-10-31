# #305: per-move iteration and time shares, by phase, with the mh_logit_p
# weight adapted (default) or pinned at its initial share.
#
# Usage (from a renamed agent build; see build-agent.sh):
#   Rscript 305_move_shares.R <dataset> <arm> <seed> <warmup> <maxTime> <out> \
#     [pkg] [lib] [pinValue] [coding]
#   dataset   a name in TreeSearch::inapplicable.phyData
#   arm       "adapted" (default scheduler), "capped" (mh_logit_p pinned
#             at pinValue) or "gkp" (gibbs_kPrime pinned at pinValue)
#   warmup    minWarmup = maxWarmup; 0 keeps the default stabilisation
#             detector
#   coding    MkPrimeModel(coding = ); "none" drops the ascertainment term
#             that makes a p move evaluate the likelihood (#305)
#
# .RunBatch is wrapped so each batch's per-move proposal, acceptance and time
# counts are kept with its phase; nothing in the sampler is changed.

args <- commandArgs(TRUE)
dataset <- args[[1]]
arm <- args[[2]]
seed <- as.integer(args[[3]])
warmup <- as.integer(args[[4]])
maxTime <- as.numeric(args[[5]])
out <- args[[6]]
pkg <- if (length(args) >= 7) args[[7]] else "MkPrime"
lib <- if (length(args) >= 8) args[[8]] else NULL
pinValue <- if (length(args) >= 9) as.numeric(args[[9]]) else NA_real_
coding <- if (length(args) >= 10) args[[10]] else "variable"

setTimeLimit(elapsed = maxTime + 120, transient = FALSE)
library(pkg, lib.loc = lib, character.only = TRUE)
ns <- asNamespace(pkg)
options(MkPrime.verbosity = 0)

logged <- new.env()
logged$batches <- list()
original <- get(".RunBatch", ns)
hook <- function(...) {
  a <- list(...)
  t0 <- proc.time()[["elapsed"]]
  res <- original(...)
  logged$batches[[length(logged$batches) + 1L]] <- list(
    weights = a[[7]], nBatch = a[[14]], batchStart = a[[15]],
    warmup = a[[16]] > 0L,
    propose = res$propose_counts[1, ], accept = res$accept_counts[1, ],
    timeNs = res$move_time_ns[1, ],
    wallSec = proc.time()[["elapsed"]] - t0
  )
  res
}
unlockBinding(".RunBatch", ns)
assign(".RunBatch", hook, envir = ns)

dat <- TreeSearch::inapplicable.phyData[[dataset]]
mkd <- suppressWarnings(MkPrimeData(dat))
set.seed(seed)

warmupArgs <- if (warmup > 0L) list(minWarmup = warmup, maxWarmup = warmup)
pinArgs <- switch(arm,
  capped = list(moveWeights = c(mh_logit_p = pinValue)),
  gkp = list(moveWeights = c(gibbs_kPrime = pinValue))
)
mcmc <- do.call(MkPrimeMCMC, c(
  list(nRuns = 1L, autoTune = FALSE, maxTime = maxTime, thin = 10L),
  warmupArgs, pinArgs
))
model <- MkPrimeModel(coding = coding)

wall <- system.time(
  post <- suppressWarnings(RunMkPrime(mkd, model = model, mcmc = mcmc))
)[["elapsed"]]

batches <- logged$batches
moveNames <- names(post$moveWeights)
for (i in seq_along(batches)) {
  for (f in c("weights", "propose", "accept", "timeNs")) {
    names(batches[[i]][[f]]) <- moveNames
  }
}

samples <- post$samples
gate <- samples[, ns$.GateCols(colnames(samples)), drop = FALSE]
kp <- samples[, grepl("^kPrime_", colnames(samples)), drop = FALSE]
ess <- ns$.EssMatrix(gate)
names(ess) <- colnames(gate)
essKp <- ns$.EssMatrix(kp)

saveRDS(list(
  dataset = dataset, arm = arm, seed = seed, warmup = warmup, coding = coding,
  pinValue = pinValue, nTips = length(dat), nTrans = sum(mkd$type ==
    "transformational"), wall = wall, batches = batches,
  finalWeights = post$moveWeights, acceptance = post$acceptance,
  stopReason = post$per_run[[1]]$stop_reason %||% NA,
  nSamples = nrow(samples),
  ess = ess, essKpMin = min(essKp), essKpMedian = stats::median(essKp),
  warmupTrace = post$warmup_trace
), out)
cat("done", dataset, arm, seed, "wall", round(wall, 1), "\n")
