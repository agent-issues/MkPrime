# Summarise 305_move_shares.R outputs: per-move shares by phase, and
# per-run ESS rates. Usage: Rscript 305_summarise.R <rds>...

PhaseTotals <- function(x, warmup) {
  b <- Filter(function(z) z$warmup == warmup, x$batches)
  if (length(b) == 0L) return(NULL)
  sumOf <- function(f) Reduce(`+`, lapply(b, `[[`, f))
  list(
    iter = sum(vapply(b, `[[`, numeric(1), "nBatch")),
    wallSec = sum(vapply(b, `[[`, numeric(1), "wallSec")),
    propose = sumOf("propose"), accept = sumOf("accept"),
    timeNs = sumOf("timeNs"), weights = b[[length(b)]]$weights
  )
}

MoveTable <- function(x, warmup) {
  p <- PhaseTotals(x, warmup)
  data.frame(
    weight = round(p$weights, 3),
    iterShare = round(p$propose / sum(p$propose), 3),
    timeShare = round(p$timeNs / sum(p$timeNs), 3),
    usPerCall = round(p$timeNs / pmax(p$propose, 1) / 1e3, 1),
    accept = round(p$accept / pmax(p$propose, 1), 2)
  )
}

RunSummary <- function(x) {
  w <- PhaseTotals(x, TRUE)
  s <- PhaseTotals(x, FALSE)
  minEss <- min(x$ess)
  data.frame(
    dataset = x$dataset, arm = x$arm, seed = x$seed,
    mhpWeightSample = round(s$weights[["mh_logit_p"]], 3),
    warmupIter = w$iter, warmupSec = round(w$wallSec, 1),
    sampleIter = s$iter, sampleSec = round(s$wallSec, 1),
    usPerIter = round(s$wallSec / s$iter * 1e6),
    moveNsShare = round(sum(s$timeNs) / 1e9 / s$wallSec, 2),
    minEss = round(minEss), minEssPar = names(which.min(x$ess)),
    essP = round(x$ess[["p"]]), essTL = round(x$ess[["tree_length"]]),
    essKpMin = round(x$essKpMin),
    minEssPerSec = round(minEss / s$wallSec, 2),
    minEssPer1kIter = round(minEss / s$iter * 1e3, 3),
    essPPerSec = round(x$ess[["p"]] / s$wallSec, 2),
    topoFloor = sum(s$weights[c("nni", "spr", "tbr")] <= 0.0105)
  )
}

# The softmax log-score each move earns from the whole run's counters, as
# .AdaptMoveWeights computes it, and the weight share it would get at the
# final temperature (0.5) if gibbs_kPrime were scored rather than pinned
# (floors ignored).
Counterfactual <- function(x, dim) {
  p <- PhaseTotals(x, TRUE)
  q <- PhaseTotals(x, FALSE)
  propose <- p$propose + q$propose
  accept <- p$accept + q$accept
  costS <- pmax((p$timeNs + q$timeNs) / (propose * 1e9), 1e-9)
  logScore <- log(pmax(accept / propose, 1e-12)) + log(pmax(dim, 1)) -
    log(costS)
  share <- exp((logScore - max(logScore)) / 0.5)
  data.frame(logScore = round(logScore, 2),
             shareAtT05 = round(share / sum(share), 3))
}

if (sys.nframe() == 0L) {
  files <- commandArgs(TRUE)
  runs <- lapply(files, readRDS)
  print(do.call(rbind, lapply(runs, RunSummary)), row.names = FALSE)
}
