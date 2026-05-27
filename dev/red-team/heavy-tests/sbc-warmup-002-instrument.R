#!/usr/bin/env Rscript
# SBC-WARMUP-002 instrumentation: log .AdaptMoveWeights inputs/outputs
# for the tree_length move during a warmup-20000 / hold-post-warm run.
#
# Writes a one-line-per-call trace to
#   dev/red-team/heavy-tests/sbc-warmup-002-trace.txt
#
# Reproduces SBC-WARMUP-002 (tree_length frozen at floor weight).

suppressPackageStartupMessages({
  pkgload::load_all(".", quiet = TRUE)
  library(ape)
  library(TreeTools)
})

logPath <- "dev/red-team/heavy-tests/sbc-warmup-002-trace.txt"
if (file.exists(logPath)) file.remove(logPath)
cat(sprintf("# .AdaptMoveWeights trace for SBC-WARMUP-002 (autoTune=FALSE, warm=20000)\n"),
    file = logPath, append = FALSE)
cat(sprintf("# columns: call,tl_in,tl_out,tl_accept,tl_prop,tl_cost_s,wMin,warmupProgress\n"),
    file = logPath, append = TRUE)

callIdx <- 0L
trace(
  MkPrime:::.AdaptMoveWeights,
  exit = quote({
    callIdx <<- callIdx + 1L
    targets <- c("tree_length", "joint_tl_rls", "branch_lengths", "rate_log_sd",
                 "kPrime", "nni", "spr")
    out <- returnValue()
    fields <- c(sprintf("call=%d wMin=%.3f prog=%.3f", callIdx, wMin, warmupProgress))
    for (nm in targets) {
      i <- match(nm, moveNames)
      if (is.na(i)) next
      a <- acceptCount[[i]]; p <- proposeCount[[i]]
      cs <- (moveTimeNs[[i]] / max(p, 1L)) / 1e9
      rate <- if (p > 0L) a / p else NA_real_
      fields <- c(fields, sprintf("%s w=%.4f acc=%.2f cost_us=%.1f",
                                    nm, out[[i]], rate, cs * 1e6))
    }
    cat(paste(fields, collapse = " | "), "\n",
        file = logPath, append = TRUE, sep = "")
  }),
  print = FALSE
)

# --- replicate sbc-warmup-trace.R sim 1 with --warm 20000 --hold-post-warm ---
N_TIP <- 8L; N_CHAR <- 30L; N_THIN <- 60L; N_WARM <- 20000L
N_ITER <- N_WARM + 6000L
EXPSTEPS_FIXED <- 50; TREE_SHAPE <- 2
seed <- 20260526L + 1000L * 2L + 1L

.simTree <- function(nTip) {
  tr <- ape::rtree(nTip, tip.label = paste0("t", seq_len(nTip)))
  nEdge <- nrow(tr$edge)
  tl <- stats::rgamma(1, shape = TREE_SHAPE, rate = TREE_SHAPE / EXPSTEPS_FIXED)
  w <- stats::rexp(nEdge, rate = 1)
  tr$edge.length <- tl * w / sum(w)
  TreeTools::Preorder(tr)
}
.simJCchar <- function(tree, kTrue) {
  nTip <- length(tree$tip.label)
  states <- integer(2L * nTip - 1L)
  rootIdx <- nTip + 1L
  states[rootIdx] <- sample.int(kTrue, 1L) - 1L
  edges <- tree$edge; el <- tree$edge.length
  for (e in rev(seq_len(nrow(edges)))) {
    pa <- edges[e, 1L]; ch <- edges[e, 2L]; t <- el[e]
    pSame <- 1 / kTrue + (1 - 1 / kTrue) * exp(-kTrue * t / (kTrue - 1))
    if (runif(1L) < pSame) {
      states[ch] <- states[pa]
    } else {
      states[ch] <- sample(setdiff(seq.int(0L, kTrue - 1L), states[pa]), 1L)
    }
  }
  states[seq_len(nTip)]
}
.canonicaliseLabels <- function(vec) {
  uvals <- sort(unique(vec))
  out <- match(vec, uvals) - 1L
  attr(out, "kObs") <- length(uvals)
  out
}

set.seed(seed)
p_true <- stats::rbeta(1, 1, 1)
rateLogSd_true <- stats::rgamma(1, shape = 1, rate = 1)
true_tree <- .simTree(N_TIP)

kSim <- rep(2L, N_CHAR)
sim_mat <- matrix(NA_integer_, N_TIP, N_CHAR,
                  dimnames = list(true_tree$tip.label, NULL))
kObs <- integer(N_CHAR)
for (j in seq_len(N_CHAR)) {
  raw <- .simJCchar(true_tree, kSim[j])
  canon <- .canonicaliseLabels(raw)
  sim_mat[, j] <- canon
  kObs[j] <- attr(canon, "kObs")
}
keep <- kObs >= 2L
sim_mat <- sim_mat[, keep, drop = FALSE]
pd <- TreeTools::MatrixToPhyDat(sim_mat)
mkd <- MkPrimeData(pd)

start_tree <- tryCatch(
  {
    st <- TreeSearch::AdditionTree(pd)
    st$edge.length <- rep_len(0.1, nrow(st$edge))
    st
  },
  error = function(e) TreeTools::NJTree(pd, edgeLengths = TRUE)
)
model <- MkPrimeModel(coding = "variable", nCat = 1L,
                       kPrimePrior = "geometric", expSteps = EXPSTEPS_FIXED,
                       kprimeHyperA = 1, kprimeHyperB = 1,
                       kprimeAlpha = 2, kprimeBeta = 2)
mcmc <- MkPrimeMCMC(nIter = N_ITER, thin = N_THIN,
                     minWarmup = N_WARM, maxWarmup = N_WARM,
                     autoTune = FALSE, nRuns = 1L, nChains = 1L)

cat(sprintf("[instrument] running RunMkPrime (nIter=%d warm=%d) ...\n",
            N_ITER, N_WARM))
res <- suppressMessages(suppressWarnings(
  RunMkPrime(mkd, start_tree, model = model, mcmc = mcmc,
              fixTopology = TRUE, overwrite = TRUE)
))
cat(sprintf("[instrument] done; %d calls logged to %s\n", callIdx, logPath))
cat(sprintf("[instrument] tl_trace unique=%d/%d (frozen if 1)\n",
            length(unique(res$samples[, "tree_length"])),
            nrow(res$samples)))
