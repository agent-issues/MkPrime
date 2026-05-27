#!/usr/bin/env Rscript
# SBC-WARMUP-001 discriminating trace check.
#
# Mirrors `sbc.R` helpers / dimensions / seed math exactly for one Mkp_geometric
# sim, then returns the full posterior samples instead of just SBC ranks. Saves
# the tree_length trace + tl_true to RDS and writes a PNG trace plot.
#
# Usage:
#   Rscript dev/red-team/heavy-tests/sbc-warmup-trace.R           # N_WARM=2000
#   Rscript dev/red-team/heavy-tests/sbc-warmup-trace.R --warm 20000
#
# Output:
#   dev/red-team/heavy-tests/sbc-warmup-trace_warm<N>.rds
#   dev/red-team/heavy-tests/sbc-warmup-trace_warm<N>.png

suppressPackageStartupMessages({
  pkgload::load_all(".", quiet = TRUE)
  library(ape)
  library(TreeTools)
})

args <- commandArgs(trailingOnly = TRUE)
nWarm <- {
  ix <- which(args == "--warm")
  if (length(ix) && length(args) >= ix + 1L) as.integer(args[ix + 1L]) else 2000L
}
simIdx <- {
  ix <- which(args == "--sim")
  if (length(ix) && length(args) >= ix + 1L) as.integer(args[ix + 1L]) else 1L
}
# When true, post-warmup iters held fixed at 6000 (so larger warmup needs
# bigger total nIter). When false, total nIter held at 6000 (v5 production).
holdPostWarm <- any(args == "--hold-post-warm")
# Override start-tree edge length (default 0.1, matching sbc.R:308).
# Use to test whether start-state distance from truth drives the bias.
startEdgeLen <- {
  ix <- which(args == "--start-edge")
  if (length(ix) && length(args) >= ix + 1L) as.numeric(args[ix + 1L]) else 0.1
}
autoTuneOn <- any(args == "--auto-tune")
useTrueTopo <- any(args == "--true-topo")
nIterOverride <- {
  ix <- which(args == "--n-iter")
  if (length(ix) && length(args) >= ix + 1L) as.integer(args[ix + 1L]) else NA_integer_
}
nThinOverride <- {
  ix <- which(args == "--n-thin")
  if (length(ix) && length(args) >= ix + 1L) as.integer(args[ix + 1L]) else NA_integer_
}

# Production dimensions (sbc.R --full)
N_TIP   <- 8L
N_CHAR  <- 30L
N_THIN  <- if (!is.na(nThinOverride)) nThinOverride else 60L
N_ITER  <- if (!is.na(nIterOverride)) nIterOverride else
  if (holdPostWarm) (nWarm + 6000L) else 6000L

EXPSTEPS_FIXED <- 50
TREE_SHAPE     <- 2

# Seed math from sbc.R:438 — arm idx of "Mkp_geometric" in ALL_ARMS is 2.
seedBase <- 20260526L
armIdx   <- 2L
seed     <- seedBase + 1000L * armIdx + simIdx

cat(sprintf("[warmup-trace] N_WARM=%d  sim=%d  seed=%d\n", nWarm, simIdx, seed))

# ---- helpers copied verbatim from sbc.R ----
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
  edges <- tree$edge
  el    <- tree$edge.length
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

# ---- one sim, matched to .runOneSim(arm = Mkp_geometric, ...) ----
set.seed(seed)

p_true <- stats::rbeta(1, 1, 1)
rateLogSd_true <- stats::rgamma(1, shape = 1, rate = 1)
true_tree <- .simTree(N_TIP)
tl_true <- sum(true_tree$edge.length)
cat(sprintf("[warmup-trace] tl_true = %.3f  (prior mean = %d)\n",
            tl_true, EXPSTEPS_FIXED))

# v5 option-α: all chars simulated as binary.
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
kObs <- kObs[keep]
n_char <- sum(keep)
cat(sprintf("[warmup-trace] n_char kept = %d/%d\n", n_char, N_CHAR))

pd <- TreeTools::MatrixToPhyDat(sim_mat)
mkd <- MkPrimeData(pd)

start_tree <- if (useTrueTopo) {
  st <- true_tree
  st$edge.length <- rep_len(startEdgeLen, nrow(st$edge))
  cat("[warmup-trace] using true_tree topology as start\n")
  st
} else {
  tryCatch(
    {
      st <- TreeSearch::AdditionTree(pd)
      st$edge.length <- rep_len(startEdgeLen, nrow(st$edge))
      st
    },
    error = function(e) TreeTools::NJTree(pd, edgeLengths = TRUE)
  )
}
# Report topology distance start_tree vs true_tree (RF, unrooted).
rf_dist <- tryCatch(
  TreeDist::RobinsonFoulds(start_tree, true_tree, normalize = FALSE),
  error = function(e) NA_real_
)
cat(sprintf("[warmup-trace] RF(start, truth) = %s\n", as.character(rf_dist)))
start_tl <- sum(start_tree$edge.length)
cat(sprintf("[warmup-trace] start_tree TL = %.3f  (n_edges = %d)\n",
            start_tl, nrow(start_tree$edge)))

model <- MkPrimeModel(
  coding = "variable",
  nCat = 1L,
  kPrimePrior = "geometric",
  expSteps = EXPSTEPS_FIXED,
  kprimeHyperA = 1, kprimeHyperB = 1,
  kprimeAlpha = 2, kprimeBeta = 2
)

mcmc <- MkPrimeMCMC(
  nIter = N_ITER, thin = N_THIN,
  minWarmup = nWarm, maxWarmup = nWarm,
  autoTune = autoTuneOn,
  nRuns = 1L, nChains = 1L
)

cat(sprintf("[warmup-trace] launching RunMkPrime (nIter=%d thin=%d warm=%d)...\n",
            N_ITER, N_THIN, nWarm))
t0 <- Sys.time()
res <- suppressMessages(suppressWarnings(
  RunMkPrime(mkd, start_tree, model = model, mcmc = mcmc,
             fixTopology = TRUE, overwrite = TRUE)
))
dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("[warmup-trace] done in %.1fs  L_samples = %d\n",
            dt, nrow(res$samples)))

tl_trace <- res$samples[, "tree_length"]
out <- list(
  seed = seed, nWarm = nWarm, simIdx = simIdx,
  tl_true = tl_true, start_tl = start_tl, n_char = n_char,
  tl_trace = tl_trace,
  rate_log_sd_trace = if ("rate_log_sd" %in% colnames(res$samples))
    res$samples[, "rate_log_sd"] else NA_real_,
  wall_sec = dt
)

tag <- sprintf("sim%d_warm%d%s_edge%g%s%s%s",
               simIdx, nWarm,
               if (holdPostWarm) "_hold" else "", startEdgeLen,
               if (autoTuneOn) "_at" else "",
               if (!is.na(nIterOverride)) sprintf("_iter%d", nIterOverride) else "",
               if (useTrueTopo) "_truetopo" else "")
rdsPath <- sprintf("dev/red-team/heavy-tests/sbc-warmup-trace_%s.rds", tag)
saveRDS(out, rdsPath)
cat(sprintf("[warmup-trace] saved %s\n", rdsPath))

# Plot
pngPath <- sprintf("dev/red-team/heavy-tests/sbc-warmup-trace_%s.png", tag)
grDevices::png(pngPath, width = 1000, height = 600)
plot(seq_along(tl_trace), tl_trace, type = "l", col = "steelblue",
     xlab = "kept sample index", ylab = "tree_length",
     main = sprintf("SBC-WARMUP-001 trace: N_WARM=%d, seed=%d, tl_true=%.2f",
                    nWarm, seed, tl_true),
     ylim = range(c(tl_trace, tl_true, start_tl, EXPSTEPS_FIXED)))
graphics::abline(h = tl_true, col = "red", lwd = 2)
graphics::abline(h = start_tl, col = "grey", lty = 3)
graphics::abline(h = EXPSTEPS_FIXED, col = "darkgreen", lty = 2)
graphics::legend("topright",
                 legend = c("tl_trace", "tl_true", "start_tl=0.1*nEdge",
                            "prior mean = 50"),
                 col = c("steelblue", "red", "grey", "darkgreen"),
                 lty = c(1, 1, 3, 2), lwd = c(1, 2, 1, 1))
grDevices::dev.off()
cat(sprintf("[warmup-trace] saved %s\n", pngPath))

# Quick numerical summary
mid <- floor(length(tl_trace) / 2)
firstHalfMean <- mean(tl_trace[1:mid])
secondHalfMean <- mean(tl_trace[(mid + 1):length(tl_trace)])
cat(sprintf("[warmup-trace] mean(first half)=%.2f mean(second half)=%.2f delta=%.2f\n",
            firstHalfMean, secondHalfMean, secondHalfMean - firstHalfMean))
cat(sprintf("[warmup-trace] sample %d = %.2f  sample %d = %.2f  rank(tl_true)=%d/%d\n",
            1, tl_trace[1], length(tl_trace), tl_trace[length(tl_trace)],
            sum(tl_trace < tl_true), length(tl_trace)))
