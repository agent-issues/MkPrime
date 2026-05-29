# eval_aware_at_reference.R
#
# Evaluate the aware logLikelihood at two "reference" topologies:
#   (A) The blind MR consensus (biologically reasonable, from converged
#       blind MCMC with 4 runs x 100k iter).
#   (B) The v6 run-3 best tree (highest logLik found so far, -5290.9,
#       but biologically nonsensical).
#
# Method: RunMkPrime(..., fixTopology = TRUE) with nChains=1 so no PT,
# only branch-length and rate moves.  500 iterations to get logLik
# distribution at each fixed topology.
#
# Interpretation:
#   logLik(A) ≈ logLik(B)  → mode-trap; aware model is fine, MCMC just
#                             cannot find the sensible-topology basin.
#   logLik(A) << logLik(B)  → aware model genuinely prefers nonsensical
#                             topologies; ecological prior or likelihood
#                             is mis-specified.

libPath <- Sys.getenv("MKP_LIB",
  file.path("/nobackup", Sys.getenv("USER"), "mkp-sim3-multirep-v3/lib"))
.libPaths(c(libPath, .libPaths()))
suppressPackageStartupMessages({
  library(ape); library(TreeTools); library(MkPrime)
})

nexFile  <- Sys.getenv("NEX_FILE",
  "/nobackup/pjjg18/mkp-rodent-blind-v2/data/mbank_X24848_2026-5-9-1135.nex")
blindDir <- Sys.getenv("BLIND_DIR",
  "/nobackup/pjjg18/mkp-rodent-MkNT-v1/blind")
v6Dir    <- Sys.getenv("V6_DIR",
  "/nobackup/pjjg18/mkp-rodent-MkNT-v6/aware")
outDir   <- Sys.getenv("OUT_DIR",
  "/nobackup/pjjg18/mkp-rodent-eval-ref")
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

cat("=== eval_aware_at_reference.R ===\n")
cat("nexFile: ", nexFile, "\n")
cat("blindDir:", blindDir, "\n")
cat("v6Dir:   ", v6Dir,    "\n")
cat("outDir:  ", outDir,   "\n\n")

# ---------------------------------------------------------------------------
# Build MkNT data (identical to run_rodent_MkNT.R)
# ---------------------------------------------------------------------------
mat       <- TreeTools::ReadCharacters(nexFile)
ecologyCol <- 220L; extantCol <- 221L
ecoVec    <- mat[, ecologyCol]
keepTaxa  <- which(mat[, extantCol] == "0")
matKeep   <- mat[keepTaxa, -extantCol, drop = FALSE]
ecoKeep   <- ecoVec[keepTaxa]
poly      <- grepl("^[(]", ecoKeep)
if (any(poly)) ecoKeep[poly] <- substr(sub("^[(]", "", ecoKeep[poly]), 1, 1)
hasEco    <- !is.na(ecoKeep) & ecoKeep != "?" &
             !is.na(suppressWarnings(as.integer(ecoKeep)))
matKeep   <- matKeep[hasEco, , drop = FALSE]
ecoKeep   <- ecoKeep[hasEco]
charMat   <- matKeep[, -ecologyCol, drop = FALSE]

pdForDetect <- MatrixToPhyDat(charMat)
neoIdx      <- AutoDetectNeomorphic(pdForDetect)
kObsRaw     <- vapply(seq_len(ncol(charMat)), function(j) {
  vals <- charMat[, j]; vals <- vals[!(vals %in% c("?", "-", NA))]
  vals <- unlist(strsplit(gsub("[()]", "", vals), ""))
  length(unique(vals))
}, integer(1L))
nonNeoVar <- setdiff(seq_len(ncol(charMat)), neoIdx)[
               kObsRaw[setdiff(seq_len(ncol(charMat)), neoIdx)] >= 2L]
knownK    <- setNames(as.integer(kObsRaw[nonNeoVar]), as.character(nonNeoVar))
mkd       <- suppressWarnings(MkPrimeData(
  pdForDetect, neomorphic = neoIdx, knownStates = knownK,
  ecology = setNames(as.integer(ecoKeep), rownames(matKeep))))
cat("MkNT: nTip=", mkd$nTip, " nChar=", mkd$nChar,
    " kEcology=", mkd$kEcology, "\n")

# ---------------------------------------------------------------------------
# Model (aware, same spec as v6)
# ---------------------------------------------------------------------------
modelMkNT <- MkPrimeModel(
  ecologyAware  = TRUE,  magnitudeMode = "global",
  rho0Alpha     = 7,     rho0Beta  = 3,
  thetaAlpha    = 2,     thetaBeta = 2,
  sigmaPhi      = 1.5,
  rateLossMeanlog = 0,   rateLossSdlog  = 2,
  rateLogSdShape  = 1,   rateLogSdRate  = 1,
  rateNeoMeanlog  = 0,   rateNeoSdlog   = 1
)

# ---------------------------------------------------------------------------
# Helper: run 1000-iter fixed-topology MCMC, return post-burnin logLik
# ---------------------------------------------------------------------------
eval_at_tree <- function(tr, label) {
  cat(sprintf("\n--- Evaluating at: %s ---\n", label))
  logf <- file.path(outDir, paste0("eval_", label, ".log"))
  ckpf <- file.path(outDir, paste0("eval_", label, ".ckp"))
  mcmc <- MkPrimeMCMC(
    nIter = 1000L, nChains = 1L, nRuns = 1L, nCore = 1L,
    thin  = 1L,    treeThin = 1L,
    minWarmup = 200L, maxWarmup = 300L,
    logFile = logf, checkpointFile = ckpf
  )
  res <- RunMkPrime(mkd, tree = tr, model = modelMkNT, mcmc = mcmc,
                    fixTopology = TRUE)
  s   <- ReadMkLog(logf)
  if (is.matrix(s)) s <- data.frame(s, check.names = FALSE)
  ll_col <- if ("log_likelihood" %in% names(s)) "log_likelihood" else "logLik"
  n   <- nrow(s)
  pb  <- s[[ll_col]][seq.int(ceiling(n * 0.25) + 1L, n)]
  cat(sprintf("  n_postburnin=%d  median=%.2f  max=%.2f  min=%.2f\n",
              length(pb), median(pb), max(pb), min(pb)))
  invisible(pb)
}

# ---------------------------------------------------------------------------
# Reference tree A: blind MR consensus (post-burnin, p=0.5)
# ---------------------------------------------------------------------------
cat("\nBuilding blind consensus...\n")
bl_files <- Sys.glob(file.path(blindDir, "rodent-MkNT-blind_trees_*.nwk"))
bl_trees <- do.call(c, lapply(bl_files, function(f) {
  tr <- read.tree(f)
  if (inherits(tr, "phylo")) list(tr) else as.list(tr)
}))
class(bl_trees) <- "multiPhylo"
n_bl  <- length(bl_trees)
bl_pb <- bl_trees[seq.int(ceiling(n_bl * 0.25) + 1L, n_bl)]
class(bl_pb) <- "multiPhylo"
cons_blind <- multi2di(consensus(bl_pb, p = 0.5, rooted = FALSE), random = FALSE)
cons_blind$edge.length <- rep(0.1, nrow(cons_blind$edge))
ll_blind <- eval_at_tree(cons_blind, "blind_consensus")

# ---------------------------------------------------------------------------
# Reference tree B: v6 run-3 best tree (highest logLik sample)
# ---------------------------------------------------------------------------
cat("\nLoading v6 run-3 best tree...\n")
s3     <- ReadMkLog(file.path(v6Dir, "rodent-MkNT-aware_3.log"))
if (is.matrix(s3)) s3 <- data.frame(s3, check.names = FALSE)
best_i <- which.max(s3$log_likelihood)
tr3_all <- read.tree(file.path(v6Dir, "rodent-MkNT-aware_trees_3.nwk"))
if (inherits(tr3_all, "phylo")) {
  tr3_best <- tr3_all
} else {
  tr3_best <- tr3_all[[min(best_i, length(tr3_all))]]
}
tr3_best$edge.length <- rep(0.1, nrow(tr3_best$edge))
ll_v6r3 <- eval_at_tree(tr3_best, "v6run3_best")

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
cat("\n=======================================================\n")
cat(sprintf("RESULT:  aware@blind_consensus  median logLik = %.2f\n",
            median(ll_blind)))
cat(sprintf("         aware@v6run3_best      median logLik = %.2f\n",
            median(ll_v6r3)))
cat(sprintf("         difference (blind - run3) = %.2f log-units\n",
            median(ll_blind) - median(ll_v6r3)))
cat("\nINTERPRETATION:\n")
diff_val <- median(ll_blind) - median(ll_v6r3)
if (abs(diff_val) < 10) {
  cat("  |diff| < 10: topologies roughly equivalent — MODE-TRAP.\n")
  cat("  Blind topology is accessible; MCMC is stuck, not misdirected.\n")
} else if (diff_val < -10) {
  cat("  blind_consensus MUCH WORSE (by", round(abs(diff_val), 1), "units).\n")
  cat("  Aware model genuinely prefers run3-type topology over realistic tree.\n")
  cat("  Possible mis-specification in ecology prior or likelihood.\n")
} else {
  cat("  blind_consensus BETTER (by", round(diff_val, 1), "units).\n")
  cat("  Sensible topology outscores run3 — run3 is itself a mode-trap.\n")
  cat("  MCMC hasn't yet found the true posterior peak.\n")
}
cat("=======================================================\n")
