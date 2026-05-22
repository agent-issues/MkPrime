# warmstart_from_blind.R -------------------------------------------------------
# Short (5 000 iter) ecology-AWARE MkNT run initialised from the blind MR
# consensus tree.  Tests whether the blind topology is accessible under the
# aware model:
#
#   If aware drifts away  -> aware genuinely prefers a different region;
#                            v5 mode-traps are within the aware posterior.
#   If aware stays near   -> v5 runs were init-artifact traps; v6 needs warm
#                            starts and/or more chains.
#
# Run from repo root:
#   Rscript --vanilla inst/ecology/scripts/rodent-comparison/warmstart_from_blind.R

suppressPackageStartupMessages({
  library(ape)
  library(TreeTools)
  library(MkPrime)
})

nexFile <- "inst/ecology/data/rodent-X24848.nex"
outDir  <- "inst/ecology/scripts/rodent-comparison"

# ---------------------------------------------------------------------------
# 1. Build MkNT data (mirrors run_rodent_MkNT.R).
# ---------------------------------------------------------------------------
mat <- TreeTools::ReadCharacters(nexFile)
cat("Raw matrix:", nrow(mat), "tips x", ncol(mat), "chars\n")

ecologyCol <- 220L; extantCol <- 221L
ecoVec <- mat[, ecologyCol]; extant <- mat[, extantCol]
keepTaxa <- which(extant == "0")
matKeep <- mat[keepTaxa, -extantCol, drop = FALSE]
ecoKeep <- ecoVec[keepTaxa]

poly <- grepl("^[(]", ecoKeep)
if (any(poly)) {
  ecoKeep[poly] <- substr(sub("^[(]", "", ecoKeep[poly]), 1, 1)
  matKeep[poly, ecologyCol] <- ecoKeep[poly]
}
hasEco <- !is.na(ecoKeep) & ecoKeep != "?" &
            !is.na(suppressWarnings(as.integer(ecoKeep)))
matKeep <- matKeep[hasEco, , drop = FALSE]
ecoKeep <- ecoKeep[hasEco]

charMat     <- matKeep[, -ecologyCol, drop = FALSE]
pdForDetect <- MatrixToPhyDat(charMat)
neoIdx      <- AutoDetectNeomorphic(pdForDetect)
cat("Neomorphic:", length(neoIdx), "  non-neomorphic:",
    ncol(charMat) - length(neoIdx), "\n")

.kObsCol <- function(col) {
  vals <- col[!(col %in% c("?", "-", NA))]
  vals <- unlist(strsplit(gsub("[()]", "", vals), ""))
  length(unique(vals))
}
kObsRaw        <- vapply(seq_len(ncol(charMat)),
                         function(j) .kObsCol(charMat[, j]), integer(1L))
nonNeoOriginal <- setdiff(seq_len(ncol(charMat)), neoIdx)
nonNeoVariable <- nonNeoOriginal[kObsRaw[nonNeoOriginal] >= 2L]
knownK         <- setNames(as.integer(kObsRaw[nonNeoVariable]),
                           as.character(nonNeoVariable))

mkd <- MkPrimeData(pdForDetect,
                   neomorphic  = neoIdx,
                   knownStates = knownK,
                   ecology     = setNames(as.integer(ecoKeep),
                                         rownames(matKeep)))
stopifnot(all(mkd$type %in% c("neomorphic", "known")))
cat("MkNT: nTip=", mkd$nTip, " nChar=", mkd$nChar,
    " kEcology=", mkd$kEcology, "\n")

# ---------------------------------------------------------------------------
# 2. Build blind consensus tree as start tree.
# ---------------------------------------------------------------------------
blind_files <- Sys.glob(file.path(outDir, "rodent-MkNT-blind_trees_*.nwk"))
blind_trees <- do.call(c, lapply(blind_files, function(f) {
  tr <- ape::read.tree(f)
  if (inherits(tr, "phylo")) list(tr) else as.list(tr)
}))
class(blind_trees) <- "multiPhylo"
n_bl   <- length(blind_trees)
bl_post <- blind_trees[seq.int(ceiling(n_bl * 0.25) + 1L, n_bl)]
cons_blind <- ape::consensus(bl_post, p = 0.5, rooted = FALSE)
cat("Blind consensus: ", length(cons_blind$tip.label), "tips\n")

# Validate tip match
blind_tips <- sort(cons_blind$tip.label)
mkd_tips   <- sort(mkd$taxon_names)
if (!identical(blind_tips, mkd_tips)) {
  only_bl <- setdiff(blind_tips, mkd_tips)
  only_mk <- setdiff(mkd_tips, blind_tips)
  stop("Tip mismatch.  Only in blind: ", paste(only_bl, collapse = ", "),
       ".  Only in mkd: ",  paste(only_mk, collapse = ", "))
}
cat("Tip sets match (", length(blind_tips), ").\n")

# Resolve polytomies (MR consensus has some; C++ kernel requires binary)
cat("Nodes before resolve:", ape::Nnode(cons_blind),
    " (binary would be", ape::Ntip(cons_blind) - 2L, ")\n")
cons_blind <- ape::multi2di(cons_blind, random = FALSE)
cat("Nodes after  resolve:", ape::Nnode(cons_blind), "\n")

# Ensure branch lengths present
if (is.null(cons_blind$edge.length))
  cons_blind$edge.length <- rep(0.1, nrow(cons_blind$edge))
startTree <- ape::unroot(TreeTools::Preorder(cons_blind))

# ---------------------------------------------------------------------------
# 3. Aware MkNT model (same spec as v5).
# ---------------------------------------------------------------------------
model <- MkPrimeModel(
  ecologyAware  = TRUE,
  magnitudeMode = "global",
  kPrimePrior   = "geometric",
  coding        = "variable",
  rho0Alpha     = 7, rho0Beta  = 3,
  thetaAlpha    = 2, thetaBeta = 2,
  sigmaPhi      = 0.5
)

# ---------------------------------------------------------------------------
# 4. Short MCMC — serial, 5 000 iter, PT with 4 chains.
#    heat = 0.35 avoids the wide-ladder collapse before adaptation can act.
# ---------------------------------------------------------------------------
nIter   <- 5000L
logFile <- file.path(outDir, "warmstart_blind.log")
ckpFile <- file.path(outDir, "warmstart_blind.ckp")
if (file.exists(logFile)) file.remove(logFile)
if (file.exists(ckpFile)) file.remove(ckpFile)

mcmc <- MkPrimeMCMC(
  nIter          = nIter,
  nChains        = 4L,
  nRuns          = 1L,
  nCore          = 1L,
  heat           = 0.35,
  thin           = 5L,
  treeThin       = 5L,
  minWarmup      = 1000L,
  maxWarmup      = 3000L,
  logFile        = logFile,
  checkpointFile = ckpFile
)
Sys.setenv(MKPRIME_ECO_RESYNC_EVERY = "200")

cat("\n=== Warm-start run (init = blind consensus) ===\n")
cat("nIter=", nIter, "  nChains=4  heat=0.35\n")
t0  <- Sys.time()
res <- RunMkPrime(mkd, tree = startTree, model = model, mcmc = mcmc)
cat("Elapsed:", format(Sys.time() - t0), "\n")

# ---------------------------------------------------------------------------
# 5. Diagnostics.
# ---------------------------------------------------------------------------
cat("\n=== PT Diagnostics ===\n")
cat("betas:            ", paste(round(res$betas, 4), collapse = " "), "\n")
cat("round_trip_count: ", res$round_trip_count, "\n")
cat("swap_rates:       ", paste(round(res$swap_rates, 4), collapse = " "), "\n")
cat("stop_reason:      ", res$stop_reason, "\n")
cat("actual_iter:      ", res$actual_iter, "\n")

s <- ReadMkLog(logFile)
if (is.matrix(s)) s <- data.frame(s, check.names = FALSE)
ll_col <- if ("log_likelihood" %in% names(s)) "log_likelihood" else "logLik"
tl_col <- if ("tree_length"    %in% names(s)) "tree_length"    else "treeLength"
cat("Samples post-MCMC:", nrow(s), "\n")
cat("logLik:    "); print(summary(s[[ll_col]]))
cat("treeLength:"); print(summary(s[[tl_col]]))

# Compare to blind consensus
treefile <- sub("\\.log$", "_trees.nwk", logFile)
if (!file.exists(treefile)) {
  treefile2 <- paste0(logFile, "_trees.nwk")
  if (file.exists(treefile2)) treefile <- treefile2
}
if (file.exists(treefile)) {
  ws_trees <- ape::read.tree(treefile)
  if (inherits(ws_trees, "phylo")) ws_trees <- list(ws_trees)
  class(ws_trees) <- "multiPhylo"
  n_ws  <- length(ws_trees)
  ws_post  <- ws_trees[seq.int(ceiling(n_ws * 0.25) + 1L, n_ws)]
  cons_ws  <- ape::consensus(ws_post, p = 0.5, rooted = FALSE)
  rf_vs_blind <- ape::dist.topo(cons_ws, cons_blind)
  max_rf      <- 2 * (length(cons_ws$tip.label) - 3)
  cat(sprintf("\nRF (warm-start MR consensus vs blind MR consensus): %d / %d (%.0f%%)\n",
              rf_vs_blind, max_rf, 100 * rf_vs_blind / max_rf))
  cat(sprintf("Trees in warm-start posterior: %d  (total: %d)\n",
              length(ws_post), n_ws))
} else {
  cat("Tree file not found:", treefile, "\n")
}
