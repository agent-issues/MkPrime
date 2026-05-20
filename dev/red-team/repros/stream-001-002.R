#!/usr/bin/env Rscript
# Reproducer for STREAM-001 (mk-arm trees not written) and STREAM-002
# (checkpoint never advances past iter=0).
#
# Run from the repo root:
#   Rscript dev/red-team/repros/stream-001-002.R [mk|mkp_eg]
#
# Expected (before fix): mk arm produces an empty `mk_trees.nwk` (1 line) and
# a `mk_checkpoint.rds` with iter=0 while the log has many Sample-phase rows;
# `mkp_eg` writes trees correctly.
#
# Expected (after fix): both arms write trees and the checkpoint iter advances.

suppressPackageStartupMessages({
  devtools::load_all(quiet = TRUE)
  library(ape)
  library(TreeTools)
})

arm <- commandArgs(trailingOnly = TRUE)[1L]
if (is.na(arm)) arm <- "mk"
stopifnot(arm %in% c("mk", "mkp_eg"))

out_dir <- tempfile("stream_repro_")
dir.create(out_dir, recursive = TRUE)
cat("out_dir:", out_dir, "\n")

# Tiny dataset: 6 taxa, 10 characters, 3 states each.
set.seed(42)
nTip <- 6L
nChar <- 10L
mat <- matrix(sample(c("0", "1", "2"), nTip * nChar, replace = TRUE),
              nrow = nTip, ncol = nChar)
rownames(mat) <- paste0("t", seq_len(nTip))
tmp_nex <- tempfile(fileext = ".nex")
ape::write.nexus.data(
  setNames(lapply(seq_len(nTip), function(i) mat[i, ]), rownames(mat)),
  file = tmp_nex, format = "standard")
pd <- TreeTools::ReadAsPhyDat(tmp_nex)
file.remove(tmp_nex)

start_tree <- TreeTools::NJTree(pd, edgeLengths = TRUE)

make_mcmc <- function(prefix) {
  MkPrimeMCMC(
    nIter      = Inf,
    thin       = 10L,
    maxWarmup  = 100L,
    minWarmup  = 50L,
    nRuns      = 2L,
    nChains    = 2L,
    heat       = 0.1,
    maxTime    = 12,
    minEss     = 1e6,  # never converge — force time-limit exit
    maxRhat    = 1.0001,
    checkEvery = 50L,
    tuningRounds = 1L,
    tuningBudget = 200L,
    checkpointFile = file.path(out_dir, paste0(prefix, "_checkpoint.rds")),
    logFile        = file.path(out_dir, paste0(prefix, "_run.log")),
    treeFile       = file.path(out_dir, paste0(prefix, "_trees.nwk"))
  )
}

if (arm == "mk") {
  kobs_raw <- apply(mat, 2L, function(col) {
    length(unique(col[!col %in% c("?", "-")]))
  })
  var_orig <- which(kobs_raw > 1L)
  kObs_for_mk <- setNames(as.integer(kobs_raw[var_orig]),
                          as.character(var_orig))
  res <- RunMkPrime(
    pd, start_tree,
    knownStates = kObs_for_mk,
    model = MkPrimeModel(coding = "variable"),
    mcmc  = make_mcmc("mk")
  )
} else {
  mkd <- MkPrimeData(pd)
  res <- RunMkPrime(
    mkd, start_tree,
    model = MkPrimeModel(coding = "variable",
                         kPrimePrior = "empirical_geometric"),
    mcmc  = make_mcmc("mkp_eg")
  )
}

# Inventory
tree_file <- file.path(out_dir, paste0(arm, "_trees.nwk"))
log_file_pat <- list.files(out_dir, pattern = paste0("^", arm, "_run.*log$"),
                           full.names = TRUE)
ckp_file <- file.path(out_dir, paste0(arm, "_checkpoint.rds"))

tree_lines <- if (file.exists(tree_file)) length(readLines(tree_file, warn = FALSE)) else NA
log_rows <- sum(vapply(log_file_pat, function(f) length(readLines(f, warn = FALSE)),
                       integer(1L)))
ckp <- if (file.exists(ckp_file)) readRDS(ckp_file) else NULL

cat("\n===\n")
cat("arm:        ", arm, "\n")
cat("tree_lines: ", tree_lines, "\n")
cat("log_rows:   ", log_rows, "\n")
cat("ckp$iter:   ", ckp$iter, "\n")
cat("ckp$phase:  ", ckp$phase %||% "NA", "\n")
cat("ckp run flushed: ",
    if (is.null(ckp$runs)) "no runs" else
      paste(vapply(ckp$runs, function(r) as.character(r$flushed %||% NA),
                   character(1L)), collapse = ","), "\n")
cat("===\n")
