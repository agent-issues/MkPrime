#!/usr/bin/env Rscript
# Run Mk and Mk' inference on one tree-inference replicate using MkPrime.
#
# Usage:
#   Rscript run_one.R <tree_idx> <rep_idx> <data_root> <out_dir>
#
# Arguments:
#   tree_idx   integer 1-26 (tree_NN in tree-inference/)
#   rep_idx    integer 1-10 (rep_MM in tree_NN/)
#   data_root  /nobackup/pjjg18/mkprime-files/tree-inference
#   out_dir    /nobackup/pjjg18/mkp-study/results

.libPaths(c("/nobackup/pjjg18/mkp-study/lib", .libPaths()))
suppressPackageStartupMessages({
  library(MkPrime)
  library(ape)
  library(TreeDist)
  library(TreeTools)
})

args      <- commandArgs(trailingOnly = TRUE)
tree_idx  <- as.integer(args[1])
rep_idx   <- as.integer(args[2])
data_root <- args[3]
out_dir   <- args[4]
test_mode <- length(args) >= 5 && args[5] == "--test"
if (test_mode) cat("*** TEST MODE: short run ***\n")

cat(sprintf("tree=%d rep=%d\n", tree_idx, rep_idx))
tag <- sprintf("t%02d_r%02d", tree_idx, rep_idx)

# ---- Locate data ------------------------------------------------------------
dataset_dir <- file.path(data_root,
  sprintf("tree_%02d/rep_%02d", tree_idx, rep_idx))
tree_file   <- file.path(data_root, sprintf("tree_%02d/tree.nwk", tree_idx))

stopifnot(dir.exists(dataset_dir), file.exists(tree_file))

# ---- Load true tree ---------------------------------------------------------
true_tree <- read.tree(tree_file)
true_tree <- UnrootTree(true_tree)
true_tree <- ape::reorder.phylo(true_tree, "cladewise")

# ---- Load character data ----------------------------------------------------
nex_files <- sort(list.files(dataset_dir, pattern = "^chr[0-9]+\\.nex$",
                              full.names = TRUE))
stopifnot(length(nex_files) > 0)

mat_list <- lapply(nex_files, TreeTools::ReadCharacters)
combined_mat <- do.call(cbind, mat_list)
n_taxa_raw <- nrow(combined_mat)
n_char_raw <- ncol(combined_mat)

# Write a temporary combined nexus file, then read as phyDat.
# This avoids phyDat constructor namespace issues (phyDat lives in phangorn,
# not directly accessible via ape::phyDat).
tmp_nex <- tempfile(fileext = ".nex")
data_list <- setNames(
  lapply(seq_len(n_taxa_raw), function(i) combined_mat[i, ]),
  rownames(combined_mat)
)
ape::write.nexus.data(data_list, file = tmp_nex, format = "standard")
pd <- TreeTools::ReadAsPhyDat(tmp_nex)
file.remove(tmp_nex)

cat(sprintf("  Loaded %d characters, %d taxa\n", n_char_raw, n_taxa_raw))

# ---- Starting tree: NJ ------------------------------------------------------
start_tree <- NJTree(pd, edgeLengths = TRUE)

# ---- Prepare data (Mk' arm) ------------------------------------------------
mkd_mkp <- MkPrimeData(pd)
n_char  <- mkd_mkp$nChar
kObs    <- mkd_mkp$kObs

# ---- knownStates for Mk arm: fix k=kObs (original column indices) ----------
# Compute kObs from combined_mat (character matrix) before phyDat conversion.
kobs_raw <- apply(combined_mat, 2L, function(col) {
  length(unique(col[!col %in% c("?", "-")]))
})
var_orig <- which(kobs_raw > 1L)
kObs_for_mk <- setNames(as.integer(kobs_raw[var_orig]), as.character(var_orig))

# ---- MCMC config (shared) ---------------------------------------------------
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ckp_dir <- file.path(out_dir, tag)
dir.create(ckp_dir, showWarnings = FALSE)

# ---- Stale-checkpoint guard --------------------------------------------------
# On Lustre/GPFS, saveRDS can return before data reaches stable storage.
# A cancelled job may leave a truncated checkpoint that passes the first readRDS
# (kernel cache hit) but fails a second read inside ResumeMkPrime.
#
# Strategy: on detecting a new SLURM_JOB_ID, VALIDATE existing checkpoints
# rather than blindly purging.
#   - CORRUPT checkpoint (readRDS fails): purge entire task directory.
#   - VALID checkpoint (from a timed-out job): preserve for resumption.
#   - No checkpoint: fresh start, nothing to do.
.validate_ckp <- function(path) {
  if (!file.exists(path)) return(TRUE)  # absent = OK (fresh start)
  tryCatch({
    ckp <- readRDS(path)
    is.list(ckp) && is.list(ckp$runs) && length(ckp$runs) > 0L &&
      is.list(ckp$mcmc) && !is.null(ckp$iter)
  }, error = function(e) FALSE)
}

.sentinel  <- file.path(ckp_dir, ".slurm_job_id")
.cur_job   <- Sys.getenv("SLURM_JOB_ID", "")
.prev_job  <- if (file.exists(.sentinel)) readLines(.sentinel, warn = FALSE)[1L] else ""

if (.cur_job != .prev_job) {
  .ckp_paths <- file.path(ckp_dir, c("mk_checkpoint.rds", "mkp_checkpoint.rds"))
  .corrupt   <- any(sapply(.ckp_paths, function(p) file.exists(p) && !.validate_ckp(p)))
  if (.corrupt) {
    message("Corrupt checkpoint from job ", .prev_job, " — purging ", ckp_dir)
    invisible(lapply(list.files(ckp_dir, full.names = TRUE), unlink))
    dir.create(ckp_dir, showWarnings = FALSE)
  } else if (.cur_job != "" && file.exists(.ckp_paths[[1L]])) {
    message("Resuming valid checkpoint (prev job ", .prev_job,
            " -> new job ", .cur_job, ")")
  }
  if (.cur_job != "") writeLines(.cur_job, .sentinel)
}

make_mcmc <- function(prefix) {
  if (test_mode) {
    MkPrimeMCMC(
      nIter      = 200L,
      thin       = 10L,
      warmup     = 100L,
      nRuns      = 2L,
      nChains    = 2L,
      heat       = 0.2,
      maxTime    = 20 * 60,  # 20 min cap
      checkEvery = 100L,
      checkpointFile = file.path(ckp_dir, paste0(prefix, "_checkpoint.rds")),
      logFile        = file.path(ckp_dir, paste0(prefix, "_run.log")),
      treeFile       = file.path(ckp_dir, paste0(prefix, "_trees.nwk"))
    )
  } else {
    MkPrimeMCMC(
      nIter      = Inf,
      thin       = 10L,
      warmup     = 5000L,
      nRuns      = 2L,
      nChains    = 4L,
      heat       = 0.2,
      maxTime    = 6 * 3600,   # 6 h safety (SLURM gives 8 h)
      minEss     = 200L,
      maxPsrf    = 1.1,
      checkEvery = 500L,
      checkpointFile = file.path(ckp_dir, paste0(prefix, "_checkpoint.rds")),
      logFile        = file.path(ckp_dir, paste0(prefix, "_run.log")),
      treeFile       = file.path(ckp_dir, paste0(prefix, "_trees.nwk"))
    )
  }
}

# Helper: run RunMkPrime and, if a checkpoint read error slips through the
# sentinel guard (e.g. under non-SLURM environments), purge and retry once.
.run_arm <- function(call_fn, label) {
  tryCatch(call_fn(), error = function(e) {
    if (grepl("reading from connection|checkpoint", conditionMessage(e),
              ignore.case = TRUE)) {
      message(label, " checkpoint error — purging and retrying fresh: ",
              conditionMessage(e))
      invisible(lapply(list.files(ckp_dir, full.names = TRUE), unlink))
      dir.create(ckp_dir, showWarnings = FALSE)
      call_fn()
    } else {
      stop(e)
    }
  })
}

# ---- Mk arm: fix k = kObs per character ------------------------------------
mk_result <- .run_arm(function() {
  RunMkPrime(
    pd,
    start_tree,
    knownStates = kObs_for_mk,
    model = MkPrimeModel(coding = "variable"),
    mcmc  = make_mcmc("mk")
  )
}, "Mk")
cat(sprintf("  Mk done: %d trees, stop=%s\n",
            length(mk_result$trees), mk_result$stop_reason))

# ---- Mk' arm: infer k -------------------------------------------------------
mkp_result <- .run_arm(function() {
  RunMkPrime(
    mkd_mkp,
    start_tree,
    model = MkPrimeModel(coding = "variable"),
    mcmc  = make_mcmc("mkp")
  )
}, "Mk'")
cat(sprintf("  Mk' done: %d trees, stop=%s\n",
            length(mkp_result$trees), mkp_result$stop_reason))

# ---- CID to true tree -------------------------------------------------------
# Reorder posterior trees so as.Splits() finds the "order" attribute it needs.
reorder_multiPhylo <- function(trees) {
  out <- lapply(trees, ape::reorder.phylo, order = "cladewise")
  class(out) <- "multiPhylo"
  out
}
cid_mk  <- as.numeric(
  ClusteringInfoDistance(reorder_multiPhylo(mk_result$trees),  true_tree, normalize = TRUE))
cid_mkp <- as.numeric(
  ClusteringInfoDistance(reorder_multiPhylo(mkp_result$trees), true_tree, normalize = TRUE))

# ---- Load samples from log if streamed to disk ------------------------------
if (is.null(mkp_result$samples) || nrow(mkp_result$samples) == 0L) {
  mkp_result$samples <- ReadMkLog(mkp_result$logFile)
}

# ---- Posterior u means (Mk' arm) -------------------------------------------
kp_cols <- grep("^kPrime_", colnames(mkp_result$samples), value = TRUE)
if (length(kp_cols) > 0) {
  k_post_means <- colMeans(mkp_result$samples[, kp_cols, drop = FALSE])
  u_post_means <- k_post_means - kObs   # u = k' - kObs
} else {
  u_post_means <- rep(NA_real_, n_char)
}

# ---- Save compact result ----------------------------------------------------
result <- list(
  tree_idx        = tree_idx,
  rep_idx         = rep_idx,
  n_char          = n_char,
  kObs            = as.integer(kObs),
  cid_mk          = cid_mk,
  cid_mkp         = cid_mkp,
  n_samples_mk    = length(mk_result$trees),
  n_samples_mkp   = length(mkp_result$trees),
  u_post_means    = u_post_means,
  stop_reason_mk  = mk_result$stop_reason,
  stop_reason_mkp = mkp_result$stop_reason,
  acceptance_mk   = mk_result$acceptance,
  acceptance_mkp  = mkp_result$acceptance
)

out_file <- file.path(out_dir, sprintf("result_%s.rds", tag))
saveRDS(result, out_file)
cat(sprintf("  Saved %s\n", out_file))
