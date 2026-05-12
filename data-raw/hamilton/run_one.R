#!/usr/bin/env Rscript
# Run Mk or Mk' inference on one tree-inference replicate, or combine results.
#
# Usage:
#   Rscript run_one.R <tree_idx> <rep_idx> <data_root> <out_dir> <arm>
#
# Arguments:
#   tree_idx   integer 1-26 (tree_NN in tree-inference/)
#   rep_idx    integer 1-10 (rep_MM in tree_NN/)
#   data_root  /nobackup/pjjg18/mkprime-files/tree-inference
#   out_dir    /nobackup/pjjg18/mkp-study/results
#   arm        "mk", "mkp", "mkp_eg" (empirical_geometric prior), or "combine"

.libPaths(c("/nobackup/pjjg18/mkp-study/lib", .libPaths()))
suppressPackageStartupMessages({
  library(MkPrime)
  library(ape)
  library(TreeTools)
})

args      <- commandArgs(trailingOnly = TRUE)
tree_idx  <- as.integer(args[1])
rep_idx   <- as.integer(args[2])
data_root <- args[3]
out_dir   <- args[4]
arm       <- match.arg(args[5], c("mk", "mkp", "mkp_eg", "combine"))

cat(sprintf("tree=%d rep=%d arm=%s\n", tree_idx, rep_idx, arm))
tag <- sprintf("t%02d_r%02d", tree_idx, rep_idx)

# ---- Output directories -----------------------------------------------------
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ckp_dir <- file.path(out_dir, tag)
dir.create(ckp_dir, showWarnings = FALSE)

# ---- Stale-checkpoint guard --------------------------------------------------
.validate_ckp <- function(path) {
  if (!file.exists(path)) return(TRUE)
  tryCatch({
    ckp <- readRDS(path)
    is.list(ckp) && is.list(ckp$runs) && length(ckp$runs) > 0L &&
      is.list(ckp$mcmc) && !is.null(ckp$iter)
  }, error = function(e) FALSE)
}

.sentinel  <- file.path(ckp_dir, ".slurm_job_id")
.cur_job   <- Sys.getenv("SLURM_JOB_ID", "")
.prev_job  <- if (file.exists(.sentinel)) readLines(.sentinel, warn = FALSE)[1L] else ""

if (.cur_job != .prev_job && arm != "combine") {
  ckp_file <- file.path(ckp_dir, paste0(arm, "_checkpoint.rds"))
  if (file.exists(ckp_file) && !.validate_ckp(ckp_file)) {
    message("Corrupt ", arm, " checkpoint from job ", .prev_job, " — purging")
    unlink(list.files(ckp_dir, pattern = paste0("^", arm, "_"),
                      full.names = TRUE))
  } else if (.cur_job != "" && file.exists(ckp_file)) {
    message("Resuming valid ", arm, " checkpoint (prev job ", .prev_job,
            " -> new job ", .cur_job, ")")
  }
  if (.cur_job != "") writeLines(.cur_job, .sentinel)
}

# ---- Combine mode: read partial results, compute CID, save ------------------
if (arm == "combine") {
  suppressPackageStartupMessages(library(TreeDist))

  mk_file    <- file.path(out_dir, sprintf("mk_%s.rds", tag))
  mkp_file   <- file.path(out_dir, sprintf("mkp_%s.rds", tag))
  mkp_eg_file <- file.path(out_dir, sprintf("mkp_eg_%s.rds", tag))
  if (!file.exists(mk_file))  stop("Missing Mk result: ", mk_file)
  if (!file.exists(mkp_file)) stop("Missing Mk' result: ", mkp_file)

  mk_res  <- readRDS(mk_file)
  mkp_res <- readRDS(mkp_file)
  mkp_eg_res <- if (file.exists(mkp_eg_file)) readRDS(mkp_eg_file) else NULL

  # True tree
  tree_file <- file.path(data_root, sprintf("tree_%02d/tree.nwk", tree_idx))
  true_tree <- read.tree(tree_file)
  true_tree <- UnrootTree(true_tree)
  true_tree <- ape::reorder.phylo(true_tree, "cladewise")

  reorder_multiPhylo <- function(trees) {
    out <- lapply(trees, ape::reorder.phylo, order = "cladewise")
    class(out) <- "multiPhylo"
    out
  }

  cid_mk  <- as.numeric(
    ClusteringInfoDistance(reorder_multiPhylo(mk_res$trees), true_tree,
                          normalize = TRUE))
  cid_mkp <- as.numeric(
    ClusteringInfoDistance(reorder_multiPhylo(mkp_res$trees), true_tree,
                          normalize = TRUE))
  cid_mkp_eg <- if (!is.null(mkp_eg_res)) as.numeric(
    ClusteringInfoDistance(reorder_multiPhylo(mkp_eg_res$trees), true_tree,
                          normalize = TRUE)) else NA_real_

  result <- list(
    tree_idx           = tree_idx,
    rep_idx            = rep_idx,
    n_char             = mkp_res$n_char,
    kObs               = mkp_res$kObs,
    cid_mk             = cid_mk,
    cid_mkp            = cid_mkp,
    cid_mkp_eg         = cid_mkp_eg,
    n_samples_mk       = length(mk_res$trees),
    n_samples_mkp      = length(mkp_res$trees),
    n_samples_mkp_eg   = if (!is.null(mkp_eg_res)) length(mkp_eg_res$trees) else 0L,
    u_post_means       = mkp_res$u_post_means,
    u_post_means_eg    = if (!is.null(mkp_eg_res)) mkp_eg_res$u_post_means else NULL,
    stop_reason_mk     = mk_res$stop_reason,
    stop_reason_mkp    = mkp_res$stop_reason,
    stop_reason_mkp_eg = if (!is.null(mkp_eg_res)) mkp_eg_res$stop_reason else NA_character_,
    acceptance_mk      = mk_res$acceptance,
    acceptance_mkp     = mkp_res$acceptance,
    acceptance_mkp_eg  = if (!is.null(mkp_eg_res)) mkp_eg_res$acceptance else NULL
  )

  out_file <- file.path(out_dir, sprintf("result_%s.rds", tag))
  saveRDS(result, out_file)
  cat(sprintf("  Combined: %s  (Mk: %d trees, Mk': %d trees)\n",
              out_file, length(mk_res$trees), length(mkp_res$trees)))
  quit(save = "no", status = 0L)
}

# ---- Locate data -------------------------------------------------------------
dataset_dir <- file.path(data_root,
  sprintf("tree_%02d/rep_%02d", tree_idx, rep_idx))
tree_file   <- file.path(data_root, sprintf("tree_%02d/tree.nwk", tree_idx))
stopifnot(dir.exists(dataset_dir), file.exists(tree_file))

# ---- Load character data -----------------------------------------------------
nex_files <- sort(list.files(dataset_dir, pattern = "^chr[0-9]+\\.nex$",
                              full.names = TRUE))
stopifnot(length(nex_files) > 0)

mat_list <- lapply(nex_files, TreeTools::ReadCharacters)
combined_mat <- do.call(cbind, mat_list)
n_taxa_raw <- nrow(combined_mat)
n_char_raw <- ncol(combined_mat)

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

# ---- MCMC config ------------------------------------------------------------
make_mcmc <- function(prefix) {
  MkPrimeMCMC(
    nIter      = Inf,
    thin       = 10L,
    warmup     = 5000L,
    nRuns      = 2L,
    nChains    = 4L,
    heat       = 0.2,
    maxTime    = 6 * 3600,
    minEss     = 200L,
    maxPsrf    = 1.1,
    checkEvery = 500L,
    checkpointFile = file.path(ckp_dir, paste0(prefix, "_checkpoint.rds")),
    logFile        = file.path(ckp_dir, paste0(prefix, "_run.log")),
    treeFile       = file.path(ckp_dir, paste0(prefix, "_trees.nwk"))
  )
}

# Retry wrapper for checkpoint read errors
.run_arm <- function(call_fn, label) {
  tryCatch(call_fn(), error = function(e) {
    if (grepl("reading from connection|checkpoint", conditionMessage(e),
              ignore.case = TRUE)) {
      message(label, " checkpoint error — purging and retrying: ",
              conditionMessage(e))
      unlink(list.files(ckp_dir, pattern = paste0("^", label, "_"),
                        full.names = TRUE))
      call_fn()
    } else {
      stop(e)
    }
  })
}

# ---- Run the requested arm ---------------------------------------------------
if (arm == "mk") {
  # kObs from character matrix (before phyDat conversion)
  kobs_raw <- apply(combined_mat, 2L, function(col) {
    length(unique(col[!col %in% c("?", "-")]))
  })
  var_orig <- which(kobs_raw > 1L)
  kObs_for_mk <- setNames(as.integer(kobs_raw[var_orig]),
                           as.character(var_orig))

  res <- .run_arm(function() {
    RunMkPrime(
      pd,
      start_tree,
      knownStates = kObs_for_mk,
      model = MkPrimeModel(coding = "variable"),
      mcmc  = make_mcmc("mk")
    )
  }, "mk")

  cat(sprintf("  Mk done: %d trees, stop=%s\n",
              length(res$trees), res$stop_reason))

  partial <- list(
    trees       = res$trees,
    stop_reason = res$stop_reason,
    acceptance  = res$acceptance
  )
  saveRDS(partial, file.path(out_dir, sprintf("mk_%s.rds", tag)))

} else if (arm == "mkp") {
  mkd_mkp <- MkPrimeData(pd)

  res <- .run_arm(function() {
    RunMkPrime(
      mkd_mkp,
      start_tree,
      model = MkPrimeModel(coding = "variable"),
      mcmc  = make_mcmc("mkp")
    )
  }, "mkp")

  cat(sprintf("  Mk' done: %d trees, stop=%s\n",
              length(res$trees), res$stop_reason))

  # Posterior u means
  if (is.null(res$samples) || nrow(res$samples) == 0L) {
    res$samples <- ReadMkLog(res$logFile)
  }
  kp_cols <- grep("^kPrime_", colnames(res$samples), value = TRUE)
  kObs <- mkd_mkp$kObs
  if (length(kp_cols) > 0) {
    k_post_means <- colMeans(res$samples[, kp_cols, drop = FALSE])
    u_post_means <- k_post_means - kObs
  } else {
    u_post_means <- rep(NA_real_, mkd_mkp$nChar)
  }

  partial <- list(
    trees        = res$trees,
    stop_reason  = res$stop_reason,
    acceptance   = res$acceptance,
    n_char       = mkd_mkp$nChar,
    kObs         = as.integer(kObs),
    u_post_means = u_post_means
  )
  saveRDS(partial, file.path(out_dir, sprintf("mkp_%s.rds", tag)))

} else if (arm == "mkp_eg") {
  # Mk' with empirical_geometric prior on k' (convolution of empirical N_obs
  # pmf with Geometric(p) prior on N_unobs).  Same likelihood as "mkp", just
  # a different prior on k'.
  mkd_mkp <- MkPrimeData(pd)

  res <- .run_arm(function() {
    RunMkPrime(
      mkd_mkp,
      start_tree,
      model = MkPrimeModel(coding = "variable",
                            kPrimePrior = "empirical_geometric"),
      mcmc  = make_mcmc("mkp_eg")
    )
  }, "mkp_eg")

  cat(sprintf("  Mk' (empirical_geometric) done: %d trees, stop=%s\n",
              length(res$trees), res$stop_reason))

  if (is.null(res$samples) || nrow(res$samples) == 0L) {
    res$samples <- ReadMkLog(res$logFile)
  }
  kp_cols <- grep("^kPrime_", colnames(res$samples), value = TRUE)
  kObs <- mkd_mkp$kObs
  if (length(kp_cols) > 0) {
    k_post_means <- colMeans(res$samples[, kp_cols, drop = FALSE])
    u_post_means <- k_post_means - kObs
  } else {
    u_post_means <- rep(NA_real_, mkd_mkp$nChar)
  }

  partial <- list(
    trees        = res$trees,
    stop_reason  = res$stop_reason,
    acceptance   = res$acceptance,
    n_char       = mkd_mkp$nChar,
    kObs         = as.integer(kObs),
    u_post_means = u_post_means
  )
  saveRDS(partial, file.path(out_dir, sprintf("mkp_eg_%s.rds", tag)))
}

cat("  Done.\n")
