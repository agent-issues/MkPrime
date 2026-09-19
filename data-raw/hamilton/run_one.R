#!/usr/bin/env Rscript
# Run Mk or Mk' inference on one tree-inference replicate.
#
# Usage:
#   Rscript run_one.R <tree_idx> <rep_idx> <data_root> <out_dir> <arm>
#
# Arguments:
#   tree_idx   integer 1-26 (tree_NN in tree-inference/)
#   rep_idx    integer 1-10 (rep_MM in tree_NN/)
#   data_root  /nobackup/pjjg18/mkprime-files/tree-inference
#   out_dir    /nobackup/pjjg18/mkp-study/results
#   arm        one of the fourteen arms listed in `.ARMS` below

.libPaths(c("/nobackup/pjjg18/mkp-study/lib", .libPaths()))
suppressPackageStartupMessages({
  library(MkPrime)
  library(ape)
  library(TreeTools)
})
options(warn = 1L)  # print warnings as they happen; avoids spurious parser-error in exit warning summary (mk_ktrue/mk_tlshrink)

args      <- commandArgs(trailingOnly = TRUE)
tree_idx  <- as.integer(args[1])
rep_idx   <- as.integer(args[2])
data_root <- args[3]
out_dir   <- args[4]
.ARMS <- c("mk", "mk_kp1", "mk_kp2", "mk_k9", "mk_k15", "mk_k24", "mk_k40",
            "mk_ktrue", "mk_tlshrink",
            "mkp", "mkp_eg", "mkp_geo", "mkp_highk", "mkp_logs")
arm       <- match.arg(args[5], .ARMS)

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

# Files an arm owns inside the shared `ckp_dir`.
#
# `ckp_dir` is per (tree, rep), NOT per arm, so all fourteen arms share one flat
# namespace. A prefix glob on `paste0("^", arm, "_")` therefore matched every
# arm whose name extends this one: purging `mk` deleted 40 foreign files across
# `mk_k9`, `mk_k15`, `mk_k24`, `mk_k40`, `mk_kp1`, `mk_kp2`, `mk_ktrue` and
# `mk_tlshrink` -- while those arms were still running and appending (#98). On
# Lustre the victim keeps its handle on the unlinked inode, so it goes on
# looking healthy and produces nothing recoverable.
#
# Anchoring on the suffixes this arm actually writes closes that: after
# `^mk_` the next characters must be one of the three literals below, which
# `k9_run_1.log` is not. The `_[0-9]+` groups are optional so the pre-2026-05-20
# single-file layout is covered too.
.arm_own_files <- function(arm) {
  list.files(
    ckp_dir,
    pattern = sprintf("^%s_(checkpoint\\.rds|run(_[0-9]+)?\\.log|trees(_[0-9]+)?\\.nwk)$",
                      arm),
    full.names = TRUE
  )
}

.sentinel  <- file.path(ckp_dir, ".slurm_job_id")
.cur_job   <- Sys.getenv("SLURM_JOB_ID", "")
# A SIGKILL during writeLines() at the foot of this block leaves a zero-byte
# sentinel, and readLines()[1L] is then NA_character_ -- so `if (.cur_job !=
# .prev_job)` raised "missing value where TRUE/FALSE needed" on exactly the
# task that most needed to resume (#103).
.prev_job  <- if (file.exists(.sentinel)) {
  .lines <- readLines(.sentinel, warn = FALSE)
  if (length(.lines) && nzchar(.lines[1L])) .lines[1L] else ""
} else ""

if (.cur_job != .prev_job) {
  ckp_file <- file.path(ckp_dir, paste0(arm, "_checkpoint.rds"))
  if (file.exists(ckp_file) && !.validate_ckp(ckp_file)) {
    message("Corrupt ", arm, " checkpoint from job ", .prev_job, " — purging")
    unlink(.arm_own_files(arm))
  } else if (.cur_job != "" && file.exists(ckp_file)) {
    message("Resuming valid ", arm, " checkpoint (prev job ", .prev_job,
            " -> new job ", .cur_job, ")")
  }
  if (.cur_job != "") writeLines(.cur_job, .sentinel)
}

# ---- Shared helpers ----------------------------------------------------------

# Observed state count per character column.
#
# A polymorphism or uncertainty token -- `{01}`, `(01)` -- is one cell, not an
# extra state, and `MkPrimeData()` resolves it to NA (R/MkPrimeData.R:94-96).
# Counting it as a state inflated `knownStates` for the fixed-k arms above the
# kObs the package itself computes, so the baseline every other arm is measured
# against was on a different convention from the package (#104).
.kObsFromMatrix <- function(mat) {
  apply(mat, 2L, function(col) {
    length(unique(col[!col %in% c("?", "-") & !grepl("^[{(]", col)]))
  })
}

# Fraction of each run discarded before any posterior mean is taken.
#
# No task converges -- every run is truncated at the wall clock (#13) -- so a
# mean over the whole stream is a mean that includes the warmup. The fraction
# is recorded alongside the means so a consumer can tell what it is holding.
.BURNIN_FRAC <- 0.25

# Posterior mean of k' - kObs, over the post-burn-in part of EVERY run.
#
# `res$logFile` is the vector of per-run paths (R/RunMkPrime.R:2555). The old
# code took `colMeans(res$samples)`, which is whichever runs the object happened
# to hold and no burn-in at all (#102). Burn-in is applied per run and then
# pooled, because the runs are separate chains: a pooled tail can otherwise
# consist entirely of one run (#99).
.UPostMeans <- function(res, mkd, burninFrac = .BURNIN_FRAC) {
  logs <- res$logFile
  mats <- if (length(logs) && all(file.exists(logs))) {
    lapply(logs, ReadMkLog)
  } else if (!is.null(res$samples) && nrow(res$samples) > 0L) {
    list(res$samples)
  } else {
    list()
  }

  keep <- lapply(mats, function(m) {
    n <- nrow(m)
    if (is.null(n) || n == 0L) return(NULL)
    m[seq.int(floor(n * burninFrac) + 1L, n), , drop = FALSE]
  })
  keep <- keep[!vapply(keep, is.null, logical(1L))]

  nChar <- mkd$nChar
  if (!length(keep)) {
    return(list(means = rep(NA_real_, nChar), n = 0L, nRuns = 0L,
                burninFrac = burninFrac))
  }

  pooled  <- do.call(rbind, keep)
  kp_cols <- grep("^kPrime_", colnames(pooled), value = TRUE)
  if (!length(kp_cols)) {
    return(list(means = rep(NA_real_, nChar), n = nrow(pooled),
                nRuns = length(keep), burninFrac = burninFrac))
  }

  list(means      = colMeans(pooled[, kp_cols, drop = FALSE]) - mkd$kObs,
       n          = nrow(pooled),
       nRuns      = length(keep),
       burninFrac = burninFrac)
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
# cbind joins positionally and keeps the first matrix's rownames, so two files
# listing the same taxa in different orders would be merged silently. The real
# data has one writer and a constant order, which is why this has never fired
# -- the guard is what keeps it that way (#104).
stopifnot(all(vapply(mat_list,
                     function(m) identical(rownames(m), rownames(mat_list[[1L]])),
                     logical(1L))))
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
# `thin_iters` used to default to 10 while every arm added since 2026-05 passed
# 500, so four arms sampled 50x more densely than the rest -- and, because
# `checkEvery` was a fixed iteration count, their convergence window was 1000
# samples against the others' 20. The arms were not being held to the same
# stopping criterion (#103). The default is now the value the majority already
# pass, and `checkEvery` is derived from `thin` so the window is the same number
# of SAMPLES whatever the thinning.
make_mcmc <- function(prefix, thin_iters = 500L) {
  checkEveryThin <- 50L                      # samples between convergence checks

  MkPrimeMCMC(
    nIter      = Inf,
    thin       = thin_iters,
    maxWarmup  = 5000L,
    nRuns      = 2L,
    nChains    = 4L,
    heat       = 0.1,  # widened from 0.2 to help mode-jumping under the
                        # empirical_geometric prior, which produces logP
                        # swings of +-50 between adjacent reports.
                        # Lower heat -> hotter hottest chain -> better
                        # discovery of distant modes.
    # HARNESS-001 fix: keep below SLURM wall (8h in *_array.slurm) so the
    # R-level graceful stop can fire and `saveRDS(partial, ...)` below runs
    # before SIGKILL. 7.5h gives ~30 min for shutdown / serialisation.
    maxTime    = 7.5 * 3600,
    minEss     = 200L,
    maxRhat    = 1.1,
    checkEvery = thin_iters * checkEveryThin,
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
      unlink(.arm_own_files(label))
      call_fn()
    } else {
      stop(e)
    }
  })
}

# ---- Run the requested arm ---------------------------------------------------
if (arm == "mk") {
  # kObs from character matrix (before phyDat conversion)
  kobs_raw <- .kObsFromMatrix(combined_mat)
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

} else if (arm == "mk_kp1") {
  # Mk with knownStates = kObs_i + 1 per character (one unobserved state
  # allowed). Reference point for EG vs mk floor: tests whether moving the
  # fixed-k floor up by one shifts tree recovery.
  kobs_raw <- .kObsFromMatrix(combined_mat)
  var_orig <- which(kobs_raw > 1L)
  kObs_for_mk <- setNames(as.integer(kobs_raw[var_orig]),
                           as.character(var_orig))

  res <- .run_arm(function() {
    RunMkPrime(
      pd,
      start_tree,
      knownStates = kObs_for_mk + 1L,
      model = MkPrimeModel(coding = "variable"),
      mcmc  = make_mcmc("mk_kp1")
    )
  }, "mk_kp1")

  cat(sprintf("  Mk(kObs+1) done: %d trees, stop=%s\n",
              length(res$trees), res$stop_reason))

  partial <- list(
    trees       = res$trees,
    stop_reason = res$stop_reason,
    acceptance  = res$acceptance
  )
  saveRDS(partial, file.path(out_dir, sprintf("mk_kp1_%s.rds", tag)))

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
  up   <- .UPostMeans(res, mkd_mkp)
  kObs <- mkd_mkp$kObs
  cat(sprintf("  u means over %d post-burn-in samples from %d run(s)\n",
              up$n, up$nRuns))

  partial <- list(
    trees        = res$trees,
    stop_reason  = res$stop_reason,
    acceptance   = res$acceptance,
    n_char       = mkd_mkp$nChar,
    kObs         = as.integer(kObs),
    u_post_means   = up$means,
    n_post_samples = up$n,
    n_runs_read    = up$nRuns,
    burnin_frac    = up$burninFrac
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

  up   <- .UPostMeans(res, mkd_mkp)
  kObs <- mkd_mkp$kObs
  cat(sprintf("  u means over %d post-burn-in samples from %d run(s)\n",
              up$n, up$nRuns))

  partial <- list(
    trees        = res$trees,
    stop_reason  = res$stop_reason,
    acceptance   = res$acceptance,
    n_char       = mkd_mkp$nChar,
    kObs         = as.integer(kObs),
    u_post_means   = up$means,
    n_post_samples = up$n,
    n_runs_read    = up$nRuns,
    burnin_frac    = up$burninFrac
  )
  saveRDS(partial, file.path(out_dir, sprintf("mkp_eg_%s.rds", tag)))

} else if (arm == "mkp_geo") {
  # Mk' with plain geometric prior on k' (EG-001 pilot arm).
  # Uses kPrimePrior = "geometric": k'_i ~ Geometric(p) shifted by kObs_i.
  # This removes the EG suspect normaliser (EG-001 HIGH) to test whether the
  # u_post≈1 anchor is caused by the missing per-character truncation
  # normaliser in the empirical_geometric prior.
  mkd_mkp <- MkPrimeData(pd)

  res <- .run_arm(function() {
    RunMkPrime(
      mkd_mkp,
      start_tree,
      model = MkPrimeModel(coding = "variable",
                            kPrimePrior = "geometric"),
      mcmc  = make_mcmc("mkp_geo", thin_iters = 500L)
    )
  }, "mkp_geo")

  cat(sprintf("  Mk' (geometric) done: %d trees, stop=%s\n",
              length(res$trees), res$stop_reason))

  up   <- .UPostMeans(res, mkd_mkp)
  kObs <- mkd_mkp$kObs
  cat(sprintf("  u means over %d post-burn-in samples from %d run(s)\n",
              up$n, up$nRuns))

  partial <- list(
    trees        = res$trees,
    stop_reason  = res$stop_reason,
    acceptance   = res$acceptance,
    n_char       = mkd_mkp$nChar,
    kObs         = as.integer(kObs),
    u_post_means   = up$means,
    n_post_samples = up$n,
    n_runs_read    = up$nRuns,
    burnin_frac    = up$burninFrac
  )
  saveRDS(partial, file.path(out_dir, sprintf("mkp_geo_%s.rds", tag)))

} else if (arm == "mk_kp2") {
  # Mk with knownStates = kObs_i + 2 per character (two unobserved states
  # allowed). Companion to mk_kp1: tests whether the +1 advantage persists
  # or reverses as we move further above the observed floor.
  kobs_raw <- .kObsFromMatrix(combined_mat)
  var_orig <- which(kobs_raw > 1L)
  kObs_for_mk <- setNames(as.integer(kobs_raw[var_orig]),
                           as.character(var_orig))

  res <- .run_arm(function() {
    RunMkPrime(
      pd, start_tree,
      knownStates = kObs_for_mk + 2L,
      model = MkPrimeModel(coding = "variable"),
      mcmc  = make_mcmc("mk_kp2", thin_iters = 500L)
    )
  }, "mk_kp2")

  cat(sprintf("  Mk(kObs+2) done: %d trees, stop=%s\n",
              length(res$trees), res$stop_reason))

  partial <- list(trees = res$trees, stop_reason = res$stop_reason,
                  acceptance = res$acceptance)
  saveRDS(partial, file.path(out_dir, sprintf("mk_kp2_%s.rds", tag)))

} else if (arm == "mk_k9") {
  # Mk with knownStates = 9 across all variable characters (fixed ceiling
  # comparator). 9 is a natural DNA-ish upper bound; well above observed
  # max kObs ~7 in this dataset.
  kobs_raw <- .kObsFromMatrix(combined_mat)
  var_orig <- which(kobs_raw > 1L)
  k9_for_mk <- setNames(rep(9L, length(var_orig)),
                         as.character(var_orig))

  res <- .run_arm(function() {
    RunMkPrime(
      pd, start_tree,
      knownStates = k9_for_mk,
      model = MkPrimeModel(coding = "variable"),
      mcmc  = make_mcmc("mk_k9", thin_iters = 500L)
    )
  }, "mk_k9")

  cat(sprintf("  Mk(9) done: %d trees, stop=%s\n",
              length(res$trees), res$stop_reason))

  partial <- list(trees = res$trees, stop_reason = res$stop_reason,
                  acceptance = res$acceptance)
  saveRDS(partial, file.path(out_dir, sprintf("mk_k9_%s.rds", tag)))

} else if (arm == "mk_k15") {
  # Mk with knownStates = 15 across all variable characters. Tests
  # whether the flexibility advantage of k=9 over kObs+2 continues to
  # climb at higher k.
  kobs_raw <- .kObsFromMatrix(combined_mat)
  var_orig <- which(kobs_raw > 1L)
  k15_for_mk <- setNames(rep(15L, length(var_orig)),
                          as.character(var_orig))

  res <- .run_arm(function() {
    RunMkPrime(
      pd, start_tree,
      knownStates = k15_for_mk,
      model = MkPrimeModel(coding = "variable"),
      mcmc  = make_mcmc("mk_k15", thin_iters = 500L)
    )
  }, "mk_k15")

  cat(sprintf("  Mk(15) done: %d trees, stop=%s\n",
              length(res$trees), res$stop_reason))

  partial <- list(trees = res$trees, stop_reason = res$stop_reason,
                  acceptance = res$acceptance)
  saveRDS(partial, file.path(out_dir, sprintf("mk_k15_%s.rds", tag)))

} else if (arm == "mk_k24") {
  # Mk with knownStates = 24 across all variable characters. High-k
  # endpoint; tests whether the trend plateaus or keeps climbing.
  kobs_raw <- .kObsFromMatrix(combined_mat)
  var_orig <- which(kobs_raw > 1L)
  k24_for_mk <- setNames(rep(24L, length(var_orig)),
                          as.character(var_orig))

  res <- .run_arm(function() {
    RunMkPrime(
      pd, start_tree,
      knownStates = k24_for_mk,
      model = MkPrimeModel(coding = "variable"),
      mcmc  = make_mcmc("mk_k24", thin_iters = 500L)
    )
  }, "mk_k24")

  cat(sprintf("  Mk(24) done: %d trees, stop=%s\n",
              length(res$trees), res$stop_reason))

  partial <- list(trees = res$trees, stop_reason = res$stop_reason,
                  acceptance = res$acceptance)
  saveRDS(partial, file.path(out_dir, sprintf("mk_k24_%s.rds", tag)))
} else if (arm == "mk_k40") {
  # Mk with knownStates = 40 across all variable characters. Extended
  # endpoint; tests whether the k-ramp continues past k=24.
  kobs_raw <- .kObsFromMatrix(combined_mat)
  var_orig <- which(kobs_raw > 1L)
  k40_for_mk <- setNames(rep(40L, length(var_orig)),
                          as.character(var_orig))

  res <- .run_arm(function() {
    RunMkPrime(
      pd, start_tree,
      knownStates = k40_for_mk,
      model = MkPrimeModel(coding = "variable"),
      mcmc  = make_mcmc("mk_k40", thin_iters = 500L)
    )
  }, "mk_k40")

  cat(sprintf("  Mk(40) done: %d trees, stop=%s\n",
              length(res$trees), res$stop_reason))

  partial <- list(trees = res$trees, stop_reason = res$stop_reason,
                  acceptance = res$acceptance)
  saveRDS(partial, file.path(out_dir, sprintf("mk_k40_%s.rds", tag)))
} else if (arm == "mkp_highk") {
  # Mk' with geometric prior on k' but a STRONG high-k Beta(1, 20) hyperprior
  # on p: E[p] = 1/21 ≈ 0.048, so E[k'] ≈ kObs + 21. Tests whether Mk' can
  # match mk_k40 if its prior is shifted to put mass on large k'. Same
  # likelihood as mkp_geo — only the hyperprior differs.
  mkd_mkp <- MkPrimeData(pd)

  res <- .run_arm(function() {
    RunMkPrime(
      mkd_mkp,
      start_tree,
      model = MkPrimeModel(coding = "variable",
                            kPrimePrior = "geometric",
                            kprimeHyperA = 1,
                            kprimeHyperB = 20),
      mcmc  = make_mcmc("mkp_highk", thin_iters = 500L)
    )
  }, "mkp_highk")

  cat(sprintf("  Mk' (geometric, high-k Beta(1,20)) done: %d trees, stop=%s\n",
              length(res$trees), res$stop_reason))

  up   <- .UPostMeans(res, mkd_mkp)
  kObs <- mkd_mkp$kObs
  cat(sprintf("  u means over %d post-burn-in samples from %d run(s)\n",
              up$n, up$nRuns))

  partial <- list(
    trees        = res$trees,
    stop_reason  = res$stop_reason,
    acceptance   = res$acceptance,
    n_char       = mkd_mkp$nChar,
    kObs         = as.integer(kObs),
    u_post_means   = up$means,
    n_post_samples = up$n,
    n_runs_read    = up$nRuns,
    burnin_frac    = up$burninFrac
  )
  saveRDS(partial, file.path(out_dir, sprintf("mkp_highk_%s.rds", tag)))

} else if (arm == "mkp_logs") {
  # Mk' with logseries prior on k': P(k) ∝ c^k / k with c = 0.95. Much heavier
  # right tail than geometric — approximates a "diffuse / 1/k-like" prior on
  # state-space cardinality. Tests whether a flatter prior over k lets Mk'
  # explore the high-k regime that mk_k40 implicitly inhabits.
  mkd_mkp <- MkPrimeData(pd)

  res <- .run_arm(function() {
    RunMkPrime(
      mkd_mkp,
      start_tree,
      model = MkPrimeModel(coding = "variable",
                            kPrimePrior = "logseries",
                            kprimeLogseriesC = 0.95),
      mcmc  = make_mcmc("mkp_logs", thin_iters = 500L)
    )
  }, "mkp_logs")

  cat(sprintf("  Mk' (logseries c=0.95) done: %d trees, stop=%s\n",
              length(res$trees), res$stop_reason))

  up   <- .UPostMeans(res, mkd_mkp)
  kObs <- mkd_mkp$kObs
  cat(sprintf("  u means over %d post-burn-in samples from %d run(s)\n",
              up$n, up$nRuns))

  partial <- list(
    trees        = res$trees,
    stop_reason  = res$stop_reason,
    acceptance   = res$acceptance,
    n_char       = mkd_mkp$nChar,
    kObs         = as.integer(kObs),
    u_post_means   = up$means,
    n_post_samples = up$n,
    n_runs_read    = up$nRuns,
    burnin_frac    = up$burninFrac
  )
  saveRDS(partial, file.path(out_dir, sprintf("mkp_logs_%s.rds", tag)))

} else if (arm == "mk_ktrue") {
  # Mk with knownStates = k_true per character (oracle / ceiling arm).
  # Reads ground_truth.csv from the rep directory and sets the state-space
  # cap to the simulator's k_true for each character. This is the
  # "best-possible" baseline since we feed inference the true generative
  # state-space cardinality. Real-world analyses cannot do this; mk_ktrue
  # exists only to bound where the kObs-based ramp asymptotes.
  gt_path <- file.path(dataset_dir, "ground_truth.csv")
  stopifnot(file.exists(gt_path))
  gt <- read.csv(gt_path)
  # ground_truth.csv has char_idx = 1..50 matching chr{N}.nex file number,
  # NOT lex sort order. combined_mat columns are in lex sort order
  # (chr1, chr10, chr11, ..., chr19, chr2, ..., chr29, chr3, ..., chr9).
  # Extract file-number from each lex-sorted file and look up k_true.
  file_nums <- as.integer(sub("^chr([0-9]+)\\.nex$", "\\1",
                              basename(nex_files)))
  # match() takes the first hit, so a duplicated char_idx would silently pair
  # every downstream k_true with the wrong character (#104).
  stopifnot(!anyDuplicated(gt$char_idx))
  k_true_lex <- gt$k_true[match(file_nums, gt$char_idx)]
  stopifnot(all(!is.na(k_true_lex)))

  kobs_raw <- .kObsFromMatrix(combined_mat)
  var_orig <- which(kobs_raw > 1L)
  # Sanity: k_true must be >= kObs for variable characters
  stopifnot(all(k_true_lex[var_orig] >= kobs_raw[var_orig]))
  k_for_mk <- setNames(as.integer(k_true_lex[var_orig]),
                       as.character(var_orig))
  cat(sprintf("  mk_ktrue: k_true range %d-%d, mean %.2f, vs kObs range %d-%d\n",
              min(k_for_mk), max(k_for_mk), mean(k_for_mk),
              min(kobs_raw[var_orig]), max(kobs_raw[var_orig])))

  res <- .run_arm(function() {
    RunMkPrime(
      pd, start_tree,
      knownStates = k_for_mk,
      model = MkPrimeModel(coding = "variable"),
      mcmc  = make_mcmc("mk_ktrue", thin_iters = 500L)
    )
  }, "mk_ktrue")

  cat(sprintf("  Mk(ktrue) done: %d trees, stop=%s\n",
              length(res$trees), res$stop_reason))

  partial <- list(trees = res$trees, stop_reason = res$stop_reason,
                  acceptance = res$acceptance,
                  k_true = k_for_mk)
  saveRDS(partial, file.path(out_dir, sprintf("mk_ktrue_%s.rds", tag)))

} else if (arm == "mk_tlshrink") {
  # Mk (kObs) with an explicit branch-length shrinkage prior.
  # Tests the regularisation-via-saturation hypothesis: if shorter posterior
  # tree length is what's driving mk_k40's CID advantage, then forcing TL
  # short via the prior (while keeping k=kObs) should close most of the gap
  # to mk_k40 without changing the state-space spec.
  #
  # Prior: Gamma(shape=20, rate=20/0.7) -> mean 0.7 (HALF of truth TL=1.4,
  # and well below mk_k40's posterior of ~1.2), sd ~0.157 — informative
  # enough to dominate the diffuse default and pull TL clearly below truth.
  # Default mk uses Gamma(2, 2/FitchScore) which has mean ~Fitch score
  # (~30-100), effectively diffuse, so the data determines TL.
  # The hypothesis: if regularisation-via-short-TL is what makes mk_k40 win,
  # then a prior pulling TL well below truth should achieve at least
  # mk_k40-level CID — without changing the state-space spec.
  kobs_raw <- .kObsFromMatrix(combined_mat)
  var_orig <- which(kobs_raw > 1L)
  kObs_for_mk <- setNames(as.integer(kobs_raw[var_orig]),
                           as.character(var_orig))

  res <- .run_arm(function() {
    RunMkPrime(
      pd, start_tree,
      knownStates = kObs_for_mk,
      model = MkPrimeModel(coding = "variable",
                            treeLengthShape = 20,
                            treeLengthRate  = 20 / 0.7),
      mcmc  = make_mcmc("mk_tlshrink", thin_iters = 500L)
    )
  }, "mk_tlshrink")

  cat(sprintf("  Mk(tlshrink) done: %d trees, stop=%s\n",
              length(res$trees), res$stop_reason))

  partial <- list(trees = res$trees, stop_reason = res$stop_reason,
                  acceptance = res$acceptance)
  saveRDS(partial, file.path(out_dir, sprintf("mk_tlshrink_%s.rds", tag)))
}

cat("  Done.\n")
