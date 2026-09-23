# dev/red-team/numerical/cache-coherence-driver.R
#
# Lane N3 — Partial-CL cache coherence stress driver
#
# INERT: src/mcmc.cpp has no mismatch counters or `diagMaxDiff`, so every
# readout below is NA. Per-move coherence is tested in
# tests/testthat/test-node-cl-cache.R.
#
# Purpose
# -------
# Force the partial-CL machinery (src/node_cl_cache.h, src/gibbs_partial_cl.h)
# through a long mixed-move MCMC trajectory and read out the existing
# per-move partial-vs-full mismatch counters (`diagDirMismatchCount`,
# `diagNniMismatchCount`, `diagBsMismatchCount`) plus the periodic
# `state->logLik` drift check (`diagDriftCount`).  The C++ inner loop already
# computes `partial` vs `cpp_log_likelihood(_partitioned)` at every
# partial-CL-eligible move (mcmc.cpp:4608-4727) and accumulates the largest
# discrepancy into `diagMaxDiff`.  This driver runs the chain twice — once
# under coding="variable" and once under coding="informative" — and emits
# a CSV plus a verdict file under
#   dev/red-team/numerical/cache-coherence-results/<coding>/
#
# Usage
# -----
#   Rscript dev/red-team/numerical/cache-coherence-driver.R --quick
#   Rscript dev/red-team/numerical/cache-coherence-driver.R --full
#
# --quick: 200 moves x both regimes, intended for sub-60s runs.
# --full : 1e4 moves x 5 random datasets x both regimes.
#
# Detection criteria
# ------------------
# A coherence violation is recorded when, for any partial-CL-eligible move:
#   abs(cached - full-eval) > 1e-9         (per-move, strict; LIKE-001 size
#                                            is several nats — well above)
# The C++ counters use a 1e-6 threshold to flag a "mismatch"; the driver
# also reads `diagMaxDiff` to characterise the maximum drift magnitude.
#
# Expected outcomes (Wave 1+2 prior):
#   coding="variable"   : all mismatch counters = 0, diagMaxDiff < 1e-9
#   coding="informative": mismatch counters >> 0 (LIKE-001 signature),
#                         diagMaxDiff in the 0.3–2 nat range per char-batch

suppressPackageStartupMessages({
  library(MkPrime)
  library(TreeTools)
})

# ---------------------------------------------------------------------------
# CLI parse
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
quick <- "--quick" %in% args || length(args) == 0
full  <- "--full" %in% args
if (full) quick <- FALSE

out_root <- file.path("dev", "red-team", "numerical",
                      "cache-coherence-results")
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

# Move budget per (dataset, coding) cell
n_moves_per_cell <- if (quick) 200L else 10000L
n_datasets       <- if (quick) 1L   else 5L

cat(sprintf("[N3] mode=%s  moves/cell=%d  datasets=%d\n",
            if (quick) "quick" else "full",
            n_moves_per_cell, n_datasets))

# ---------------------------------------------------------------------------
# Helper: build a small all-transformational dataset
# ---------------------------------------------------------------------------
make_setup <- function(seed) {
  set.seed(seed)
  nTip  <- 8L
  nChar <- 10L
  tree  <- ape::rtree(nTip, rooted = FALSE)
  tree  <- Preorder(tree)
  mat   <- matrix(sample(0:2, nTip * nChar, replace = TRUE),
                  nrow = nTip,
                  dimnames = list(tree$tip.label,
                                  paste0("c", seq_len(nChar))))
  pd  <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)  # no neomorphic → all transformational
  list(tree = tree, mkd = mkd)
}

# ---------------------------------------------------------------------------
# Run one MCMC and read diag counters
# ---------------------------------------------------------------------------
# We deliberately enable NNI + beta-simplex + Dirichlet (all partial-CL),
# disable Gibbs SPR / subtree-swap (these are always-accept and would mask
# acceptance-bias effects; LIKE-001 still applies to them but they are
# probed under a separate lane).  Slice samplers stay on by default and
# exercise the M-145 invalidation regression.
run_one <- function(setup, coding, n_moves, seed) {
  set.seed(seed)
  model <- MkPrimeModel(coding = coding)
  # nIter is the inner-iter count of run_mcmc_batch_cpp; total moves
  # ~= nIter * nMovesPerIter — but we set thin=1, autoTune=FALSE and
  # use nIter as the move budget driver.
  mcmc <- MkPrimeMCMC(
    nIter     = n_moves,
    minWarmup = max(50L, as.integer(n_moves %/% 4L)),
    maxWarmup = max(50L, as.integer(n_moves %/% 4L)),
    thin      = 1L,
    autoTune  = FALSE,
    nRuns     = 1L,
    gibbsSpr        = FALSE,
    gibbsSubtreeSwap = FALSE
  )
  res <- suppressWarnings(RunMkPrime(
    data  = setup$mkd,
    tree  = setup$tree,
    model = model,
    mcmc  = mcmc
  ))
  attr(res, "coding") <- coding
  res
}

# ---------------------------------------------------------------------------
# Extract a single-row summary from a run
# ---------------------------------------------------------------------------
summarise_run <- function(res, coding, dataset_id, n_moves) {
  # diag_counters lives on the run's runs[[1]] object, populated from C++
  # via run_mcmc_batch_cpp.  RunMkPrime stashes it into res$diag (per
  # existing convention) or res$runs[[1]]$diag.  We check both.
  diag <- res$diag
  if (is.null(diag)) diag <- res$runs[[1]]$diag
  if (is.null(diag)) diag <- attr(res, "diag")
  # Fallback: pcl_diag.txt is appended after each batch by mcmc.cpp:5274.
  # Parse the most recent line if direct R-side diag is missing.
  if (is.null(diag)) {
    if (file.exists("pcl_diag.txt")) {
      ln <- tail(readLines("pcl_diag.txt"), 1L)
      mm <- regmatches(ln, regexec(
        "nni=(\\d+)\\(mm=(\\d+)\\) bs=(\\d+)\\(mm=(\\d+)\\) dir=(\\d+)\\(mm=(\\d+),fb=(\\d+)\\) drift=(\\d+) maxD=([0-9.eE+-]+)",
        ln
      ))[[1]]
      if (length(mm) == 10L) {
        diag <- list(
          nni_partial    = as.integer(mm[2]),
          nni_mismatch   = as.integer(mm[3]),
          bs_partial     = as.integer(mm[4]),
          bs_mismatch    = as.integer(mm[5]),
          dir_partial    = as.integer(mm[6]),
          dir_mismatch   = as.integer(mm[7]),
          dir_fullback   = as.integer(mm[8]),
          drift          = as.integer(mm[9]),
          max_diff       = as.numeric(mm[10])
        )
      }
    }
  }
  if (is.null(diag)) {
    diag <- list(nni_partial = NA, nni_mismatch = NA,
                 bs_partial = NA,  bs_mismatch = NA,
                 dir_partial = NA, dir_mismatch = NA, dir_fullback = NA,
                 drift = NA, max_diff = NA)
  }
  data.frame(
    dataset      = dataset_id,
    coding       = coding,
    n_moves      = n_moves,
    nni_partial  = diag$nni_partial  %||% NA,
    nni_mismatch = diag$nni_mismatch %||% NA,
    bs_partial   = diag$bs_partial   %||% NA,
    bs_mismatch  = diag$bs_mismatch  %||% NA,
    dir_partial  = diag$dir_partial  %||% NA,
    dir_mismatch = diag$dir_mismatch %||% NA,
    dir_fullback = diag$dir_fullback %||% NA,
    drift        = diag$drift        %||% NA,
    max_diff     = diag$max_diff     %||% NA,
    stringsAsFactors = FALSE
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x

# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------
rows <- list()

for (ds_id in seq_len(n_datasets)) {
  setup <- make_setup(seed = 1000L + ds_id)
  for (coding in c("variable", "informative")) {
    cat(sprintf("[N3] dataset=%d coding=%-11s ... ", ds_id, coding))
    # Clear stale pcl_diag.txt before each run so fallback parse picks up
    # the right line
    if (file.exists("pcl_diag.txt")) {
      file.remove("pcl_diag.txt")
    }
    t0 <- Sys.time()
    res <- tryCatch(
      run_one(setup, coding, n_moves_per_cell, seed = 7000L + ds_id),
      error = function(e) {
        cat(sprintf("ERROR: %s\n", conditionMessage(e)))
        NULL
      }
    )
    dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (is.null(res)) {
      rows[[length(rows) + 1L]] <- data.frame(
        dataset = ds_id, coding = coding, n_moves = n_moves_per_cell,
        nni_partial = NA, nni_mismatch = NA,
        bs_partial = NA,  bs_mismatch = NA,
        dir_partial = NA, dir_mismatch = NA, dir_fullback = NA,
        drift = NA, max_diff = NA, stringsAsFactors = FALSE
      )
      next
    }
    row <- summarise_run(res, coding, ds_id, n_moves_per_cell)
    rows[[length(rows) + 1L]] <- row
    cat(sprintf("done (%.1fs)  dir_mm=%s/%s  drift=%s  maxD=%.3e\n",
                dt,
                format(row$dir_mismatch), format(row$dir_partial),
                format(row$drift),
                row$max_diff))
  }
}

drift_log <- do.call(rbind, rows)

# ---------------------------------------------------------------------------
# Emit per-coding CSV and verdict
# ---------------------------------------------------------------------------
for (coding in c("variable", "informative")) {
  out_dir <- file.path(out_root, coding)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  sub <- drift_log[drift_log$coding == coding, , drop = FALSE]
  write.csv(sub,
            file = file.path(out_dir, "drift-log.csv"),
            row.names = FALSE)

  # Verdict
  total_partial <- sum(c(sub$nni_partial, sub$bs_partial, sub$dir_partial),
                       na.rm = TRUE)
  total_mism    <- sum(c(sub$nni_mismatch, sub$bs_mismatch,
                         sub$dir_mismatch), na.rm = TRUE)
  max_drift     <- suppressWarnings(max(sub$max_diff, na.rm = TRUE))
  if (!is.finite(max_drift)) max_drift <- NA_real_
  drift_events  <- sum(sub$drift, na.rm = TRUE)

  verdict <- if (coding == "variable") {
    if (isTRUE(total_mism == 0L) && isTRUE(max_drift < 1e-9)) {
      "STABLE: no partial-vs-full mismatches detected under coding=\"variable\"."
    } else {
      sprintf("UNEXPECTED DRIFT under coding=\"variable\": %d mismatches across %d partial-CL moves, maxDiff=%.3e",
              total_mism, total_partial, max_drift)
    }
  } else {
    if (total_mism > 0L && isTRUE(max_drift > 0.1)) {
      sprintf("EXPECTED DRIFT (LIKE-001 signature reproduced) under coding=\"informative\": %d mismatches across %d partial-CL moves, maxDiff=%.3f nats, periodic drift events=%d",
              total_mism, total_partial, max_drift, drift_events)
    } else if (total_mism == 0L) {
      "NO DRIFT under coding=\"informative\" — LIKE-001 may already be patched. Confirm by inspecting src/node_cl_cache.h:697-720 for singleton_site_prob_jc calls."
    } else {
      sprintf("PARTIAL DRIFT under coding=\"informative\": %d mismatches across %d partial-CL moves, maxDiff=%.3e",
              total_mism, total_partial, max_drift)
    }
  }

  writeLines(c(
    sprintf("Lane N3 — partial-CL cache coherence — coding=\"%s\"", coding),
    sprintf("mode=%s  moves/cell=%d  datasets=%d",
            if (quick) "quick" else "full",
            n_moves_per_cell, n_datasets),
    "",
    sprintf("partial-CL moves total : %d", total_partial),
    sprintf("mismatch count (>1e-6) : %d", total_mism),
    sprintf("max |partial - full|   : %.6e", max_drift),
    sprintf("periodic drift events  : %d", drift_events),
    "",
    "Verdict:",
    verdict
  ), con = file.path(out_dir, "verdict.txt"))

  cat(sprintf("[N3] coding=%-11s -> %s\n", coding,
              file.path(out_dir, "verdict.txt")))
}

# Also emit a combined CSV at the root for convenience
write.csv(drift_log,
          file = file.path(out_root, "drift-log-combined.csv"),
          row.names = FALSE)

cat("[N3] done.\n")
