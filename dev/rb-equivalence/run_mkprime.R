#!/usr/bin/env Rscript
# MkPrime driver for the RB-oracle benchmark.
#
# Usage:
#   Rscript dev/rb-equivalence/run_mkprime.R <pid> <model> \
#     [--rhat=1.025] [--ess=128] [--max-time=3600] [--out-dir=dev/rb-equivalence/out]
#
# <model> ∈ {by_nt_9v, by_nt_kv}.
#
# Output:
#   dev/rb-equivalence/out/mkprime_<pid>_<model>.rds containing the posterior,
#   timing, host, and diagnostics. Schema documented at the bottom of the file.

# Print warnings immediately instead of queueing for R's exit-time formatter,
# which has been observed to trip an "object 'seconds' not found" error on
# multi-warning streamed runs (see MEMORY feedback_resumable_runs.md and
# mkp commit 75dd825 for the analogous fix in run_one.R).
options(warn = 1L)

suppressPackageStartupMessages({
  library(MkPrime)
  library(TreeTools)
})

# --- Locate ourselves ------------------------------------------------------

script_dir <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  fileArg <- args[grep("^--file=", args)]
  if (length(fileArg)) dirname(normalizePath(sub("^--file=", "", fileArg[1])))
  else getwd()
})()
source(file.path(script_dir, "R", "utils.R"))

# --- Parse CLI -------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("Usage: run_mkprime.R <pid> <model> [--rhat=1.025] [--ess=128] ",
       "[--max-time=3600] [--out-dir=dev/rb-equivalence/out]")
}
pid <- args[[1]]
model <- args[[2]]

opt <- list(rhat = 1.025, ess = 128, max_time = 3600,
            out_dir = file.path(script_dir, "out"))
for (a in args[-(1:2)]) {
  kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
  if (length(kv) != 2L) stop("Bad option: ", a)
  opt[[gsub("-", "_", kv[[1]])]] <- kv[[2]]
}
opt$rhat <- as.numeric(opt$rhat)
opt$ess <- as.numeric(opt$ess)
opt$max_time <- as.numeric(opt$max_time)
dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("[run_mkprime] pid=%s model=%s target=(rhat<%.3f, ess>%d) max_time=%ds\n",
            pid, model, opt$rhat, opt$ess, opt$max_time))

# --- Load frozen matrix selection ------------------------------------------

matrices_path <- file.path(script_dir, "matrices.rds")
if (!file.exists(matrices_path)) {
  stop("matrices.rds missing; run select_matrices.R first")
}
matrices <- readRDS(matrices_path)
info <- matrices[[pid]]
if (is.null(info)) stop("pid '", pid, "' not in matrices.rds")

# --- Read splits, combine, build phyDat ------------------------------------

neotrans_proj <- if (nzchar(Sys.getenv("NEOTRANS"))) {
  file.path(Sys.getenv("NEOTRANS"), "inst", "projects")
} else {
  normalizePath(file.path(script_dir, "..", "..", "..", "neotrans",
                          "inst", "projects"),
                mustWork = TRUE)
}

trans_mat <- ReadSplit(file.path(neotrans_proj,
                                 sprintf("project%s.trans.nex", pid)))
neo_mat <- ReadSplit(file.path(neotrans_proj,
                               sprintf("project%s.neo.nex", pid)))
comb <- CombineSplits(trans_mat, neo_mat)
stopifnot(comb$nTrans == info$nTrans, comb$nNeo == info$nNeo)

phy <- MatrixToCombinedPhyDat(comb$matrix)

# --- Build MkPrimeData -----------------------------------------------------

neoIdx <- BuildNeoIdx(info, model)
knownStates <- BuildKnownStates(info, model)

mkd <- MkPrimeData(phy, neomorphic = neoIdx, knownStates = knownStates)

cat(sprintf("[run_mkprime] data: %d tips, %d chars (after invariant drop), \\
%d neo, %d trans-known\n",
            mkd$nTip, mkd$nChar,
            sum(mkd$type == "neomorphic"),
            sum(mkd$type == "known")))

# --- Build model -----------------------------------------------------------

mkm <- RBMatchedModel()

# --- Build MCMC ------------------------------------------------------------

run_id <- sprintf("mkprime_%s_%s", pid, model)
log_file <- file.path(opt$out_dir, sprintf("%s.log", run_id))
tree_file <- file.path(opt$out_dir, sprintf("%s.trees", run_id))
checkpoint_file <- file.path(opt$out_dir, sprintf("%s.ckp", run_id))

mcmc <- MkPrimeMCMC(
  nRuns = 2L,
  nChains = 4L,        # PT minimum per the plan
  heat = 0.2,
  minEss = opt$ess,
  maxRhat = opt$rhat,
  # minTreeEss intentionally NOT set -- tree convergence is reported, not gating
  maxTime = opt$max_time,
  checkpointFile = checkpoint_file,
  treeFile = tree_file,
  logFile = log_file,
  checkEvery = 500L
)

# --- Run -------------------------------------------------------------------

t0 <- Sys.time()
posterior <- RunMkPrime(data = mkd, model = mkm, mcmc = mcmc)
t1 <- Sys.time()
wall_total <- as.numeric(difftime(t1, t0, units = "secs"))

cat(sprintf("[run_mkprime] wall_total = %.1f s\n", wall_total))

# --- Diagnostics (scalar + tree) -------------------------------------------

diag_scalar <- ConvergenceDiagnostics(posterior, trees = FALSE)
diag_tree <- tryCatch(
  ConvergenceDiagnostics(posterior, trees = TRUE, frechetESS = TRUE),
  error = function(e) { cli::cli_warn("Tree ESS failed: {e$message}"); NULL }
)

cat(sprintf("[run_mkprime] scalar: minEss=%.0f  maxRhat=%.3f\n",
            diag_scalar$minEss, diag_scalar$maxRhat %||% NA))
if (!is.null(diag_tree)) {
  cat(sprintf("[run_mkprime] tree:   medianPseudoESS=%.0f  frechetESS=%s\n",
              diag_tree$treeEss[["medianPseudoESS"]] %||% NA,
              format(diag_tree$treeEss[["frechetCorrelationESS"]] %||% NA)))
}

# --- Pro-rated wall-to-tree-target estimate --------------------------------

target_ess <- opt$ess
treeEssNow <- diag_tree$treeEss[["medianPseudoESS"]] %||% NA_real_
wall_to_tree_target_est <- if (is.finite(treeEssNow) && treeEssNow > 0) {
  wall_total * (target_ess / treeEssNow)
} else NA_real_

# --- Materialise per-run scalar traces from streamed logs ------------------
# In streaming mode the per-run scalar samples live on disk only; load them
# now so the rds is self-contained for compare.R.

per_run_log_files <- sprintf("%s_%d.log", tools::file_path_sans_ext(log_file),
                             seq_len(2L))
per_run_scalars <- lapply(per_run_log_files, function(f) {
  if (!file.exists(f)) {
    cli::cli_warn("Missing per-run log: {f}"); return(NULL)
  }
  m <- ReadMkLog(f)
  as.data.frame(m, stringsAsFactors = FALSE)
})
# Also collect per-run trees in the same shape RB uses
per_run_trees <- lapply(posterior$per_run, function(pr) pr$trees)

# --- Save ------------------------------------------------------------------

out_path <- file.path(opt$out_dir, sprintf("%s.rds", run_id))

# Wall-to-target: posterior may carry per-iter timing; fall back to total
wall_to_target <- wall_total  # MkPrime auto-stops on rhat/ess hit, so total ≈ target

result <- list(
  pid = pid,
  model = model,
  wall_to_target = wall_to_target,
  wall_total = wall_total,
  wall_to_tree_target_estimated = wall_to_tree_target_est,
  diag_scalar = diag_scalar,
  diag_tree = diag_tree,
  per_run_scalars = per_run_scalars,
  per_run_trees = per_run_trees,
  posterior = posterior,
  host = Sys.info()[["nodename"]],
  sysinfo = list(
    sysname = Sys.info()[["sysname"]],
    release = Sys.info()[["release"]],
    machine = Sys.info()[["machine"]],
    rversion = R.version.string,
    mkprime_version = as.character(packageVersion("MkPrime")),
    started = format(t0, "%Y-%m-%dT%H:%M:%S%z"),
    finished = format(t1, "%Y-%m-%dT%H:%M:%S%z")
  ),
  coding = "variable",
  nCat = 6L,
  prior_spec = list(  # for assertion in compare.R
    treeLengthShape = 2, treeLengthRate = 2,
    rateLogSdShape = 1, rateLogSdRate = 1,
    rateLossMeanlog = 0, rateLossSdlog = 2,
    rateNeoMeanlog = 0, rateNeoSdlog = 2
  )
)
saveRDS(result, out_path)
cat(sprintf("[run_mkprime] Saved: %s\n", out_path))

# === rds schema ============================================================
# $pid                                 character
# $model                               "by_nt_9v" | "by_nt_kv"
# $wall_to_target                      numeric (seconds)
# $wall_total                          numeric (seconds)
# $wall_to_tree_target_estimated       numeric or NA
# $diag_scalar                         MkpDiagnostics (trees = FALSE)
# $diag_tree                           MkpDiagnostics (trees = TRUE) or NULL
# $posterior                           MkPosterior (full chains)
# $host                                character (nodename)
# $sysinfo                             list (OS, R, MkPrime version, timestamps)
# $coding                              "variable"
# $nCat                                6L
# $prior_spec                          list of prior hyperparams (for assert)
# ===========================================================================
