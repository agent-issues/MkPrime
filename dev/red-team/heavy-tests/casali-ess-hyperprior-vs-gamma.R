#!/usr/bin/env Rscript
#
# casali-ess-hyperprior-vs-gamma.R — measure the ESS gain delivered by
# the pooled half-normal hyperprior on per-class `class_rate_log_sd`
# (commit 35f65fd) against the legacy independent-gamma prior, on the
# AutoPart Casali production cells where the gamma prior leaves the
# small-class σ_c diffusing.
#
# Per invocation: ONE (matrix, treatment, prior) MCMC at the production
# config (nGen, thin, nChains from auto-part config.yml). Designed to be
# fanned out as a SLURM array — one task per (cell × prior) — with seeds
# paired across the prior pair so the only difference between the two
# runs is the prior on σ_c.
#
# CLI:
#   --matrix NAME        (required) AutoPart casali_matrices entry
#   --treatment T0|T1|T2a|T2b|T4  (required) — T3 not supported (needs Layer 2)
#   --prior  hyperprior_pooled|gamma_independent  (required)
#   --partitions-dir DIR (required) cache/partitions/ from auto-part repo
#   --out FILE           (required) RDS to write
#   --n-gen N            default 5000000     (auto-part config.yml nGen)
#   --thin  N            default 1000        (auto-part config.yml thin)
#   --burnin-frac F      default 0.25        (auto-part config.yml burninFrac)
#   --n-chains N         default 4           (auto-part config.yml chains)
#   --n-cat N            default 4           (auto-part config.yml nCat)
#   --cores N            default 1           (within-run cores)
#   --seed N             default deterministic from (matrix, treatment)
#                                       — IDENTICAL for both priors of the
#                                       same cell, see .CellSeed below
#   --autopart-lib DIR   library path that has installed AutoPart (so
#                        casali_matrices loads via data())
#
# Output RDS (a list):
#   cell           : list(matrix, treatment, nClasses, classSizes, partition)
#   prior          : "hyperprior_pooled" | "gamma_independent"
#   ess_sigma      : numeric, per-class σ_c ESS (names class<c>_rate_log_sd)
#   ess_tree_len   : tree_length ESS
#   ess_log_post   : log_posterior ESS
#   ess_hyper_tau  : hyper_tau ESS (NA under gamma_independent)
#   nSamples       : nrow(result$samples) after burn-in trimming
#   wallSecs       : MCMC wall-clock seconds
#   config         : echoed parsed CLI options
#   mkpHead        : git short SHA of the active MkPrime build
#   sessionInfo    : utils::sessionInfo()
#
# Modelled on dev/red-team/heavy-tests/funnel-stress-hyperprior-sigma.R
# but loads real AutoPart Casali cells instead of synthetic fixtures.

suppressPackageStartupMessages({
  library(optparse)
  library(coda)
})

`%||%` <- function(a, b) if (is.null(a)) b else a


# ---------- CLI ----------------------------------------------------------

.parseArgs <- function(argv = commandArgs(trailingOnly = TRUE)) {
  optList <- list(
    make_option("--matrix",          type = "character", default = NULL),
    make_option("--treatment",       type = "character", default = NULL),
    make_option("--prior",           type = "character", default = NULL),
    make_option("--partitions-dir",  type = "character", default = NULL),
    make_option("--out",             type = "character", default = NULL),
    make_option("--n-gen",           type = "double",    default = 5e6),
    make_option("--thin",            type = "integer",   default = 1000L),
    make_option("--burnin-frac",     type = "double",    default = 0.25),
    make_option("--n-chains",        type = "integer",   default = 4L),
    make_option("--n-cat",           type = "integer",   default = 4L),
    make_option("--cores",           type = "integer",   default = 1L),
    make_option("--seed",            type = "integer",   default = NA_integer_),
    make_option("--autopart-lib",    type = "character", default = NULL)
  )
  opt <- parse_args(OptionParser(option_list = optList,
                                  usage = "Rscript casali-ess-hyperprior-vs-gamma.R [options]"),
                    args = argv)

  miss <- c(if (is.null(opt$matrix))             "--matrix",
            if (is.null(opt$treatment))          "--treatment",
            if (is.null(opt$prior))              "--prior",
            if (is.null(opt[["partitions-dir"]])) "--partitions-dir",
            if (is.null(opt$out))                "--out")
  if (length(miss)) {
    stop("Missing required argument(s): ", paste(miss, collapse = ", "),
         call. = FALSE)
  }
  if (!opt$prior %in% c("hyperprior_pooled", "gamma_independent")) {
    stop("--prior must be hyperprior_pooled or gamma_independent",
         call. = FALSE)
  }
  if (!opt$treatment %in% c("T0", "T1", "T2a", "T2b", "T4")) {
    stop("--treatment must be one of T0, T1, T2a, T2b, T4 (T3 unsupported)",
         call. = FALSE)
  }
  opt$nGen        <- as.integer(opt[["n-gen"]])
  opt$thin        <- opt$thin
  opt$burninFrac  <- opt[["burnin-frac"]]
  opt$nChains     <- opt[["n-chains"]]
  opt$nCat        <- opt[["n-cat"]]
  opt$cores       <- opt$cores
  opt$partsDir    <- opt[["partitions-dir"]]
  opt$autopartLib <- opt[["autopart-lib"]]
  opt
}


# ---------- Seed pairing across priors -----------------------------------
#
# Same (matrix, treatment) -> same seed under BOTH priors, so the only
# difference between paired runs is the σ_c prior structure.

.CellSeed <- function(matrixName, treatment) {
  key   <- list(matrixName, treatment)
  hex   <- substr(digest::digest(key), 1L, 8L)
  bytes <- as.integer(as.hexmode(strsplit(hex, "")[[1L]]))
  v     <- sum(bytes * 16 ^ (7L:0L))
  as.integer(v %% .Machine$integer.max)
}


# ---------- Cell loader --------------------------------------------------

.LoadCell <- function(matrixName, treatment, partsDir) {
  if (!requireNamespace("AutoPart", quietly = TRUE)) {
    stop("Package 'AutoPart' must be installed (provides casali_matrices). ",
         "Set --autopart-lib if installed to a non-default library.",
         call. = FALSE)
  }
  env <- new.env(parent = emptyenv())
  utils::data("casali_matrices", package = "AutoPart", envir = env)
  cm <- env$casali_matrices
  if (!matrixName %in% names(cm)) {
    stop("Matrix '", matrixName, "' not in casali_matrices.", call. = FALSE)
  }

  partFile <- file.path(partsDir, sprintf("%s__%s.rds", matrixName, treatment))
  if (!file.exists(partFile)) {
    stop("Partition file not found: ", partFile, call. = FALSE)
  }
  partRec <- readRDS(partFile)
  if (!identical(partRec$matrix, matrixName) ||
      !identical(partRec$treatment, treatment)) {
    stop("Partition file disagrees with --matrix/--treatment: ", partFile,
         call. = FALSE)
  }

  rawMat  <- cm[[matrixName]]$matrix
  keep    <- AutoPart::FilterCharacters(rawMat, return = "keep")
  filtMat <- rawMat[, keep, drop = FALSE]
  if (ncol(filtMat) != length(partRec$partition) &&
      !(treatment == "T0" && length(partRec$partition) == ncol(filtMat))) {
    stop(sprintf(
      "Character-count mismatch: NVI-filtered matrix has %d chars, partition %d.",
      ncol(filtMat), length(partRec$partition)), call. = FALSE)
  }
  filtMat[is.na(filtMat)] <- "?"
  matrixData <- TreeTools::MatrixToPhyDat(filtMat)

  list(matrixData = matrixData,
       partition  = as.integer(partRec$partition),
       nClasses   = as.integer(partRec$nClasses %||% length(unique(partRec$partition))),
       classSizes = as.integer(partRec$classSizes %||% table(partRec$partition)),
       branchModel = partRec$branchModel %||% "linked")
}


# ---------- Single-prior MCMC run ----------------------------------------

.RunOne <- function(cell, opt) {
  isT0 <- identical(opt$treatment, "T0") || identical(cell$nClasses, 1L)
  burninIter <- as.integer(opt$nGen * opt$burninFrac)

  model <- MkPrime::MkPrimeModel(
    kPrimePrior           = "geometric",
    nCat                  = opt$nCat,
    priorOnClassRateLogSd = opt$prior
  )
  mcmc <- MkPrime::MkPrimeMCMC(
    nIter      = opt$nGen,
    thin       = opt$thin,
    maxWarmup  = burninIter,
    nRuns      = 1L,
    nChains    = opt$nChains,
    nCore      = opt$cores,
    autoTune   = TRUE
  )

  callArgs <- list(data = cell$matrixData, model = model, mcmc = mcmc)
  if (isT0) {
    callArgs$partition <- NULL
    callArgs$unlink    <- character(0)
  } else {
    callArgs$partition <- cell$partition
    callArgs$unlink    <- c("shape", "ratemultiplier")
  }

  set.seed(opt$seed)
  t0 <- Sys.time()
  res <- do.call(MkPrime::RunMkPrime, callArgs)
  list(result   = res,
       wallSecs = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}


# ---------- ESS extractors -----------------------------------------------

.SafeEss <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 10L) return(NA_real_)
  as.numeric(coda::effectiveSize(coda::as.mcmc(x)))
}

.PerClassSigmaEss <- function(samples, nClasses) {
  cols <- paste0("class", seq_len(nClasses), "_rate_log_sd")
  have <- intersect(cols, colnames(samples))
  out  <- setNames(rep_len(NA_real_, length(cols)), cols)
  for (col in have) out[[col]] <- .SafeEss(samples[, col])
  out
}

.NamedEss <- function(samples, col) {
  if (!col %in% colnames(samples)) return(NA_real_)
  .SafeEss(samples[, col])
}


# ---------- Main ---------------------------------------------------------

.Main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  opt <- .parseArgs(argv)
  if (!is.null(opt$autopartLib) && nzchar(opt$autopartLib)) {
    .libPaths(c(opt$autopartLib, .libPaths()))
  }
  if (!requireNamespace("MkPrime", quietly = TRUE)) {
    stop("Package 'MkPrime' must be installed (R CMD INSTALL the worktree).",
         call. = FALSE)
  }
  if (is.na(opt$seed)) opt$seed <- .CellSeed(opt$matrix, opt$treatment)

  message("== casali-ess-hyperprior-vs-gamma.R ==")
  message("  matrix    : ", opt$matrix)
  message("  treatment : ", opt$treatment)
  message("  prior     : ", opt$prior)
  message("  seed      : ", opt$seed)
  message("  nGen      : ", opt$nGen, "  thin: ", opt$thin,
          "  burninFrac: ", opt$burninFrac,
          "  nChains: ", opt$nChains)

  cell <- .LoadCell(opt$matrix, opt$treatment, opt$partsDir)
  message("  nClasses  : ", cell$nClasses,
          "  classSizes: ", paste(cell$classSizes, collapse = ","))

  run <- .RunOne(cell, opt)

  samples  <- run$result$samples
  ess_sig  <- .PerClassSigmaEss(samples, cell$nClasses)
  ess_tl   <- .NamedEss(samples, "tree_length")
  ess_lp   <- .NamedEss(samples, "log_posterior") |>
                suppressWarnings() # may be log_likelihood depending on cols
  if (is.na(ess_lp)) ess_lp <- .NamedEss(samples, "log_likelihood")
  ess_tau  <- .NamedEss(samples, "hyper_tau")

  mkpHead <- tryCatch(
    utils::packageDescription("MkPrime")$GithubSHA1 %||%
    utils::packageDescription("MkPrime")$RemoteSha   %||% NA_character_,
    error = function(e) NA_character_
  )

  out <- list(
    cell = list(
      matrix      = opt$matrix,
      treatment   = opt$treatment,
      nClasses    = cell$nClasses,
      classSizes  = cell$classSizes,
      partition   = cell$partition,
      branchModel = cell$branchModel
    ),
    prior         = opt$prior,
    ess_sigma     = ess_sig,
    ess_tree_len  = ess_tl,
    ess_log_post  = ess_lp,
    ess_hyper_tau = ess_tau,
    nSamples      = NROW(samples),
    wallSecs      = run$wallSecs,
    config        = opt,
    mkpHead       = mkpHead,
    sessionInfo   = utils::sessionInfo(),
    writtenAt     = Sys.time()
  )
  dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
  saveRDS(out, opt$out)

  cat("\n== ESS summary ==\n")
  cat(sprintf("  tree_length:  %s\n", format(round(ess_tl, 1))))
  cat(sprintf("  log_post:     %s\n", format(round(ess_lp, 1))))
  cat(sprintf("  hyper_tau:    %s\n", format(round(ess_tau, 1))))
  cat("  σ_c per class:\n")
  for (nm in names(ess_sig)) {
    cat(sprintf("    %-30s %s\n", nm, format(round(ess_sig[[nm]], 1))))
  }
  cat(sprintf("\nWrote %s (%.1f s MCMC wall)\n", opt$out, run$wallSecs))
  invisible(0L)
}

if (sys.nframe() == 0L) {
  status <- tryCatch(.Main(), error = function(e) {
    message("ERROR: ", conditionMessage(e))
    1L
  })
  quit(status = if (is.numeric(status)) as.integer(status) else 0L)
}
