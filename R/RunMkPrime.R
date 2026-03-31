# Main MCMC entry point for MkPrime
#
# Phase 3: single chain, fixed topology, R-side loop.
# Phase 4: topology moves (NNI, SPR) via mutable state$tree.
# Phase 5: parallel tempering, independent runs, convergence, stopping.

#' Run Bayesian MCMC under the MkPrime model
#'
#' Metropolis-Hastings MCMC sampling tree topology, branch lengths, model
#' parameters, and per-character k' (for transformational characters).
#' Supports parallel tempering with a geometric temperature ladder,
#' multiple independent runs, convergence monitoring, and early stopping.
#'
#' @param data A `phyDat` object or `MkPrimeData` object.
#' @param tree A `phylo` object (starting topology), or `NULL` (default)
#'   to start from a neighbour-joining tree built from the data.
#' @param neomorphic,knownStates Passed to [MkPrimeData()] if `data` is
#'   a `phyDat` object.
#' @param model An `MkPrimeModel` object, or `NULL` for defaults.
#' @param mcmc An `MkPrimeMCMC` object, or `NULL` for defaults.
#'   When `NULL`, any additional arguments in `...` are forwarded to
#'   [MkPrimeMCMC()], so MCMC options can be passed inline without
#'   constructing a separate object.
#' @param fixTopology Logical. If `TRUE`, tree topology is fixed (Phase 3
#'   behaviour). Default `FALSE` enables NNI and SPR topology proposals.
#' @param overwrite Logical. If `FALSE` (the default) and
#'   `mcmc$checkpointFile` points to an existing file, the run is
#'   automatically resumed from that checkpoint.
#'   Set to `TRUE` to discard the existing checkpoint and start fresh.
#' @param ... Additional arguments forwarded to [MkPrimeMCMC()]. Allows
#'   passing MCMC configuration inline (e.g. `nIter`, `logFile`, `nChains`)
#'   without constructing a separate object. Cannot be combined with an
#'   explicit `mcmc` argument.
#'
#' @return An `MkPosterior` object.
#'
#' @section Inline MCMC options:
#'
#' For quick runs, pass MCMC options directly:
#'
#' ```r
#' result <- RunMkPrime(data, nIter = 10000, nRuns = 1L)
#' ```
#'
#' For complex configurations, construct the object explicitly:
#'
#' ```r
#' mcmc <- MkPrimeMCMC(nRuns = 2, nChains = 4, minEss = 200)
#' result <- RunMkPrime(data, mcmc = mcmc)
#' ```
#'
#' @section Parallel independent runs (HPC usage):
#'
#' Set `parallel = TRUE` in [MkPrimeMCMC()] to run independent chains
#' concurrently. The \pkg{future} package (in `Suggests`) provides the
#' backend-agnostic parallelism. Call `future::plan()` **before**
#' `RunMkPrime()`:
#'
#' ```r
#' # Workstation — use local cores
#' future::plan("multisession", workers = 4)
#' result <- RunMkPrime(data, tree,
#'   mcmc = MkPrimeMCMC(nRuns = 4, parallel = TRUE))
#'
#' # HPC (SLURM) — requires the future.batchtools package
#' future::plan(future.batchtools::batchtools_slurm(
#'   resources = list(ncpus = 1, memory = "4gb", walltime = "24:00:00")
#' ))
#' result <- RunMkPrime(data, tree,
#'   mcmc = MkPrimeMCMC(nRuns = 4, parallel = TRUE,
#'                       logFile  = "/scratch/myrun/run.log",
#'                       pollInterval = 60L))
#' ```
#'
#' Notes for HPC:
#' - Set `logFile` explicitly to a path on a shared filesystem (tempdir is
#'   node-local and workers on different nodes cannot read it).
#' - Use `pollInterval = 30` -- `60` seconds; job startup latency makes
#'   frequent polling wasteful.
#' - Set `checkpointFile` so runs can be resumed if the master job times out.
#' - Each `future` worker becomes a separate job submission on SLURM/PBS/LSF.
#'
#' @export
RunMkPrime <- function(data, tree = NULL,
                       neomorphic = integer(0),
                       knownStates = integer(0),
                       model = NULL,
                       mcmc = NULL,
                       fixTopology = FALSE,
                       overwrite = FALSE,
                       ...) {

  # --- Build or validate MCMC config ---
  dots <- list(...)
  if (length(dots) > 0L && !is.null(mcmc)) {
    cli::cli_abort(c(
      "Supply MCMC options via {.arg mcmc} or {.code ...}, not both.",
      "i" = "Either pass {.code mcmc = MkPrimeMCMC(...)}, or pass \\
            MCMC arguments directly (e.g. {.code nIter = 50000})."
    ))
  }
  if (is.null(mcmc)) {
    mcmc <- do.call(MkPrimeMCMC, dots)
  }

  # --- Auto-resume from checkpoint ---
  cpFile <- mcmc$checkpointFile
  if (!overwrite && !is.null(cpFile) && file.exists(cpFile)) {
    cli::cli_alert_info("Resuming from checkpoint {.file {cpFile}}.")
    return(ResumeMkPrime(
      checkpointFile = cpFile,
      data = data,
      tree = tree,
      neomorphic = neomorphic,
      knownStates = knownStates,
      model = model
    ))
  }

  # Discard any temp logs from a previously interrupted run that the

  # user chose not to recover.
  .CleanupStaleTempLogs()

  # --- Input processing ---
  if (inherits(data, "MkPrimeData")) {
    mkd <- data
  } else {
    mkd <- MkPrimeData(data, neomorphic = neomorphic,
                       knownStates = knownStates)
  }

  if (is.null(tree)) {
    njInput <- if (inherits(data, "phyDat")) data else mkd$phyDat
    tree <- TreeTools::NJTree(njInput, edgeLengths = TRUE)
    cli::cli_alert_info("No starting tree supplied; using neighbour-joining tree.")
  } else if (!inherits(tree, "phylo")) {
    cli::cli_abort("{.arg tree} must be a {.cls phylo} object, or {.val NULL} to use a neighbour-joining tree.")
  }
  if (is.null(tree$edge.length)) {
    cli::cli_abort(c(
      "{.arg tree} has no branch lengths.",
      "i" = "Supply a tree with edge lengths, e.g. \\
             {.code TreeTools::NJTree(data)}."
    ))
  }
  nNeg <- sum(tree$edge.length <= 0)
  if (nNeg > 0L) {
    cli::cli_warn(c(
      "{nNeg} non-positive branch length{?s} clamped to 1e-8.",
      "i" = "Zero or negative lengths arise in NJ trees when taxa are very similar. \\
             They are invalid for likelihood computation."
    ))
    tree$edge.length[tree$edge.length <= 0] <- 1e-8
  }

  # Tip labels must match data taxon names
  treeTips <- tree$tip.label
  dataTaxa <- rownames(mkd$matrix)
  missing  <- setdiff(dataTaxa, treeTips)
  extra    <- setdiff(treeTips, dataTaxa)
  if (length(missing) > 0L || length(extra) > 0L) {
    msgs <- character(0)
    if (length(missing) > 0L) {
      msgs <- c(msgs,
        "x" = "{length(missing)} taxon{?/a} in data but not in tree: {.val {missing}}.")
    }
    if (length(extra) > 0L) {
      msgs <- c(msgs,
        "x" = "{length(extra)} tip{?s} in tree but not in data: {.val {extra}}.")
    }
    cli::cli_abort(c(
      "Tip labels in {.arg tree} do not match taxa in {.arg data}.",
      msgs,
      "i" = "Every taxon in the data must appear as a tip label in the tree, and vice versa."
    ))
  }

  if (is.null(model)) model <- MkPrimeModel()
  # mcmc already defaulted above (before auto-resume check)

  model <- .FinalizeModel(model, tree, mkd)

  # PREORDER INVARIANT: all topology proposals maintain canonical preorder.
  tree <- TreeTools::Preorder(tree)
  nEdge <- nrow(tree$edge)
  tipLabels <- tree$tip.label

  hasNeo <- any(mkd$type == "neomorphic")
  transIdx <- which(mkd$type == "transformational")
  nTrans <- length(transIdx)

  qHet <- isTRUE(model$qHeterogeneity)
  moves <- .BuildMoves(nEdge, nTrans, hasNeo, mcmc,
                       fixTopology = fixTopology,
                       kPrimePrior = model$kPrimePrior %||% "geometric",
                       qHeterogeneity = qHet,
                       joint2d = isTRUE(mcmc$joint2d))

  mcmc$thinWasAuto <- identical(mcmc$thin, "auto")
  if (mcmc$thinWasAuto) {
    mcmc$thin <- length(moves)
  }
  # Resolve treeThin: NULL → same as thin; validate multiple-of-thin
  mcmc$treeThinWasAuto <- is.null(mcmc$treeThin)
  if (mcmc$treeThinWasAuto) {
    mcmc$treeThin <- mcmc$thin
  }
  if (mcmc$treeThin %% mcmc$thin != 0L) {
    cli::cli_abort(
      "{.arg treeThin} ({mcmc$treeThin}) must be a multiple of \\
       {.arg thin} ({mcmc$thin})."
    )
  }
  treeEvery <- as.integer(mcmc$treeThin / mcmc$thin)

  nRuns <- mcmc$nRuns

  # --- Initialize per-run state ---
  runs <- vector("list", nRuns)
  for (run in seq_len(nRuns)) {
    startTree <- if (run == 1L) tree else .PerturbStart(tree)
    runs[[run]] <- .InitRun(startTree, mkd, model, mcmc, moves)
  }

  paramNames  <- .ParamNames(mkd, nEdge,
                             kPrimePrior = model$kPrimePrior %||% "geometric",
                             qHeterogeneity = qHet)

  # --- Log file setup ---
  # Always stream to a log file for interrupt recovery.  When the user
  # didn't supply logFile, write to a temp file and load samples into
  # memory on clean completion.
  userLogFile <- mcmc$logFile
  # Capture user-supplied checkpointFile before auto-derivation from temp log.
  # If the user (or MkPrimeMCMC auto-derive from logFile) set it, preserve it.
  userCkpFile <- mcmc$checkpointFile
  isTempLog   <- is.null(userLogFile)
  if (isTempLog) {
    mcmc$logFile <- tempfile("mkp_run_", fileext = ".log")
  }
  # Always checkpoint — derive from log path if not already set
  if (is.null(mcmc$checkpointFile)) {
    mcmc$checkpointFile <- sub("\\.[^.]+$", ".ckp", mcmc$logFile)
  }
  isStreaming     <- TRUE
  convWindowSize  <- .ComputeConvWindowSize(mcmc)
  logFilePaths    <- .OpenLogFiles(mcmc$logFile, paramNames, nRuns)

  # Register temp files so cleanup can find them (crash, new run, etc.)
  # Only treat the checkpoint as temp if it was auto-derived from the temp log.
  # A user-supplied checkpointFile should survive cleanup.
  isTempCkp <- isTempLog && is.null(userCkpFile)
  tempFiles <- if (isTempLog) {
    c(logFilePaths, if (isTempCkp) mcmc$checkpointFile)
  } else {
    NULL
  }
  .mkp_env$active_temp_logs <- tempFiles

  # Clean up temp files on normal exit or error — but NOT on interrupt,
  # where we want MkPrimeRecover() to find them.
  if (isTempLog) {
    on.exit(.CleanupTempLogs(tempFiles), add = TRUE)
  }

  treeFile <- mcmc$treeFile
  if (!is.null(treeFile)) writeLines("", treeFile)

  # Column indices for tree reconstruction in scalar_samples (1-based R).
  # Layout: log_post, log_lik, tree_length, [rate_loss — if hasNeo],
  #         rate_log_sd, [p — geometric only], [rate_neo — if hasNeo],
  #         [beta_scale — if qHet], kPrime_i..., br_j...
  isLogseries <- identical(model$kPrimePrior, "logseries")
  pCols       <- if (isLogseries) 0L else 1L
  neoCols     <- if (hasNeo) 2L else 0L   # rate_loss + rate_neo
  brColStart  <- 4L + neoCols + pCols + qHet + nTrans + 1L

  # --- Execute MCMC (with interrupt recovery) ---
  execResult <- .RunWithRecovery(
    mkd, model, mcmc, runs, moves, tipLabels, paramNames, nEdge, brColStart,
    logFilePaths, convWindowSize, treeFile, isTempLog
  )

  # If interrupted, on.exit cleanup is cancelled and we return early
  if (identical(execResult, "interrupted")) {
    # Cancel the on.exit cleanup — temp logs must survive for recovery
    on.exit(NULL, add = FALSE)
    return(invisible(NULL))
  }

  runs       <- execResult$runs
  stopReason <- execResult$stopReason
  actualIter <- execResult$actualIter

  # --- Build result ---
  result <- .BuildResult(runs, model, mkd, mcmc, paramNames, logFilePaths,
                         actualIter, stopReason, isTempLog = isTempLog)

  # For temp-log runs, load samples into memory so the result is
  # self-contained (temp files will be deleted by on.exit).
  if (isTempLog && result$nSamples > 0L) {
    result$samples <- ReadMkLog(logFilePaths)
    # Load per-run samples too (streaming sets them to NULL)
    if (!is.null(result$per_run)) {
      for (i in seq_along(result$per_run)) {
        if (is.null(result$per_run[[i]]$samples) &&
            i <= length(logFilePaths)) {
          result$per_run[[i]]$samples <- ReadMkLog(logFilePaths[i])
        }
      }
    }
    result$logFile <- NULL
  }

  result
}


# --- Interrupt recovery helpers ---

#' Execute MCMC with interrupt recovery
#'
#' Wraps the parallel/sequential execution block in a tryCatch so that
#' Ctrl-C interrupts are caught and partial results are preserved in
#' the temp log files for later retrieval via [MkPrimeRecover()].
#'
#' @return On success, a list with `runs`, `stopReason`, `actualIter`.
#'   On interrupt, the string `"interrupted"` (after storing recovery
#'   metadata in `.mkp_env`).
#' @keywords internal
.RunWithRecovery <- function(mkd, model, mcmc, runs, moves, tipLabels,
                              paramNames, nEdge, brColStart, logFilePaths,
                              convWindowSize, treeFile, isTempLog) {
  nRuns <- mcmc$nRuns
  stopReason <- "max_iter"
  actualIter <- if (is.finite(mcmc$nIter)) mcmc$nIter else mcmc$warmup

  # M-149: Shared mutable state for interrupt-safe checkpointing.
  # Inner functions (.RunMkPrimeSingleRun, .RunSerialRuns) update this
  # environment at each batch boundary.  The interrupt handler reads it
  # to save a checkpoint even if the run hasn't completed.
  shared <- new.env(parent = emptyenv())
  shared$runs        <- runs      # initial R-serializable states
  shared$actualIter  <- 0L
  shared$phase       <- "Warmup"
  shared$moveWeights <- NULL

  # Save an initial checkpoint so the .ckp file always exists, even if

  # the interrupt fires before the first batch boundary update.
  if (!is.null(mcmc$checkpointFile)) {
    .SaveCheckpoint(runs, mcmc, 0L, paramNames, mcmc$checkpointFile,
                    model = model)
  }

  tryCatch({
    if (isTRUE(mcmc$parallel) && nRuns > 1L) {
      # Parallel path: shared env is not accessible from future workers.
      # Checkpointing handled by the orchestrator after completion.
      parResult    <- .RunParallelRuns(mkd, model, mcmc, runs, moves,
                                        tipLabels, paramNames, nEdge,
                                        brColStart, treeFile,
                                        TRUE, logFilePaths, convWindowSize)
      runs         <- parResult$runs
      logFilePaths <- parResult$logFilePaths
      stopReason   <- parResult$stopReason
      actualIter   <- parResult$actualIter

      if (!is.null(treeFile)) {
        for (r in runs) {
          for (tr in r$tree_samples) {
            if (!is.null(tr)) cat(ape::write.tree(tr), "\n",
                                  file = treeFile, append = TRUE)
          }
        }
      }

      if (!is.null(mcmc$checkpointFile)) {
        .SaveCheckpoint(runs, mcmc, actualIter, paramNames,
                        mcmc$checkpointFile, model = model)
      }
    } else if (nRuns >= 2L && !is.null(mcmc$maxRhat)) {
      # M-146: cross-run R-hat convergence orchestrator
      serialResult <- .RunSerialRuns(mkd, model, mcmc, runs, moves,
                                      tipLabels, paramNames, nEdge,
                                      brColStart, logFilePaths,
                                      convWindowSize, treeFile,
                                      shared = shared)
      runs       <- serialResult$runs
      stopReason <- serialResult$stopReason
      actualIter <- serialResult$actualIter
    } else {
      # Simple sequential path (1 run or no maxRhat)
      for (run in seq_len(nRuns)) {
        runs[[run]] <- .RunMkPrimeSingleRun(
          mkd, model, mcmc, runs[[run]], moves, tipLabels, run,
          paramNames, nEdge, brColStart,
          logFilePath    = logFilePaths[run],
          cancelFile     = mcmc$cancelFile,
          checkpointFile = if (nRuns == 1L) mcmc$checkpointFile else NULL,
          startIter      = 1L,
          isStreaming     = TRUE,
          convWindowSize = convWindowSize,
          treeFile       = treeFile,
          shared         = shared
        )
        shared$runs[[run]] <- runs[[run]]
        shared$actualIter  <- runs[[run]]$actual_iter
        stopReason <- runs[[run]]$stop_reason
        actualIter <- runs[[run]]$actual_iter
        if (stopReason == "cancelled") break
      }

      if (nRuns > 1L && !is.null(mcmc$checkpointFile)) {
        .SaveCheckpoint(runs, mcmc, actualIter, paramNames,
                        mcmc$checkpointFile, model = model)
      }
    }

    list(runs = runs, stopReason = stopReason, actualIter = actualIter)
  },
  interrupt = function(cond) {
    # M-149: Best-effort checkpoint from shared state.
    ckpSaved <- FALSE
    if (!is.null(mcmc$checkpointFile) && shared$actualIter > 0L) {
      tryCatch({
        .SaveCheckpoint(shared$runs, mcmc, shared$actualIter,
                        paramNames, mcmc$checkpointFile,
                        moveWeights = shared$moveWeights,
                        phase = shared$phase, model = model)
        ckpSaved <- TRUE
      }, error = function(e) NULL)
    }

    # Flush any buffered samples to disk (best-effort).
    # Use shared$runs when available (latest state); fall back to local runs.
    flushRuns <- if (shared$actualIter > 0L) shared$runs else runs
    for (i in seq_along(logFilePaths)) {
      tryCatch({
        if (!is.null(flushRuns[[i]]$flush_idx) &&
            flushRuns[[i]]$flush_idx > 0L) {
          .FlushBuffer(flushRuns[[i]]$flush_buf, flushRuns[[i]]$flush_idx,
                       flushRuns[[i]]$flush_iter, logFilePaths[i])
        }
      }, error = function(e) NULL)
    }

    # Store recovery metadata so MkPrimeRecover() can reconstruct results.
    .mkp_env$recovery <- list(
      logFiles   = logFilePaths,
      paramNames = paramNames,
      model      = model,
      data       = mkd,
      mcmc       = mcmc,
      isTempLog  = isTempLog,
      time       = Sys.time()
    )
    .mkp_env$active_temp_logs <- NULL  # prevent on.exit cleanup

    nSaved <- sum(vapply(logFilePaths, function(f) {
      tryCatch(length(readLines(f, warn = FALSE)) - 1L,
               error = function(e) 0L)
    }, integer(1L)))

    if (ckpSaved) {
      cli::cli_alert_warning(c(
        "Run interrupted at iteration {shared$actualIter}. \\
         {nSaved} sample{?s} saved to log file{?s}.",
        "i" = "Checkpoint saved to {.file {mcmc$checkpointFile}}.",
        "i" = "Re-run the same {.fn RunMkPrime} call to resume."
      ))
    } else {
      cli::cli_alert_warning(c(
        "Run interrupted. {nSaved} sample{?s} saved to temporary log file{?s}.",
        "i" = "Retrieve partial results: {.code posterior <- MkPrimeRecover()}"
      ))
    }

    "interrupted"
  })
}


#' Delete temporary files (logs + checkpoint)
#' @keywords internal
.CleanupTempLogs <- function(tempFiles) {
  for (f in tempFiles) {
    tryCatch(unlink(f), error = function(e) NULL)
  }
  .mkp_env$active_temp_logs <- NULL
  .mkp_env$recovery <- NULL
}


#' Discard stale temp logs from a previously interrupted run
#'
#' Called at the start of RunMkPrime() so that starting a new run
#' cleans up any leftover temp files the user chose not to recover.
#' @keywords internal
.CleanupStaleTempLogs <- function() {
  # Clean up active temp logs from a previous run that was interrupted
  # and never recovered.
  stale <- .mkp_env$active_temp_logs
  if (!is.null(stale)) {
    for (f in stale) tryCatch(unlink(f), error = function(e) NULL)
  }
  .mkp_env$active_temp_logs <- NULL
  .mkp_env$recovery <- NULL
}


# --- Run initialization ---

#' Initialize state for a single run
#'
#' Returns a fully R-serializable run state: `chains` holds per-chain
#' parameter lists in checkpoint-compatible format; no XPtrs are created here.
#' XPtrs are built inside [.RunMkPrimeSingleRun()] so the state can be sent
#' to `future` workers without serialisation errors.
#' @keywords internal
.InitRun <- function(tree, mkd, model, mcmc, moves) {
  nChains <- mcmc$nChains
  betas <- .BuildTemperatureLadder(nChains, mcmc$heat)

  # Build per-chain state as R lists (checkpoint-compatible format).
  chains <- vector("list", nChains)
  for (ch in seq_len(nChains)) {
    s <- .InitState(tree, mkd, model)
    chains[[ch]] <- list(
      edge           = s$tree$edge,
      rel_br_lengths = s$rel_br_lengths,
      tree_length    = s$tree_length,
      rate_loss      = s$rate_loss,
      rate_log_sd    = s$rate_log_sd,
      rate_neo       = s$rate_neo %||% 1.0,
      p              = s$p %||% 0.5,
      kPrime         = as.integer(s$kPrime),
      log_lik        = s$log_lik,
      log_prior      = s$log_prior,
      log_post       = s$log_post
    )
  }

  moveNames <- vapply(moves, `[[`, character(1), "name")
  chainAccept <- chainPropose <- chainTuning <- vector("list", nChains)
  chainTimeNs <- chainSliceExp <- vector("list", nChains)
  for (ch in seq_len(nChains)) {
    chainAccept[[ch]] <- integer(length(moves))
    chainPropose[[ch]] <- integer(length(moves))
    names(chainAccept[[ch]]) <- names(chainPropose[[ch]]) <- moveNames
    chainTimeNs[[ch]] <- numeric(length(moves))
    names(chainTimeNs[[ch]]) <- moveNames
    chainSliceExp[[ch]] <- numeric(length(moves))
    names(chainSliceExp[[ch]]) <- moveNames
    chainTuning[[ch]] <- mcmc$tuning
  }

  swapAccept <- swapPropose <- if (nChains > 1L) {
    integer(nChains - 1L)
  } else {
    integer(0)
  }

  # M-120: Per-chain rho estimates for 2D joint Bactrian (start at 0)
  chainRhos <- vector("list", nChains)
  for (ch in seq_len(nChains)) {
    chainRhos[[ch]] <- list(rho_tl_rls = 0.0, rho_tl_rl = 0.0)
  }

  list(
    chains        = chains,
    betas         = betas,
    chain_accept    = chainAccept,
    chain_propose   = chainPropose,
    chain_time_ns   = chainTimeNs,
    chain_slice_exp = chainSliceExp,
    chain_tuning    = chainTuning,
    chain_rhos      = chainRhos,
    swap_accept     = swapAccept,
    swap_propose    = swapPropose
  )
}


# --- Single-run MCMC engine ---

#' Run the complete MCMC batch loop for one independent run
#'
#' Self-contained: accepts only R-serializable inputs (no XPtrs), reconstructs
#' C++ state internally, runs the full `repeat` loop, and returns a serializable
#' updated run state.  Used by the sequential `for (run)` loop in
#' [RunMkPrime()] / [ResumeMkPrime()] and (Phase 10b) by `future` workers.
#'
#' @param initialState Run list from [.InitRun()] or a checkpoint, with
#'   `chains` in checkpoint-compatible format (no `chainStates`).
#' @param runIdx Integer index of this run (1-based), used for display only.
#' @param startIter First iteration to run (1L for fresh, >1 when resuming).
#' @keywords internal
.RunMkPrimeSingleRun <- function(mkd, model, mcmc, initialState, moves,
                                  tipLabels, runIdx, paramNames, nEdge,
                                  brColStart, logFilePath, cancelFile,
                                  checkpointFile, startIter = 1L,
                                  isStreaming = FALSE, convWindowSize = 0L,
                                  treeFile = NULL, shared = NULL,
                                  resumeMoveWeights = NULL) {
  nChains <- mcmc$nChains
  treeEvery <- as.integer(mcmc$treeThin / mcmc$thin)

  # --- Reconstruct C++ XPtrs from R-serializable chain lists ---
  mcmcData <- .InitMcmcData(mkd, model)
  set_branch_bins(mcmcData, mcmc$nBranchBins)
  r <- initialState
  r$chainStates <- vector("list", nChains)
  for (ch in seq_len(nChains)) {
    ch_r <- r$chains[[ch]]
    r$chainStates[[ch]] <- init_mcmc_state(
      ch_r$edge[, 1L], ch_r$edge[, 2L],
      ch_r$rel_br_lengths, ch_r$tree_length,
      ch_r$rate_loss, ch_r$rate_log_sd,
      ch_r$rate_neo %||% 1.0, ch_r$p %||% 0.5,
      as.integer(ch_r$kPrime),
      ch_r$log_lik, ch_r$log_prior,
      ch_r$beta_scale %||% 1.0
    )
  }
  for (ch in seq_len(nChains)) {
    fill_partition_cache(mcmcData, r$chainStates[[ch]])
    allocate_cl_workspace(mcmcData, r$chainStates[[ch]])
  }

  # --- Sample storage ---
  savedIdx     <- r$saved_idx %||% 0L
  treeSavedIdx <- r$tree_saved_idx %||% 0L
  nSavedPerRun <- if (is.finite(mcmc$nIter)) {
    as.integer((mcmc$nIter - mcmc$warmup) / mcmc$thin)
  } else {
    2000L
  }
  nTreePerRun <- as.integer(ceiling(nSavedPerRun / treeEvery))

  if (isStreaming) {
    # Clear stale streaming fields before merging fresh buffers (resume path).
    r$flush_idx <- NULL; r$flushed     <- NULL
    r$conv_head <- NULL; r$conv_filled <- NULL
    r$flush_buf <- NULL; r$flush_iter  <- NULL
    r$conv_window <- NULL
    bufs <- .InitStreamBuffers(length(paramNames), paramNames,
                               mcmc$bufferSize, convWindowSize)
    r <- c(r, bufs)
    r$saved_idx      <- savedIdx
    r$tree_saved_idx <- treeSavedIdx
    if (is.null(r$tree_samples)) r$tree_samples <- vector("list", 0L)
  } else {
    if (is.null(r$samples)) {
      # Fresh run
      r$samples <- matrix(NA_real_, nrow = nSavedPerRun,
                          ncol = length(paramNames),
                          dimnames = list(NULL, paramNames))
      r$tree_samples <- vector("list", nTreePerRun)
    } else if (is.finite(mcmc$nIter)) {
      # Resuming: extend matrix if needed
      currentRows <- nrow(r$samples)
      if (currentRows < nSavedPerRun) {
        extra <- matrix(NA_real_, nrow = nSavedPerRun - currentRows,
                        ncol = ncol(r$samples),
                        dimnames = list(NULL, colnames(r$samples)))
        r$samples <- rbind(r$samples, extra)
      }
      if (length(r$tree_samples) < nTreePerRun) {
        r$tree_samples <- c(r$tree_samples,
                            vector("list", nTreePerRun - length(r$tree_samples)))
      }
    }
    r$saved_idx      <- savedIdx
    r$tree_saved_idx <- treeSavedIdx
  }

  # M-126: Rho estimation buffer for 2D joint Bactrian moves.
  # During warmup, C++ saves no samples (nSaved=0), so we accumulate

  # cold-chain state snapshots after each batch instead.
  rhoSampleBuf <- NULL

  # --- Batch loop constants ---
  # Adaptive batch size (M-106): fewer iterations per batch during warmup
  # (for tuning adaptation responsiveness), larger batches during sampling
  # to minimize R<->C++ round-trip overhead (90%+ of CPU at batchSize=200).
  warmupBatch   <- 500L
  tuningBatch   <- 500L
  samplingBatch <- 5000L
  moveNames     <- vapply(moves, `[[`, character(1), "name")
  moveWeights   <- vapply(moves, `[[`, numeric(1), "weight")
  moveWeights   <- moveWeights / sum(moveWeights)
  names(moveWeights) <- moveNames

  # M-149 #2: restore adapted weights from checkpoint / previous run
  if (!is.null(resumeMoveWeights)) {
    common <- intersect(names(resumeMoveWeights), moveNames)
    if (length(common) > 0L) {
      moveWeights[common] <- resumeMoveWeights[common]
      moveWeights <- moveWeights / sum(moveWeights)
    }
  }

  moveDim       <- vapply(moves, function(m) m$dim %||% 1L, integer(1L))
  names(moveDim) <- moveNames
  moveTypeCodes <- vapply(moves, function(m) .kMoveTypes[[m$name]], integer(1L))
  # Slice param index: 0=treeLength, 1=rateLoss, 2=rateLogSd, 3=rateNeo, 4=betaScale
  sliceParamCodes <- vapply(moves, function(m) m$sliceParamIdx %||% 0L, integer(1L))
  # Per-move integer parameter (e.g. nCats for dirichlet_branch; 0 = use chain default)
  moveIntParams <- vapply(moves, function(m) {
    if (!is.null(m$nCats)) as.integer(m$nCats) else 0L
  }, integer(1L))
  transIdx      <- which(mkd$type == "transformational")
  transIdx0     <- if (length(transIdx) > 0L) transIdx - 1L else integer(0L)
  hasNeo        <- any(mkd$type == "neomorphic")

  # Auto-pin always-accept moves (Gibbs, slice) at initial weights.
  # The warmup scheduler's score (accept_rate × dim / cost) gives these
  # astronomical scores because acceptance = 1.0 and cost ≈ 0; this inflates
  # their weight and starves bottleneck MH moves.  One Gibbs draw or slice
  # sample per cycle is already optimal, so freeze them.
  # Also pin hyperparameter moves (kprime_alpha/beta): they are cheap but

  # their score gets inflated relative to expensive topology moves.
  alwaysAcceptTypes <- c("gibbs_p", "slice", "gibbs_kprime_sweep",
                         "kprime_alpha", "kprime_beta")
  moveTypes <- vapply(moves, `[[`, character(1), "type")
  autoPin <- moveWeights[moveTypes %in% alwaysAcceptTypes]

  # Merge with user-specified pins (user takes precedence)
  userPins <- .ResolvePinnedWeights(mcmc$moveWeights, moveNames)
  if (length(autoPin) > 0L || !is.null(userPins)) {
    allPins <- autoPin
    if (!is.null(userPins)) allPins[names(userPins)] <- userPins
    pinnedWeights <- allPins
    moveWeights <- .NormalizeMoveWeights(moveWeights, pinnedWeights)
  } else {
    pinnedWeights <- NULL
  }

  # Ensure chain_time_ns exists (may be absent in older checkpoints)
  if (is.null(r$chain_time_ns)) {
    r$chain_time_ns <- vector("list", nChains)
    for (ch in seq_len(nChains)) {
      r$chain_time_ns[[ch]] <- numeric(length(moves))
      names(r$chain_time_ns[[ch]]) <- moveNames
    }
  }

  startTime     <- proc.time()["elapsed"]
  hasProgressFn <- !is.null(mcmc$progressFn) && !is.null(mcmc$plotEvery) &&
    mcmc$plotEvery > 0L

  # --- Three-phase state machine ---
  # Determine initial phase from checkpoint or fresh start.
  # Phases: "Warmup" → "Tuning" → "Sample"
  phase <- r$phase %||% (if (startIter <= mcmc$warmup) "Warmup" else "Sample")
  # M-141: wall-clock start of sample phase (for ETA estimation)
  sampleWallStart <- if (phase == "Sample") startTime else NA_real_
  etaStr          <- NULL

  # Stabilisation detector state (Warmup phase)
  logPostHistory     <- r$logPostHistory %||% numeric(0)
  nStableConsecutive <- r$nStableConsecutive %||% 0L

  # samplePhaseStart: iteration at which sampling began (for iterNum in log files).
  # For legacy checkpoints without this field, default to mcmc$warmup.
  if (is.null(r$samplePhaseStart)) {
    r$samplePhaseStart <- mcmc$warmup
  }

  # Tuning phase state
  tuningIterUsed   <- r$tuningIterUsed %||% 0L
  tuningRoundsDone <- r$tuningRoundsDone %||% 0L
  tuningBuf        <- NULL
  tuningBufIdx     <- 0L
  tuningTreeBuf    <- list()
  tuneWithTreeEss  <- !is.null(mcmc$minTreeEss) &&
    requireNamespace("TreeDist", quietly = TRUE)
  tuningWindowStart <- NULL
  bestMinEssPerSec <- -Inf
  bestWeights      <- moveWeights
  tuningCandidates <- list()
  tuningCandIdx    <- 0L
  effectiveTuningBudget <- r$effectiveTuningBudget %||% mcmc$tuningBudget

  # M-149: Re-allocate tuning infrastructure on Tuning-phase resume.
  # The tuningBuf is normally allocated at the Warmup→Tuning transition,
  # but that code doesn't re-run on resume.  Without this, the first saved
  # sample during Tuning hits nrow(NULL) → crash.
  if (phase == "Tuning") {
    tuningBufSize <- as.integer(effectiveTuningBudget / mcmc$thin) + 100L
    tuningBuf <- matrix(NA_real_, nrow = tuningBufSize,
                        ncol = length(paramNames),
                        dimnames = list(NULL, paramNames))
    tuningWindowStart <- proc.time()["elapsed"]
    tuningCandidates <- .PerturbMoveWeights(
      moveWeights, pinnedWeights, moveNames, nPerturbations = 3L
    )
  }

  # The C++ warmup parameter controls when samples are saved.
  # During Warmup: set to Inf so no samples saved.
  # During Tuning: set to 0 so all samples saved (collected into tuning buffer).
  # During Sample: set to 0 so all samples saved (collected into main storage).
  cppWarmup <- if (phase == "Warmup") mcmc$warmup else 0L

  weightsLogged <- phase == "Sample"
  thinAdapted   <- FALSE

  # --- Progress bar (M-097 rotating ticker) ---
  coldLogpost    <- {s <- get_mcmc_state(r$chainStates[[1]]); s$logPost}
  recentAcc      <- 0
  batchEnd       <- startIter - 1L
  phaseLabel     <- phase
  progressLabel  <- if (startIter == 1L) "MCMC" else "Resuming MCMC"
  progressTotal  <- if (is.finite(mcmc$nIter)) {
    if (startIter == 1L) mcmc$nIter else mcmc$nIter - startIter + 1L
  } else NA

  tickerPage  <- ""
  tickerPages <- if (phase == "Tuning") {
    "minESS/s: ?"
  } else {
    "minESS: ?"
  }
  logPWidth   <- 5L

  cli::cli_progress_bar(
    progressLabel,
    total  = progressTotal,
    format = "{tickerPage}",
    format_done = "{tickerPage}",
    clear  = FALSE
  )

  stopReason <- "max_iter"
  actualIter <- if (is.finite(mcmc$nIter)) mcmc$nIter else startIter - 1L

  # --- Main batch loop ---
  batchStart <- startIter
  repeat {
    # Batch size depends on phase
    batchSize <- switch(phase,
      Warmup  = warmupBatch,
      Tuning  = tuningBatch,
      Sample  = samplingBatch
    )
    batchEnd <- batchStart + batchSize - 1L
    if (is.finite(mcmc$nIter)) batchEnd <- min(batchEnd, mcmc$nIter)
    # Don't straddle the warmup boundary: end at maxWarmup so adaptation fires
    if (phase == "Warmup" && batchEnd > mcmc$warmup)
      batchEnd <- mcmc$warmup
    nBatch   <- batchEnd - batchStart + 1L

    scaleTunings <- .BuildScaleTuningMatrix(r$chain_tuning, moves)
    sliceWidths  <- .BuildSliceWidthMatrix(r$chain_tuning, moves)
    jointRhos    <- .BuildJointRhoMatrix(r$chain_rhos, moves, nChains)
    bsTunings    <- vapply(r$chain_tuning,
                           function(t) t$beta_simplex, numeric(1L))
    iwWins       <- vapply(r$chain_tuning,
                           function(t) as.integer(t$int_walk_window), integer(1L))

    result <- run_mcmc_batch_cpp(
      mcmcData, r$chainStates, r$betas,
      moveTypeCodes, transIdx0, sliceParamCodes, moveWeights,
      scaleTunings, bsTunings, iwWins, moveIntParams,
      sliceWidths, jointRhos,
      nBatch, batchStart, cppWarmup, mcmc$thin,
      hasNeo, nEdge,
      mcmc$cacheBonus
    )

    # Accept/propose counts, timing (M-092), and slice expansion counts
    for (ch in seq_len(nChains)) {
      r$chain_accept[[ch]]  <- r$chain_accept[[ch]]  +
        as.integer(result$accept_counts[ch, ])
      r$chain_propose[[ch]] <- r$chain_propose[[ch]] +
        as.integer(result$propose_counts[ch, ])
      r$chain_time_ns[[ch]] <- r$chain_time_ns[[ch]] +
        as.numeric(result$move_time_ns[ch, ])
      r$chain_slice_exp[[ch]] <- r$chain_slice_exp[[ch]] +
        as.numeric(result$slice_expansions[ch, ])
    }
    if (nChains > 1L) {
      r$swap_accept  <- r$swap_accept  + result$swap_accept
      r$swap_propose <- r$swap_propose + result$swap_propose
    }

    # --- Sample handling depends on phase ---
    nSaved <- result$n_saved
    if (phase == "Sample" && nSaved > 0L) {
      # Save to main storage (posterior samples)
      for (i in seq_len(nSaved)) {
        r$saved_idx <- r$saved_idx + 1L
        row <- result$scalar_samples[i, ]

        if (isStreaming) {
          iterNum <- r$samplePhaseStart + r$saved_idx * mcmc$thin
          r <- .AddToStreamBuffer(r, row, iterNum, logFilePath,
                                  mcmc$bufferSize, convWindowSize)
        } else {
          if (r$saved_idx > nrow(r$samples)) {
            n <- nrow(r$samples)
            extra <- matrix(NA_real_, nrow = n, ncol = ncol(r$samples),
                            dimnames = list(NULL, colnames(r$samples)))
            r$samples <- rbind(r$samples, extra)
          }
          r$samples[r$saved_idx, ] <- row
        }

        # Tree storage: only on treeThin boundary
        if (r$saved_idx %% treeEvery == 0L) {
          r$tree_saved_idx <- r$tree_saved_idx + 1L
          if (r$tree_saved_idx > length(r$tree_samples)) {
            n <- max(length(r$tree_samples), 1L)
            r$tree_samples <- c(r$tree_samples, vector("list", n))
          }
          tl    <- row[3L]
          relBr <- row[brColStart:(brColStart + nEdge - 1L)]
          curTree <- TreeTools::Preorder(structure(
            list(edge        = result$edge_samples[[i]],
                 edge.length = tl * relBr,
                 Nnode       = length(tipLabels) - 1L,
                 tip.label   = tipLabels),
            class = "phylo"
          ))
          r$tree_samples[[r$tree_saved_idx]] <- curTree
          if (!is.null(treeFile))
            cat(ape::write.tree(curTree), "\n", file = treeFile, append = TRUE)
        }
      }
    } else if (phase == "Tuning" && nSaved > 0L) {
      # Collect into tuning buffer (discarded after tuning)
      for (i in seq_len(nSaved)) {
        tuningBufIdx <- tuningBufIdx + 1L
        if (tuningBufIdx > nrow(tuningBuf)) {
          extra <- matrix(NA_real_, nrow = nrow(tuningBuf),
                          ncol = ncol(tuningBuf),
                          dimnames = list(NULL, colnames(tuningBuf)))
          tuningBuf <- rbind(tuningBuf, extra)
        }
        tuningBuf[tuningBufIdx, ] <- result$scalar_samples[i, ]

        # Store trees during tuning for tree-ESS-aware bandit
        if (tuneWithTreeEss) {
          row <- result$scalar_samples[i, ]
          tl    <- row[3L]
          relBr <- row[brColStart:(brColStart + nEdge - 1L)]
          tuningTreeBuf[[tuningBufIdx]] <- TreeTools::Preorder(structure(
            list(edge        = result$edge_samples[[i]],
                 edge.length = tl * relBr,
                 Nnode       = length(tipLabels) - 1L,
                 tip.label   = tipLabels),
            class = "phylo"
          ))
        }
      }
    }

    # ===== PHASE-SPECIFIC LOGIC =====

    # Cold-chain state (reused for rho estimation + progress bar)
    s <- get_mcmc_state(r$chainStates[[1]])
    coldLogpost <- s$logPost

    if (phase == "Warmup") {
      # --- Warmup: adapt tuning, temperatures, and move weights ---
      for (ch in seq_len(nChains)) {
        r$chain_tuning[[ch]] <- .AdaptTuning(
          r$chain_tuning[[ch]], r$chain_accept[[ch]],
          r$chain_propose[[ch]], moves
        )
        r$chain_tuning[[ch]] <- .AdaptSliceWidths(
          r$chain_tuning[[ch]], r$chain_propose[[ch]],
          r$chain_slice_exp[[ch]], moves
        )
      }
      if (nChains > 1L)
        r$betas <- .AdaptTemperatures(r$betas, r$swap_accept, r$swap_propose)
      # Adaptive move weight scheduling (M-092): acceptance-rate heuristic
      moveWeights <- .AdaptMoveWeights(
        moveWeights, r$chain_accept[[1L]], r$chain_propose[[1L]],
        r$chain_time_ns[[1L]], moveNames, moveDim = moveDim,
        pinnedWeights = pinnedWeights,
        warmupProgress = min(1, batchEnd / mcmc$warmup)
      )

      # M-126: Accumulate cold-chain state snapshots for rho estimation.
      # C++ saves no samples during warmup, so we use the chain state
      # after each batch (already queried for the progress bar above).
      rhoSampleBuf <- .AccumulateRhoSnapshot(
        rhoSampleBuf, s, hasNeo, paramNames
      )
      newRhos <- .EstimateJointRhos(rhoSampleBuf, hasNeo)
      for (ch in seq_len(nChains)) {
        r$chain_rhos[[ch]] <- newRhos
      }

      # Stabilisation detection
      logPostHistory <- c(logPostHistory, coldLogpost)
      r$logPostHistory <- logPostHistory
      if (batchEnd >= mcmc$minWarmup) {
        stabResult <- .CheckStabilisation(
          logPostHistory, nStableConsecutive
        )
        nStableConsecutive <- stabResult$nStableConsecutive
        r$nStableConsecutive <- nStableConsecutive

        if (stabResult$stable || batchEnd >= mcmc$warmup) {
          # Transition: Warmup → Tuning (or Sample if autoTune = FALSE)
          if (batchEnd >= mcmc$warmup && !stabResult$stable) {
            cli::cli_warn(
              "Warmup reached {.arg maxWarmup} ({mcmc$warmup}) without stabilisation."
            )
          } else {
            cli::cli_alert_success(
              "Chain stabilised at iteration {batchEnd}."
            )
          }

          # Check if there's enough remaining budget for tuning + sampling
          remainingIter <- if (is.finite(mcmc$nIter)) {
            mcmc$nIter - batchEnd
          } else {
            Inf
          }
          # Need at least tuningBatch * 2 for tuning + some for sampling
          canTune <- mcmc$autoTune && remainingIter > tuningBatch * 4L

          if (canTune) {
            phase      <- "Tuning"
            r$phase    <- phase
            phaseLabel <- "Tuning"
            cppWarmup  <- 0L
            # Scale tuning budget: ensure each candidate window has enough
            # samples for meaningful ESS. 4 windows per round (current + 3
            # perturbations), ~100 samples each.
            nCandidates <- 4L
            minSamplesPerWindow <- 100L
            scaledBudget <- as.integer(
              mcmc$thin * minSamplesPerWindow * nCandidates
            )
            baseBudget <- max(mcmc$tuningBudget, scaledBudget)
            effectiveTuningBudget <- if (is.finite(remainingIter)) {
              min(baseBudget, as.integer(remainingIter / 2))
            } else {
              baseBudget
            }
            r$effectiveTuningBudget <- effectiveTuningBudget
            # Allocate tuning buffer
            tuningBufSize <- as.integer(effectiveTuningBudget / mcmc$thin) + 100L
            tuningBuf <- matrix(NA_real_, nrow = tuningBufSize,
                                ncol = length(paramNames),
                                dimnames = list(NULL, paramNames))
            tuningBufIdx     <- 0L
            tuningTreeBuf    <- list()
            tuningWindowStart <- proc.time()["elapsed"]
            # Reset acceptance/timing counters for clean tuning measurement
            for (ch in seq_len(nChains)) {
              r$chain_accept[[ch]][]    <- 0L
              r$chain_propose[[ch]][]   <- 0L
              r$chain_time_ns[[ch]][]   <- 0
              r$chain_slice_exp[[ch]][] <- 0
            }
            bestMinEssPerSec <- -Inf
            bestWeights      <- moveWeights
            tuningCandidates <- .PerturbMoveWeights(
              moveWeights, pinnedWeights, moveNames,
              nPerturbations = 3L
            )
            tuningCandIdx    <- 0L
            tickerPages      <- "minESS/s: ?"
          } else {
            # Skip tuning, go straight to Sample
            phase      <- "Sample"
            r$phase    <- phase
            phaseLabel <- "Sample"
            cppWarmup  <- 0L
            r$samplePhaseStart <- batchEnd
            sampleWallStart <- proc.time()["elapsed"]
            if (isStreaming && !is.null(logFilePath))
              .LogMoveWeights(moveWeights, moveNames, logFilePath)
            .PrintMoveWeights(moveWeights, moveNames)
            weightsLogged <- TRUE
          }
        }
      }
    } else if (phase == "Tuning") {
      # --- Tuning: min-ESS/s perturbation bandit ---
      tuningIterUsed <- tuningIterUsed + nBatch
      r$tuningIterUsed <- tuningIterUsed

      # M-126: Continue rho estimation from tuning buffer samples
      if (tuningBufIdx >= 50L) {
        newRhos <- .EstimateJointRhos(
          tuningBuf[seq_len(tuningBufIdx), , drop = FALSE], hasNeo
        )
        for (ch in seq_len(nChains))
          r$chain_rhos[[ch]] <- newRhos
      }

      # Evaluate current weight vector after each tuning window
      if (tuningBufIdx >= 10L) {
        windowTime <- proc.time()["elapsed"] - tuningWindowStart
        currentEssPerSec <- .MinEssPerSec(
          tuningBuf[seq_len(tuningBufIdx), , drop = FALSE],
          windowTime,
          tuningTrees = if (tuneWithTreeEss && tuningBufIdx >= 20L) {
            tuningTreeBuf[seq_len(tuningBufIdx)]
          }
        )

        if (!is.na(currentEssPerSec)) {
          tickerPages <- sprintf("minESS/s: %.2f", currentEssPerSec)
          if (currentEssPerSec > bestMinEssPerSec) {
            bestMinEssPerSec <- currentEssPerSec
            bestWeights      <- moveWeights
          }
        }

        # Move to next candidate or next round
        tuningCandIdx <- tuningCandIdx + 1L
        if (tuningCandIdx <= length(tuningCandidates)) {
          # Try next perturbation candidate
          moveWeights <- tuningCandidates[[tuningCandIdx]]
          tuningBufIdx      <- 0L
          tuningTreeBuf     <- list()
          tuningWindowStart <- proc.time()["elapsed"]
          # Reset counters for clean measurement
          for (ch in seq_len(nChains)) {
            r$chain_accept[[ch]][]    <- 0L
            r$chain_propose[[ch]][]   <- 0L
            r$chain_time_ns[[ch]][]   <- 0
            r$chain_slice_exp[[ch]][] <- 0
          }
        } else {
          # End of round: adopt best weights, start new round
          tuningRoundsDone <- tuningRoundsDone + 1L
          r$tuningRoundsDone <- tuningRoundsDone
          moveWeights <- bestWeights

          if (tuningRoundsDone >= mcmc$tuningRounds ||
              tuningIterUsed >= effectiveTuningBudget) {
            # Transition: Tuning → Sample
            phase      <- "Sample"
            r$phase    <- phase
            phaseLabel <- "Sample"
            r$samplePhaseStart <- batchEnd
            sampleWallStart <- proc.time()["elapsed"]
            if (isStreaming && !is.null(logFilePath))
              .LogMoveWeights(moveWeights, moveNames, logFilePath)
            .PrintMoveWeights(moveWeights, moveNames)
            if (bestMinEssPerSec > 0) {
              cli::cli_alert_info(
                "Tuning complete ({tuningRoundsDone} round{?s}). Best minESS/s: {sprintf('%.2f', bestMinEssPerSec)}"
              )
            }
            weightsLogged <- TRUE
            tickerPages   <- "minESS: ?"
          } else {
            # Start new round with fresh perturbations
            tuningCandidates <- .PerturbMoveWeights(
              moveWeights, pinnedWeights, moveNames,
              nPerturbations = 3L
            )
            tuningCandIdx     <- 0L
            tuningBufIdx      <- 0L
            tuningTreeBuf     <- list()
            tuningWindowStart <- proc.time()["elapsed"]
            bestMinEssPerSec  <- -Inf
            for (ch in seq_len(nChains)) {
              r$chain_accept[[ch]][]    <- 0L
              r$chain_propose[[ch]][]   <- 0L
              r$chain_time_ns[[ch]][]   <- 0
              r$chain_slice_exp[[ch]][] <- 0
            }
          }
        }
      }
    }
    # Sample phase: no adaptation needed (weights frozen)

    # Streaming checkpoint: fire when buffer was flushed this batch (Sample)
    if (phase == "Sample" && isStreaming &&
        !is.null(checkpointFile) && isTRUE(r$flushed)) {
      if (r$flush_idx > 0L) {
        .FlushBuffer(r$flush_buf, r$flush_idx, r$flush_iter, logFilePath)
        r$flush_idx <- 0L
      }
      r$flushed <- FALSE
      .SaveCheckpoint(list(r), mcmc, batchEnd, paramNames, checkpointFile,
                      moveWeights = moveWeights, phase = phase,
                      model = model)
    }

    # Warmup/tuning checkpoint at checkEvery intervals so long warmups
    # are recoverable.  No samples to flush — just save chain state.
    if (phase != "Sample" && !is.null(checkpointFile) &&
        !is.null(mcmc$checkEvery) && mcmc$checkEvery > 0L &&
        (batchEnd %/% mcmc$checkEvery) >
          ((batchStart - 1L) %/% mcmc$checkEvery)) {
      .SaveCheckpoint(list(r), mcmc, batchEnd, paramNames, checkpointFile,
                      moveWeights = moveWeights, phase = phase,
                      model = model)
    }

    # Persist move weights in run state so checkpoints capture them (M-149 #2).
    r$moveWeights <- moveWeights

    # Update shared state for interrupt-safe checkpointing (M-149).
    # The interrupt handler in .RunWithRecovery() reads from this env.
    if (!is.null(shared)) {
      shared$runs[[runIdx]] <- r
      shared$actualIter     <- max(shared$actualIter, batchEnd)
      shared$phase          <- phase
      shared$moveWeights    <- moveWeights
    }

    # Progress update (M-097 rotating ticker)
    batchAcc  <- sum(result$accept_counts[1L, ])
    batchProp <- sum(result$propose_counts[1L, ])
    if (batchProp > 0L) recentAcc <- batchAcc / batchProp

    # Fixed-width logP: ratchet width up as magnitude grows, never shrink
    logPRaw   <- format(round(coldLogpost, 1), nsmall = 1)
    logPWidth <- max(logPWidth, nchar(logPRaw))
    logPStr   <- formatC(round(coldLogpost, 1), width = logPWidth,
                         format = "f", digits = 1)
    # Dim separators; phase prefix silver for visual separation
    sep <- cli::col_silver("\u2502")
    tickerPage <- paste(
      cli::col_silver(paste(phaseLabel, batchEnd)),
      sep, sprintf("logP:%s", logPStr),
      sep, tickerPages
    )

    cli::cli_progress_update(
      set = if (startIter == 1L) batchEnd else batchEnd - startIter + 1L
    )
    if (hasProgressFn &&
        (batchEnd %/% mcmc$plotEvery) > ((batchStart - 1L) %/% mcmc$plotEvery)) {
      info <- .BuildProgressInfo(list(r), batchEnd, mcmc, startTime,
                                 recentAcc, paramNames, phase = phase)
      mcmc$progressFn(info)
    }

    # Stopping: max wall-clock time (M-150: flush + checkpoint before break)
    if (!is.null(mcmc$maxTime) &&
        proc.time()["elapsed"] - startTime >= mcmc$maxTime) {
      if (isStreaming && r$flush_idx > 0L) {
        .FlushBuffer(r$flush_buf, r$flush_idx, r$flush_iter, logFilePath)
        r$flush_idx <- 0L
        r$flushed   <- FALSE
      }
      if (!is.null(checkpointFile)) {
        .SaveCheckpoint(list(r), mcmc, batchEnd, paramNames, checkpointFile,
                        moveWeights = moveWeights, phase = phase,
                        model = model)
      }
      stopReason <- "max_time"
      actualIter <- batchEnd
      break
    }

    # Stopping: cancel file
    if (!is.null(cancelFile) && file.exists(cancelFile)) {
      if (isStreaming && r$flush_idx > 0L) {
        .FlushBuffer(r$flush_buf, r$flush_idx, r$flush_iter, logFilePath)
        r$flush_idx <- 0L
        r$flushed   <- FALSE
      }
      if (!is.null(checkpointFile)) {
        .SaveCheckpoint(list(r), mcmc, batchEnd, paramNames, checkpointFile,
                        moveWeights = moveWeights, phase = phase,
                        model = model)
      }
      stopReason <- "cancelled"
      actualIter <- batchEnd
      break
    }

    # Convergence check + checkpoint at checkEvery intervals (Sample phase only)
    doCheck <- phase == "Sample" && !is.null(mcmc$checkEvery) &&
      mcmc$checkEvery > 0L &&
      (batchEnd %/% mcmc$checkEvery) > ((batchStart - 1L) %/% mcmc$checkEvery)

    if (doCheck) {
      if (!is.null(checkpointFile)) {
        if (isStreaming && r$flush_idx > 0L) {
          .FlushBuffer(r$flush_buf, r$flush_idx, r$flush_iter, logFilePath)
          r$flush_idx <- 0L
          r$flushed   <- FALSE
        }
        .SaveCheckpoint(list(r), mcmc, batchEnd, paramNames, checkpointFile,
                        moveWeights = moveWeights, phase = phase,
                        model = model)
      }

      diagCheck <- .CheckConvergence(list(r), paramNames, mcmc, isStreaming)
      if (!is.null(diagCheck)) {
        # M-141: ETA from worst-case ESS accumulation rate.
        # Use whichever criterion (scalar ESS or tree ESS) has the
        # worst current/target ratio — that's the binding constraint.
        elapsedSample <- proc.time()["elapsed"] - sampleWallStart
        etaCurrent <- diagCheck$minEss
        etaTarget  <- mcmc$minEss
        if (!is.null(mcmc$minTreeEss) && !is.na(diagCheck$treeEss) &&
            is.finite(diagCheck$treeEss) && !is.null(mcmc$minEss) &&
            is.finite(diagCheck$minEss)) {
          scalarRatio <- diagCheck$minEss / mcmc$minEss
          treeRatio   <- diagCheck$treeEss / mcmc$minTreeEss
          if (treeRatio < scalarRatio) {
            etaCurrent <- diagCheck$treeEss
            etaTarget  <- mcmc$minTreeEss
          }
        }
        etaStr <- .EstimateEta(etaCurrent, etaTarget, elapsedSample)
        # Refresh ticker pages from latest diagnostics (M-097)
        tickerPages <- .BuildTickerPages(diagCheck, etaStr)
        if (diagCheck$converged) {
          stopReason <- "converged"
          actualIter <- batchEnd
          break
        }
      }

      # M-135: adapt thin from observed ACT (fires once, first check only)
      if (!thinAdapted && isTRUE(mcmc$thinWasAuto)) {
        thinAdapted <- TRUE
        mat <- if (isStreaming) {
          .ConvWindowRows(r, minRows = 50L)
        } else if (r$saved_idx >= 50L) {
          r$samples[seq_len(r$saved_idx), , drop = FALSE]
        } else {
          NULL
        }
        if (!is.null(mat) && nrow(mat) >= 50L) {
          newThin <- .AdaptThinning(mat, mcmc$thin, length(moves))
          if (newThin != mcmc$thin) {
            oldThin <- mcmc$thin
            mcmc$thin <- newThin
            if (isTRUE(mcmc$treeThinWasAuto)) {
              mcmc$treeThin <- newThin
              treeEvery <- 1L
            } else {
              mcmc$treeThin <- max(mcmc$treeThin, newThin)
              if (mcmc$treeThin %% newThin != 0L)
                mcmc$treeThin <- newThin * ceiling(mcmc$treeThin / newThin)
              treeEvery <- as.integer(mcmc$treeThin / newThin)
            }
            cli::cli_alert_info(
              "Adapted thin: {oldThin} \u2192 {newThin} (max ACT \u2248 {round(newThin / log(2))} iter)"
            )
          }
        }
      }
    }

    actualIter <- batchEnd
    batchStart <- batchEnd + 1L
    if (is.finite(mcmc$nIter) && batchStart > mcmc$nIter) break
  }
  tickerPage <- paste(
    cli::col_silver(paste(phaseLabel, batchEnd)),
    cli::col_silver("\u2502"), cli::col_green("done \u2714")
  )
  cli::cli_progress_done()

  # Flush any remaining streaming buffer (belt-and-suspenders; .BuildResult

  # also flushes, but doing it here means the log file is up-to-date
  # immediately on return)
  if (isStreaming && r$flush_idx > 0L) {
    .FlushBuffer(r$flush_buf, r$flush_idx, r$flush_iter, logFilePath)
    r$flush_idx <- 0L
  }

  # --- Serialize and return ---
  r$chains <- lapply(r$chainStates, function(ptr) {
    s <- get_mcmc_state(ptr)
    list(
      log_lik        = s$logLik,
      log_prior      = s$logPrior,
      log_post       = s$logPost,
      tree_length    = s$treeLength,
      rel_br_lengths = s$relBrLengths,
      rate_loss      = s$rateLoss,
      rate_log_sd    = s$rateLogSd,
      rate_neo       = s$rateNeo,
      p              = s$p,
      kPrime         = s$kPrime,
      edge           = s$edge
    )
  })
  r$chainStates  <- NULL
  r$stop_reason  <- stopReason
  r$actual_iter  <- actualIter
  r$phase        <- phase
  r$moveWeights  <- moveWeights   # M-149 #2: persist for resume / Phase 2
  r
}


# --- Parallel run orchestration ---

# --- Serial multi-run orchestration (M-146) ---

#' Run multiple serial MCMC runs with cross-run R-hat convergence
#'
#' Called by [.RunWithRecovery()] when `parallel = FALSE`, `nRuns >= 2`, and
#' `maxRhat` is set.  Phase 1 runs each run sequentially until per-run ESS
#' convergence (or nIter / maxTime / cancel).  Phase 2 checks cross-run R-hat
#' from log files; if not met and iteration headroom remains, resumes each run
#' for one epoch (`checkEvery` iterations) and re-checks.
#'
#' @return Named list: `runs`, `stopReason`, `actualIter`.
#' @keywords internal
.RunSerialRuns <- function(mkd, model, mcmc, runs, moves, tipLabels,
                            paramNames, nEdge, brColStart, logFilePaths,
                            convWindowSize, treeFile, startIters = NULL,
                            shared = NULL, startPhase = 1L) {
  nRuns     <- length(runs)
  epochSize <- max(mcmc$checkEvery %||% 1000L, 1000L)
  startTime <- proc.time()["elapsed"]

  if (is.null(startIters)) startIters <- rep(1L, nRuns)

  # --- Phase 1: first pass (ESS-based stopping per run) ---
  if (startPhase <= 1L) {
    # Strip maxRhat so it doesn't block ESS-based stopping inside each run.
    innerMcmc <- mcmc
    innerMcmc$maxRhat <- NULL

    for (run in seq_len(nRuns)) {
      runs[[run]] <- .RunMkPrimeSingleRun(
        mkd, model, innerMcmc, runs[[run]], moves, tipLabels, run,
        paramNames, nEdge, brColStart,
        logFilePath    = logFilePaths[run],
        cancelFile     = mcmc$cancelFile,
        checkpointFile = NULL,
        startIter      = startIters[run],
        isStreaming     = TRUE,
        convWindowSize = convWindowSize,
        treeFile       = treeFile,
        shared         = shared,
        resumeMoveWeights = runs[[run]]$moveWeights
      )
      if (runs[[run]]$stop_reason == "cancelled") {
        return(list(runs = runs,
                    stopReason = "cancelled",
                    actualIter = runs[[run]]$actual_iter))
      }
      startIters[run] <- runs[[run]]$actual_iter + 1L
    }

    # No cross-run convergence needed when nRuns < 2 or no maxRhat
    # (defensive — caller should not route here in those cases).
    if (nRuns < 2L || is.null(mcmc$maxRhat)) {
      return(list(runs = runs,
                  stopReason = runs[[nRuns]]$stop_reason,
                  actualIter = runs[[nRuns]]$actual_iter))
    }

    # Checkpoint after first pass (M-149 #6: record serialPhase)
    if (!is.null(mcmc$checkpointFile)) {
      maxActual <- max(vapply(runs, `[[`, 0, "actual_iter"))
      .SaveCheckpoint(runs, mcmc, maxActual, paramNames, mcmc$checkpointFile,
                      moveWeights = runs[[nRuns]]$moveWeights,
                      model = model, serialPhase = 2L)
    }
  }

  # --- Phase 2: cross-run R-hat loop ---
  repeat {
    diagCheck <- .CheckConvergenceFromLogs(logFilePaths, paramNames, mcmc)
    if (!is.null(diagCheck) && diagCheck$converged) {
      return(list(runs = runs,
                  stopReason = "converged",
                  actualIter = max(vapply(runs, `[[`, 0, "actual_iter"))))
    }

    # Report R-hat status
    if (!is.null(diagCheck) && !is.na(diagCheck$maxRhat)) {
      cli::cli_alert_info(
        "Cross-run max R-hat = {round(diagCheck$maxRhat, 3)} \\
         (target: {mcmc$maxRhat}). Extending runs\u2026"
      )
    }

    # Check hard limits before extending
    maxActual <- max(vapply(runs, `[[`, 0, "actual_iter"))
    if (is.finite(mcmc$nIter) && maxActual >= mcmc$nIter) {
      return(list(runs = runs,
                  stopReason = "max_iter",
                  actualIter = maxActual))
    }
    elapsed <- proc.time()["elapsed"] - startTime
    if (!is.null(mcmc$maxTime) && elapsed >= mcmc$maxTime) {
      # M-150: save checkpoint before returning
      if (!is.null(mcmc$checkpointFile)) {
        .SaveCheckpoint(runs, mcmc, maxActual, paramNames,
                        mcmc$checkpointFile,
                        moveWeights = runs[[nRuns]]$moveWeights,
                        model = model, serialPhase = 2L)
      }
      return(list(runs = runs,
                  stopReason = "max_time",
                  actualIter = maxActual))
    }

    # Resumption epoch: strip convergence criteria, use finite nIter cap.
    epochMcmc <- mcmc
    epochMcmc$maxRhat <- NULL
    epochMcmc$minEss  <- NULL

    for (run in seq_len(nRuns)) {
      epochEnd <- startIters[run] + epochSize - 1L
      if (is.finite(mcmc$nIter)) epochEnd <- min(epochEnd, mcmc$nIter)
      epochMcmc$nIter <- epochEnd

      runs[[run]] <- .RunMkPrimeSingleRun(
        mkd, model, epochMcmc, runs[[run]], moves, tipLabels, run,
        paramNames, nEdge, brColStart,
        logFilePath    = logFilePaths[run],
        cancelFile     = mcmc$cancelFile,
        checkpointFile = NULL,
        startIter      = startIters[run],
        isStreaming     = TRUE,
        convWindowSize = convWindowSize,
        treeFile       = treeFile,
        shared         = shared,
        resumeMoveWeights = runs[[run]]$moveWeights
      )
      if (runs[[run]]$stop_reason == "cancelled") {
        return(list(runs = runs,
                    stopReason = "cancelled",
                    actualIter = runs[[run]]$actual_iter))
      }
      startIters[run] <- runs[[run]]$actual_iter + 1L
    }

    # Checkpoint after each epoch (M-149 #6: record serialPhase = 2)
    if (!is.null(mcmc$checkpointFile)) {
      maxActual <- max(vapply(runs, `[[`, 0, "actual_iter"))
      .SaveCheckpoint(runs, mcmc, maxActual, paramNames, mcmc$checkpointFile,
                      moveWeights = runs[[nRuns]]$moveWeights,
                      model = model, serialPhase = 2L)
    }
  }
}


# --- Parallel run orchestration ---

#' Launch and manage parallel independent MCMC runs via `future`
#'
#' Called by [RunMkPrime()] when `mcmc$parallel = TRUE` and `nRuns > 1`.
#' Spawns `nRuns` non-blocking futures each calling [.RunMkPrimeSingleRun()],
#' then polls for convergence / time limits / user cancel. When a stopping
#' criterion fires, writes per-run cancel files so workers exit cleanly.
#'
#' @return Named list: `runs`, `logFilePaths`, `stopReason`, `actualIter`.
#' @keywords internal
.RunParallelRuns <- function(mkd, model, mcmc, runs, moves, tipLabels,
                              paramNames, nEdge, brColStart, treeFile,
                              isStreaming, logFilePaths, convWindowSize) {
  nRuns <- mcmc$nRuns

  if (!requireNamespace("future", quietly = TRUE)) {
    cli::cli_abort(c(
      "Package {.pkg future} is required for parallel runs.",
      "i" = "Install it with: {.code install.packages(\"future\")}",
      "i" = "Then set a plan before calling RunMkPrime(): \\
             {.code future::plan(\"multisession\", workers = {nRuns})}"
    ))
  }

  # Parallel mode requires streaming so the orchestrator can read samples.
  # Since RunMkPrime() now always streams (temp log), this branch is a
  # safety net for any future direct callers.
  if (!isStreaming) {
    tmpLog <- tempfile(fileext = ".log")
    cli::cli_alert_info(c(
      "Parallel mode requires {.arg logFile} (workers share samples via disk).",
      "i" = "Auto-assigning: {.file {tmpLog}}"
    ))
    mcmc$logFile   <- tmpLog
    convWindowSize <- .ComputeConvWindowSize(mcmc)
    logFilePaths   <- .OpenLogFiles(mcmc$logFile, paramNames, nRuns)
  }

  # Per-run cancel files (orchestrator signals each worker individually).
  cancelFiles <- vapply(seq_len(nRuns), function(i) tempfile(), character(1L))

  # Launch nRuns persistent workers — each runs its full batch loop.
  fList <- vector("list", nRuns)
  for (run in seq_len(nRuns)) {
    runState <- runs[[run]]
    logPath  <- logFilePaths[run]
    cfPath   <- cancelFiles[run]
    fList[[run]] <- future::future(
      {
        .RunMkPrimeSingleRun(
          mkd, model, mcmc, runState, moves, tipLabels, run,
          paramNames, nEdge, brColStart,
          logFilePath    = logPath,
          cancelFile     = cfPath,
          checkpointFile = NULL,
          startIter      = 1L,
          isStreaming    = TRUE,
          convWindowSize = convWindowSize,
          treeFile       = NULL
        )
      },
      seed = TRUE
    )
  }

  # Polling loop: sleep → check stopping criteria → signal workers if needed.
  startTime    <- proc.time()["elapsed"]
  pollInterval <- mcmc$pollInterval %||% 10L
  stopReason   <- "max_iter"
  actualIter   <- if (is.finite(mcmc$nIter)) mcmc$nIter else mcmc$warmup

  # Progress display and live trace plot
  hasProgressFn <- !is.null(mcmc$progressFn) && is.function(mcmc$progressFn)
  pollStatus <- "Waiting for workers\u2026"
  cli::cli_progress_bar(
    "Parallel MCMC ({nRuns} runs)",
    format       = "{cli::pb_spin} {pollStatus}",
    format_done  = "{pollStatus}",
    clear        = FALSE
  )

  repeat {
    Sys.sleep(pollInterval)

    # User cancel file (shared across runs)
    if (!is.null(mcmc$cancelFile) && file.exists(mcmc$cancelFile)) {
      for (cf in cancelFiles) file.create(cf)
      stopReason <- "cancelled"
      break
    }

    # Wall-clock time limit
    elapsed <- proc.time()["elapsed"] - startTime
    if (!is.null(mcmc$maxTime) && elapsed >= mcmc$maxTime) {
      for (cf in cancelFiles) file.create(cf)
      stopReason <- "max_time"
      break
    }

    # Convergence (reads log files from disk)
    diagCheck <- .CheckConvergenceFromLogs(logFilePaths, paramNames, mcmc)
    if (!is.null(diagCheck)) {
      elStr  <- .FormatElapsed(elapsed)
      essStr <- round(diagCheck$minEss)
      etaStr <- .EstimateEta(diagCheck$minEss, mcmc$minEss, elapsed)
      pollStatus <- paste0(
        elStr, " | min ESS = ", essStr,
        if (!is.null(mcmc$minEss)) paste0(" / ", mcmc$minEss) else "",
        if (!is.na(diagCheck$maxRhat))
          paste0(" | max Rhat = ", round(diagCheck$maxRhat, 3),
                 if (!is.null(mcmc$maxRhat))
                   paste0(" / ", mcmc$maxRhat))
        else "",
        if (!is.null(etaStr)) paste0(" | ETA: ", etaStr) else ""
      )
      cli::cli_progress_update()

      # Live trace plot from log-file samples
      if (hasProgressFn) {
        nSamp <- nrow(diagCheck$perRunSamples[[1]])
        info <- list(
          iter           = nSamp * if (is.numeric(mcmc$thin)) mcmc$thin else 1L,
          nIter          = mcmc$nIter,
          warmup         = mcmc$warmup,
          inWarmup       = FALSE,
          phase          = "Sample",
          nRuns          = nRuns,
          nChains        = mcmc$nChains,
          runSamples     = diagCheck$perRunSamples,
          currentState   = NULL,
          recentAcceptance = NA_real_,
          elapsed        = elapsed,
          paramNames     = paramNames
        )
        tryCatch(mcmc$progressFn(info), error = function(e) NULL)
      }

      if (diagCheck$converged) {
        for (cf in cancelFiles) file.create(cf)
        stopReason <- "converged"
        break
      }
    }

    # All workers finished naturally
    if (all(vapply(fList, future::resolved, logical(1L)))) break
  }

  pollStatus <- paste0(
    "Parallel MCMC (", nRuns, " runs) \u2014 ",
    stopReason, " [", .FormatElapsed(proc.time()["elapsed"] - startTime), "]"
  )
  cli::cli_progress_done()

  # Collect results (blocks until each worker is done)
  completedRuns <- lapply(fList, future::value)

  # Take actualIter from the first completed run
  actualIter <- completedRuns[[1L]]$actual_iter %||% actualIter

  list(
    runs         = completedRuns,
    logFilePaths = logFilePaths,
    stopReason   = stopReason,
    actualIter   = actualIter
  )
}


# --- Convergence check during MCMC ---

#' Check convergence criteria (called during the loop)
#'
#' Works for any number of runs. ESS is always computed on combined samples;
#' R-hat is computed only when `nRuns >= 2`. Returns full per-parameter `ess`
#' and `rhat` vectors so the caller can display a progress table.
#' @keywords internal
.CheckConvergence <- function(runs, paramNames, mcmc, isStreaming = FALSE) {
  nRuns <- length(runs)
  keyCols <- .KeyParamCols(
    matrix(0, 1, length(paramNames), dimnames = list(NULL, paramNames))
  )

  # Gather saved samples from each run (convergence window in streaming mode)
  perRunSamples <- lapply(runs, function(r) {
    if (isStreaming) {
      rows <- .ConvWindowRows(r, minRows = 10L)
      if (is.null(rows)) return(NULL)
      rows[, keyCols, drop = FALSE]
    } else {
      idx <- r$saved_idx
      if (idx < 10L) return(NULL)
      r$samples[seq_len(idx), keyCols, drop = FALSE]
    }
  })

  if (any(vapply(perRunSamples, is.null, logical(1)))) return(NULL)

  # ESS on combined samples (works for any nRuns)
  combined <- do.call(rbind, perRunSamples)
  ess <- .EssMatrix(combined)

  # kPrime are discrete nuisance parameters — exclude from convergence criteria
  # (M-098). They remain in the `ess` vector for display in .PrintProgressTable.
  isConvParam <- !grepl("^kPrime_", names(ess)) & names(ess) != "log_likelihood"
  minEss <- min(ess[isConvParam], na.rm = TRUE)

  # R-hat (requires >= 2 runs)
  rhat    <- NULL
  maxRhat <- NA_real_
  if (nRuns >= 2L) {
    paramNms <- colnames(perRunSamples[[1]])
    rhat <- vapply(seq_along(paramNms), function(j) {
      chainMat <- do.call(cbind, lapply(perRunSamples, function(s) s[, j]))
      .Rhat(chainMat)
    }, numeric(1))
    names(rhat) <- paramNms
    maxRhat <- max(rhat[isConvParam[names(rhat) %in% names(ess)]],
                   na.rm = TRUE)
  }

  # --- Adaptive tree ESS ---
  # Three tiers: skip (scalars far off), coarse (500 trees), fine (1000 trees).
  # Avoids expensive RF distance computation when it can't affect the stopping

  # decision, and upgrades to full precision when tree ESS is the binding
  # constraint.
  treeEss <- NA_real_
  treeEssPrecision <- "skip"
  if (!is.null(mcmc$minTreeEss) &&
      requireNamespace("TreeDist", quietly = TRUE)) {

    scalarsFarOff <- !is.null(mcmc$minEss) && minEss < 0.5 * mcmc$minEss
    scalarsConverged <-
      (is.null(mcmc$minEss)  || minEss >= mcmc$minEss) &&
      (is.null(mcmc$maxRhat) || (nRuns >= 2L && !is.na(maxRhat) &&
                                  maxRhat <= mcmc$maxRhat))

    treeEssPrecision <- if (scalarsFarOff) "skip"
                        else if (scalarsConverged) "fine"
                        else "coarse"

    if (treeEssPrecision != "skip") {
      maxPerRun <- if (treeEssPrecision == "fine") 1000L else 500L
      treeEss <- .ComputeTreeEssInLoop(runs, maxPerRun, isStreaming)

      # Upgrade coarse -> fine if estimate is close to threshold
      if (treeEssPrecision == "coarse" && !is.na(treeEss) &&
          treeEss >= 0.8 * mcmc$minTreeEss) {
        treeEss <- .ComputeTreeEssInLoop(runs, 1000L, isStreaming)
        treeEssPrecision <- "fine"
      }
    }
  }

  # Converged only when at least one criterion is set AND all set criteria pass.
  # (Avoids spurious early stopping when no criteria are configured.)
  hasCriteria <- !is.null(mcmc$minEss) || !is.null(mcmc$maxRhat) ||
                 !is.null(mcmc$minTreeEss)
  converged   <- hasCriteria &&
    (is.null(mcmc$minEss)     || minEss >= mcmc$minEss) &&
    (is.null(mcmc$maxRhat)    || (nRuns >= 2L && !is.na(maxRhat) &&
                                   maxRhat <= mcmc$maxRhat)) &&
    (is.null(mcmc$minTreeEss) || (!is.na(treeEss) &&
                                   treeEss >= mcmc$minTreeEss))

  list(converged = converged, minEss = minEss, maxRhat = maxRhat,
       treeEss = treeEss, treeEssPrecision = treeEssPrecision,
       ess = ess, rhat = rhat)
}


#' Check convergence by reading log files from disk (parallel mode)
#'
#' Reads each run's log file via [ReadMkLog()], extracts key parameters,
#' and computes ESS (all runs combined) and R-hat (when `nRuns >= 2`).
#' Returns `NULL` if any log is missing or has fewer than 10 rows.
#' @keywords internal
.CheckConvergenceFromLogs <- function(logFilePaths, paramNames, mcmc) {
  nRuns   <- length(logFilePaths)
  keyCols <- .KeyParamCols(
    matrix(0, 1, length(paramNames), dimnames = list(NULL, paramNames))
  )

  perRunSamples <- lapply(logFilePaths, function(p) {
    if (!file.exists(p)) return(NULL)
    m <- tryCatch(ReadMkLog(p), error = function(e) NULL)
    if (is.null(m) || nrow(m) < 10L) return(NULL)
    m[, keyCols, drop = FALSE]
  })

  if (any(vapply(perRunSamples, is.null, logical(1L)))) return(NULL)

  # Equalise chain lengths — serial runs may produce different sample counts.
  # Keep the most recent (tail) samples to avoid penalising early convergers.
  nRows   <- vapply(perRunSamples, nrow, integer(1L))
  minRows <- min(nRows)
  if (any(nRows != minRows)) {
    perRunSamples <- lapply(perRunSamples, function(s) {
      tail(s, minRows)
    })
  }

  combined <- do.call(rbind, perRunSamples)
  ess <- .EssMatrix(combined)

  # Exclude kPrime nuisance parameters from convergence criteria (M-098)
  isConvParam <- !grepl("^kPrime_", names(ess)) & names(ess) != "log_likelihood"
  minEss <- min(ess[isConvParam], na.rm = TRUE)

  rhat    <- NULL
  maxRhat <- NA_real_
  if (nRuns >= 2L) {
    paramNms <- colnames(perRunSamples[[1]])
    rhat <- vapply(seq_along(paramNms), function(j) {
      chainMat <- do.call(cbind, lapply(perRunSamples, function(s) s[, j]))
      .Rhat(chainMat)
    }, numeric(1))
    names(rhat) <- paramNms
    maxRhat <- max(rhat[isConvParam[names(rhat) %in% names(ess)]],
                   na.rm = TRUE)
  }

  # Tree ESS not available in log-based mode (scalar logs don't contain trees).
  # minTreeEss is only enforced by .CheckConvergence() which has in-memory trees.
  hasCriteria <- !is.null(mcmc$minEss) || !is.null(mcmc$maxRhat)
  converged   <- hasCriteria &&
    (is.null(mcmc$minEss)  || minEss >= mcmc$minEss) &&
    (is.null(mcmc$maxRhat) || (nRuns >= 2L && !is.na(maxRhat) &&
                                maxRhat <= mcmc$maxRhat))

  list(converged = converged, minEss = minEss, maxRhat = maxRhat,
       treeEss = NA_real_, treeEssPrecision = "skip",
       ess = ess, rhat = rhat, perRunSamples = perRunSamples)
}


#' Compute tree ESS from in-memory MCMC runs
#'
#' Extracts tree samples from each run, subsamples to `maxPerRun`,
#' and returns the **minimum** median pseudo-ESS across runs (conservative).
#' Used during convergence checks when `minTreeEss` is set.
#'
#' @param runs List of run state objects (each with `$tree_samples` and
#'   `$tree_saved_idx`).
#' @param maxPerRun Maximum trees per run to use (controls coarse vs fine).
#' @param isStreaming Logical; when `TRUE`, trees are in a flat list
#'   (streaming mode stores all trees, no `saved_idx` subsetting needed).
#' @return Scalar minimum median pseudo-ESS, or `NA_real_` on failure.
#' @keywords internal
.ComputeTreeEssInLoop <- function(runs, maxPerRun, isStreaming) {
  perRunTrees <- lapply(runs, function(r) {
    if (isStreaming) {
      ts <- r$tree_samples
      if (is.null(ts)) return(NULL)
      # Filter out NULL slots (pre-allocated but unused)
      ts <- Filter(Negate(is.null), ts)
      n <- length(ts)
      if (n < 5L) return(NULL)
      ts
    } else {
      idx <- r$tree_saved_idx %||% 0L
      if (idx < 5L) return(NULL)
      r$tree_samples[seq_len(idx)]
    }
  })

  perRunTrees <- Filter(Negate(is.null), perRunTrees)
  if (length(perRunTrees) == 0L) return(NA_real_)

  # Subsample and convert to multiPhylo
  perRunTrees <- lapply(perRunTrees, function(ts) {
    n <- length(ts)
    if (n > maxPerRun) {
      ts <- ts[round(seq(1, n, length.out = maxPerRun))]
    }
    structure(ts, class = "multiPhylo")
  })

  tryCatch({
    essVals <- vapply(perRunTrees, function(chain) {
      .TreeESS(chain, dist_fn = TreeDist::RobinsonFoulds,
               frechet = FALSE)[["medianPseudoESS"]]
    }, double(1))
    essVals <- essVals[is.finite(essVals)]
    if (length(essVals) == 0L) return(NA_real_)
    min(essVals)
  }, error = function(e) NA_real_)
}


# --- Build final result ---

#' Build MkPosterior from all runs
#' @keywords internal
.BuildResult <- function(runs, model, mkd, mcmc, paramNames, logFilePaths,
                         actualIter, stopReason, isTempLog = FALSE) {
  nRuns       <- length(runs)
  # Use logFilePaths (not mcmc$logFile) to determine streaming mode: in
  # parallel runs, .RunParallelRuns() may auto-assign a tempfile log even
  # when mcmc$logFile is NULL, so mcmc$logFile would be stale here.
  isStreaming <- !is.null(logFilePaths)

  # Flush any remaining streaming buffer rows, then trim to actual save count
  for (run in seq_len(nRuns)) {
    if (isStreaming && runs[[run]]$flush_idx > 0L) {
      .FlushBuffer(runs[[run]]$flush_buf, runs[[run]]$flush_idx,
                   runs[[run]]$flush_iter, logFilePaths[run])
      runs[[run]]$flush_idx <- 0L
    }
    idx <- runs[[run]]$saved_idx
    treeIdx <- runs[[run]]$tree_saved_idx %||% idx
    runs[[run]]$tree_samples <- runs[[run]]$tree_samples[seq_len(max(treeIdx, 0L))]
    if (!isStreaming) {
      if (idx > 0L) {
        runs[[run]]$samples <- runs[[run]]$samples[seq_len(idx), , drop = FALSE]
      } else {
        runs[[run]]$samples <- runs[[run]]$samples[integer(0), , drop = FALSE]
      }
    }
  }

  # Per-run summaries
  perRunSummaries <- lapply(runs, function(r) {
    coldAcc <- r$chain_accept[[1]] / pmax(r$chain_propose[[1]], 1L)
    result <- list(
      samples    = if (isStreaming) NULL else r$samples,
      trees      = r$tree_samples,
      acceptance = coldAcc,
      saved_idx  = r$saved_idx
    )
    if (mcmc$nChains > 1L) {
      result$betas      <- r$betas
      result$swap_rates <- r$swap_accept / pmax(r$swap_propose, 1L)
    }
    result
  })

  totalSaved <- sum(vapply(perRunSummaries, `[[`, integer(1), "saved_idx"))

  if (isStreaming) {
    # Streaming mode: return empty sample matrix with correct columns.
    # Caller loads samples via ReadMkLog(result$logFile).
    emptySamples <- matrix(numeric(0), nrow = 0L, ncol = length(paramNames),
                           dimnames = list(NULL, paramNames))
    allTrees <- do.call(c, lapply(perRunSummaries, `[[`, "trees"))
    avgAcceptance <- Reduce(`+`, lapply(perRunSummaries, `[[`,
                                        "acceptance")) / nRuns

    result <- MkPosterior(
      samples = emptySamples,
      trees   = allTrees,
      acceptance = avgAcceptance,
      model = model, data = mkd, mcmc = mcmc,
      warmup = mcmc$warmup, tuning = runs[[1]]$chain_tuning[[1]]
    )
    result$logFile  <- logFilePaths
    result$nSamples <- totalSaved

    if (nRuns > 1L) {
      result$nRuns   <- nRuns
      result$per_run <- perRunSummaries
    }
    if (mcmc$nChains > 1L) {
      result$betas <- runs[[1]]$betas
      result$swap_rates <- perRunSummaries[[1]]$swap_rates
      result$chain_acceptance <- lapply(seq_len(mcmc$nChains), function(ch) {
        runs[[1]]$chain_accept[[ch]] / pmax(runs[[1]]$chain_propose[[ch]], 1L)
      })
    }
    if (!isTempLog) {
      cli::cli_alert_info(c(
        "Streaming mode: {totalSaved} sample{?s} written to \\
         {.file {logFilePaths}}.",
        "i" = "Load with: {.code result$samples <- ReadMkLog(result$logFile)}"
      ))
    }

  } else {
    # In-memory mode (unchanged)
    if (nRuns == 1L) {
      r <- perRunSummaries[[1]]
      result <- MkPosterior(
        samples = r$samples, trees = r$trees,
        acceptance = r$acceptance,
        model = model, data = mkd, mcmc = mcmc,
        warmup = mcmc$warmup, tuning = runs[[1]]$chain_tuning[[1]]
      )
      if (!is.null(r$betas)) {
        result$betas <- r$betas
        result$swap_rates <- r$swap_rates
        result$chain_acceptance <- lapply(seq_len(mcmc$nChains), function(ch) {
          runs[[1]]$chain_accept[[ch]] /
            pmax(runs[[1]]$chain_propose[[ch]], 1L)
        })
      }
    } else {
      allSamples <- do.call(rbind, lapply(perRunSummaries, `[[`, "samples"))
      allTrees   <- do.call(c,     lapply(perRunSummaries, `[[`, "trees"))
      avgAcceptance <- Reduce(`+`, lapply(perRunSummaries, `[[`,
                                          "acceptance")) / nRuns

      result <- MkPosterior(
        samples = allSamples, trees = allTrees,
        acceptance = avgAcceptance,
        model = model, data = mkd, mcmc = mcmc,
        warmup = mcmc$warmup, tuning = runs[[1]]$chain_tuning[[1]]
      )
      result$nRuns   <- nRuns
      result$per_run <- perRunSummaries

      if (mcmc$nChains > 1L) {
        result$betas <- runs[[1]]$betas
        result$swap_rates <- perRunSummaries[[1]]$swap_rates
        result$chain_acceptance <- lapply(seq_len(mcmc$nChains), function(ch) {
          runs[[1]]$chain_accept[[ch]] /
            pmax(runs[[1]]$chain_propose[[ch]], 1L)
        })
      }
    }
  }

  result$stop_reason <- stopReason
  result$actual_iter <- actualIter
  result$treeThin    <- mcmc$treeThin
  result
}


# --- Checkpointing ---

#' Save MCMC checkpoint to RDS
#'
#' Version 1 (in-memory mode): stores full run history (samples, trees).
#' Version 2 (streaming mode): stores chain state only; samples live in the
#' log file.  The large flush_buf and conv_window matrices are excluded.
#'
#' @keywords internal
.SaveCheckpoint <- function(runs, mcmc, iter, paramNames, file,
                            moveWeights = NULL, phase = NULL,
                            model = NULL, serialPhase = NULL) {
  isStreaming <- !is.null(mcmc$logFile)

  serialRuns <- lapply(runs, function(r) {
    # Serialize C++ chain state (XPtrs cannot cross session boundaries).
    # Guard: .RunMkPrimeSingleRun() already serializes on return (chainStates
    # is NULL, chains is populated), so skip when chains are already R lists.
    if (!is.null(r$chainStates)) {
      r$chains <- lapply(r$chainStates, function(ptr) {
        s <- get_mcmc_state(ptr)
        list(
          log_lik        = s$logLik,
          log_prior      = s$logPrior,
          log_post       = s$logPost,
          tree_length    = s$treeLength,
          rel_br_lengths = s$relBrLengths,
          rate_loss      = s$rateLoss,
          rate_log_sd    = s$rateLogSd,
          rate_neo       = s$rateNeo,
          p              = s$p,
          kPrime         = s$kPrime,
          edge           = s$edge,
          beta_scale     = s$betaScale
        )
      })
      r$chainStates <- NULL
    }  # else: chains already serialized, chainStates already NULL
    if (isStreaming) {
      # Omit the large buffer matrices and transient flags; recreated on resume.
      # flush_idx should be 0 (caller flushes before checkpointing).
      r$flush_buf   <- NULL
      r$flush_iter  <- NULL
      r$conv_window <- NULL
      r$flushed     <- NULL
    }
    r
  })

  version <- if (isStreaming) 2L else 1L
  payload <- list(runs = serialRuns, mcmc = mcmc, iter = iter,
                  timestamp = Sys.time(), version = version)
  if (isStreaming) {
    payload$logFilePaths <- .LogFilePaths(mcmc$logFile, length(runs))
    payload$paramNames   <- paramNames
  }
  if (!is.null(moveWeights))  payload$moveWeights  <- moveWeights
  if (!is.null(phase))        payload$phase        <- phase
  if (!is.null(model))        payload$model        <- model
  if (!is.null(serialPhase))  payload$serialPhase  <- serialPhase
  saveRDS(payload, file)
}


# Flush any pending streaming buffers then save a checkpoint, if configured.
# Returns the (possibly modified) runs list so flush_idx resets propagate.
#
# @keywords internal
.FlushAndSaveCheckpoint <- function(runs, nRuns, mcmc, batchEnd,
                                    paramNames, isStreaming, logFilePaths,
                                    moveWeights = NULL, model = NULL) {
  if (!is.null(mcmc$checkpointFile)) {
    if (isStreaming) {
      for (run in seq_len(nRuns)) {
        if (runs[[run]]$flush_idx > 0L) {
          .FlushBuffer(runs[[run]]$flush_buf, runs[[run]]$flush_idx,
                       runs[[run]]$flush_iter, logFilePaths[run])
          runs[[run]]$flush_idx <- 0L
        }
      }
    }
    .SaveCheckpoint(runs, mcmc, batchEnd, paramNames, mcmc$checkpointFile,
                    moveWeights = moveWeights, model = model)
  }
  runs
}


#' Resume MCMC from a checkpoint
#'
#' Loads a checkpoint file and continues the MCMC from where it left off.
#' The remaining iterations will be appended to the existing samples.
#' Called automatically by [RunMkPrime()] when a checkpoint file exists
#' and `overwrite = FALSE`.
#'
#' @param checkpointFile Path to the checkpoint RDS file.
#' @param data A `phyDat` or `MkPrimeData` object (must match original).
#' @param model An `MkPrimeModel` object.  If `NULL` (the default),
#'   the model stored in the checkpoint is used; otherwise it is finalized
#'   using `tree` or re-derived from data.
#' @param tree A `phylo` object (used only for model finalization via
#'   Fitch parsimony).  Can be `NULL` when the checkpoint already contains
#'   a finalized model (all checkpoints since M-149).
#' @param neomorphic,knownStates Passed to [MkPrimeData()] if `data`
#'   is a `phyDat` object.
#'
#' @return An `MkPosterior` object with combined samples.
#' @export
ResumeMkPrime <- function(checkpointFile, data, tree = NULL,
                           neomorphic = integer(0),
                           knownStates = integer(0),
                           model = NULL) {
  checkpoint <- readRDS(checkpointFile)

  version <- checkpoint$version %||% 1L
  if (!version %in% c(1L, 2L)) {
    cli::cli_abort("Unsupported checkpoint version: {version}.")
  }
  isStreaming <- version == 2L

  if (inherits(data, "MkPrimeData")) {
    mkd <- data
  } else {
    mkd <- MkPrimeData(data, neomorphic = neomorphic,
                       knownStates = knownStates)
  }

  # Model resolution: prefer checkpoint model (already finalized), then
  # user-supplied model + tree, then fresh MkPrimeModel + NJ tree.
  if (is.null(model) && !is.null(checkpoint$model)) {
    model <- checkpoint$model
  } else {
    if (is.null(model)) model <- MkPrimeModel()
    if (is.null(tree)) {
      # Need tree only for .FinalizeModel (Fitch parsimony score)
      njInput <- if (inherits(data, "phyDat")) data else mkd$phyDat
      tree <- TreeTools::NJTree(njInput, edgeLengths = TRUE)
    }
    tree <- TreeTools::Preorder(tree)
    model <- .FinalizeModel(model, tree, mkd)
  }

  # Validate tip labels when a tree is available
  if (!is.null(tree)) {
    treeTips <- tree$tip.label
    dataTaxa <- rownames(mkd$matrix)
    missing  <- setdiff(dataTaxa, treeTips)
    extra    <- setdiff(treeTips, dataTaxa)
    if (length(missing) > 0L || length(extra) > 0L) {
      msgs <- character(0)
      if (length(missing) > 0L) {
        msgs <- c(msgs,
          "x" = "{length(missing)} taxon{?/a} in data but not in tree: {.val {missing}}.")
      }
      if (length(extra) > 0L) {
        msgs <- c(msgs,
          "x" = "{length(extra)} tip{?s} in tree but not in data: {.val {extra}}.")
      }
      cli::cli_abort(c(
        "Tip labels in {.arg tree} do not match taxa in {.arg data}.",
        msgs,
        "i" = "Every taxon in the data must appear as a tip label in the tree, and vice versa."
      ))
    }
  }

  runs <- checkpoint$runs
  mcmc <- checkpoint$mcmc
  startIter <- checkpoint$iter + 1L
  nRuns <- mcmc$nRuns

  # Derive nEdge and paramNames: version 2 stores them directly;
  # version 1 derives from the sample matrix column names.
  if (isStreaming) {
    nEdge      <- nrow(runs[[1]]$chains[[1]]$edge)
    paramNames <- checkpoint$paramNames
    logFilePaths   <- checkpoint$logFilePaths

    # Check whether log files exist.  If they were temp files (deleted on
    # clean exit), recreate them so the chain can resume from checkpoint
    # state.  Previous samples are lost (they were loaded into memory in
    # the original result).
    logsRecreated <- FALSE
    for (i in seq_along(logFilePaths)) {
      if (!file.exists(logFilePaths[i])) {
        logFilePaths[i] <- tempfile("mkp_resume_", fileext = ".log")
        writeLines(paste(c("Sample", paramNames), collapse = "\t"),
                   logFilePaths[i])
        runs[[i]]$saved_idx <- 0L
        logsRecreated <- TRUE
      }
    }
    if (logsRecreated) {
      cli::cli_alert_info(
        "Original log files not found (temp files cleaned up). \\
         Resuming from checkpoint state; previous samples unavailable."
      )
    }
    convWindowSize <- .ComputeConvWindowSize(mcmc)
  } else {
    nEdge      <- sum(grepl("^br_", colnames(runs[[1]]$samples)))
    paramNames <- colnames(runs[[1]]$samples)
    logFilePaths   <- NULL
    convWindowSize <- 0L
  }

  hasNeo <- any(mkd$type == "neomorphic")
  nTrans <- sum(mkd$type == "transformational")

  qHet <- isTRUE(model$qHeterogeneity)
  moves <- .BuildMoves(nEdge, nTrans, hasNeo, mcmc, fixTopology = FALSE,
                       kPrimePrior = model$kPrimePrior %||% "geometric",
                       qHeterogeneity = qHet,
                       joint2d = isTRUE(mcmc$joint2d))

  if (identical(mcmc$thin, "auto")) {
    mcmc$thin <- length(moves)
  }
  # On resume, don't re-adapt thin (checkpoint has the adapted value)
  mcmc$thinWasAuto <- FALSE
  mcmc$treeThinWasAuto <- FALSE
  if (is.null(mcmc$treeThin)) mcmc$treeThin <- mcmc$thin
  treeEvery <- as.integer(mcmc$treeThin / mcmc$thin)

  if (isStreaming) {
    # Rewind each log file to the checkpoint's saved_idx.  Any samples
    # flushed after the last checkpoint are discarded — the chain state
    # doesn't cover them.  Buffer reinit happens inside .RunMkPrimeSingleRun.
    for (run in seq_len(nRuns)) {
      .TruncateLogToN(logFilePaths[run],
                      as.integer(runs[[run]]$saved_idx %||% 0L))
    }
  }

  tipLabels   <- tree$tip.label %||% rownames(mkd$matrix)
  isLogseries <- identical(model$kPrimePrior, "logseries")
  pCols       <- if (isLogseries) 0L else 1L
  brColStart  <- 5L + pCols + (any(mkd$type == "neomorphic")) + qHet + nTrans + 1L

  # --- Sequential per-run execution ---
  stopReason <- "max_iter"
  actualIter <- if (is.finite(mcmc$nIter)) mcmc$nIter else startIter - 1L

  # M-149 #7: per-run startIters from individual actual_iter
  perRunStarts <- vapply(
    runs,
    function(r) as.integer((r$actual_iter %||% (startIter - 1L)) + 1L),
    integer(1L)
  )

  if (isStreaming && nRuns >= 2L && !is.null(mcmc$maxRhat)) {
    # M-146: cross-run R-hat convergence orchestrator
    # M-149 #6: skip Phase 1 if checkpoint was during Phase 2
    serialResult <- .RunSerialRuns(mkd, model, mcmc, runs, moves,
                                    tipLabels, paramNames, nEdge,
                                    brColStart, logFilePaths,
                                    convWindowSize, treeFile = NULL,
                                    startIters = perRunStarts,
                                    startPhase = checkpoint$serialPhase %||% 1L)
    runs       <- serialResult$runs
    stopReason <- serialResult$stopReason
    actualIter <- serialResult$actualIter
  } else {
    for (run in seq_len(nRuns)) {
      runs[[run]] <- .RunMkPrimeSingleRun(
        mkd, model, mcmc, runs[[run]], moves, tipLabels, run,
        paramNames, nEdge, brColStart,
        logFilePath    = if (isStreaming) logFilePaths[run] else NULL,
        cancelFile     = mcmc$cancelFile,
        checkpointFile = if (nRuns == 1L) mcmc$checkpointFile else NULL,
        startIter      = perRunStarts[run],
        isStreaming    = isStreaming,
        convWindowSize = convWindowSize,
        treeFile       = NULL,
        resumeMoveWeights = runs[[run]]$moveWeights %||% checkpoint$moveWeights
      )
      stopReason <- runs[[run]]$stop_reason
      actualIter <- runs[[run]]$actual_iter
      if (stopReason == "cancelled") break
    }
  }

  .BuildResult(runs, model, mkd, mcmc, paramNames, logFilePaths,
               actualIter, stopReason)
}


# --- Start perturbation ---

#' Generate a perturbed starting tree for independent runs
#' @keywords internal
.PerturbStart <- function(tree) {
  nTip <- length(tree$tip.label)
  if (nTip < 4L) return(tree)

  nNni <- sample(2:5, 1)
  for (i in seq_len(nNni)) {
    treeLength <- sum(tree$edge.length)
    relBr <- tree$edge.length / treeLength
    prop <- ProposeNni(tree, treeLength, relBr)
    tree <- prop$tree
    tree$edge.length <- treeLength * prop$rel_br_lengths
  }

  tree$edge.length <- tree$edge.length *
    exp(rnorm(length(tree$edge.length), sd = 0.1))
  tree
}


# --- Internal helpers ---

#' Build the progress info list for callbacks
#' @keywords internal
.BuildProgressInfo <- function(runs, iter, mcmc, startTime,
                               recentAcc, paramNames,
                               phase = "Sample") {
  nRuns       <- length(runs)
  isStreaming <- !is.null(mcmc$logFile)
  runSamples <- lapply(runs, function(r) {
    if (isStreaming) {
      .ConvWindowRows(r, minRows = 1L)
    } else {
      idx <- r$saved_idx
      if (idx > 0L) r$samples[seq_len(idx), , drop = FALSE] else NULL
    }
  })

  hasNeo <- "rate_loss" %in% paramNames
  currentState <- lapply(runs, function(r) {
    s <- get_mcmc_state(r$chainStates[[1]])
    st <- list(log_lik = s$logLik, log_prior = s$logPrior,
               tree_length = s$treeLength,
               rate_log_sd = s$rateLogSd, p = s$p)
    if (hasNeo) st$rate_loss <- s$rateLoss
    st
  })

  list(
    iter = iter,
    nIter = mcmc$nIter,
    warmup = mcmc$warmup,
    inWarmup = phase == "Warmup",
    phase = phase,
    nRuns = nRuns,
    nChains = mcmc$nChains,
    runSamples = runSamples,
    currentState = currentState,
    recentAcceptance = recentAcc,
    elapsed = as.numeric(proc.time()["elapsed"] - startTime),
    paramNames = paramNames
  )
}


#' Build scale tuning matrix for run_mcmc_batch_cpp (nChains × nMoves)
#' @keywords internal
.BuildScaleTuningMatrix <- function(chainTuning, moves) {
  nChains <- length(chainTuning)
  nMoves <- length(moves)
  mat <- matrix(0.5, nChains, nMoves)
  for (ch in seq_len(nChains)) {
    tun <- chainTuning[[ch]]
    for (m in seq_along(moves)) {
      mat[ch, m] <- switch(moves[[m]]$name,
        tree_length = tun$scale_tree_length,
        rate_loss   = tun$scale_rate_loss,
        rate_log_sd = tun$scale_rate_log_sd,
        rate_neo    = tun$scale_rate_neo %||% 0.5,
        neo_joint   = tun$scale_neo_joint %||% tun$scale_rate_loss,
        beta_scale  = tun$scale_beta_scale %||% 0.5,
        kprime_alpha = tun$scale_kprime_alpha %||% 0.05,
        kprime_beta  = tun$scale_kprime_beta %||% 0.05,
        joint_tl_rls = tun$scale_joint_tl_rls %||% 0.5,
        joint_tl_rl  = tun$scale_joint_tl_rl %||% 0.5,
        dirichlet_branch = tun$dirichlet_alpha %||% 10,
        local_dirichlet  = tun$local_dirichlet_alpha %||% 10,
        p           = 0.5,  # Gibbs move: scale ignored by C++; placeholder
        0.5
      )
    }
  }
  mat
}


#' Build slice width matrix for run_mcmc_batch_cpp (nChains × nMoves)
#' @keywords internal
.BuildSliceWidthMatrix <- function(chainTuning, moves) {
  nChains <- length(chainTuning)
  nMoves <- length(moves)
  mat <- matrix(1.0, nChains, nMoves)  # Default width 1.0 (ignored for non-slice)
  for (ch in seq_len(nChains)) {
    tun <- chainTuning[[ch]]
    for (m in seq_along(moves)) {
      mat[ch, m] <- switch(moves[[m]]$name,
        slice_rate_loss   = tun$slice_width_rate_loss %||% 1.0,
        slice_rate_neo    = tun$slice_width_rate_neo %||% 1.0,
        slice_rate_log_sd = tun$slice_width_rate_log_sd %||% 1.0,
        slice_tree_length = tun$slice_width_tree_length %||% 1.0,
        slice_beta_scale  = tun$slice_width_beta_scale %||% 1.0,
        1.0
      )
    }
  }
  mat
}


#' Build geometric temperature ladder
#' @keywords internal
.BuildTemperatureLadder <- function(nChains, heat) {
  if (nChains == 1L) return(1.0)
  heat^(seq(0, 1, length.out = nChains))
}


#' Build joint-rho matrix (nChains × nMoves) for 2D joint Bactrian moves
#' @keywords internal
.BuildJointRhoMatrix <- function(chainRhos, moves, nChains) {
  nMoves <- length(moves)
  mat <- matrix(0.0, nChains, nMoves)
  for (ch in seq_len(nChains)) {
    rhos <- chainRhos[[ch]]
    for (m in seq_along(moves)) {
      mat[ch, m] <- switch(moves[[m]]$name,
        joint_tl_rls = rhos$rho_tl_rls %||% 0.0,
        joint_tl_rl  = rhos$rho_tl_rl  %||% 0.0,
        0.0
      )
    }
  }
  mat
}


#' Accumulate a cold-chain state snapshot into the rho sample buffer
#'
#' During warmup, C++ saves no thinned samples, so we extract scalar
#' parameters from the chain state after each batch.  The buffer is
#' capped at 500 rows (rolling window).
#' @keywords internal
.AccumulateRhoSnapshot <- function(buf, state, hasNeo, paramNames) {
  row <- numeric(length(paramNames))
  names(row) <- paramNames
  row["tree_length"] <- state$treeLength
  row["rate_log_sd"] <- state$rateLogSd
  if (hasNeo) row["rate_loss"] <- state$rateLoss
  buf <- if (is.null(buf)) {
    matrix(row, nrow = 1, dimnames = list(NULL, paramNames))
  } else {
    rbind(buf, row)
  }
  if (nrow(buf) > 500L) buf <- buf[(nrow(buf) - 499L):nrow(buf), , drop = FALSE]
  buf
}


#' Estimate posterior correlations for 2D joint proposals from recent samples
#' @keywords internal
.EstimateJointRhos <- function(samples, hasNeo) {
  rhos <- list(rho_tl_rls = 0.0, rho_tl_rl = 0.0)
  if (is.null(samples) || nrow(samples) < 50) return(rhos)

  # tree_length × rate_log_sd
  if (all(c("tree_length", "rate_log_sd") %in% colnames(samples))) {
    tl <- samples[, "tree_length"]
    rls <- samples[, "rate_log_sd"]
    ok <- tl > 0 & rls > 0
    if (sum(ok) >= 30) {
      # u.123: cor() returns NaN/NA for constant input; guard before clamping
      rho <- suppressWarnings(cor(log(tl[ok]), log(rls[ok])))
      if (is.finite(rho)) {
        rhos$rho_tl_rls <- max(-0.95, min(0.95, rho))
      }
    }
  }

  # tree_length × rate_loss
  if (hasNeo && all(c("tree_length", "rate_loss") %in% colnames(samples))) {
    tl <- samples[, "tree_length"]
    rl <- samples[, "rate_loss"]
    ok <- tl > 0 & rl > 0
    if (sum(ok) >= 30) {
      rho <- suppressWarnings(cor(log(tl[ok]), log(rl[ok])))
      if (is.finite(rho)) {
        rhos$rho_tl_rl <- max(-0.95, min(0.95, rho))
      }
    }
  }

  rhos
}


#' Initialize MCMC state from tree and data
#' @keywords internal
.InitState <- function(tree, mkd, model) {
  treeLength <- sum(tree$edge.length)
  relBr <- tree$edge.length / treeLength

  kPrime <- mkd$kObs
  knownIdx <- which(mkd$type == "known")
  if (length(knownIdx)) kPrime[knownIdx] <- mkd$known_k[knownIdx]

  hasNeo <- any(mkd$type == "neomorphic")

  state <- list(
    tree = tree,
    tree_length = treeLength,
    rel_br_lengths = relBr,
    rate_loss = 1.0,
    rate_log_sd = 0.5,
    kPrime = as.integer(kPrime)
  )

  # kPrime hyperparameters depend on prior choice
  if (identical(model$kPrimePrior, "beta_geometric")) {
    state$kprime_alpha <- model$kprimeAlpha
    state$kprime_beta  <- model$kprimeBeta
  } else if (!identical(model$kPrimePrior, "logseries")) {
    # Hierarchical geometric: shared p
    state$p <- 0.5
  }

  # Partition rate scalar for neomorphic characters
  if (hasNeo) {
    state$rate_neo <- 1.0
  }

  # M-052: beta_scale for Q-matrix heterogeneity
  if (isTRUE(model$qHeterogeneity)) {
    state$beta_scale <- 1.0
  }

  # Tree is already preorder (reordered at init); use internal fast-path
  state$log_lik <- .MkpLogLikelihood(
    tree, mkd,
    kPrime = state$kPrime,
    rate_loss = state$rate_loss,
    rate_log_sd = state$rate_log_sd,
    nCat = model$nCat,
    coding = model$coding,
    rate_neo = state$rate_neo %||% 1.0,
    relabel = model$relabel
  )
  state$log_prior <- LogPrior(state, model, mkd)
  state$log_post <- state$log_lik + state$log_prior

  state
}


#' Build move schedule
#' @keywords internal
.BuildMoves <- function(nEdge, nTrans, hasNeo, mcmc,
                        fixTopology = FALSE,
                        kPrimePrior = "geometric",
                        qHeterogeneity = FALSE,
                        joint2d = TRUE) {
  moves <- list(
    list(name = "tree_length", type = "scale", target = "tree_length",
         weight = 1, dim = 1L),
    list(name = "branch_lengths", type = "beta_simplex",
         target = "rel_br_lengths", weight = max(1, nEdge / 3), dim = 1L)
  )

  if (!fixTopology && nEdge >= 5L) {
    moves <- c(moves, list(
      list(name = "nni", type = "nni", target = NULL,
           weight = max(1, nEdge / 2), dim = 1L),
      list(name = "spr", type = "spr", target = NULL,
           weight = max(1, nEdge / 4), dim = 1L)
    ))
    # Gibbs topology moves (M-090).
    # Per-invocation cost is O(nEdge * depth * nChar) because Gibbs evaluates
    # all candidate regraft/swap positions.  M-113 showed this outweighs the
    # mixing gain on trees with > ~20 tips, so cap the initial weight (M-115).
    # The adaptive scheduler refines from here during warmup.
    gibbsCap <- 10L
    if (isTRUE(mcmc$gibbsSpr)) {
      moves <- c(moves, list(
        list(name = "gibbs_spr", type = "gibbs_spr", target = NULL,
             weight = max(1L, min(nEdge / 4, gibbsCap)), dim = 1L)
      ))
    }
    if (isTRUE(mcmc$gibbsSubtreeSwap)) {
      moves <- c(moves, list(
        list(name = "gibbs_subtree_swap", type = "gibbs_subtree_swap",
             target = NULL,
             weight = max(1L, min(nEdge / 6, gibbsCap)), dim = 1L)
      ))
    }
    # Weighted moves (M-090)
    if (isTRUE(mcmc$weightedBranchScale)) {
      moves <- c(moves, list(
        list(name = "weighted_branch_lengths", type = "weighted_branch_scale",
             target = "rel_br_lengths", weight = max(1, nEdge / 6), dim = 1L)
      ))
    }
    if (isTRUE(mcmc$weightedSpr)) {
      moves <- c(moves, list(
        list(name = "weighted_spr", type = "weighted_spr", target = NULL,
             weight = max(1, nEdge / 8), dim = 1L)
      ))
    }
    if (isTRUE(mcmc$weightedSubtreeSwap)) {
      moves <- c(moves, list(
        list(name = "weighted_subtree_swap", type = "weighted_subtree_swap",
             target = NULL, weight = max(1, nEdge / 8), dim = 1L)
      ))
    }
    # TBR topology move (M-053)
    if (isTRUE(mcmc$tbr)) {
      moves <- c(moves, list(
        list(name = "tbr", type = "tbr", target = NULL,
             weight = max(1, nEdge / 4), dim = 1L)
      ))
    }
    # Parsimony-guided SPR (M-119)
    if (isTRUE(mcmc$pSpr)) {
      moves <- c(moves, list(
        list(name = "pspr", type = "pspr", target = NULL,
             weight = max(1, nEdge / 4), dim = 1L)
      ))
    }
    # Block Gibbs branch-length sweep (M-054 reframed)
    if (isTRUE(mcmc$blockGibbsBranch)) {
      moves <- c(moves, list(
        list(name = "block_gibbs_branch", type = "block_gibbs_branch",
             target = "rel_br_lengths",
             weight = max(1L, min(nEdge / 4, gibbsCap)),
             dim = as.integer(nEdge))
      ))
    }
  }

  # Block Dirichlet simplex branch-length move (M-125)
  if (isTRUE(mcmc$dirichletBranch) && nEdge >= 4L) {
    # M-127: K=5 empirically optimal (K sweep: 22× baseline at K=5 vs 2× at K=10)
    nCatsDirichlet <- as.integer(mcmc$dirichletK %||% min(nEdge, 5L))
    moves <- c(moves, list(
      list(name = "dirichlet_branch", type = "dirichlet_simplex",
           target = "rel_br_lengths",
           weight = max(1L, nEdge %/% 4L),
           dim = nCatsDirichlet,
           nCats = nCatsDirichlet)
    ))
  }

  # M-127: localized Dirichlet — connected edges for compact partial eval
  if (isTRUE(mcmc$localDirichlet) && nEdge >= 4L) {
    nCatsLocal <- as.integer(mcmc$localDirichletK %||% min(nEdge, 6L))
    moves <- c(moves, list(
      list(name = "local_dirichlet", type = "local_dirichlet_simplex",
           target = "rel_br_lengths",
           weight = max(1L, nEdge %/% 4L),
           dim = nCatsLocal,
           nCats = nCatsLocal)
    ))
  }

  if (nTrans > 0) {
    kPrimeMoves <- list(
      # Univariate integer walk (reduced weight — Gibbs sweep does heavy lifting)
      list(name = "kPrime", type = "int_walk", target = "kPrime",
           weight = max(1, nTrans), dim = 1L),
      # Gibbs kPrime sweep: sample all k'_i from full conditionals
      list(name = "gibbs_kPrime", type = "gibbs_kprime_sweep",
           target = "kPrime", weight = max(1, nTrans), dim = as.integer(nTrans)),
      # Block kPrime shift: shift all trans chars by same delta
      list(name = "block_kPrime", type = "block_kprime_shift",
           target = "kPrime", weight = 2, dim = 1L)
    )
    if (identical(kPrimePrior, "beta_geometric")) {
      # Scale proposals for shared (α, β) hyperparameters
      kPrimeMoves <- c(kPrimeMoves, list(
        list(name = "kprime_alpha", type = "kprime_alpha",
             target = "kprime_alpha", weight = 0.05, dim = 1L),
        list(name = "kprime_beta", type = "kprime_beta",
             target = "kprime_beta", weight = 0.05, dim = 1L),
        # Prior-only slice samplers — robust, tuning-free exploration
        list(name = "slice_kprime_alpha", type = "slice_kprime_hyper",
             target = "kprime_alpha", weight = 1, dim = 1L,
             sliceParamIdx = 0L),
        list(name = "slice_kprime_beta", type = "slice_kprime_hyper",
             target = "kprime_beta", weight = 1, dim = 1L,
             sliceParamIdx = 1L)
      ))
    } else if (!identical(kPrimePrior, "logseries")) {
      # Conjugate Gibbs draw: p | k' ~ Beta(a + nTrans, b + sum(k' - kObs))
      kPrimeMoves <- c(kPrimeMoves, list(
        list(name = "p", type = "gibbs_p", target = "p", weight = 1,
             dim = 1L)
      ))
    }
    moves <- c(moves, kPrimeMoves)
  }

  if (hasNeo) {
    moves <- c(moves, list(
      list(name = "rate_loss", type = "scale", target = "rate_loss",
           weight = 1.5, dim = 1L),
      list(name = "rate_neo", type = "scale", target = "rate_neo",
           weight = 1, dim = 1L),
      list(name = "neo_joint", type = "neo_joint", target = "rate_loss",
           weight = 1.5, dim = 2L)
    ))
  }

  moves <- c(moves, list(
    list(name = "rate_log_sd", type = "scale", target = "rate_log_sd",
         weight = 1.5, dim = 1L)
  ))

  # M-052: beta_scale for Q-matrix heterogeneity
  if (isTRUE(qHeterogeneity)) {
    moves <- c(moves, list(
      list(name = "beta_scale", type = "scale", target = "beta_scale",
           weight = 1, dim = 1L)
    ))
  }

  # --- Slice sampling moves for scalar parameters ---
  # paramIdx: 0=treeLength, 1=rateLoss, 2=rateLogSd, 3=rateNeo, 4=betaScale
  if (hasNeo) {
    moves <- c(moves, list(
      list(name = "slice_rate_loss", type = "slice", target = "rate_loss",
           weight = 1.5, dim = 1L, sliceParamIdx = 1L),
      list(name = "slice_rate_neo", type = "slice", target = "rate_neo",
           weight = 1, dim = 1L, sliceParamIdx = 3L)
    ))
  }
  moves <- c(moves, list(
    list(name = "slice_rate_log_sd", type = "slice", target = "rate_log_sd",
         weight = 1.5, dim = 1L, sliceParamIdx = 2L)
  ))
  if (isTRUE(qHeterogeneity)) {
    moves <- c(moves, list(
      list(name = "slice_beta_scale", type = "slice", target = "beta_scale",
           weight = 1, dim = 1L, sliceParamIdx = 4L)
    ))
  }

  # --- M-120: 2D joint Bactrian proposals ---
  if (isTRUE(joint2d)) {
    moves <- c(moves, list(
      list(name = "joint_tl_rls", type = "joint_2d",
           target = "tree_length", weight = 1, dim = 2L)
    ))
    if (hasNeo) {
      moves <- c(moves, list(
        list(name = "joint_tl_rl", type = "joint_2d",
             target = "tree_length", weight = 1, dim = 2L)
      ))
    }
  }

  # --- Scalar weight floor ---
  # Scalar model-parameter moves (dim=1, non-topology) can be starved when
  # kPrime and branch_lengths dominate the weight budget. Guarantee each
  # scalar move gets at least 2% of the pre-floor total weight.
  # Joint 2D moves also get the floor so they're comparable to individual
  # scalar moves they complement.
  scalarTypes <- c("scale", "int_walk", "gibbs_p", "scale_p", "slice",
                    "kprime_alpha", "kprime_beta")
  totalWeight <- sum(vapply(moves, `[[`, numeric(1), "weight"))
  floorVal <- totalWeight * 0.02
  for (i in seq_along(moves)) {
    m <- moves[[i]]
    if ((m$dim == 1L && m$type %in% scalarTypes) || m$type == "joint_2d") {
      moves[[i]]$weight <- max(m$weight, floorVal)
    }
  }

  moves
}


# --- Move type integer codes (must match src/mcmc.cpp) ---
# 0=scale_tl, 1=scale_rl, 2=scale_rls, 3=scale_rn,
# 4=beta_simplex, 5=nni, 6=spr, 7=int_walk, 8=scale_p (legacy),
# 9=gibbs_p, 10=gibbs_spr, 11=gibbs_subtree_swap,
# 12=weighted_br_scale, 13=weighted_spr, 14=weighted_subtree_swap,
# 15=block_gibbs_branch, 16=beta_scale (M-052), 17=tbr (M-053),
# 25=gibbs_kprime_sweep, 26=block_kprime_shift,
# 27=scale_kprime_alpha, 28=scale_kprime_beta
.kMoveTypes <- c(
  tree_length = 0L, rate_loss = 1L, rate_log_sd = 2L,
  rate_neo = 3L, branch_lengths = 4L,
  nni = 5L, spr = 6L, kPrime = 7L, p = 9L,
  gibbs_spr = 10L, gibbs_subtree_swap = 11L,
  weighted_branch_lengths = 12L,
  weighted_spr = 13L,
  weighted_subtree_swap = 14L,
  block_gibbs_branch = 15L,
  beta_scale = 16L,
  tbr = 17L,
  neo_joint = 18L,
  slice_rate_loss = 19L,
  slice_rate_neo = 19L,
  slice_rate_log_sd = 19L,
  slice_tree_length = 19L,
  slice_beta_scale = 19L,
  pspr = 20L,
  joint_tl_rls = 21L,
  joint_tl_rl = 22L,
  dirichlet_branch = 23L,
  local_dirichlet = 24L,
  gibbs_kPrime = 25L,
  block_kPrime = 26L,
  kprime_alpha = 27L,
  kprime_beta = 28L
)

#' Initialize the C++ MCMC data structure (call once before loop)
#' @keywords internal
.InitMcmcData <- function(mkd, model) {
  # Replace NA with -1 in tip states for C++
  parts <- lapply(mkd$partitions, function(p) {
    ts <- p$tip_states
    ts[is.na(ts)] <- -1L
    storage.mode(ts) <- "integer"
    p$tip_states <- ts
    # Ensure k is integer (NA_integer_ for non-known)
    if (is.na(p$k)) p$k <- 0L
    p
  })
  prepare_mcmc_data(
    parts, as.integer(mkd$kObs), mkd$type,
    any(mkd$type == "neomorphic"),
    model$nCat, model$coding, model$relabel,
    model$treeLengthShape, model$treeLengthRate,
    model$rateLossMeanlog, model$rateLossSdlog,
    model$rateLogSdShape, model$rateLogSdRate,
    model$rateNeoMeanlog, model$rateNeoSdlog,
    model$kprimeHyperA, model$kprimeHyperB,
    identical(model$kPrimePrior, "logseries"),
    model$kprimeLogseriesC %||% 0.7,
    identical(model$kPrimePrior, "beta_geometric"),
    isTRUE(model$qHeterogeneity),
    model$nBetaCat %||% 4L,
    model$betaScaleShape %||% 1.0,
    model$betaScaleRate %||% 1.0
  )
}

#' Convert an R state to a C++ XPtr<McmcState>
#' @keywords internal
.InitMcmcChain <- function(state) {
  init_mcmc_state(
    state$tree$edge[, 1], state$tree$edge[, 2],
    state$rel_br_lengths, state$tree_length,
    state$rate_loss, state$rate_log_sd,
    state$rate_neo %||% 1.0, state$p %||% 0.5,
    as.integer(state$kPrime),
    state$log_lik, state$log_prior,
    state$beta_scale %||% 1.0,
    state$kprime_alpha %||% 1.0,
    state$kprime_beta %||% 1.0
  )
}

#' Execute a move via C++ XPtr engine
#'
#' @return list(accept, statePtr) where statePtr is the (possibly updated) XPtr.
#' @keywords internal
#' Execute a move via C++ XPtr engine (or R fallback for tests)
#'
#' @keywords internal
.DoMove <- function(move, stateOrPtr, mkdOrData = NULL, modelOrTuning = NULL,
                    tuning = NULL, beta = 1.0,
                    transIdx = integer(0), mcmcData = NULL) {
  # Detect XPtr mode: mcmcData is provided and stateOrPtr is externalptr
  if (!is.null(mcmcData) && inherits(stateOrPtr, "externalptr")) {
    moveCode <- .kMoveTypes[[move$name]]
    charIdx <- if (moveCode == 7L && length(transIdx) > 0L) {
      sample(transIdx, 1L) - 1L
    } else {
      0L
    }
    scaleTun <- switch(move$name,
      tree_length = tuning$scale_tree_length,
      rate_loss   = tuning$scale_rate_loss,
      rate_log_sd = tuning$scale_rate_log_sd,
      rate_neo    = tuning$scale_rate_neo,
      neo_joint   = tuning$scale_neo_joint %||% tuning$scale_rate_loss,
      beta_scale  = tuning$scale_beta_scale,
      kprime_alpha = tuning$scale_kprime_alpha %||% 0.05,
      kprime_beta  = tuning$scale_kprime_beta %||% 0.05,
      dirichlet_branch = tuning$dirichlet_alpha %||% 0.1,
      local_dirichlet = tuning$local_dirichlet_alpha %||% 0.1,
      0.5  # default; gibbs_p ignores scaleTun (returns before using it)
    )
    # For dirichlet_branch / local_dirichlet, intWalkWindow carries nCats
    iww <- if (!is.null(move$nCats)) as.integer(move$nCats)
           else tuning$int_walk_window
    accepted <- do_move_cpp(
      mcmcData, stateOrPtr, moveCode, charIdx,
      scaleTun, tuning$beta_simplex, iww, beta
    )
    return(list(accept = accepted, statePtr = stateOrPtr))
  }

  # Fallback: old R interface for tests
  # stateOrPtr = R state list, mkdOrData = mkd, modelOrTuning = model
  state <- stateOrPtr
  mkd <- mkdOrData
  model <- modelOrTuning
  proposed <- state

  switch(move$type,
    scale = {
      prop <- ProposeScale(state[[move$target]],
                           tuning = tuning[[paste0("scale_", move$target)]])
      proposed[[move$target]] <- prop$value
      logHastings <- prop$logHastings
    },
    beta_simplex = {
      prop <- ProposeBetaSimplex(state$rel_br_lengths,
                                 tuning = tuning$beta_simplex)
      proposed$rel_br_lengths <- prop$value
      logHastings <- prop$logHastings
    },
    nni = {
      prop <- ProposeNni(state$tree, state$tree_length,
                         state$rel_br_lengths)
      proposed$tree <- prop$tree
      proposed$rel_br_lengths <- prop$rel_br_lengths
      logHastings <- prop$logHastings
    },
    spr = {
      prop <- ProposeSpr(state$tree, state$tree_length,
                         state$rel_br_lengths)
      proposed$tree <- prop$tree
      proposed$rel_br_lengths <- prop$rel_br_lengths
      logHastings <- prop$logHastings
    },
    int_walk = {
      transIdx <- which(mkd$type == "transformational")
      charI <- sample(transIdx, 1L)
      prop <- ProposeBoundedIntWalk(
        state$kPrime[charI],
        lower = mkd$kObs[charI],
        window = tuning$int_walk_window
      )
      proposed$kPrime[charI] <- prop$value
      logHastings <- prop$logHastings
    },
    gibbs_p = {
      # Conjugate Beta draw; acceptance = 1, no MH step needed
      transIdx <- which(mkd$type == "transformational")
      sumU <- sum(state$kPrime[transIdx] - mkd$kObs[transIdx])
      proposed$p <- rbeta(1L,
        shape1 = model$kprimeHyperA + length(transIdx),
        shape2 = model$kprimeHyperB + sumU)
      proposed$log_prior <- LogPrior(proposed, model, mkd)
      return(list(accept = TRUE, state = proposed))
    }
  )

  if (!is.finite(logHastings)) {
    return(list(accept = FALSE, state = state))
  }

  proposed$log_prior <- LogPrior(proposed, model, mkd)
  if (!is.finite(proposed$log_prior)) {
    return(list(accept = FALSE, state = state))
  }

  tmpTree <- proposed$tree
  tmpTree$edge.length <- proposed$tree_length * proposed$rel_br_lengths
  proposed$log_lik <- .MkpLogLikelihood(
    tmpTree, mkd,
    kPrime = proposed$kPrime,
    rate_loss = proposed$rate_loss,
    rate_log_sd = proposed$rate_log_sd,
    nCat = model$nCat,
    coding = model$coding,
    rate_neo = proposed$rate_neo %||% 1.0,
    relabel = model$relabel
  )
  proposed$log_post <- proposed$log_lik + proposed$log_prior

  logAlpha <- beta * (proposed$log_lik - state$log_lik) +
              (proposed$log_prior - state$log_prior) +
              logHastings

  if (is.finite(logAlpha) && log(runif(1)) < logAlpha) {
    list(accept = TRUE, state = proposed)
  } else {
    list(accept = FALSE, state = state)
  }
}
#' Propose a swap between two adjacent chains
#' @keywords internal
.ProposeChainSwap <- function(chains, betas) {
  nChains <- length(chains)
  if (nChains < 2L) {
    return(list(chains = chains, pair = NULL, accepted = FALSE))
  }

  i <- sample.int(nChains - 1L, 1L)
  j <- i + 1L

  # Handle both XPtr and R-list states
  lik_j <- if (inherits(chains[[j]], "externalptr")) {
    get_state_log_lik(chains[[j]])
  } else {
    chains[[j]]$log_lik
  }
  lik_i <- if (inherits(chains[[i]], "externalptr")) {
    get_state_log_lik(chains[[i]])
  } else {
    chains[[i]]$log_lik
  }

  logAlpha <- (betas[i] - betas[j]) * (lik_j - lik_i)

  accepted <- is.finite(logAlpha) && log(runif(1)) < logAlpha
  if (accepted) {
    tmp <- chains[[i]]
    chains[[i]] <- chains[[j]]
    chains[[j]] <- tmp
  }

  list(chains = chains, pair = c(i, j), accepted = accepted)
}


#' Adapt temperature ladder based on swap acceptance rates
#' @keywords internal
.AdaptTemperatures <- function(betas, swapAccept, swapPropose,
                               target = 0.25) {
  nChains <- length(betas)
  if (nChains < 2L) return(betas)

  totalPropose <- sum(swapPropose)
  if (totalPropose < 20L) return(betas)

  overallRate <- sum(swapAccept) / totalPropose

  heat <- betas[nChains]
  adj <- exp(0.5 * (overallRate - target))
  heatNew <- heat^adj

  heatNew <- max(0.01, min(0.95, heatNew))

  .BuildTemperatureLadder(nChains, heatNew)
}


#' Parameter names for the sample matrix
#' @keywords internal
.ParamNames <- function(mkd, nEdge, kPrimePrior = "geometric",
                        qHeterogeneity = FALSE) {
  hasNeo <- any(mkd$type == "neomorphic")

  nms <- c("log_posterior", "log_likelihood", "tree_length")
  if (hasNeo) nms <- c(nms, "rate_loss")
  nms <- c(nms, "rate_log_sd")

  # kPrime hyperparameter columns depend on prior choice
  if (identical(kPrimePrior, "beta_geometric")) {
    nms <- c(nms, "kprime_alpha", "kprime_beta")
  } else if (!identical(kPrimePrior, "logseries")) {
    nms <- c(nms, "p")
  }

  if (hasNeo) {
    nms <- c(nms, "rate_neo")
  }

  # M-052: beta_scale column when Het is enabled
  if (isTRUE(qHeterogeneity)) {
    nms <- c(nms, "beta_scale")
  }

  # Diagnostic columns (always present from C++ batch)
  nms <- c(nms, "swap_cold", "topo_hash")

  transIdx <- which(mkd$type == "transformational")
  if (length(transIdx)) {
    nms <- c(nms, paste0("kPrime_", transIdx))
  }

  nms <- c(nms, paste0("br_", seq_len(nEdge)))

  nms
}


#' Extract state values to a row vector for storage
#' @keywords internal
#' Extract state row for sample storage
#' @param statePtr XPtr<McmcState> or R list (for backward compat)
#' @keywords internal
.StateToRow <- function(statePtr, mkd, nEdge, tipLabels = NULL,
                        kPrimePrior = "geometric",
                        qHeterogeneity = FALSE) {
  state <- get_mcmc_state(statePtr)

  hasNeo <- any(mkd$type == "neomorphic")
  rateLossVal <- if (hasNeo) state$rateLoss else numeric(0)
  rateNeoVal  <- if (hasNeo) state$rateNeo else numeric(0)

  # kPrime hyperparameter columns depend on prior choice
  kpHyperVal <- if (identical(kPrimePrior, "beta_geometric")) {
    c(state$kprimeAlpha, state$kprimeBeta)
  } else if (!identical(kPrimePrior, "logseries")) {
    state$p
  } else {
    numeric(0)
  }

  # M-052: beta_scale
  bsVal <- if (isTRUE(qHeterogeneity)) state$betaScale else numeric(0)

  transIdx <- which(mkd$type == "transformational")
  kp <- if (length(transIdx)) as.numeric(state$kPrime[transIdx]) else numeric(0)


  # Topology hash: FNV-1a of parent vector (shared C++ implementation)
  topoHash <- compute_topo_hash(state$edge[, 1])

  c(state$logPost, state$logLik, state$treeLength,
    rateLossVal, state$rateLogSd, kpHyperVal,
    rateNeoVal,
    bsVal,
    0,          # swap_cold: not applicable for R-side row extraction
    topoHash,
    kp,
    state$relBrLengths)
}



#' Reconstruct a phylo object from XPtr state
#' @keywords internal
.StateToTree <- function(statePtr, tipLabels) {
  state <- get_mcmc_state(statePtr)
  TreeTools::Preorder(structure(
    list(
      edge = state$edge,
      edge.length = state$treeLength * state$relBrLengths,
      Nnode = length(tipLabels) - 1L,
      tip.label = tipLabels
    ),
    class = "phylo"
  ))
}


# --- Adaptive move weight scheduler (M-092) ---

#' Resolve user-pinned weights for moves in the current pool.
#'
#' Returns a named numeric vector with entries only for moves present in
#' the pool, or NULL if no pinned weights apply.
#' @keywords internal
.ResolvePinnedWeights <- function(userMoveWeights, moveNames) {
  if (is.null(userMoveWeights)) return(NULL)
  keep <- intersect(names(userMoveWeights), moveNames)
  if (length(keep) == 0L) return(NULL)
  dropped <- setdiff(names(userMoveWeights), moveNames)
  if (length(dropped) > 0L) {
    cli::cli_warn(
      "Pinned move weight{?s} ignored (not in move pool): {.val {dropped}}."
    )
  }
  userMoveWeights[keep]
}


#' Apply pinned weights and renormalize free moves.
#' @keywords internal
.NormalizeMoveWeights <- function(weights, pinnedWeights) {
  moveNames <- names(weights)
  pinnedIdx <- match(names(pinnedWeights), moveNames)
  pinnedIdx <- pinnedIdx[!is.na(pinnedIdx)]
  if (length(pinnedIdx) == 0L) return(weights / sum(weights))

  budget <- 1.0 - sum(pinnedWeights)
  freeIdx <- setdiff(seq_along(weights), pinnedIdx)
  if (length(freeIdx) == 0L || budget < 1e-12) {
    # All pinned: normalize pinned to sum to 1
    weights[pinnedIdx] <- pinnedWeights / sum(pinnedWeights)
    if (length(freeIdx) > 0L) weights[freeIdx] <- 0
    return(weights)
  }
  freeSum <- sum(weights[freeIdx])
  if (freeSum > 0) {
    weights[freeIdx] <- weights[freeIdx] / freeSum * budget
  } else {
    weights[freeIdx] <- budget / length(freeIdx)
  }
  weights[pinnedIdx] <- pinnedWeights[names(weights)[pinnedIdx]]
  weights
}


#' Adapt move weights based on per-move acceptance rates and wall-time.
#'
#' Uses softmax reweighting with temperature annealing. Pinned moves
#' are excluded from adaptation. Moves with fewer than `minProposals`
#' proposals keep their current weight.
#'
#' The score formula accounts for multi-dimensional moves via the
#' `moveDim` parameter: `score = accept_rate * dim / cost`. For
#' single-parameter moves (`dim = 1`), this reduces to the original
#' `accept_rate / cost`. Block moves (e.g. block Gibbs branch sweep)
#' set `dim = nEdge` so the scheduler values them proportionally to
#' the number of parameters they update per call.
#'
#' @param currentWeights Numeric vector (current move probabilities,
#'   sums to 1).
#' @param acceptCount Named integer vector (cold chain, cumulative).
#' @param proposeCount Named integer vector (cold chain, cumulative).
#' @param moveTimeNs Named numeric vector (cold chain, cumulative ns).
#' @param moveNames Character vector of move names.
#' @param moveDim Integer vector of per-move dimensionality (number of
#'   parameters updated per call). Default: all 1s.
#' @param pinnedWeights Named numeric vector or NULL.
#' @param warmupProgress Fraction of warmup completed (0 to 1).
#' @param tStart Starting softmax temperature (default 2.0).
#' @param tEnd Ending softmax temperature (default 0.5).
#' @param wMin Floor per free move as fraction of 1 (default 0.01).
#' @param minProposals Minimum proposals before adapting (default 20).
#'
#' @return Updated weight vector (sums to 1).
#' @keywords internal
.AdaptMoveWeights <- function(currentWeights, acceptCount, proposeCount,
                               moveTimeNs, moveNames,
                               moveDim = rep(1L, length(currentWeights)),
                               pinnedWeights,
                               warmupProgress, tStart = 2.0, tEnd = 0.5,
                               wMin = 0.01, minProposals = 20L) {
  nMoves <- length(currentWeights)
  stopifnot(length(acceptCount) == nMoves,
            length(proposeCount) == nMoves,
            length(moveTimeNs) == nMoves)

  # Identify pinned and free moves
  pinnedIdx <- integer(0)
  if (!is.null(pinnedWeights)) {
    pinnedIdx <- match(names(pinnedWeights), moveNames)
    pinnedIdx <- pinnedIdx[!is.na(pinnedIdx)]
  }
  freeIdx <- setdiff(seq_len(nMoves), pinnedIdx)
  if (length(freeIdx) == 0L) return(currentWeights)

  budget <- if (length(pinnedIdx) > 0L) {
    1.0 - sum(pinnedWeights)
  } else {
    1.0
  }
  if (budget < 1e-12) return(currentWeights)

  # Identify scoreable free moves (enough proposals)
  scoreableIdx <- freeIdx[proposeCount[freeIdx] >= minProposals]

  if (length(scoreableIdx) == 0L) return(currentWeights)

  # Compute scores: dim-adjusted acceptances per second. log-scores avoid
  # overflow when accept_rate*dim/cost_s spans many orders of magnitude.
  acceptRate <- acceptCount[scoreableIdx] / proposeCount[scoreableIdx]
  dimAdj <- pmax(moveDim[scoreableIdx], 1L)
  meanCostS <- moveTimeNs[scoreableIdx] /
    (proposeCount[scoreableIdx] * 1e9)
  meanCostS <- pmax(meanCostS, 1e-9)
  logScores <- log(pmax(acceptRate, 1e-12)) + log(dimAdj) - log(meanCostS)

  # Softmax with annealed temperature
  temp <- tStart + (tEnd - tStart) * min(1, warmupProgress)
  scaledLogScores <- logScores / temp
  # Overflow guard
  scaledLogScores <- scaledLogScores - max(scaledLogScores)
  rawWeights <- exp(scaledLogScores)

  # Allocate budget: scoreable get softmax, non-scoreable keep current
  nonScoreableIdx <- setdiff(freeIdx, scoreableIdx)
  nonScoreableShare <- sum(currentWeights[nonScoreableIdx])

  scoreableBudget <- budget - nonScoreableShare
  if (scoreableBudget < 1e-12) return(currentWeights)

  softmaxWeights <- rawWeights / sum(rawWeights) * scoreableBudget

  # Apply floor: each free move gets at least wMin of the total budget
  nFree <- length(freeIdx)
  floorVal <- wMin
  nScoreable <- length(scoreableIdx)
  floored <- softmaxWeights < floorVal
  if (any(floored) && !all(floored)) {
    nFloored <- sum(floored)
    floorTotal <- nFloored * floorVal
    freeTotal <- scoreableBudget - floorTotal
    if (freeTotal > 0) {
      softmaxWeights[floored] <- floorVal
      softmaxWeights[!floored] <- softmaxWeights[!floored] /
        sum(softmaxWeights[!floored]) * freeTotal
    }
  }

  newWeights <- currentWeights
  newWeights[scoreableIdx] <- softmaxWeights
  # Enforce pinned values
  if (length(pinnedIdx) > 0L) {
    for (nm in names(pinnedWeights)) {
      idx <- match(nm, moveNames)
      if (!is.na(idx)) newWeights[idx] <- pinnedWeights[nm]
    }
  }
  newWeights
}


#' Format move weights as a compact string for display.
#'
#' Colour-coded via cli: high-weight moves are bright, low-weight moves
#' are dim, and names/equals are silver (M-140).
#' @keywords internal
#' Category definitions for move types
#' @keywords internal
.moveCategoryMap <- c(
  nni = "Topology", spr = "Topology", tbr = "Topology", pspr = "Topology",
  gibbs_spr = "Topology", gibbs_subtree_swap = "Topology",
  weighted_spr = "Topology", weighted_subtree_swap = "Topology",

  tree_length = "Branches", branch_lengths = "Branches",
  dirichlet_branch = "Branches", local_dirichlet = "Branches",
  block_gibbs_branch = "Branches", weighted_branch_lengths = "Branches",

  kPrime = "Characters", gibbs_kPrime = "Characters",
  block_kPrime = "Characters", p = "Characters",

  rate_loss = "Rates", rate_neo = "Rates", rate_log_sd = "Rates",
  beta_scale = "Rates", neo_joint = "Rates",
  slice_rate_loss = "Rates", slice_rate_neo = "Rates",
  slice_rate_log_sd = "Rates", slice_beta_scale = "Rates",
  joint_tl_rls = "Rates", joint_tl_rl = "Rates"
)

.moveCategoryOrder <- c("Topology", "Branches", "Characters", "Rates")

#' Format move weights as styled, categorized lines
#'
#' Returns a character vector (one element per category line).
#' Within each category, moves are sorted from highest to lowest weight.
#' Names are silver; values are green (>=10%), yellow (5-10%), white (<5%).
#' @keywords internal
.FormatMoveWeights <- function(weights, moveNames) {
  pct <- weights * 100
  cats <- .moveCategoryMap[moveNames]
  cats[is.na(cats)] <- "Other"

  presentCats <- intersect(.moveCategoryOrder, unique(cats))
  if (any(cats == "Other")) presentCats <- c(presentCats, "Other")

  vapply(presentCats, function(cat) {
    idx <- which(cats == cat)
    idx <- idx[order(pct[idx], decreasing = TRUE)]
    parts <- vapply(idx, function(i) {
      name <- cli::col_silver(paste0(moveNames[i], "="))
      val  <- sprintf("%.1f%%", pct[i])
      val  <- if (pct[i] >= 10) {
        cli::col_green(val)
      } else if (pct[i] >= 5) {
        cli::col_yellow(val)
      } else {
        cli::col_white(val)
      }
      paste0(name, val)
    }, character(1))
    paste0(cli::col_silver(paste0(cat, ": ")), paste(parts, collapse = " "))
  }, character(1), USE.NAMES = FALSE)
}

#' Format move weights as plain text (for log files).
#'
#' Categorized and sorted to match the styled version.
#' @keywords internal
.FormatMoveWeightsPlain <- function(weights, moveNames) {
  pct <- weights * 100
  cats <- .moveCategoryMap[moveNames]
  cats[is.na(cats)] <- "Other"

  presentCats <- intersect(.moveCategoryOrder, unique(cats))
  if (any(cats == "Other")) presentCats <- c(presentCats, "Other")

  lines <- vapply(presentCats, function(cat) {
    idx <- which(cats == cat)
    idx <- idx[order(pct[idx], decreasing = TRUE)]
    entries <- paste0(moveNames[idx], "=", sprintf("%.1f%%", pct[idx]))
    paste0(cat, ": ", paste(entries, collapse = " "))
  }, character(1), USE.NAMES = FALSE)
  paste(lines, collapse = "\n")
}

#' Print styled move weights to console (multi-line)
#' @keywords internal
.PrintMoveWeights <- function(weights, moveNames) {
  lines <- .FormatMoveWeights(weights, moveNames)
  cli::cli_alert_info("Move weights frozen:")
  for (line in lines) cli::cli_text("
 {line}")
}


#' Write adapted move weights as a comment in the log file.
#' @keywords internal
.LogMoveWeights <- function(weights, moveNames, logFilePaths) {
  plain <- .FormatMoveWeightsPlain(weights, moveNames)
  lines <- paste0("# ", strsplit(plain, "\n", fixed = TRUE)[[1]])
  block <- paste0(paste(lines, collapse = "\n"), "\n")
  for (p in logFilePaths) {
    cat(block, file = p, append = TRUE, sep = "")
  }
}


#' Adapt tuning parameters based on acceptance rates
#' @keywords internal
.AdaptTuning <- function(tuning, acceptCount, proposeCount, moves) {
  targets <- c(
    tree_length = 0.35, branch_lengths = 0.23,
    nni = 0.23, spr = 0.10,
    kPrime = 0.35,
    p = 0.35, rate_loss = 0.35, rate_log_sd = 0.35,
    rate_neo = 0.35, neo_joint = 0.35,
    beta_scale = 0.35,
    pspr = 0.10,
    dirichlet_branch = 0.234,
    local_dirichlet = 0.234,
    joint_tl_rls = 0.25, joint_tl_rl = 0.25,
    kprime_alpha = 0.35, kprime_beta = 0.35,
    # Gibbs/weighted/block/kPrime/slice moves: no MH tuning to adapt
    gibbs_kPrime = NA_real_, block_kPrime = 0.234,
    gibbs_spr = NA_real_, gibbs_subtree_swap = NA_real_,
    weighted_branch_lengths = NA_real_,
    weighted_spr = NA_real_, weighted_subtree_swap = NA_real_,
    block_gibbs_branch = NA_real_,
    slice_rate_loss = NA_real_, slice_rate_neo = NA_real_,
    slice_rate_log_sd = NA_real_, slice_tree_length = NA_real_,
    slice_beta_scale = NA_real_,
    slice_kprime_alpha = NA_real_, slice_kprime_beta = NA_real_
  )

  tuningKeys <- c(
    tree_length = "scale_tree_length",
    branch_lengths = "beta_simplex",
    nni = NA_character_,
    spr = NA_character_,
    kPrime = "int_walk_window",
    p = NA_character_,       # Gibbs move: no tuning needed
    rate_loss = "scale_rate_loss",
    rate_log_sd = "scale_rate_log_sd",
    rate_neo = "scale_rate_neo",
    neo_joint = "scale_neo_joint",
    pspr = NA_character_,
    joint_tl_rls = "scale_joint_tl_rls",
    joint_tl_rl = "scale_joint_tl_rl",
    dirichlet_branch = "dirichlet_alpha",
    local_dirichlet = "local_dirichlet_alpha",
    # Gibbs/weighted/block/kPrime/slice moves: no tuning to adapt
    gibbs_kPrime = NA_character_,
    block_kPrime = "int_walk_window",  # uses intWalkWindow for shift range
    gibbs_spr = NA_character_, gibbs_subtree_swap = NA_character_,
    weighted_branch_lengths = NA_character_,
    weighted_spr = NA_character_, weighted_subtree_swap = NA_character_,
    block_gibbs_branch = NA_character_,
    beta_scale = "scale_beta_scale",
    kprime_alpha = "scale_kprime_alpha",
    kprime_beta = "scale_kprime_beta",
    slice_rate_loss = NA_character_, slice_rate_neo = NA_character_,
    slice_rate_log_sd = NA_character_, slice_tree_length = NA_character_,
    slice_beta_scale = NA_character_
  )

  for (move in moves) {
    nm <- move$name
    if (proposeCount[nm] < 20) next
    rate <- acceptCount[nm] / proposeCount[nm]
    target <- targets[nm]
    tk <- tuningKeys[nm]

    if (!is.na(tk) && !is.null(tuning[[tk]])) {
      adj <- exp(0.5 * (rate - target))
      if (nm == "kPrime") {
        tuning[[tk]] <- max(1L, as.integer(round(tuning[[tk]] * adj)))
      } else if (nm == "branch_lengths") {
        tuning[[tk]] <- tuning[[tk]] / adj
        tuning[[tk]] <- max(2, tuning[[tk]])
      } else if (nm %in% c("dirichlet_branch", "local_dirichlet")) {
        # Inverted: higher alpha = tighter concentration = more conservative
        tuning[[tk]] <- tuning[[tk]] / adj
        tuning[[tk]] <- max(1.0, min(tuning[[tk]], 1000))
      } else {
        tuning[[tk]] <- tuning[[tk]] * adj
        tuning[[tk]] <- max(0.01, tuning[[tk]])
      }
    }
  }

  tuning
}


#' Adapt slice sampler widths based on stepping-out expansion counts
#'
#' Targets ~3 total expansions (left + right) per slice call. Fewer
#' expansions means the width is too wide; more means too narrow.
#'
#' @param tuning Per-chain tuning list (modified in place conceptually).
#' @param proposeCount Named integer vector of cumulative proposal counts.
#' @param sliceExpCount Named numeric vector of cumulative expansion counts.
#' @param moves List of move specifications.
#' @param target Target average expansions per slice call.
#' @return Updated tuning list.
#' @keywords internal
.AdaptSliceWidths <- function(tuning, proposeCount, sliceExpCount,
                               moves, target = 3.0) {
  sliceKeys <- c(
    slice_rate_loss   = "slice_width_rate_loss",
    slice_rate_neo    = "slice_width_rate_neo",
    slice_rate_log_sd = "slice_width_rate_log_sd",
    slice_tree_length = "slice_width_tree_length",
    slice_beta_scale  = "slice_width_beta_scale"
  )
  for (move in moves) {
    nm <- move$name
    tk <- sliceKeys[nm]
    if (is.na(tk) || is.null(tuning[[tk]])) next
    nProp <- proposeCount[nm]
    if (nProp < 10) next
    avgExp <- sliceExpCount[nm] / nProp
    ratio <- avgExp / target
    ratio <- max(0.25, min(ratio, 4.0))
    tuning[[tk]] <- tuning[[tk]] * ratio
    tuning[[tk]] <- max(0.05, min(tuning[[tk]], 10.0))
  }
  tuning
}


# --- Stabilisation detector ---

#' Check whether the MCMC chain has reached stationarity.
#'
#' Uses a Geweke-style z-score comparing recent and previous windows of
#' cold-chain log-posterior values. Returns `TRUE` when `|z| < zThreshold`
#' for `nStableRequired` consecutive checks.
#'
#' @param logPostHistory Numeric vector of log-posterior snapshots
#'   (one per warmup batch endpoint, chronological order).
#' @param nStableConsecutive Integer counter of consecutive stable checks
#'   so far (carried across calls).
#' @param windowSize Number of snapshots per comparison window.
#'   Default 10 (= 10 × 500 = 5000 iterations at default batch size).
#' @param zThreshold Absolute z-score threshold for declaring stability.
#'   Default 1.5.
#' @param nStableRequired Number of consecutive stable checks required.
#'   Default 3.
#'
#' @return A list with `stable` (logical) and `nStableConsecutive`
#'   (updated counter).
#' @keywords internal
.CheckStabilisation <- function(logPostHistory, nStableConsecutive,
                                 windowSize = 10L, zThreshold = 1.5,
                                 nStableRequired = 3L) {
  n <- length(logPostHistory)
  # Need at least 2 full windows

if (n < 2L * windowSize) {
    return(list(stable = FALSE, nStableConsecutive = 0L))
  }

  recent <- logPostHistory[(n - windowSize + 1L):n]
  prev   <- logPostHistory[(n - 2L * windowSize + 1L):(n - windowSize)]

  meanR <- mean(recent)
  meanP <- mean(prev)
  varR  <- var(recent)
  varP  <- var(prev)
  nR    <- length(recent)
  nP    <- length(prev)

  denom <- sqrt(varR / nR + varP / nP)
  # If both windows have zero variance, chain is flat → stable
  if (denom < .Machine$double.eps) {
    nStableConsecutive <- nStableConsecutive + 1L
  } else {
    z <- (meanR - meanP) / denom
    if (abs(z) < zThreshold) {
      nStableConsecutive <- nStableConsecutive + 1L
    } else {
      nStableConsecutive <- 0L
    }
  }

  list(stable = nStableConsecutive >= nStableRequired,
       nStableConsecutive = nStableConsecutive)
}


# --- Tuning-phase bandit (min-ESS/s optimisation) ---

#' Compute min-ESS/s for a set of samples collected over a known wall-time.
#'
#' @param sampleMatrix Numeric matrix (rows = samples, columns = parameters).
#' @param wallTimeSec Wall-clock seconds for the evaluation window.
#' @param excludePattern Regex pattern for column names to exclude from
#'   the min-ESS calculation (e.g. `"^kPrime_"`).
#'
#' @return Numeric scalar: min(ESS) / wallTimeSec, or `NA` if ESS
#'   cannot be computed.
#' @keywords internal
.MinEssPerSec <- function(sampleMatrix, wallTimeSec,
                           excludePattern = "^(kPrime_|br_|log_likelihood)",
                           tuningTrees = NULL) {
  if (nrow(sampleMatrix) < 10L || wallTimeSec < 1e-6) return(NA_real_)

  keyCols <- grep(excludePattern, colnames(sampleMatrix), invert = TRUE)
  if (length(keyCols) == 0L) return(NA_real_)

  ess <- .EssMatrix(sampleMatrix[, keyCols, drop = FALSE])

  minEss <- min(ess, na.rm = TRUE)
  if (!is.finite(minEss)) return(NA_real_)

  # Include tree ESS in the minimum when topology trees are available.
  # This gives topology moves credit in the bandit, preventing the
  # starvation that M-152 described.
  if (!is.null(tuningTrees) && length(tuningTrees) >= 20L) {
    trees <- structure(tuningTrees, class = "multiPhylo")
    treeEss <- tryCatch(
      .TreeESS(trees, dist_fn = TreeDist::RobinsonFoulds,
               frechet = FALSE)[["medianPseudoESS"]],
      error = function(e) NA_real_
    )
    if (!is.na(treeEss) && is.finite(treeEss)) {
      minEss <- min(minEss, treeEss)
    }
  }

  minEss / wallTimeSec
}


#' Estimate optimal thinning interval from observed autocorrelation
#'
#' Computes per-parameter integrated autocorrelation time (ACT) for key
#' scalar parameters and returns `max(nMoves, round(maxACT * log(2)))`.
#' The `log(2)` factor targets ~50% correlation between consecutive stored
#' samples.
#'
#' @param sampleMatrix Matrix of posterior samples (rows = draws, cols =
#'   parameters).
#' @param currentThin Current thinning interval (iterations per stored
#'   sample).
#' @param nMoves Number of active MCMC moves (floor for thinning).
#' @param excludePattern Regex for columns to exclude from ACT estimation.
#' @return Integer thinning interval.
#' @keywords internal
.AdaptThinning <- function(sampleMatrix, currentThin, nMoves,
                           excludePattern = "^(kPrime_|br_|log_likelihood)") {
  n <- nrow(sampleMatrix)
  if (n < 50L) return(currentThin)

  keyCols <- grep(excludePattern, colnames(sampleMatrix), invert = TRUE)
  if (length(keyCols) == 0L) return(currentThin)

  ess <- .EssMatrix(sampleMatrix[, keyCols, drop = FALSE])
  ess <- ess[is.finite(ess) & ess > 0]
  if (length(ess) == 0L) return(currentThin)

  # ACT in iterations for the worst-mixing parameter
  maxAct <- max(n / ess) * currentThin
  newThin <- as.integer(max(nMoves, round(maxAct * log(2))))

  # Cap: never more than 50x the move count
  newThin <- min(newThin, 50L * as.integer(nMoves))
  newThin
}


#' Generate perturbed move weight vectors.
#'
#' Produces `nPerturbations` candidate weight vectors by randomly
#' shifting one free (unpinned) move's weight by a small delta and
#' renormalising.
#'
#' @param currentWeights Named numeric vector (sums to 1).
#' @param pinnedWeights Named numeric vector or `NULL`.
#' @param moveNames Character vector of move names.
#' @param nPerturbations Number of candidates to generate.
#' @param deltaRange Numeric length-2 vector: range of absolute
#'   perturbation magnitude. Default `c(0.02, 0.10)`.
#' @param wMin Floor per free move. Default 0.01.
#'
#' @return A list of `nPerturbations` named numeric vectors.
#' @keywords internal
.PerturbMoveWeights <- function(currentWeights, pinnedWeights, moveNames,
                                 nPerturbations = 3L,
                                 deltaRange = c(0.02, 0.10),
                                 wMin = 0.01) {
  pinnedIdx <- integer(0)
  if (!is.null(pinnedWeights)) {
    pinnedIdx <- match(names(pinnedWeights), moveNames)
    pinnedIdx <- pinnedIdx[!is.na(pinnedIdx)]
  }
  freeIdx <- setdiff(seq_along(currentWeights), pinnedIdx)
  if (length(freeIdx) < 2L) {
    # Can't meaningfully perturb with fewer than 2 free moves
    return(list())
  }

  candidates <- vector("list", nPerturbations)
  for (i in seq_len(nPerturbations)) {
    w <- currentWeights
    # Pick a random free move to perturb
    target <- sample(freeIdx, 1L)
    delta <- runif(1, deltaRange[1], deltaRange[2]) * sample(c(-1, 1), 1)
    w[target] <- w[target] + delta

    # Enforce floor on free moves
    w[freeIdx] <- pmax(w[freeIdx], wMin)

    # Renormalise free moves to their budget
    budget <- if (length(pinnedIdx) > 0L) {
      1.0 - sum(pinnedWeights)
    } else {
      1.0
    }
    freeSum <- sum(w[freeIdx])
    if (freeSum > 0) {
      w[freeIdx] <- w[freeIdx] / freeSum * budget
    }
    # Restore pinned
    if (length(pinnedIdx) > 0L) {
      for (nm in names(pinnedWeights)) {
        idx <- match(nm, moveNames)
        if (!is.na(idx)) w[idx] <- pinnedWeights[nm]
      }
    }
    names(w) <- moveNames
    candidates[[i]] <- w
  }
  candidates
}
