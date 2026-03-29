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
#' @param fixTopology Logical. If `TRUE`, tree topology is fixed (Phase 3
#'   behaviour). Default `FALSE` enables NNI and SPR topology proposals.
#' @param overwrite Logical. If `FALSE` (the default) and
#'   `mcmc$checkpointFile` points to an existing file, the run is
#'   automatically resumed from that checkpoint.
#'   Set to `TRUE` to discard the existing checkpoint and start fresh.
#'
#' @return An `MkPosterior` object.
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
                       overwrite = FALSE) {

  # --- Auto-resume from checkpoint ---
  if (is.null(mcmc)) mcmc <- MkPrimeMCMC()
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
                       qHeterogeneity = qHet)

  if (identical(mcmc$thin, "auto")) {
    mcmc$thin <- length(moves)
  }

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
  isStreaming <- !is.null(mcmc$logFile)

  if (isStreaming) {
    convWindowSize <- .ComputeConvWindowSize(mcmc)
    logFilePaths   <- .OpenLogFiles(mcmc$logFile, paramNames, nRuns)
  } else {
    logFilePaths   <- NULL
    convWindowSize <- 0L
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

  # --- Parallel or sequential execution ---
  stopReason <- "max_iter"
  actualIter <- if (is.finite(mcmc$nIter)) mcmc$nIter else mcmc$warmup

  if (isTRUE(mcmc$parallel) && nRuns > 1L) {
    # Parallel: launch future workers; orchestrator polls for convergence.
    parResult    <- .RunParallelRuns(mkd, model, mcmc, runs, moves, tipLabels,
                                      paramNames, nEdge, brColStart, treeFile,
                                      isStreaming, logFilePaths, convWindowSize)
    runs         <- parResult$runs
    logFilePaths <- parResult$logFilePaths
    stopReason   <- parResult$stopReason
    actualIter   <- parResult$actualIter
    isStreaming   <- !is.null(logFilePaths)

    # Write tree samples to treeFile if set (workers passed treeFile = NULL
    # to avoid concurrent file writes; we flush from the collected states).
    if (!is.null(treeFile)) {
      for (r in runs) {
        for (tr in r$tree_samples) {
          if (!is.null(tr)) cat(ape::write.tree(tr), "\n",
                                file = treeFile, append = TRUE)
        }
      }
    }

    # Save combined checkpoint after parallel runs.  Workers pass
    # checkpointFile = NULL (no per-batch checkpointing from workers), so
    # this is the only checkpoint written for the parallel path.  Useful
    # for resuming a cancelled run (convergence, maxTime, or cancel file).
    if (!is.null(mcmc$checkpointFile)) {
      .SaveCheckpoint(runs, mcmc, actualIter, paramNames, mcmc$checkpointFile)
    }
  } else {
    # Sequential: run each run to completion before starting the next.
    # .RunMkPrimeSingleRun() accepts R-serializable state, reconstructs
    # C++ XPtrs internally, and returns serialized state — making each
    # call safe to replace with a future::future() worker (M-095 / M-096).
    for (run in seq_len(nRuns)) {
      runs[[run]] <- .RunMkPrimeSingleRun(
        mkd, model, mcmc, runs[[run]], moves, tipLabels, run,
        paramNames, nEdge, brColStart,
        logFilePath    = if (isStreaming) logFilePaths[run] else NULL,
        cancelFile     = mcmc$cancelFile,
        # Per-run checkpointing only for nRuns = 1; for nRuns > 1 the
        # combined checkpoint is saved below after all runs complete.
        checkpointFile = if (nRuns == 1L) mcmc$checkpointFile else NULL,
        startIter      = 1L,
        isStreaming    = isStreaming,
        convWindowSize = convWindowSize,
        treeFile       = treeFile
      )
      stopReason <- runs[[run]]$stop_reason
      actualIter <- runs[[run]]$actual_iter
      if (stopReason == "cancelled") break
    }

    # For nRuns > 1, save combined checkpoint at run granularity.
    if (nRuns > 1L && !is.null(mcmc$checkpointFile)) {
      .SaveCheckpoint(runs, mcmc, actualIter, paramNames, mcmc$checkpointFile)
    }
  }

  # --- Build result ---
  .BuildResult(runs, model, mkd, mcmc, paramNames, logFilePaths,
               actualIter, stopReason)
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
  chainTimeNs <- vector("list", nChains)
  for (ch in seq_len(nChains)) {
    chainAccept[[ch]] <- integer(length(moves))
    chainPropose[[ch]] <- integer(length(moves))
    names(chainAccept[[ch]]) <- names(chainPropose[[ch]]) <- moveNames
    chainTimeNs[[ch]] <- numeric(length(moves))
    names(chainTimeNs[[ch]]) <- moveNames
    chainTuning[[ch]] <- mcmc$tuning
  }

  swapAccept <- swapPropose <- if (nChains > 1L) {
    integer(nChains - 1L)
  } else {
    integer(0)
  }

  list(
    chains        = chains,
    betas         = betas,
    chain_accept  = chainAccept,
    chain_propose = chainPropose,
    chain_time_ns = chainTimeNs,
    chain_tuning  = chainTuning,
    swap_accept   = swapAccept,
    swap_propose  = swapPropose
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
                                  treeFile = NULL) {
  nChains <- mcmc$nChains

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
  savedIdx <- r$saved_idx %||% 0L
  nSavedPerRun <- if (is.finite(mcmc$nIter)) {
    as.integer((mcmc$nIter - mcmc$warmup) / mcmc$thin)
  } else {
    2000L
  }

  if (isStreaming) {
    # Clear stale streaming fields before merging fresh buffers (resume path).
    r$flush_idx <- NULL; r$flushed     <- NULL
    r$conv_head <- NULL; r$conv_filled <- NULL
    r$flush_buf <- NULL; r$flush_iter  <- NULL
    r$conv_window <- NULL
    bufs <- .InitStreamBuffers(length(paramNames), paramNames,
                               mcmc$bufferSize, convWindowSize)
    r <- c(r, bufs)
    r$saved_idx <- savedIdx
    if (is.null(r$tree_samples)) r$tree_samples <- vector("list", 0L)
  } else {
    if (is.null(r$samples)) {
      # Fresh run
      r$samples <- matrix(NA_real_, nrow = nSavedPerRun,
                          ncol = length(paramNames),
                          dimnames = list(NULL, paramNames))
      r$tree_samples <- vector("list", nSavedPerRun)
    } else if (is.finite(mcmc$nIter)) {
      # Resuming: extend matrix if needed
      currentRows <- nrow(r$samples)
      if (currentRows < nSavedPerRun) {
        extra <- matrix(NA_real_, nrow = nSavedPerRun - currentRows,
                        ncol = ncol(r$samples),
                        dimnames = list(NULL, colnames(r$samples)))
        r$samples <- rbind(r$samples, extra)
        r$tree_samples <- c(r$tree_samples,
                            vector("list", nSavedPerRun - length(r$tree_samples)))
      }
    }
    r$saved_idx <- savedIdx
  }

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
  moveDim       <- vapply(moves, function(m) m$dim %||% 1L, integer(1L))
  names(moveDim) <- moveNames
  moveTypeCodes <- vapply(moves, function(m) .kMoveTypes[[m$name]], integer(1L))
  # Slice param index: 0=treeLength, 1=rateLoss, 2=rateLogSd, 3=rateNeo, 4=betaScale
  sliceParamCodes <- vapply(moves, function(m) m$sliceParamIdx %||% 0L, integer(1L))
  transIdx      <- which(mkd$type == "transformational")
  transIdx0     <- if (length(transIdx) > 0L) transIdx - 1L else integer(0L)
  hasNeo        <- any(mkd$type == "neomorphic")

  # Auto-pin always-accept moves (Gibbs, slice) at initial weights.
  # The warmup scheduler's score (accept_rate × dim / cost) gives these
  # astronomical scores because acceptance = 1.0 and cost ≈ 0; this inflates
  # their weight and starves bottleneck MH moves.  One Gibbs draw or slice
  # sample per cycle is already optimal, so freeze them.
  alwaysAcceptTypes <- c("gibbs_p", "slice")
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
  tuningWindowStart <- NULL
  bestMinEssPerSec <- -Inf
  bestWeights      <- moveWeights
  tuningCandidates <- list()
  tuningCandIdx    <- 0L
  effectiveTuningBudget <- mcmc$tuningBudget

  # The C++ warmup parameter controls when samples are saved.
  # During Warmup: set to Inf so no samples saved.
  # During Tuning: set to 0 so all samples saved (collected into tuning buffer).
  # During Sample: set to 0 so all samples saved (collected into main storage).
  cppWarmup <- if (phase == "Warmup") mcmc$warmup else 0L

  weightsLogged <- phase == "Sample"

  # --- Progress bar (M-097 rotating ticker) ---
  coldLogpost    <- {s <- get_mcmc_state(r$chainStates[[1]]); s$logPost}
  recentAcc      <- 0
  batchEnd       <- startIter - 1L
  phaseLabel     <- phase
  progressLabel  <- if (startIter == 1L) "MCMC" else "Resuming MCMC"
  progressTotal  <- if (is.finite(mcmc$nIter)) {
    if (startIter == 1L) mcmc$nIter else mcmc$nIter - startIter + 1L
  } else NA

  # Ticker state: pages rotate every ~1.5 s wall-clock time.
  # Summary (minESS/PSRF) interleaved with detail (2 params each).
  tickerStart <- proc.time()["elapsed"]
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
    bsTunings    <- vapply(r$chain_tuning,
                           function(t) t$beta_simplex, numeric(1L))
    iwWins       <- vapply(r$chain_tuning,
                           function(t) as.integer(t$int_walk_window), integer(1L))

    result <- run_mcmc_batch_cpp(
      mcmcData, r$chainStates, r$betas,
      moveTypeCodes, transIdx0, sliceParamCodes, moveWeights,
      scaleTunings, bsTunings, iwWins, sliceWidths,
      nBatch, batchStart, cppWarmup, mcmc$thin,
      hasNeo, nEdge
    )

    # Accept/propose counts and timing (M-092)
    for (ch in seq_len(nChains)) {
      r$chain_accept[[ch]]  <- r$chain_accept[[ch]]  +
        as.integer(result$accept_counts[ch, ])
      r$chain_propose[[ch]] <- r$chain_propose[[ch]] +
        as.integer(result$propose_counts[ch, ])
      r$chain_time_ns[[ch]] <- r$chain_time_ns[[ch]] +
        as.numeric(result$move_time_ns[ch, ])
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
          if (r$saved_idx > length(r$tree_samples)) {
            n <- max(length(r$tree_samples), 1L)
            r$tree_samples <- c(r$tree_samples, vector("list", n))
          }
        } else {
          if (r$saved_idx > nrow(r$samples)) {
            n <- nrow(r$samples)
            extra <- matrix(NA_real_, nrow = n, ncol = ncol(r$samples),
                            dimnames = list(NULL, colnames(r$samples)))
            r$samples <- rbind(r$samples, extra)
            r$tree_samples <- c(r$tree_samples, vector("list", n))
          }
          r$samples[r$saved_idx, ] <- row
        }

        tl    <- row[3L]
        relBr <- row[brColStart:(brColStart + nEdge - 1L)]
        curTree <- structure(
          list(edge        = result$edge_samples[[i]],
               edge.length = tl * relBr,
               Nnode       = length(tipLabels) - 2L,
               tip.label   = tipLabels),
          class = "phylo"
        )
        r$tree_samples[[r$saved_idx]] <- curTree
        if (!is.null(treeFile))
          cat(ape::write.tree(curTree), "\n", file = treeFile, append = TRUE)
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
      }
    }

    # ===== PHASE-SPECIFIC LOGIC =====

    if (phase == "Warmup") {
      # --- Warmup: adapt tuning, temperatures, and move weights ---
      for (ch in seq_len(nChains)) {
        r$chain_tuning[[ch]] <- .AdaptTuning(
          r$chain_tuning[[ch]], r$chain_accept[[ch]],
          r$chain_propose[[ch]], moves
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

      # Stabilisation detection
      logPostHistory <- c(logPostHistory, coldLogpost)
      if (batchEnd >= mcmc$minWarmup) {
        stabResult <- .CheckStabilisation(
          logPostHistory, nStableConsecutive
        )
        nStableConsecutive <- stabResult$nStableConsecutive

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
            # Allocate tuning buffer
            tuningBufSize <- as.integer(effectiveTuningBudget / mcmc$thin) + 100L
            tuningBuf <- matrix(NA_real_, nrow = tuningBufSize,
                                ncol = length(paramNames),
                                dimnames = list(NULL, paramNames))
            tuningBufIdx     <- 0L
            tuningWindowStart <- proc.time()["elapsed"]
            # Reset acceptance/timing counters for clean tuning measurement
            for (ch in seq_len(nChains)) {
              r$chain_accept[[ch]][]  <- 0L
              r$chain_propose[[ch]][] <- 0L
              r$chain_time_ns[[ch]][] <- 0
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
            if (isStreaming && !is.null(logFilePath))
              .LogMoveWeights(moveWeights, moveNames, logFilePath)
            cli::cli_alert_info(
              "Move weights frozen: {(.FormatMoveWeights(moveWeights, moveNames))}"
            )
            weightsLogged <- TRUE
          }
        }
      }
    } else if (phase == "Tuning") {
      # --- Tuning: min-ESS/s perturbation bandit ---
      tuningIterUsed <- tuningIterUsed + nBatch

      # Evaluate current weight vector after each tuning window
      if (tuningBufIdx >= 10L) {
        windowTime <- proc.time()["elapsed"] - tuningWindowStart
        currentEssPerSec <- .MinEssPerSec(
          tuningBuf[seq_len(tuningBufIdx), , drop = FALSE],
          windowTime
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
          tuningWindowStart <- proc.time()["elapsed"]
          # Reset counters for clean measurement
          for (ch in seq_len(nChains)) {
            r$chain_accept[[ch]][]  <- 0L
            r$chain_propose[[ch]][] <- 0L
            r$chain_time_ns[[ch]][] <- 0
          }
        } else {
          # End of round: adopt best weights, start new round
          tuningRoundsDone <- tuningRoundsDone + 1L
          moveWeights <- bestWeights

          if (tuningRoundsDone >= mcmc$tuningRounds ||
              tuningIterUsed >= effectiveTuningBudget) {
            # Transition: Tuning → Sample
            phase      <- "Sample"
            r$phase    <- phase
            phaseLabel <- "Sample"
            r$samplePhaseStart <- batchEnd
            if (isStreaming && !is.null(logFilePath))
              .LogMoveWeights(moveWeights, moveNames, logFilePath)
            cli::cli_alert_info(
              "Move weights frozen: {(.FormatMoveWeights(moveWeights, moveNames))}"
            )
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
            tuningWindowStart <- proc.time()["elapsed"]
            bestMinEssPerSec  <- -Inf
            for (ch in seq_len(nChains)) {
              r$chain_accept[[ch]][]  <- 0L
              r$chain_propose[[ch]][] <- 0L
              r$chain_time_ns[[ch]][] <- 0
            }
          }
        }
      }
    }
    # Sample phase: no adaptation needed (weights frozen)

    # Streaming checkpoint: fire when buffer was flushed this batch
    if (phase == "Sample" && isStreaming &&
        !is.null(checkpointFile) && isTRUE(r$flushed)) {
      if (r$flush_idx > 0L) {
        .FlushBuffer(r$flush_buf, r$flush_idx, r$flush_iter, logFilePath)
        r$flush_idx <- 0L
      }
      r$flushed <- FALSE
      .SaveCheckpoint(list(r), mcmc, batchEnd, paramNames, checkpointFile,
                      moveWeights = moveWeights, phase = phase)
    }

    # Progress update (M-097 rotating ticker)
    coldLogpost <- {s <- get_mcmc_state(r$chainStates[[1]]); s$logPost}
    batchAcc  <- sum(result$accept_counts[1L, ])
    batchProp <- sum(result$propose_counts[1L, ])
    if (batchProp > 0L) recentAcc <- batchAcc / batchProp

    # Fixed-width logP: ratchet width up as magnitude grows, never shrink
    logPRaw   <- format(round(coldLogpost, 1), nsmall = 1)
    logPWidth <- max(logPWidth, nchar(logPRaw))
    logPStr   <- formatC(round(coldLogpost, 1), width = logPWidth,
                         format = "f", digits = 1)
    pageIdx   <- floor((proc.time()["elapsed"] - tickerStart) / 1.5) %%
                   length(tickerPages)
    # Dim separators; phase prefix silver for visual separation
    sep <- cli::col_silver("\u2502")
    tickerPage <- paste(
      cli::col_silver(paste(phaseLabel, batchEnd)),
      sep, sprintf("logP:%s", logPStr),
      sep, tickerPages[pageIdx + 1L]
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

    # Stopping: max wall-clock time
    if (!is.null(mcmc$maxTime) &&
        proc.time()["elapsed"] - startTime >= mcmc$maxTime) {
      stopReason <- "max_time"
      actualIter <- batchEnd
      break
    }

    # Stopping: cancel file
    if (!is.null(cancelFile) && file.exists(cancelFile)) {
      if (!is.null(checkpointFile)) {
        if (isStreaming && r$flush_idx > 0L) {
          .FlushBuffer(r$flush_buf, r$flush_idx, r$flush_iter, logFilePath)
          r$flush_idx <- 0L
          r$flushed   <- FALSE
        }
        .SaveCheckpoint(list(r), mcmc, batchEnd, paramNames, checkpointFile,
                        moveWeights = moveWeights, phase = phase)
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
                        moveWeights = moveWeights, phase = phase)
      }

      diagCheck <- .CheckConvergence(list(r), paramNames, mcmc, isStreaming)
      if (!is.null(diagCheck)) {
        # Refresh ticker pages from latest diagnostics (M-097)
        tickerPages <- .BuildTickerPages(diagCheck)
        if (diagCheck$converged) {
          stopReason <- "converged"
          actualIter <- batchEnd
          break
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
  r$chainStates <- NULL
  r$stop_reason <- stopReason
  r$actual_iter <- actualIter
  r$phase       <- phase
  r
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
      cli::cli_alert_info(
        "Parallel poll ({.field {.FormatElapsed(elapsed)}}): \\
         min ESS = {round(diagCheck$minEss)}"
      )
      if (diagCheck$converged) {
        for (cf in cancelFiles) file.create(cf)
        stopReason <- "converged"
        break
      }
    }

    # All workers finished naturally
    if (all(vapply(fList, future::resolved, logical(1L)))) break
  }

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
#' PSRF is computed only when `nRuns >= 2`. Returns full per-parameter `ess`
#' and `psrf` vectors so the caller can display a progress table.
#' @keywords internal
.CheckConvergence <- function(runs, paramNames, mcmc, isStreaming = FALSE) {
  if (!requireNamespace("coda", quietly = TRUE)) return(NULL)

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
  ess <- apply(combined, 2, function(col) {
    s <- sd(col, na.rm = TRUE)
    # Treat near-constant columns (FP noise only) as NA to avoid ESS = 0
    if (is.na(s) || s < sqrt(.Machine$double.eps) * (max(abs(col), na.rm = TRUE) + 1))
      return(NA_real_)
    coda::effectiveSize(coda::mcmc(col))
  })

  # kPrime are discrete nuisance parameters — exclude from convergence criteria
  # (M-098). They remain in the `ess` vector for display in .PrintProgressTable.
  isConvParam <- !grepl("^kPrime_", names(ess)) & names(ess) != "log_likelihood"
  minEss <- min(ess[isConvParam], na.rm = TRUE)

  # PSRF (requires >= 2 runs)
  psrf    <- NULL
  maxPsrf <- NA_real_
  if (nRuns >= 2L) {
    chainList <- lapply(perRunSamples, function(s) coda::mcmc(s))
    mcmcList  <- coda::mcmc.list(chainList)
    gd <- tryCatch(
      coda::gelman.diag(mcmcList, multivariate = FALSE),
      error = function(e) NULL
    )
    if (!is.null(gd)) {
      psrf    <- gd$psrf[, 1]
      maxPsrf <- max(psrf[isConvParam[names(psrf) %in% names(ess)]],
                     na.rm = TRUE)
    }
  }

  # Converged only when at least one criterion is set AND all set criteria pass.
  # (Avoids spurious early stopping when no criteria are configured.)
  hasCriteria <- !is.null(mcmc$minEss) || !is.null(mcmc$maxPsrf)
  converged   <- hasCriteria &&
    (is.null(mcmc$minEss)  || minEss >= mcmc$minEss) &&
    (is.null(mcmc$maxPsrf) || (nRuns >= 2L && !is.na(maxPsrf) && maxPsrf <= mcmc$maxPsrf))

  list(converged = converged, minEss = minEss, maxPsrf = maxPsrf,
       ess = ess, psrf = psrf)
}


#' Check convergence by reading log files from disk (parallel mode)
#'
#' Reads each run's log file via [ReadMkLog()], extracts key parameters,
#' and computes ESS (all runs combined) and PSRF (when `nRuns >= 2`).
#' Returns `NULL` if any log is missing or has fewer than 10 rows.
#' @keywords internal
.CheckConvergenceFromLogs <- function(logFilePaths, paramNames, mcmc) {
  if (!requireNamespace("coda", quietly = TRUE)) return(NULL)

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

  combined <- do.call(rbind, perRunSamples)
  ess <- apply(combined, 2, function(col) {
    s <- sd(col, na.rm = TRUE)
    if (is.na(s) || s < sqrt(.Machine$double.eps) * (max(abs(col), na.rm = TRUE) + 1))
      return(NA_real_)
    coda::effectiveSize(coda::mcmc(col))
  })

  # Exclude kPrime nuisance parameters from convergence criteria (M-098)
  isConvParam <- !grepl("^kPrime_", names(ess)) & names(ess) != "log_likelihood"
  minEss <- min(ess[isConvParam], na.rm = TRUE)

  psrf    <- NULL
  maxPsrf <- NA_real_
  if (nRuns >= 2L) {
    chainList <- lapply(perRunSamples, function(s) coda::mcmc(s))
    mcmcList  <- coda::mcmc.list(chainList)
    gd <- tryCatch(
      coda::gelman.diag(mcmcList, multivariate = FALSE),
      error = function(e) NULL
    )
    if (!is.null(gd)) {
      psrf    <- gd$psrf[, 1]
      maxPsrf <- max(psrf[isConvParam[names(psrf) %in% names(ess)]],
                     na.rm = TRUE)
    }
  }

  hasCriteria <- !is.null(mcmc$minEss) || !is.null(mcmc$maxPsrf)
  converged   <- hasCriteria &&
    (is.null(mcmc$minEss)  || minEss >= mcmc$minEss) &&
    (is.null(mcmc$maxPsrf) || (nRuns >= 2L && !is.na(maxPsrf) &&
                                maxPsrf <= mcmc$maxPsrf))

  list(converged = converged, minEss = minEss, maxPsrf = maxPsrf,
       ess = ess, psrf = psrf)
}


# --- Build final result ---

#' Build MkPosterior from all runs
#' @keywords internal
.BuildResult <- function(runs, model, mkd, mcmc, paramNames, logFilePaths,
                         actualIter, stopReason) {
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
    runs[[run]]$tree_samples <- runs[[run]]$tree_samples[seq_len(max(idx, 0L))]
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
    cli::cli_alert_info(c(
      "Streaming mode: {totalSaved} sample{?s} written to {.file {logFilePaths}}.",
      "i" = "Load with: {.code result$samples <- ReadMkLog(result$logFile)}"
    ))

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
                            moveWeights = NULL, phase = NULL) {
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
  if (!is.null(moveWeights)) payload$moveWeights <- moveWeights
  if (!is.null(phase))       payload$phase       <- phase
  saveRDS(payload, file)
}


# Flush any pending streaming buffers then save a checkpoint, if configured.
# Returns the (possibly modified) runs list so flush_idx resets propagate.
#
# @keywords internal
.FlushAndSaveCheckpoint <- function(runs, nRuns, mcmc, batchEnd,
                                    paramNames, isStreaming, logFilePaths,
                                    moveWeights = NULL) {
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
                    moveWeights = moveWeights)
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
#' @param model An `MkPrimeModel` object (must match original).
#' @param tree A `phylo` object (used only for model finalization).
#' @param neomorphic,knownStates Passed to [MkPrimeData()] if `data`
#'   is a `phyDat` object.
#'
#' @return An `MkPosterior` object with combined samples.
#' @export
ResumeMkPrime <- function(checkpointFile, data, tree,
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

  if (is.null(model)) model <- MkPrimeModel()
  tree <- TreeTools::Preorder(tree)

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

  model <- .FinalizeModel(model, tree, mkd)

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
    # Validate log files exist before committing to resume
    for (p in logFilePaths) {
      if (!file.exists(p)) {
        cli::cli_abort(c(
          "Log file not found: {.file {p}}.",
          "i" = "Cannot resume streaming run without its log file."
        ))
      }
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
                       qHeterogeneity = qHet)

  if (identical(mcmc$thin, "auto")) {
    mcmc$thin <- length(moves)
  }

  if (isStreaming) {
    # Rewind each log file to the checkpoint's saved_idx.  Any samples
    # flushed after the last checkpoint are discarded — the chain state
    # doesn't cover them.  Buffer reinit happens inside .RunMkPrimeSingleRun.
    for (run in seq_len(nRuns)) {
      .TruncateLogToN(logFilePaths[run],
                      as.integer(runs[[run]]$saved_idx %||% 0L))
    }
  }

  tipLabels   <- tree$tip.label
  isLogseries <- identical(model$kPrimePrior, "logseries")
  pCols       <- if (isLogseries) 0L else 1L
  brColStart  <- 5L + pCols + (any(mkd$type == "neomorphic")) + qHet + nTrans + 1L

  # --- Sequential per-run loop (same structure as RunMkPrime) ---
  stopReason <- "max_iter"
  actualIter <- if (is.finite(mcmc$nIter)) mcmc$nIter else startIter - 1L

  for (run in seq_len(nRuns)) {
    runs[[run]] <- .RunMkPrimeSingleRun(
      mkd, model, mcmc, runs[[run]], moves, tipLabels, run,
      paramNames, nEdge, brColStart,
      logFilePath    = if (isStreaming) logFilePaths[run] else NULL,
      cancelFile     = mcmc$cancelFile,
      checkpointFile = if (nRuns == 1L) mcmc$checkpointFile else NULL,
      startIter      = startIter,
      isStreaming    = isStreaming,
      convWindowSize = convWindowSize,
      treeFile       = NULL
    )
    stopReason <- runs[[run]]$stop_reason
    actualIter <- runs[[run]]$actual_iter
    if (stopReason == "cancelled") break
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

  # p hyperparameter only exists for the hierarchical geometric prior
  if (!identical(model$kPrimePrior, "logseries")) {
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
                        qHeterogeneity = FALSE) {
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

  if (nTrans > 0) {
    kPrimeMoves <- list(
      list(name = "kPrime", type = "int_walk", target = "kPrime",
           weight = max(1, 2 * nTrans), dim = 1L)
    )
    # p hyperparameter only exists for hierarchical geometric prior
    if (!identical(kPrimePrior, "logseries")) {
      kPrimeMoves <- c(kPrimeMoves, list(
        # Conjugate Gibbs draw: p | k' ~ Beta(a + nTrans, b + sum(k' - kObs))
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

  # --- Scalar weight floor ---
  # Scalar model-parameter moves (dim=1, non-topology) can be starved when
  # kPrime and branch_lengths dominate the weight budget. Guarantee each
  # scalar move gets at least 2% of the pre-floor total weight.
  scalarTypes <- c("scale", "int_walk", "gibbs_p", "scale_p", "slice")
  totalWeight <- sum(vapply(moves, `[[`, numeric(1), "weight"))
  floorVal <- totalWeight * 0.02
  for (i in seq_along(moves)) {
    m <- moves[[i]]
    if (m$dim == 1L && m$type %in% scalarTypes) {
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
# 15=block_gibbs_branch, 16=beta_scale (M-052), 17=tbr (M-053)
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
  pspr = 20L
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
    state$beta_scale %||% 1.0
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
      0.5  # default; gibbs_p ignores scaleTun (returns before using it)
    )
    accepted <- do_move_cpp(
      mcmcData, stateOrPtr, moveCode, charIdx,
      scaleTun, tuning$beta_simplex, tuning$int_walk_window, beta
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

  # p hyperparameter column only exists for hierarchical geometric prior
  if (!identical(kPrimePrior, "logseries")) {
    nms <- c(nms, "p")
  }

  if (hasNeo) {
    nms <- c(nms, "rate_neo")
  }

  # M-052: beta_scale column when Het is enabled
  if (isTRUE(qHeterogeneity)) {
    nms <- c(nms, "beta_scale")
  }

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

  # p only included in row when using hierarchical geometric prior
  pVal <- if (!identical(kPrimePrior, "logseries")) state$p else numeric(0)

  # M-052: beta_scale
  bsVal <- if (isTRUE(qHeterogeneity)) state$betaScale else numeric(0)

  transIdx <- which(mkd$type == "transformational")
  kp <- if (length(transIdx)) as.numeric(state$kPrime[transIdx]) else numeric(0)

  c(state$logPost, state$logLik, state$treeLength,
    rateLossVal, state$rateLogSd, pVal,
    rateNeoVal,
    bsVal,
    kp,
    state$relBrLengths)
}


#' Reconstruct a phylo object from XPtr state
#' @keywords internal
.StateToTree <- function(statePtr, tipLabels) {
  state <- get_mcmc_state(statePtr)
  structure(
    list(
      edge = state$edge,
      edge.length = state$treeLength * state$relBrLengths,
      Nnode = length(tipLabels) - 2L,
      tip.label = tipLabels
    ),
    class = "phylo"
  )
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
#' @keywords internal
.FormatMoveWeights <- function(weights, moveNames) {
  pct <- sprintf("%.1f%%", weights * 100)
  paste(paste0(moveNames, "=", pct), collapse = " ")
}


#' Write adapted move weights as a comment in the log file.
#' @keywords internal
.LogMoveWeights <- function(weights, moveNames, logFilePaths) {
  line <- paste0("# Adapted move weights: ",
                 .FormatMoveWeights(weights, moveNames))
  for (p in logFilePaths) {
    cat(line, "\n", file = p, append = TRUE, sep = "")
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
    # Gibbs/weighted/block/slice moves: no MH tuning to adapt
    gibbs_spr = NA_real_, gibbs_subtree_swap = NA_real_,
    weighted_branch_lengths = NA_real_,
    weighted_spr = NA_real_, weighted_subtree_swap = NA_real_,
    block_gibbs_branch = NA_real_,
    slice_rate_loss = NA_real_, slice_rate_neo = NA_real_,
    slice_rate_log_sd = NA_real_, slice_tree_length = NA_real_,
    slice_beta_scale = NA_real_
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
    # Gibbs/weighted/block/slice moves: no tuning to adapt
    gibbs_spr = NA_character_, gibbs_subtree_swap = NA_character_,
    weighted_branch_lengths = NA_character_,
    weighted_spr = NA_character_, weighted_subtree_swap = NA_character_,
    block_gibbs_branch = NA_character_,
    beta_scale = "scale_beta_scale",
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
      } else {
        tuning[[tk]] <- tuning[[tk]] * adj
        tuning[[tk]] <- max(0.01, tuning[[tk]])
      }
    }
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
                           excludePattern = "^(kPrime_|br_|log_likelihood)") {
  if (!requireNamespace("coda", quietly = TRUE)) return(NA_real_)
  if (nrow(sampleMatrix) < 10L || wallTimeSec < 1e-6) return(NA_real_)

  keyCols <- grep(excludePattern, colnames(sampleMatrix), invert = TRUE)
  if (length(keyCols) == 0L) return(NA_real_)

  ess <- apply(sampleMatrix[, keyCols, drop = FALSE], 2, function(col) {
    s <- sd(col, na.rm = TRUE)
    if (is.na(s) || s == 0) return(NA_real_)
    as.numeric(coda::effectiveSize(coda::mcmc(col)))
  })

  minEss <- min(ess, na.rm = TRUE)
  if (!is.finite(minEss)) return(NA_real_)
  minEss / wallTimeSec
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
