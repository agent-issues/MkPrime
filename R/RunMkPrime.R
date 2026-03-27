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
#' @param tree A `phylo` object (starting topology).
#' @param neomorphic,knownStates Passed to [MkPrimeData()] if `data` is
#'   a `phyDat` object.
#' @param model An `MkPrimeModel` object, or `NULL` for defaults.
#' @param mcmc An `MkPrimeMCMC` object, or `NULL` for defaults.
#' @param fixTopology Logical. If `TRUE`, tree topology is fixed (Phase 3
#'   behaviour). Default `FALSE` enables NNI and SPR topology proposals.
#'
#' @return An `MkPosterior` object.
#' @export
RunMkPrime <- function(data, tree,
                       neomorphic = integer(0),
                       knownStates = integer(0),
                       model = NULL,
                       mcmc = NULL,
                       fixTopology = FALSE) {

  # --- Input processing ---
  if (inherits(data, "MkPrimeData")) {
    mkd <- data
  } else {
    mkd <- MkPrimeData(data, neomorphic = neomorphic,
                       knownStates = knownStates)
  }

  if (!inherits(tree, "phylo")) {
    cli::cli_abort("{.arg tree} must be a {.cls phylo} object.")
  }
  if (is.null(tree$edge.length)) {
    cli::cli_abort(c(
      "{.arg tree} has no branch lengths.",
      "i" = "Supply a tree with edge lengths, e.g. \\
             {.code TreeTools::NJTree(data, edgeLengths = TRUE)}."
    ))
  }
  nNeg <- sum(tree$edge.length < 0)
  if (nNeg > 0L) {
    cli::cli_warn(c(
      "{nNeg} negative branch length{?s} clamped to 1e-8.",
      "i" = "Negative lengths arise in NJ trees when taxa are very similar. \\
             They are invalid for likelihood computation."
    ))
    tree$edge.length[tree$edge.length < 0] <- 1e-8
  }

  if (is.null(model)) model <- MkPrimeModel()
  if (is.null(mcmc)) mcmc <- MkPrimeMCMC()

  model <- .FinalizeModel(model, tree, mkd)

  # POSTORDER INVARIANT: all topology proposals maintain this ordering.
  tree <- TreeTools::Postorder(tree)
  nEdge <- nrow(tree$edge)
  tipLabels <- tree$tip.label

  hasNeo <- any(mkd$type == "neomorphic")
  transIdx <- which(mkd$type == "transformational")
  nTrans <- length(transIdx)

  moves <- .BuildMoves(nEdge, nTrans, hasNeo, mcmc,
                       fixTopology = fixTopology)

  nRuns <- mcmc$nRuns

  # --- Initialize per-run state ---
  runs <- vector("list", nRuns)
  for (run in seq_len(nRuns)) {
    startTree <- if (run == 1L) tree else .PerturbStart(tree)
    runs[[run]] <- .InitRun(startTree, mkd, model, mcmc, moves)
  }

  # --- Interleaved MCMC loop ---
  nSavedPerRun <- as.integer((mcmc$nIter - mcmc$warmup) / mcmc$thin)
  paramNames <- .ParamNames(mkd, nEdge)

  # Pre-allocate storage per run
  for (run in seq_len(nRuns)) {
    runs[[run]]$samples <- matrix(NA_real_, nrow = nSavedPerRun,
                                  ncol = length(paramNames),
                                  dimnames = list(NULL, paramNames))
    runs[[run]]$tree_samples <- vector("list", nSavedPerRun)
    runs[[run]]$saved_idx <- 0L
  }

  # Tree file
  treeFile <- mcmc$treeFile
  if (!is.null(treeFile)) writeLines("", treeFile)

  startTime <- proc.time()["elapsed"]
  hasProgressFn <- !is.null(mcmc$progressFn) && !is.null(mcmc$plotEvery)

  # Progress bar state
  coldLogpost <- {s <- get_mcmc_state(runs[[1]]$chainStates[[1]]); s$logPost}
  recentAcc <- 0
  convergeSummary <- ""

  # Constant inputs for the C++ batch function
  batchSize   <- 200L
  moveWeights <- vapply(moves, `[[`, numeric(1), "weight")
  moveTypeCodes <- vapply(moves, function(m) .kMoveTypes[[m$name]], integer(1L))
  transIdx0 <- if (nTrans > 0L) transIdx - 1L else integer(0L)

  # Initialize C++ data pointer (created once)
  mcmcData <- .InitMcmcData(mkd, model)

  # Populate per-partition log-likelihood cache (M-064) and pre-allocate
  # flat CL workspace to eliminate per-call heap allocations (M-063).
  for (run in seq_len(nRuns)) {
    for (ch in seq_len(mcmc$nChains)) {
      fill_partition_cache(mcmcData, runs[[run]]$chainStates[[ch]])
      allocate_cl_workspace(mcmcData, runs[[run]]$chainStates[[ch]])
    }
  }

  # Column indices in scalar_samples for tree reconstruction (1-based R)
  brColStart <- 6L + hasNeo + nTrans + 1L  # first br_ column

  cli::cli_progress_bar(
    "MCMC", total = mcmc$nIter,
    format = paste0(
      "{cli::pb_bar} {cli::pb_current}/{cli::pb_total}",
      " | logP: {format(round(coldLogpost, 1), nsmall = 1)}",
      " | accept: {format(round(recentAcc * 100, 0))}%",
      "{convergeSummary}"
    )
  )

  stopReason <- "max_iter"
  actualIter <- mcmc$nIter

  for (batchStart in seq(1L, mcmc$nIter, by = batchSize)) {
    batchEnd <- min(batchStart + batchSize - 1L, mcmc$nIter)
    nBatch   <- batchEnd - batchStart + 1L
    r1Result <- NULL

    for (run in seq_len(nRuns)) {
      r <- runs[[run]]

      scaleTunings <- .BuildScaleTuningMatrix(r$chain_tuning, moves)
      bsTunings <- vapply(r$chain_tuning,
                          function(t) t$beta_simplex, numeric(1L))
      iwWins <- vapply(r$chain_tuning,
                       function(t) as.integer(t$int_walk_window), integer(1L))

      result <- run_mcmc_batch_cpp(
        mcmcData, r$chainStates, r$betas,
        moveTypeCodes, transIdx0, moveWeights,
        scaleTunings, bsTunings, iwWins,
        nBatch, batchStart, mcmc$warmup, mcmc$thin,
        hasNeo, nEdge
      )
      if (run == 1L) r1Result <- result

      # Accumulate accept/propose counts (cumulative for adaptation)
      for (ch in seq_len(mcmc$nChains)) {
        r$chain_accept[[ch]]  <- r$chain_accept[[ch]]  +
          as.integer(result$accept_counts[ch, ])
        r$chain_propose[[ch]] <- r$chain_propose[[ch]] +
          as.integer(result$propose_counts[ch, ])
      }
      if (mcmc$nChains > 1L) {
        r$swap_accept  <- r$swap_accept  + result$swap_accept
        r$swap_propose <- r$swap_propose + result$swap_propose
      }

      # Store samples and reconstruct trees
      nSaved <- result$n_saved
      if (nSaved > 0L) {
        for (i in seq_len(nSaved)) {
          r$saved_idx <- r$saved_idx + 1L
          r$samples[r$saved_idx, ] <- result$scalar_samples[i, ]
          tl    <- result$scalar_samples[i, 3L]
          relBr <- result$scalar_samples[i, brColStart:(brColStart + nEdge - 1L)]
          curTree <- structure(
            list(edge        = result$edge_samples[[i]],
                 edge.length = tl * relBr,
                 Nnode       = as.integer(nEdge / 2L + 1L),
                 tip.label   = tipLabels),
            class = "phylo"
          )
          r$tree_samples[[r$saved_idx]] <- curTree
          if (!is.null(treeFile))
            cat(ape::write.tree(curTree), "\n", file = treeFile, append = TRUE)
        }
      }

      # Adapt tuning during warmup (at every batch boundary within warmup)
      if (batchEnd <= mcmc$warmup) {
        for (ch in seq_len(mcmc$nChains)) {
          r$chain_tuning[[ch]] <- .AdaptTuning(
            r$chain_tuning[[ch]], r$chain_accept[[ch]],
            r$chain_propose[[ch]], moves
          )
        }
        if (mcmc$nChains > 1L)
          r$betas <- .AdaptTemperatures(r$betas, r$swap_accept,
                                        r$swap_propose)
      }

      runs[[run]] <- r
    }

    # Progress update after each batch
    coldLogpost <- {s <- get_mcmc_state(runs[[1]]$chainStates[[1]]); s$logPost}
    if (!is.null(r1Result)) {
      batchAcc  <- sum(r1Result$accept_counts[1L, ])
      batchProp <- sum(r1Result$propose_counts[1L, ])
      if (batchProp > 0L) recentAcc <- batchAcc / batchProp
    }
    cli::cli_progress_update(set = batchEnd)

    # Progress callback (trace plots); trigger if batchEnd crosses a plotEvery multiple
    if (hasProgressFn &&
        (batchEnd %/% mcmc$plotEvery) > ((batchStart - 1L) %/% mcmc$plotEvery)) {
      info <- .BuildProgressInfo(runs, batchEnd, mcmc, startTime,
                                 recentAcc, paramNames)
      mcmc$progressFn(info)
    }

    # Stopping rule: max wall-clock time
    if (!is.null(mcmc$maxTime)) {
      if (proc.time()["elapsed"] - startTime >= mcmc$maxTime) {
        stopReason <- "max_time"
        actualIter <- batchEnd
        break
      }
    }

    # Convergence check + checkpoint; trigger if batchEnd crosses a checkEvery multiple
    doCheck <- batchEnd > mcmc$warmup && !is.null(mcmc$checkEvery) &&
               (batchEnd %/% mcmc$checkEvery) > ((batchStart - 1L) %/% mcmc$checkEvery)

    if (doCheck) {
      if (!is.null(mcmc$checkpointFile))
        .SaveCheckpoint(runs, mcmc, batchEnd, mcmc$checkpointFile)

      if (nRuns >= 2L) {
        diagCheck <- .CheckConvergence(runs, paramNames, mcmc)
        if (!is.null(diagCheck)) {
          convergeSummary <- sprintf(
            " | ESS: %d | PSRF: %.3f",
            as.integer(diagCheck$minEss), diagCheck$maxPsrf
          )
          if (diagCheck$converged) {
            stopReason <- "converged"
            actualIter <- batchEnd
            break
          }
        }
      }
    }
  }
  cli::cli_progress_done()

  # --- Build result ---
  .BuildResult(runs, model, mkd, mcmc, actualIter, stopReason)
}


# --- Run initialization ---

#' Initialize state for a single run
#' @keywords internal
.InitRun <- function(tree, mkd, model, mcmc, moves) {
  nChains <- mcmc$nChains
  betas <- .BuildTemperatureLadder(nChains, mcmc$heat)

  chains <- vector("list", nChains)
  for (ch in seq_len(nChains)) {
    chains[[ch]] <- .InitState(tree, mkd, model)
  }

  # Convert R states to C++ XPtr states
  chainStates <- vector("list", nChains)
  for (ch in seq_len(nChains)) {
    chainStates[[ch]] <- .InitMcmcChain(chains[[ch]])
  }

  moveNames <- vapply(moves, `[[`, character(1), "name")
  chainAccept <- chainPropose <- chainTuning <- vector("list", nChains)
  for (ch in seq_len(nChains)) {
    chainAccept[[ch]] <- integer(length(moves))
    chainPropose[[ch]] <- integer(length(moves))
    names(chainAccept[[ch]]) <- names(chainPropose[[ch]]) <- moveNames
    chainTuning[[ch]] <- mcmc$tuning
  }

  swapAccept <- swapPropose <- if (nChains > 1L) {
    integer(nChains - 1L)
  } else {
    integer(0)
  }

  list(
    chainStates = chainStates,
    betas = betas,
    chain_accept = chainAccept,
    chain_propose = chainPropose,
    chain_tuning = chainTuning,
    swap_accept = swapAccept,
    swap_propose = swapPropose
  )
}


# --- Convergence check during MCMC ---

#' Check convergence criteria (called during the loop)
#' @keywords internal
.CheckConvergence <- function(runs, paramNames, mcmc) {
  if (!requireNamespace("coda", quietly = TRUE)) return(NULL)

  keyCols <- .KeyParamCols(
    matrix(0, 1, length(paramNames), dimnames = list(NULL, paramNames))
  )

  # Gather saved samples from each run
  perRunSamples <- lapply(runs, function(r) {
    idx <- r$saved_idx
    if (idx < 10L) return(NULL)
    r$samples[seq_len(idx), keyCols, drop = FALSE]
  })

  if (any(vapply(perRunSamples, is.null, logical(1)))) return(NULL)

  # Compute PSRF
  chainList <- lapply(perRunSamples, function(s) coda::mcmc(s))
  mcmcList <- coda::mcmc.list(chainList)

  gd <- tryCatch(
    coda::gelman.diag(mcmcList, multivariate = FALSE),
    error = function(e) NULL
  )
  if (is.null(gd)) return(NULL)

  maxPsrf <- max(gd$psrf[, 1], na.rm = TRUE)

  # Compute min ESS across all runs combined
  combined <- do.call(rbind, perRunSamples)
  ess <- apply(combined, 2, function(col) {
    s <- sd(col, na.rm = TRUE); if (is.na(s) || s == 0) return(NA_real_)
    coda::effectiveSize(coda::mcmc(col))
  })
  minEss <- min(ess, na.rm = TRUE)

  converged <- TRUE
  if (!is.null(mcmc$minEss) && minEss < mcmc$minEss) {
    converged <- FALSE
  }
  if (!is.null(mcmc$maxPsrf) && maxPsrf > mcmc$maxPsrf) {
    converged <- FALSE
  }

  list(converged = converged, minEss = minEss, maxPsrf = maxPsrf)
}


# --- Build final result ---

#' Build MkPosterior from all runs
#' @keywords internal
.BuildResult <- function(runs, model, mkd, mcmc, actualIter, stopReason) {
  nRuns <- length(runs)

  # Trim samples to actual saved count
  for (run in seq_len(nRuns)) {
    idx <- runs[[run]]$saved_idx
    if (idx > 0L) {
      runs[[run]]$samples <- runs[[run]]$samples[seq_len(idx), , drop = FALSE]
      runs[[run]]$tree_samples <- runs[[run]]$tree_samples[seq_len(idx)]
    } else {
      runs[[run]]$samples <- runs[[run]]$samples[integer(0), , drop = FALSE]
      runs[[run]]$tree_samples <- list()
    }
  }

  # Per-run summaries
  perRunSummaries <- lapply(runs, function(r) {
    coldAcc <- r$chain_accept[[1]] / pmax(r$chain_propose[[1]], 1L)
    result <- list(
      samples = r$samples,
      trees = r$tree_samples,
      acceptance = coldAcc
    )
    if (mcmc$nChains > 1L) {
      result$betas <- r$betas
      result$swap_rates <- r$swap_accept / pmax(r$swap_propose, 1L)
    }
    result
  })

  if (nRuns == 1L) {
    r <- perRunSummaries[[1]]
    result <- MkPosterior(
      samples = r$samples,
      trees = r$trees,
      acceptance = r$acceptance,
      model = model,
      data = mkd,
      mcmc = mcmc,
      warmup = mcmc$warmup,
      tuning = runs[[1]]$chain_tuning[[1]]
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
    allTrees <- do.call(c, lapply(perRunSummaries, `[[`, "trees"))
    avgAcceptance <- Reduce(`+`, lapply(perRunSummaries, `[[`,
                                        "acceptance")) / nRuns

    result <- MkPosterior(
      samples = allSamples,
      trees = allTrees,
      acceptance = avgAcceptance,
      model = model,
      data = mkd,
      mcmc = mcmc,
      warmup = mcmc$warmup,
      tuning = runs[[1]]$chain_tuning[[1]]
    )

    result$nRuns <- nRuns
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

  result$stop_reason <- stopReason
  result$actual_iter <- actualIter
  result
}


# --- Checkpointing ---

#' Save MCMC checkpoint to RDS
#' @keywords internal
.SaveCheckpoint <- function(runs, mcmc, iter, file) {
  # XPtr<McmcState> objects cannot be serialized across R sessions.
  # Convert chain states to R lists via get_mcmc_state() before saving.
  serialRuns <- lapply(runs, function(r) {
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
    r$chainStates <- NULL  # never serialize raw XPtrs
    r
  })
  saveRDS(list(runs = serialRuns, mcmc = mcmc, iter = iter,
               timestamp = Sys.time(), version = 1L), file)
}


#' Resume MCMC from a checkpoint
#'
#' Loads a checkpoint file and continues the MCMC from where it left off.
#' The remaining iterations will be appended to the existing samples.
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

  if (is.null(checkpoint$version) || checkpoint$version != 1L) {
    cli::cli_abort("Unsupported checkpoint version.")
  }

  if (inherits(data, "MkPrimeData")) {
    mkd <- data
  } else {
    mkd <- MkPrimeData(data, neomorphic = neomorphic,
                       knownStates = knownStates)
  }

  if (is.null(model)) model <- MkPrimeModel()
  tree <- TreeTools::Postorder(tree)
  model <- .FinalizeModel(model, tree, mkd)

  runs <- checkpoint$runs
  mcmc <- checkpoint$mcmc
  startIter <- checkpoint$iter + 1L

  nRuns <- mcmc$nRuns
  # Infer nEdge from br_ columns
  nEdge <- sum(grepl("^br_", colnames(runs[[1]]$samples)))

  hasNeo <- any(mkd$type == "neomorphic")
  nTrans <- sum(mkd$type == "transformational")

  moves <- .BuildMoves(nEdge, nTrans, hasNeo, mcmc,
                       fixTopology = FALSE)
  paramNames <- colnames(runs[[1]]$samples)

  # Extend sample storage if needed
  totalSaved <- as.integer((mcmc$nIter - mcmc$warmup) / mcmc$thin)
  for (run in seq_len(nRuns)) {
    currentRows <- nrow(runs[[run]]$samples)
    if (currentRows < totalSaved) {
      extra <- matrix(NA_real_, nrow = totalSaved - currentRows,
                      ncol = ncol(runs[[run]]$samples),
                      dimnames = list(NULL, paramNames))
      runs[[run]]$samples <- rbind(runs[[run]]$samples, extra)
      runs[[run]]$tree_samples <- c(
        runs[[run]]$tree_samples,
        vector("list", totalSaved - length(runs[[run]]$tree_samples))
      )
    }
  }

  # Resume loop
  stopReason <- "max_iter"
  startTime <- proc.time()["elapsed"]
  hasProgressFn <- !is.null(mcmc$progressFn) && !is.null(mcmc$plotEvery)

  # Pre-compute move weight vector (constant across all iterations).
  moveWeights <- vapply(moves, `[[`, numeric(1), "weight")
  tipLabels <- tree$tip.label
  transIdx <- which(mkd$type == "transformational")
  mcmcData <- .InitMcmcData(mkd, model)

  # Rebuild XPtr<McmcState> from serialized R chain state lists.
  # Checkpoints store r$chains as plain R lists (see .SaveCheckpoint).
  for (run in seq_len(nRuns)) {
    nChains <- mcmc$nChains
    runs[[run]]$chainStates <- vector("list", nChains)
    for (ch in seq_len(nChains)) {
      ch_r <- runs[[run]]$chains[[ch]]
      runs[[run]]$chainStates[[ch]] <- init_mcmc_state(
        ch_r$edge[, 1], ch_r$edge[, 2],
        ch_r$rel_br_lengths, ch_r$tree_length,
        ch_r$rate_loss, ch_r$rate_log_sd,
        ch_r$rate_neo %||% 1.0, ch_r$p %||% 0.5,
        as.integer(ch_r$kPrime),
        ch_r$log_lik, ch_r$log_prior
      )
    }
  }

  # Populate per-partition log-likelihood cache (M-064) and pre-allocate
  # flat CL workspace (M-063) — must run after chainStates are rebuilt.
  for (run in seq_len(nRuns)) {
    for (ch in seq_len(mcmc$nChains)) {
      fill_partition_cache(mcmcData, runs[[run]]$chainStates[[ch]])
      allocate_cl_workspace(mcmcData, runs[[run]]$chainStates[[ch]])
    }
  }

  # Batch-based resume loop (mirrors RunMkPrime batch approach)
  batchSize     <- 200L
  moveTypeCodes <- vapply(moves, function(m) .kMoveTypes[[m$name]], integer(1L))
  transIdx0     <- if (length(transIdx) > 0L) transIdx - 1L else integer(0L)
  hasNeo        <- any(mkd$type == "neomorphic")
  brColStart    <- 6L + hasNeo + length(transIdx) + 1L

  cli::cli_progress_bar(
    "Resuming MCMC", total = mcmc$nIter - startIter + 1L,
    format = "{cli::pb_bar} {cli::pb_current}/{cli::pb_total} (from iter {startIter})"
  )

  stopReason <- "max_iter"
  actualIter <- mcmc$nIter

  for (batchStart in seq(startIter, mcmc$nIter, by = batchSize)) {
    batchEnd <- min(batchStart + batchSize - 1L, mcmc$nIter)
    nBatch   <- batchEnd - batchStart + 1L

    for (run in seq_len(nRuns)) {
      r <- runs[[run]]

      scaleTunings <- .BuildScaleTuningMatrix(r$chain_tuning, moves)
      bsTunings <- vapply(r$chain_tuning,
                          function(t) t$beta_simplex, numeric(1L))
      iwWins <- vapply(r$chain_tuning,
                       function(t) as.integer(t$int_walk_window), integer(1L))

      result <- run_mcmc_batch_cpp(
        mcmcData, r$chainStates, r$betas,
        moveTypeCodes, transIdx0, moveWeights,
        scaleTunings, bsTunings, iwWins,
        nBatch, batchStart, mcmc$warmup, mcmc$thin,
        hasNeo, nEdge
      )

      for (ch in seq_len(mcmc$nChains)) {
        r$chain_accept[[ch]]  <- r$chain_accept[[ch]]  +
          as.integer(result$accept_counts[ch, ])
        r$chain_propose[[ch]] <- r$chain_propose[[ch]] +
          as.integer(result$propose_counts[ch, ])
      }
      if (mcmc$nChains > 1L) {
        r$swap_accept  <- r$swap_accept  + result$swap_accept
        r$swap_propose <- r$swap_propose + result$swap_propose
      }

      nSaved <- result$n_saved
      if (nSaved > 0L) {
        for (i in seq_len(nSaved)) {
          r$saved_idx <- r$saved_idx + 1L
          r$samples[r$saved_idx, ] <- result$scalar_samples[i, ]
          tl    <- result$scalar_samples[i, 3L]
          relBr <- result$scalar_samples[i, brColStart:(brColStart + nEdge - 1L)]
          r$tree_samples[[r$saved_idx]] <- structure(
            list(edge        = result$edge_samples[[i]],
                 edge.length = tl * relBr,
                 Nnode       = as.integer(nEdge / 2L + 1L),
                 tip.label   = tipLabels),
            class = "phylo"
          )
        }
      }

      # Adapt tuning if still in warmup
      if (batchEnd <= mcmc$warmup) {
        for (ch in seq_len(mcmc$nChains)) {
          r$chain_tuning[[ch]] <- .AdaptTuning(
            r$chain_tuning[[ch]], r$chain_accept[[ch]],
            r$chain_propose[[ch]], moves
          )
        }
        if (mcmc$nChains > 1L)
          r$betas <- .AdaptTemperatures(r$betas, r$swap_accept,
                                        r$swap_propose)
      }

      runs[[run]] <- r
    }

    cli::cli_progress_update(set = batchEnd - startIter + 1L)

    if (hasProgressFn &&
        !is.null(mcmc$plotEvery) &&
        (batchEnd %/% mcmc$plotEvery) > ((batchStart - 1L) %/% mcmc$plotEvery)) {
      info <- .BuildProgressInfo(runs, batchEnd, mcmc, startTime,
                                 0, paramNames)
      mcmc$progressFn(info)
    }

    if (!is.null(mcmc$maxTime)) {
      if (proc.time()["elapsed"] - startTime >= mcmc$maxTime) {
        stopReason <- "max_time"
        actualIter <- batchEnd
        break
      }
    }
  }
  cli::cli_progress_done()

  .BuildResult(runs, model, mkd, mcmc, actualIter, stopReason)
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
                               recentAcc, paramNames) {
  nRuns <- length(runs)
  runSamples <- lapply(runs, function(r) {
    idx <- r$saved_idx
    if (idx > 0L) r$samples[seq_len(idx), , drop = FALSE] else NULL
  })

  currentState <- lapply(runs, function(r) {
    s <- get_mcmc_state(r$chainStates[[1]])
    list(log_lik = s$logLik, log_prior = s$logPrior,
         tree_length = s$treeLength, rate_loss = s$rateLoss,
         rate_log_sd = s$rateLogSd, p = s$p)
  })

  list(
    iter = iter,
    nIter = mcmc$nIter,
    warmup = mcmc$warmup,
    inWarmup = iter <= mcmc$warmup,
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
        p           = tun$scale_p,
        0.5
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
    kPrime = as.integer(kPrime),
    p = 0.5
  )

  # Partition rate scalar for neomorphic characters
  if (hasNeo) {
    state$rate_neo <- 1.0
  }

  # Tree is already postorder (reordered at init); use internal fast-path
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
                        fixTopology = FALSE) {
  moves <- list(
    list(name = "tree_length", type = "scale", target = "tree_length",
         weight = 1),
    list(name = "branch_lengths", type = "beta_simplex",
         target = "rel_br_lengths", weight = max(1, nEdge / 3))
  )

  if (!fixTopology && nEdge >= 5L) {
    moves <- c(moves, list(
      list(name = "nni", type = "nni", target = NULL,
           weight = max(1, nEdge / 2)),
      list(name = "spr", type = "spr", target = NULL,
           weight = max(1, nEdge / 4))
    ))
  }

  if (nTrans > 0) {
    moves <- c(moves, list(
      list(name = "kPrime", type = "int_walk", target = "kPrime",
           weight = max(1, 2 * nTrans)),
      list(name = "p", type = "scale", target = "p",
           weight = 1)
    ))
  }

  if (hasNeo) {
    moves <- c(moves, list(
      list(name = "rate_loss", type = "scale", target = "rate_loss",
           weight = 1.5),
      list(name = "rate_neo", type = "scale", target = "rate_neo",
           weight = 1)
    ))
  }

  moves <- c(moves, list(
    list(name = "rate_log_sd", type = "scale", target = "rate_log_sd",
         weight = 1.5)
  ))

  moves
}


# --- Move type integer codes (must match src/mcmc.cpp) ---
# 0=scale_tl, 1=scale_rl, 2=scale_rls, 3=scale_rn,
# 4=beta_simplex, 5=nni, 6=spr, 7=int_walk, 8=scale_p
.kMoveTypes <- c(
  tree_length = 0L, rate_loss = 1L, rate_log_sd = 2L,
  rate_neo = 3L, branch_lengths = 4L,
  nni = 5L, spr = 6L, kPrime = 7L, p = 8L
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
    model$kprimeHyperA, model$kprimeHyperB
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
    state$log_lik, state$log_prior
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
      p           = tuning$scale_p,
      0.5
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
.ParamNames <- function(mkd, nEdge) {
  nms <- c("log_posterior", "log_likelihood", "tree_length",
           "rate_loss", "rate_log_sd", "p")

  if (any(mkd$type == "neomorphic")) {
    nms <- c(nms, "rate_neo")
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
.StateToRow <- function(statePtr, mkd, nEdge, tipLabels = NULL) {
  state <- get_mcmc_state(statePtr)

  rateNeoVal <- if (!is.null(state$rateNeo) && state$rateNeo != 1.0) {
    state$rateNeo
  } else {
    numeric(0)
  }

  transIdx <- which(mkd$type == "transformational")
  kp <- if (length(transIdx)) as.numeric(state$kPrime[transIdx]) else numeric(0)

  c(state$logPost, state$logLik, state$treeLength,
    state$rateLoss, state$rateLogSd, state$p,
    rateNeoVal,
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
      Nnode = nrow(state$edge) / 2L + 1L,
      tip.label = tipLabels
    ),
    class = "phylo"
  )
}


#' Adapt tuning parameters based on acceptance rates
#' @keywords internal
.AdaptTuning <- function(tuning, acceptCount, proposeCount, moves) {
  targets <- c(
    tree_length = 0.35, branch_lengths = 0.23,
    nni = 0.23, spr = 0.10,
    kPrime = 0.35,
    p = 0.35, rate_loss = 0.35, rate_log_sd = 0.35,
    rate_neo = 0.35
  )

  tuningKeys <- c(
    tree_length = "scale_tree_length",
    branch_lengths = "beta_simplex",
    nni = NA_character_,
    spr = NA_character_,
    kPrime = "int_walk_window",
    p = "scale_p",
    rate_loss = "scale_rate_loss",
    rate_log_sd = "scale_rate_log_sd",
    rate_neo = "scale_rate_neo"
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
