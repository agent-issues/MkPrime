# MCMC configuration for MkPrime

#' Configure MCMC settings
#'
#' @param nIter Total iterations (including warmup). Default `Inf`, which
#'   relies on stopping criteria (`minEss`, `maxPsrf`, `maxTime`) to
#'   terminate the run. Pass a finite integer to cap the number of iterations
#'   regardless of convergence.
#' @param thin Thinning interval (save every `thin`-th iteration).
#'   Default 10.
#' @param warmup Number of warmup (adaptation) iterations. Default
#'   `nIter / 2` for finite `nIter`, or 5,000 when `nIter = Inf`.
#' @param nRuns Number of independent runs. Default 2. Each run has its
#'   own set of `nChains` chains. Convergence diagnostics (PSRF) require
#'   `nRuns >= 2`.
#' @param nChains Number of chains in the temperature ladder (parallel
#'   tempering). Default 1 (no tempering). Set to 4 for typical analyses.
#'   Chain 1 is the cold chain (beta = 1).
#' @param heat Temperature of the hottest chain (0 < heat < 1). Default
#'   0.2. The temperature ladder uses geometric spacing:
#'   `beta_i = heat^((i-1)/(nChains-1))` for i = 1, ..., nChains.
#'   Ignored when `nChains = 1`.
#' @param maxTime Maximum wall-clock time in seconds. `NULL` (default)
#'   means no time limit.
#' @param minEss Minimum effective sample size for early stopping.
#'   `NULL` (default) disables ESS-based stopping.
#' @param maxPsrf Maximum PSRF for early stopping. `NULL` (default)
#'   disables PSRF-based stopping. Requires `nRuns >= 2`.
#' @param checkEvery Check convergence every this many iterations
#'   (default 1000). Only used when stopping criteria are set.
#' @param cancelFile Path to a cancel-signal file. `NULL` (default) disables
#'   cancel-file checking. When set, [RunMkPrime()] checks every 200
#'   iterations whether this file exists. If it does, the run flushes any
#'   buffered samples, saves a checkpoint (if `checkpointFile` is set), and
#'   exits with `stop_reason = "cancelled"`. Create the file to request a
#'   clean stop: `file.create(cancelFile)`. See also [MkCancelPath()].
#' @param checkpointFile Path to write checkpoint RDS files. `NULL`
#'   (default) disables checkpointing. Checkpoints are saved at each
#'   convergence check interval and on cancel.
#' @param treeFile Path to write sampled trees in Newick format.
#'   `NULL` (default) disables file logging. Trees are always stored
#'   in the returned `MkPosterior` object regardless.
#' @param logFile Path to write a tab-separated scalar parameter log
#'   (Tracer-compatible). `NULL` (default) keeps all samples in memory only.
#'   When set, samples are flushed to disk in batches of `bufferSize`, so the
#'   file can be opened in Tracer during the run.
#'   For multiple runs (`nRuns > 1`), separate files are created automatically
#'   by appending `_1`, `_2`, … before the file extension
#'   (e.g. `"run.log"` → `"run_1.log"`, `"run_2.log"`).
#'   See [ReadMkLog()] to load the log back into R after the run.
#' @param bufferSize Integer. Number of thinned samples to accumulate per run
#'   before flushing to `logFile`. Ignored when `logFile = NULL`. Default 500.
#'   Smaller values give more frequent flushes (lower Tracer latency, slightly
#'   more I/O); larger values amortise write overhead at the cost of slightly
#'   higher peak memory.
#' @param plotEvery Integer; invoke the progress callback every this many
#'   iterations. `NULL` (default) disables progress plotting (the `cli`
#'   progress bar still runs). Typical values: 100--500.
#' @param progressFn Progress callback function, the string `"default"`,
#'   or `NULL`. When `"default"`, uses [MkpTracePlot()] for live
#'   base-R trace plots. A custom function must accept a single argument:
#'   a named list with fields `iter`, `nIter`, `warmup`, `inWarmup`,
#'   `nRuns`, `nChains`, `runSamples`, `currentState`,
#'   `recentAcceptance`, and `elapsed`. See [MkpTracePlot()] for
#'   details.
#' @param gibbsSpr Logical; include the Gibbs SPR move (default `TRUE`).
#'   Cost: O(N) likelihood evaluations per proposal, where N is the number
#'   of candidate reattachment edges.
#' @param gibbsSubtreeSwap Logical; include the Gibbs subtree-swap move
#'   (default `TRUE`). Cost: O(N) likelihood evaluations.
#' @param weightedBranchScale Logical; include the weighted branch-length
#'   scale move (default `FALSE`). Cost: O(B) likelihood evaluations, where
#'   B = `nBranchBins`.
#' @param weightedSpr Logical; include the weighted SPR move (default
#'   `FALSE`). Cost: O(N * B) likelihood evaluations. Enable for improved
#'   mixing on difficult tree spaces at the cost of slower iterations.
#' @param weightedSubtreeSwap Logical; include the weighted subtree-swap
#'   move (default `FALSE`). Cost: O(N * B) likelihood evaluations.
#' @param nBranchBins Integer; number of branch-fraction bins for weighted
#'   moves (default `10L`). Used by `weightedBranchScale`, `weightedSpr`,
#'   and `weightedSubtreeSwap`. Ignored when all weighted moves are
#'   disabled. Higher values increase accuracy of the Gibbs approximation
#'   but cost more likelihood evaluations.
#' @param moveWeights Named numeric vector of user-pinned move weights,
#'   or `NULL` (default). When non-NULL, each named entry fixes the
#'   probability of proposing that move type. Names must match valid move
#'   names (e.g., `"nni"`, `"spr"`, `"gibbs_spr"`, `"tree_length"`, etc.).
#'   Values must be positive and sum to at most 1. Remaining probability
#'   is distributed among un-pinned moves by the adaptive scheduler during
#'   warmup. Example: `moveWeights = c(nni = 0.3, spr = 0.2)` fixes NNI
#'   at 30% and SPR at 20%, with the remaining 50% allocated adaptively
#'   among other moves. To disable adaptive scheduling entirely, pin all
#'   moves (sum to 1). See section **Adaptive move scheduling** below.
#' @param parallel Logical. If `TRUE` and `nRuns > 1`, independent runs are
#'   launched as non-blocking `future::future()` workers and the main process
#'   polls for convergence. Requires the \pkg{future} package (in `Suggests`).
#'   Set a parallel plan before calling [RunMkPrime()]:
#'   `future::plan("multisession", workers = nRuns)`. Default `FALSE` (sequential).
#' @param pollInterval Integer. Seconds between convergence polls in parallel
#'   mode. Ignored when `parallel = FALSE`. Default `10L`.
#' @param tuning Named list of initial tuning parameters for each move
#'   type. See Details.
#'
#' @details
#' Default tuning values:
#' - `scale_tree_length`: 0.5
#' - `beta_simplex`: 10
#' - `scale_rate_loss`: 0.5
#' - `scale_rate_log_sd`: 0.5
#' - `scale_p`: 0.5
#' - `int_walk_window`: 1
#'
#' ## Parallel tempering
#'
#' When `nChains > 1`, chains run at different temperatures. The cold
#' chain (beta = 1) targets the true posterior. Heated chains (beta < 1)
#' flatten the likelihood surface, aiding exploration. Chain swap proposals
#' (see M-030) exchange states between adjacent temperatures.
#'
#' ## Gibbs and weighted moves
#'
#' The Gibbs moves (`gibbsSpr`, `gibbsSubtreeSwap`) evaluate all candidate
#' topologies and sample proportional to their posterior weight. They cost
#' O(N) likelihood evaluations per proposal (where N is the number of
#' candidates, roughly the number of edges) but often achieve much better
#' mixing than standard NNI/SPR.
#'
#' The weighted moves (`weightedBranchScale`, `weightedSpr`,
#' `weightedSubtreeSwap`) additionally marginalise over a discrete grid
#' of branch-fraction placements. This produces proposals that are
#' approximately independent of the current state but costs
#' O(B) or O(N * B) likelihood evaluations. These are off by default
#' and recommended only when standard + Gibbs moves show poor mixing.
#'
#' ## Adaptive move scheduling
#'
#' During warmup, the MCMC engine tracks per-move acceptance rates and
#' wall-clock cost, then reweights the move pool every 200 iterations
#' to favor moves with high "acceptances per second" (a proxy for
#' ESS/wall-time efficiency). The reweighting uses softmax with a
#' temperature that anneals from 2.0 (near-uniform) to 0.5 (more
#' peaked) over warmup. At the end of warmup, weights are frozen to
#' preserve detailed balance. Use `moveWeights` to pin specific move
#' frequencies and exclude them from adaptation.
#'
#' @return An S3 object of class `MkPrimeMCMC`.
#' @export
MkPrimeMCMC <- function(
    nIter = Inf,
    thin = 10L,
    warmup = NULL,
    nRuns = 2L,
    nChains = 1L,
    heat = 0.2,
    maxTime = NULL,
    minEss = NULL,
    maxPsrf = NULL,
    checkEvery = 1000L,
    cancelFile = NULL,
    checkpointFile = NULL,
    treeFile = NULL,
    logFile = NULL,
    bufferSize = 500L,
    plotEvery = NULL,
    progressFn = NULL,
    gibbsSpr = TRUE,
    gibbsSubtreeSwap = TRUE,
    weightedBranchScale = FALSE,
    weightedSpr = FALSE,
    weightedSubtreeSwap = FALSE,
    nBranchBins = 10L,
    moveWeights = NULL,
    tuning = list(),
    parallel = FALSE,
    pollInterval = 10L
) {
  nIter <- if (is.infinite(nIter)) Inf else as.integer(nIter)
  thin <- as.integer(thin)
  if (is.null(warmup)) {
    warmup <- if (is.finite(nIter)) as.integer(nIter / 2L) else 5000L
  }
  warmup <- as.integer(warmup)
  nRuns <- as.integer(nRuns)
  nChains <- as.integer(nChains)

  if (nRuns < 1L) {
    cli::cli_abort("{.arg nRuns} must be at least 1.")
  }
  if (nChains < 1L) {
    cli::cli_abort("{.arg nChains} must be at least 1.")
  }
  if (nChains > 1L) {
    if (heat <= 0 || heat >= 1) {
      cli::cli_abort("{.arg heat} must be in (0, 1), got {heat}.")
    }
  }
  if (!is.null(cancelFile)) {
    if (!is.character(cancelFile) || length(cancelFile) != 1L) {
      cli::cli_abort("{.arg cancelFile} must be a length-1 character string or NULL.")
    }
  }
  if (!is.null(logFile)) {
    if (!is.character(logFile) || length(logFile) != 1L) {
      cli::cli_abort("{.arg logFile} must be a length-1 character string or NULL.")
    }
  }
  bufferSize <- as.integer(bufferSize)
  if (bufferSize < 1L) {
    cli::cli_abort("{.arg bufferSize} must be a positive integer.")
  }

  # Validate move toggles
  gibbsSpr <- as.logical(gibbsSpr)
  gibbsSubtreeSwap <- as.logical(gibbsSubtreeSwap)
  weightedBranchScale <- as.logical(weightedBranchScale)
  weightedSpr <- as.logical(weightedSpr)
  weightedSubtreeSwap <- as.logical(weightedSubtreeSwap)
  nBranchBins <- as.integer(nBranchBins)
  if (nBranchBins < 2L) {
    cli::cli_abort("{.arg nBranchBins} must be at least 2, got {nBranchBins}.")
  }

  # Validate moveWeights (M-092: adaptive scheduler)
  if (!is.null(moveWeights)) {
    if (is.list(moveWeights)) moveWeights <- unlist(moveWeights)
    if (!is.numeric(moveWeights) || is.null(names(moveWeights))) {
      cli::cli_abort(
        "{.arg moveWeights} must be a named numeric vector or NULL."
      )
    }
    validNames <- c(
      "tree_length", "branch_lengths", "nni", "spr", "kPrime", "p",
      "rate_loss", "rate_log_sd", "rate_neo",
      "gibbs_spr", "gibbs_subtree_swap",
      "weighted_branch_lengths", "weighted_spr", "weighted_subtree_swap"
    )
    bad <- setdiff(names(moveWeights), validNames)
    if (length(bad) > 0L) {
      cli::cli_abort(
        "{.arg moveWeights} contains unknown move name{?s}: {.val {bad}}."
      )
    }
    if (any(moveWeights <= 0)) {
      cli::cli_abort("All {.arg moveWeights} values must be positive.")
    }
    if (sum(moveWeights) > 1.0 + 1e-8) {
      cli::cli_abort(
        "{.arg moveWeights} sum to {sum(moveWeights)}, which exceeds 1.0."
      )
    }
  }

  defaults <- list(
    scale_tree_length = 0.5,
    beta_simplex = 10,
    scale_rate_loss = 0.5,
    scale_rate_log_sd = 0.5,
    scale_p = 0.5,
    scale_rate_neo = 0.5,
    int_walk_window = 1L
  )
  tuning <- modifyList(defaults, tuning)

  if (!is.null(checkEvery)) checkEvery <- as.integer(checkEvery)
  if (!is.null(plotEvery)) plotEvery <- as.integer(plotEvery)

  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel)) {
    cli::cli_abort("{.arg parallel} must be a length-1 logical (TRUE or FALSE).")
  }
  pollInterval <- as.integer(pollInterval)
  if (pollInterval < 1L) {
    cli::cli_abort("{.arg pollInterval} must be a positive integer.")
  }

  # Resolve progressFn
  if (identical(progressFn, "default")) {
    progressFn <- MkpTracePlot
  }
  if (!is.null(progressFn) && !is.function(progressFn)) {
    cli::cli_abort(
      "{.arg progressFn} must be a function, {.val default}, or {.val NULL}."
    )
  }

  structure(
    list(nIter = nIter, thin = thin, warmup = warmup,
         nRuns = nRuns, nChains = nChains, heat = heat,
         maxTime = maxTime, minEss = minEss, maxPsrf = maxPsrf,
         checkEvery = checkEvery, cancelFile = cancelFile,
         checkpointFile = checkpointFile,
         treeFile = treeFile, logFile = logFile, bufferSize = bufferSize,
         plotEvery = plotEvery, progressFn = progressFn,
         gibbsSpr = gibbsSpr, gibbsSubtreeSwap = gibbsSubtreeSwap,
         weightedBranchScale = weightedBranchScale,
         weightedSpr = weightedSpr,
         weightedSubtreeSwap = weightedSubtreeSwap,
         nBranchBins = nBranchBins,
         moveWeights = moveWeights,
         tuning = tuning,
         parallel = parallel, pollInterval = pollInterval),
    class = "MkPrimeMCMC"
  )
}
