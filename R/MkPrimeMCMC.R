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
    tuning = list()
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
         tuning = tuning),
    class = "MkPrimeMCMC"
  )
}
