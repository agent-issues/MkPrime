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
#' @param checkpointFile Path to write checkpoint RDS files. `NULL`
#'   (default) disables checkpointing. Checkpoints are saved at each
#'   convergence check interval.
#' @param treeFile Path to write sampled trees in Newick format.
#'   `NULL` (default) disables file logging. Trees are always stored
#'   in the returned `MkPosterior` object regardless.
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
    checkpointFile = NULL,
    treeFile = NULL,
    plotEvery = NULL,
    progressFn = NULL,
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
         checkEvery = checkEvery, checkpointFile = checkpointFile,
         treeFile = treeFile, plotEvery = plotEvery,
         progressFn = progressFn, tuning = tuning),
    class = "MkPrimeMCMC"
  )
}
