# MCMC configuration for MkPrime

#' Configure MCMC settings
#'
#' @param nIter Total iterations (including warmup). Default 10,000.
#' @param thin Thinning interval (save every `thin`-th iteration).
#'   Default 10.
#' @param warmup Number of warmup (adaptation) iterations. Default
#'   `nIter / 2`.
#' @param nChains Number of chains in the temperature ladder (parallel
#'   tempering). Default 1 (no tempering). Set to 4 for typical analyses.
#'   Chain 1 is the cold chain (beta = 1).
#' @param heat Temperature of the hottest chain (0 < heat < 1). Default
#'   0.2. The temperature ladder uses geometric spacing:
#'   `beta_i = heat^((i-1)/(nChains-1))` for i = 1, ..., nChains.
#'   Ignored when `nChains = 1`.
#' @param tree_file Path to write sampled trees in Newick format.
#'   `NULL` (default) disables file logging. Trees are always stored
#'   in the returned `MkPosterior` object regardless.
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
    nIter = 10000L,
    thin = 10L,
    warmup = NULL,
    nChains = 1L,
    heat = 0.2,
    tree_file = NULL,
    tuning = list()
) {
  nIter <- as.integer(nIter)
  thin <- as.integer(thin)
  if (is.null(warmup)) warmup <- as.integer(nIter / 2)
  warmup <- as.integer(warmup)
  nChains <- as.integer(nChains)

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
    int_walk_window = 1L
  )
  tuning <- modifyList(defaults, tuning)

  structure(
    list(nIter = nIter, thin = thin, warmup = warmup,
         nChains = nChains, heat = heat,
         tree_file = tree_file, tuning = tuning),
    class = "MkPrimeMCMC"
  )
}
