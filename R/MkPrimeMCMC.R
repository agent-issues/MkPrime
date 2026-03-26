# MCMC configuration for MkPrime

#' Configure MCMC settings
#'
#' @param nIter Total iterations (including warmup). Default 10,000.
#' @param thin Thinning interval (save every `thin`-th iteration).
#'   Default 10.
#' @param warmup Number of warmup (adaptation) iterations. Default
#'   `nIter / 2`.
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
#' @return An S3 object of class `MkPrimeMCMC`.
#' @export
MkPrimeMCMC <- function(
    nIter = 10000L,
    thin = 10L,
    warmup = NULL,
    tuning = list()
) {
  nIter <- as.integer(nIter)
  thin <- as.integer(thin)
  if (is.null(warmup)) warmup <- as.integer(nIter / 2)
  warmup <- as.integer(warmup)

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
    list(nIter = nIter, thin = thin, warmup = warmup, tuning = tuning),
    class = "MkPrimeMCMC"
  )
}
