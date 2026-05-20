# Tree topology ESS methods
#
# frechetCorrelationESS() and medianPseudoESS() are adapted from the
# treess package by Andrew F. Magee (GPL-3).
#   Source: https://github.com/afmagee/treess
#   Paper:  Magee et al. (2021), arXiv:2109.07629
#
# Modifications from the original:
# - frechetCorrelationESS inner loop reimplemented in C++ for performance
# - Thin wrapper replaces treess::treess() dispatcher

#' Compute tree-topology ESS for one MCMC chain
#'
#' Computes tree-topology effective sample size for a sample of
#' phylogenetic trees using Robinson–Foulds distances.
#'
#' By default only the median pseudo-ESS is returned.  Set
#' `frechet = TRUE` to also compute the Fréchet correlation ESS
#' (Magee et al. 2021); this adds negligible cost since the full
#' distance matrix is computed either way.
#'
#' @param trees A `multiPhylo` list of trees from a single MCMC chain.
#' @param dist_fn Distance function applied to `trees`; must return a
#'   `dist` object.  Default: [TreeDist::RobinsonFoulds].
#' @param min_nsamples Integer; minimum number of samples used when
#'   computing lag-k statistics (default 5).
#' @param frechet Logical; if `TRUE`, also compute the Fréchet
#'   correlation ESS.  Default: `FALSE`.
#' @param maxRows Maximum rows for median pseudo-ESS subsampling
#'   (default 200).
#' @return Named numeric vector with elements `frechetCorrelationESS`
#'   and `medianPseudoESS`.  `frechetCorrelationESS` is `NA` when
#'   `frechet = FALSE`.
#' @references
#' Magee AF, Karcher MD, Matsen IV FA, Minin VN (2021).
#' "How trustworthy is your tree? Bayesian phylogenetic effective sample
#' size through the lens of Monte Carlo error."
#' \emph{arXiv:2109.07629}.
#' \doi{10.48550/arXiv.2109.07629}
#'
#' Lanfear R, Hua X, Warren DL (2016).
#' "Estimating the effective sample size of tree topologies from
#' Bayesian phylogenetic analyses."
#' \emph{Genome Biology and Evolution}, 8(8), 2319--2332.
#' @examples
#' if (requireNamespace("ape", quietly = TRUE) &&
#'     requireNamespace("TreeDist", quietly = TRUE)) {
#'   set.seed(1)
#'   trees <- lapply(seq_len(10), function(i) ape::rtree(8))
#'   class(trees) <- "multiPhylo"
#'   TreeESS(trees)
#' }
#' @export
TreeESS <- function(trees, dist_fn = TreeDist::RobinsonFoulds,
                     min_nsamples = 5L, frechet = FALSE,
                     maxRows = 200L) {
  dmat <- as.matrix(dist_fn(trees))
  c(
    frechetCorrelationESS = if (frechet) {
      .FrechetCorrelationESS(dmat, min_nsamples)
    } else {
      NA_real_
    },
    medianPseudoESS = .MedianPseudoESS(dmat, min_nsamples, maxRows)
  )
}

#' Fréchet correlation ESS
#'
#' Generalises the univariate autocorrelation ESS to non-Euclidean
#' spaces via Fréchet variance.  The inner loop is implemented in C++.
#'
#' @param dmat Numeric square distance matrix (not squared).
#' @param min_nsamples Minimum samples for lag computation.
#' @return Scalar ESS estimate.
#' @references Magee et al. (2021), arXiv:2109.07629.
#' @keywords internal
.FrechetCorrelationESS <- function(dmat, min_nsamples = 5L) {
  n <- nrow(dmat)
  if (n < min_nsamples + 2L) return(NA_real_)

  # All-zero matrix means a single unique topology; ESS = 1

if (all(dmat == 0)) return(1)

  dmat_sq <- dmat * dmat
  frechet_correlation_ess_cpp(dmat_sq, as.integer(min_nsamples))
}

#' Median pseudo-ESS (Lanfear et al. 2016)
#'
#' Treats each row of the distance matrix as a univariate time series
#' and returns the median of the per-row ESS estimates.
#'
#' Uses a C++ implementation of the Geyer (1992) initial-monotone-sequence
#' estimator instead of [coda::effectiveSize()] for performance.  For
#' large matrices, a deterministic subsample of rows is used (the median
#' stabilises well before all rows are evaluated).
#'
#' @param dmat Numeric square distance matrix.
#' @param min_nsamples Minimum samples for lag computation (default 5).
#' @param maxRows Maximum rows to evaluate; 0 = all rows.  Default 200.
#' @return Scalar ESS estimate.
#' @references
#' Lanfear R, Hua X, Warren DL (2016). \emph{Genome Biology
#'   and Evolution}, 8(8), 2319--2332.
#'
#' Geyer CJ (1992). "Practical Markov Chain Monte Carlo."
#'   \emph{Statistical Science}, 7(4), 473--483.
#' @keywords internal
.MedianPseudoESS <- function(dmat, min_nsamples = 5L, maxRows = 200L) {
  median_pseudo_ess_cpp(dmat, as.integer(min_nsamples), as.integer(maxRows))
}
