# Native ESS and R-hat convergence diagnostics
#
# Implements rank-normalized split-chain R-hat and FFT-based effective
# sample size following:
#
#   Vehtari A, Gelman A, Simpson D, Carpenter B, Bürkner P-C (2021).
#   "Rank-Normalization, Folding, and Localization: An Improved R-hat
#   for Assessing Convergence of MCMC."
#   Bayesian Analysis, 16(2), 667-718. doi:10.1214/20-BA1221
#
#   Geyer CJ (1992). "Practical Markov Chain Monte Carlo."
#   Statistical Science, 7(4), 473-483.
#
# Written from scratch from the published algorithm descriptions.


# --- R-hat ----------------------------------------------------------------

#' Rank-normalized split-chain R-hat
#'
#' Computes the improved R-hat convergence diagnostic (Vehtari et al. 2021)
#' as `max(bulk R-hat, tail R-hat)`.  Supersedes the classical Gelman-Rubin
#' PSRF, which assumes approximate normality and ignores tail convergence.
#'
#' @param x Numeric matrix with dimensions `nIter x nChains`, or a numeric
#'   vector (treated as a single chain).
#' @return Scalar R-hat value (>= 1.0 if chains have converged).
#'   Returns `NA_real_` if fewer than 4 draws per split-half (< 8 per chain).
#' @references
#' \insertRef{Vehtari2021}{MkPrime}
#' @keywords internal
.Rhat <- function(x) {
  x <- .AsChainMatrix(x)
  splits <- .SplitChains(x)

  # Need >= 4 draws per split-chain for meaningful variance estimates
  if (nrow(splits) < 4L) return(NA_real_)

  bulkRhat <- .RhatClassical(.ZScale(splits))
  tailRhat <- .RhatClassical(.ZScale(.FoldDraws(splits)))
  max(bulkRhat, tailRhat)
}


#' Classical R-hat on a pre-transformed matrix
#'
#' Computes the standard between-chain / within-chain variance ratio.
#' Intended to be called on rank-normalized (z-scaled) draws.
#'
#' @param x Numeric matrix (nIter x nChains).
#' @return Scalar R-hat.
#' @noRd
.RhatClassical <- function(x) {
  if (.IsConstant(x)) return(NA_real_)
  nIter <- nrow(x)
  chainMeans <- colMeans(x)
  chainVars <- apply(x, 2, var)
  varBetween <- nIter * var(chainMeans)
  varWithin <- mean(chainVars)
  if (varWithin == 0) return(NA_real_)
  sqrt((varBetween / varWithin + nIter - 1) / nIter)
}


# --- Effective sample size ------------------------------------------------

#' Effective sample size (single chain)
#'
#' FFT-based autocovariance with Geyer's (1992) initial positive sequence
#' and initial monotone sequence estimators, plus the improved truncated
#' estimate from Vehtari et al. (2021).
#'
#' @param x Numeric vector of MCMC draws from one chain.
#' @return Scalar ESS estimate. Returns `NA_real_` for < 3 draws or
#'   constant input.
#' @references
#' \insertRef{Vehtari2021}{MkPrime}
#'
#' Geyer CJ (1992). "Practical Markov Chain Monte Carlo."
#' \emph{Statistical Science}, 7(4), 473--483.
#' @keywords internal
.Ess <- function(x) {
  x <- as.numeric(x)
  n <- length(x)
  if (n < 3L) return(NA_real_)
  if (.IsConstant(x)) return(NA_real_)

  acov <- .Autocovariance(x)
  meanVar <- acov[1] * n / (n - 1)
  varPlus <- meanVar * (n - 1) / n

  # Geyer's initial positive sequence
  rhoHat <- numeric(n)
  rhoHat[1] <- 1
  t <- 0L

  rhoEven <- 1
  rhoOdd <- 1 - (meanVar - acov[2]) / varPlus
  rhoHat[1] <- rhoEven
  rhoHat[2] <- rhoOdd

  while (t < n - 5L &&
         !is.nan(rhoEven + rhoOdd) &&
         (rhoEven + rhoOdd) > 0) {
    t <- t + 2L
    rhoEven <- 1 - (meanVar - acov[t + 1L]) / varPlus
    rhoOdd  <- 1 - (meanVar - acov[t + 2L]) / varPlus
    if ((rhoEven + rhoOdd) >= 0) {
      rhoHat[t + 1L] <- rhoEven
      rhoHat[t + 2L] <- rhoOdd
    }
  }
  maxT <- t

  # Save for improved estimate
  if (rhoEven > 0) rhoHat[maxT + 1L] <- rhoEven

  # Geyer's initial monotone sequence
  t <- 0L
  while (t <= maxT - 4L) {
    t <- t + 2L
    if (rhoHat[t + 1L] + rhoHat[t + 2L] >
        rhoHat[t - 1L] + rhoHat[t]) {
      rhoHat[t + 1L] <- (rhoHat[t - 1L] + rhoHat[t]) / 2
      rhoHat[t + 2L] <- rhoHat[t + 1L]
    }
  }

  # Improved truncated estimate (Vehtari et al. 2021)
  tauHat <- -1 + 2 * sum(rhoHat[seq_len(maxT)]) + rhoHat[maxT + 1L]

  ess <- n / tauHat

  # Safety: cap to avoid unstable estimates
  tauBound <- 1 / log10(n)
  if (tauHat < tauBound) {
    ess <- n * log10(n)
  }

  max(ess, 1)
}


#' Bulk effective sample size
#'
#' ESS for the bulk (central tendency) of the distribution.
#' Rank-normalizes pooled draws before computing ESS.
#'
#' @param x Numeric matrix (nIter x nChains) or vector.
#' @return Scalar bulk ESS.
#' @references
#' \insertRef{Vehtari2021}{MkPrime}
#' @keywords internal
.EssBulk <- function(x) {
  x <- .AsChainMatrix(x)
  splits <- .SplitChains(x)
  .EssMultiChain(.ZScale(splits))
}


#' Tail effective sample size
#'
#' ESS for the tails of the distribution.  Defined as the minimum of the
#' ESS for the 5% and 95% quantiles.
#'
#' @param x Numeric matrix (nIter x nChains) or vector.
#' @return Scalar tail ESS.
#' @references
#' \insertRef{Vehtari2021}{MkPrime}
#' @keywords internal
.EssTail <- function(x) {
  x <- .AsChainMatrix(x)
  splits <- .SplitChains(x)
  q05Ess <- .EssQuantile(splits, 0.05)
  q95Ess <- .EssQuantile(splits, 0.95)
  min(q05Ess, q95Ess)
}


#' ESS for a quantile
#'
#' Binarizes draws at the given quantile and computes ESS on the indicator.
#'
#' @param x Numeric matrix (nIter x nChains).
#' @param prob Quantile probability.
#' @return Scalar ESS.
#' @noRd
.EssQuantile <- function(x, prob) {
  if (.IsConstant(x)) return(NA_real_)
  indicator <- x <= quantile(x, prob)
  dim(indicator) <- dim(x)
  .EssMultiChain(indicator)
}


#' Multi-chain ESS
#'
#' Combines autocovariances from multiple chains using the split-chain
#' approach from Vehtari et al. (2021).
#'
#' @param x Numeric matrix (nIter x nChains).
#' @return Scalar ESS.
#' @noRd
.EssMultiChain <- function(x) {
  x <- as.matrix(x)
  nChains <- ncol(x)
  nIter <- nrow(x)

  if (nIter < 3L || .IsConstant(x)) return(NA_real_)

  # Per-chain autocovariances, averaged
  acovList <- apply(x, 2, .Autocovariance, simplify = FALSE)
  acovMeans <- Reduce(`+`, acovList) / nChains

  meanVar <- acovMeans[1] * nIter / (nIter - 1)
  varPlus <- meanVar * (nIter - 1) / nIter
  if (nChains > 1L) {
    varPlus <- varPlus + var(colMeans(x))
  }

  if (varPlus <= 0 || !is.finite(varPlus)) return(NA_real_)

  # Geyer's initial positive sequence
  rhoHat <- numeric(nIter)
  t <- 0L

  rhoEven <- 1
  rhoOdd <- 1 - (meanVar - acovMeans[2]) / varPlus
  rhoHat[1] <- rhoEven
  rhoHat[2] <- rhoOdd

  while (t < nIter - 5L &&
         !is.nan(rhoEven + rhoOdd) &&
         (rhoEven + rhoOdd) > 0) {
    t <- t + 2L
    rhoEven <- 1 - (meanVar - acovMeans[t + 1L]) / varPlus
    rhoOdd  <- 1 - (meanVar - acovMeans[t + 2L]) / varPlus
    if ((rhoEven + rhoOdd) >= 0) {
      rhoHat[t + 1L] <- rhoEven
      rhoHat[t + 2L] <- rhoOdd
    }
  }
  maxT <- t

  if (rhoEven > 0) rhoHat[maxT + 1L] <- rhoEven

  # Geyer's initial monotone sequence
  t <- 0L
  while (t <= maxT - 4L) {
    t <- t + 2L
    if (rhoHat[t + 1L] + rhoHat[t + 2L] >
        rhoHat[t - 1L] + rhoHat[t]) {
      rhoHat[t + 1L] <- (rhoHat[t - 1L] + rhoHat[t]) / 2
      rhoHat[t + 2L] <- rhoHat[t + 1L]
    }
  }

  # Improved truncated estimate
  tauHat <- -1 + 2 * sum(rhoHat[seq_len(maxT)]) + rhoHat[maxT + 1L]

  ess <- nChains * nIter / tauHat

  # Safety cap: consistent with .Ess() single-chain version
  totalN <- nChains * nIter
  tauBound <- 1 / log10(totalN)
  if (tauHat < tauBound) {
    ess <- totalN * log10(totalN)
  }

  max(ess, 1)
}


# --- Drop-in replacements for coda:: call sites --------------------------

#' Effective sample size for a single numeric vector
#'
#' Drop-in replacement for `coda::effectiveSize(coda::mcmc(x))`.
#' Guards against constant or degenerate input.
#'
#' @param x Numeric vector of MCMC draws.
#' @return Scalar ESS (at least 1).  Returns `NA_real_` for constant or
#'   very short input.
#' @keywords internal
.EssVector <- function(x) {
  n <- length(x)
  if (n < 2L) return(NA_real_)
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(NA_real_)
  .Ess(x)
}


#' Column-wise effective sample size for a matrix
#'
#' Drop-in replacement for `coda::effectiveSize(mat)`.
#'
#' @param mat Numeric matrix (rows = draws, columns = parameters).
#' @return Named numeric vector of ESS values.
#' @keywords internal
.EssMatrix <- function(mat) {
  mat <- as.matrix(mat)
  ess <- apply(mat, 2, .EssVector)
  names(ess) <- colnames(mat)
  ess
}


# --- Helpers --------------------------------------------------------------

#' FFT-based autocovariance
#'
#' Compute autocovariance at every lag using a zero-padded FFT.
#' Uses the "biased" (divide by N) normalization as recommended by
#' Geyer (1992).
#'
#' @param x Numeric vector.
#' @return Numeric vector of autocovariances (length `length(x)`).
#' @noRd
.Autocovariance <- function(x) {
  n <- length(x)
  vx <- var(x)
  if (vx == 0) return(rep(0, n))

  # Zero-pad to next efficient FFT size (powers of 2/3/5)
  m <- nextn(n)
  m2 <- 2L * m
  yc <- x - mean(x)
  yc <- c(yc, rep.int(0, m2 - n))

  # FFT-based unnormalized autocovariance
  ac <- Re(fft(Mod(fft(yc))^2, inverse = TRUE))[seq_len(n)]

  # Biased estimate (Geyer 1992 recommendation)
  ac <- ac / ac[1] * vx * (n - 1) / n
  ac
}


#' Rank-normalize draws via Blom transform
#'
#' Replaces each value with its rank, then applies the inverse normal
#' transform using Blom's (1958) fractional offset: `qnorm((r - 3/8) / (N + 1/4))`.
#'
#' @param x Numeric matrix or vector.
#' @return Numeric array of the same dimensions with rank-normalized values.
#' @noRd
.ZScale <- function(x) {
  r <- rank(x, ties.method = "average")
  n <- length(x)
  z <- qnorm((r - 3 / 8) / (n + 1 / 4))
  z[is.na(x)] <- NA
  if (!is.null(dim(x))) {
    dim(z) <- dim(x)
  }
  z
}


#' Fold draws around their median
#' @noRd
.FoldDraws <- function(x) {
  abs(x - median(x))
}


#' Split each chain in half
#'
#' Takes an nIter x nChains matrix and returns an (nIter/2) x (2*nChains)
#' matrix where each original chain contributes two split-halves.
#'
#' @param x Numeric matrix (nIter x nChains).
#' @return Numeric matrix with twice as many columns and half the rows.
#' @noRd
.SplitChains <- function(x) {
  x <- as.matrix(x)
  nIter <- nrow(x)
  half <- nIter %/% 2L
  if (half < 1L) return(x)
  # Odd nIter: drop the middle draw (Vehtari et al. 2021; `posterior`),
  # not the last one -- the two halves stay adjacent to the discarded draw.
  offset <- if (nIter %% 2L == 1L) 1L else 0L
  cbind(x[seq_len(half), , drop = FALSE],
        x[half + offset + seq_len(half), , drop = FALSE])
}


#' Coerce input to nIter x nChains matrix
#'
#' Vectors become single-column matrices.
#' @noRd
.AsChainMatrix <- function(x) {
  if (is.null(dim(x))) {
    matrix(x, ncol = 1L)
  } else {
    as.matrix(x)
  }
}


#' Check if input is effectively constant
#' @noRd
.IsConstant <- function(x, tol = .Machine$double.eps) {
  if (anyNA(x) || any(is.infinite(x))) return(TRUE)
  rng <- range(x)
  abs(rng[2] - rng[1]) < tol
}
