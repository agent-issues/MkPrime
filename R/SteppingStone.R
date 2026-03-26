#' Marginal likelihood via stepping-stone sampling
#'
#' Estimates the marginal likelihood of the data under the Mk' model using
#' stepping-stone sampling (Xie et al., 2011). This runs MCMC at a series of
#' "power posteriors" interpolating between the prior (beta = 0) and the full
#' posterior (beta = 1), then combines the results into a marginal likelihood
#' estimate.
#'
#' @param data A `phyDat` object or `MkPrimeData` object.
#' @param tree A starting tree (`phylo` object).
#' @param neomorphic Integer vector of neomorphic character indices (if `data`
#'   is `phyDat`). Ignored if `data` is `MkPrimeData`.
#' @param model An `MkPrimeModel` object. Default: `MkPrimeModel()`.
#' @param nStones Number of stepping stones (power posterior levels). Default 50.
#' @param nIter Number of MCMC iterations per stone. Default 1000.
#' @param warmup Warmup iterations per stone (discarded). Default 200.
#' @param alpha Shape parameter for the Beta(alpha, 1) quantile schedule.
#'   Smaller values concentrate more stones near beta = 0 where the integrand
#'   changes fastest. Default 0.3 (Xie et al. recommendation).
#' @param fixTopology If `TRUE`, fix the tree topology and only sample
#'   continuous parameters + branch lengths. Default `FALSE`.
#' @param verbose Print progress? Default `TRUE`.
#'
#' @return A list with components:
#'   \describe{
#'     \item{log_marginal}{Log marginal likelihood estimate.}
#'     \item{se}{Standard error of the estimate (delta-method approximation;
#'       see Details).}
#'     \item{log_ratios}{Per-stone log ratios (length `nStones`).}
#'     \item{betas}{The beta schedule used.}
#'   }
#'
#' @details
#' The standard error is computed via the delta method applied to each
#' stepping-stone ratio, following Xie et al. (2011) and the implementation in
#' the \pkg{mcmc3r} package (dos Reis). For each stone \eqn{k},
#' the importance weights \eqn{L_{k,i} = \exp(\Delta\beta_k \cdot \ell_i - C_k)}
#' are computed (where \eqn{C_k} is a centering constant for numerical
#' stability). The per-stone ratio is \eqn{r_k = \bar{L}_k} and its log
#' contributes to the overall log-marginal likelihood. By the delta method,
#' \deqn{\mathrm{Var}(\log r_k) \approx \frac{\mathrm{Var}(L_k)}{N_{\mathrm{eff},k}
#'   \cdot r_k^2}}
#' where \eqn{N_{\mathrm{eff},k}} is the effective sample size (accounting
#' for MCMC autocorrelation). The total SE is
#' \eqn{\sqrt{\sum_k \mathrm{Var}(\log r_k)}}. A warning is issued for any
#' stone where \eqn{\mathrm{Var}(L_k) / (N_{\mathrm{eff},k} \cdot r_k^2) > 0.1},
#' indicating the delta approximation may be unreliable; consider increasing
#' \code{nIter} or \code{nStones}.
#'
#' If the \pkg{coda} package is not available, a fallback ESS estimate is used
#' (initial positive sequence estimator on the autocorrelation function).
#'
#' @references
#' Xie, W., Lewis, P.O., Fan, Y., Kuo, L., & Chen, M.-H. (2011).
#' Improving marginal likelihood estimation for Bayesian phylogenetic model
#' selection. *Systematic Biology*, 60(2), 150--160.
#'
#' @export
mkp_stepping_stone <- function(data, tree,
                              neomorphic = integer(0),
                              model = NULL,
                              nStones = 50L,
                              nIter = 1000L,
                              warmup = 200L,
                              alpha = 0.3,
                              fixTopology = FALSE,
                              verbose = TRUE) {

  # --- Input processing ---
  if (inherits(data, "MkPrimeData")) {
    mkd <- data
  } else {
    mkd <- MkPrimeData(data, neomorphic = neomorphic)
  }

  if (is.null(model)) model <- MkPrimeModel()
  model <- .FinalizeModel(model, tree, mkd)

  nStones <- as.integer(nStones)
  nIter <- as.integer(nIter)
  warmup <- as.integer(warmup)

  # --- Beta schedule: quantiles of Beta(alpha, 1) ---
  # beta_k = (k/K)^(1/alpha), k = 0, ..., K
  # This concentrates stones near beta = 0
  betas <- ((seq_len(nStones + 1L) - 1L) / nStones)^(1 / alpha)
  # betas[1] = 0 (prior only), betas[nStones + 1] = 1 (full posterior)

  # --- Initialize MCMC state ---
  # POSTORDER INVARIANT: .MkpLogLikelihood() (internal fast-path) requires
  # edges in postorder. Reorder here so .InitState() and all subsequent
  # topology proposals maintain the invariant.
  tree <- TreeTools::Postorder(tree)
  state <- .InitState(tree, mkd, model)

  nEdge <- nrow(tree$edge)
  nTrans <- sum(mkd$type == "transformational")
  hasNeo <- any(mkd$type == "neomorphic")
  moves <- .BuildMoves(nEdge, nTrans, hasNeo, NULL,
                       fixTopology = fixTopology)
  moveWeights <- vapply(moves, `[[`, numeric(1), "weight")

  mcmcTuning <- MkPrimeMCMC(nIter = 100L)$tuning

  # Build C++ data struct and XPtr state (shared across all stones)
  mcmcData <- .InitMcmcData(mkd, model)
  transIdx  <- which(mkd$type == "transformational")
  statePtr  <- .InitMcmcChain(state)

  logRatios <- numeric(nStones)

  # Per-stone importance weights (centered) and their means,
  # retained for delta-method SE calculation
  stoneWeights <- vector("list", nStones)
  stoneMeans <- numeric(nStones)

  for (stone in seq_len(nStones)) {
    betaLo <- betas[stone]
    betaHi <- betas[stone + 1L]

    if (verbose) {
      cli::cli_progress_step(
        "Stone {stone}/{nStones}: beta = [{round(betaLo, 4)}, {round(betaHi, 4)}]"
      )
    }

    # Run MCMC at power posterior betaLo
    logLiks <- numeric(nIter)

    for (iter in seq_len(warmup + nIter)) {
      moveIdx <- sample.int(length(moves), 1L, prob = moveWeights)
      move <- moves[[moveIdx]]
      .DoMove(move, statePtr, tuning = mcmcTuning,
              beta = betaLo, transIdx = transIdx, mcmcData = mcmcData)

      if (iter > warmup) {
        logLiks[iter - warmup] <- get_state_log_lik(statePtr)
      }
    }

    # Stepping-stone ratio for this interval:
    # log r_k = log( (1/n) * sum_i exp((betaHi - betaLo) * logLik_i) )
    # Use log-sum-exp for numerical stability
    deltaBeta <- betaHi - betaLo
    # Remove -Inf logLiks (can occur when sampling near the prior; those
    # samples contribute zero importance weight and are numerically harmless
    # if at least one finite sample exists, but cause NaN via -Inf - (-Inf)
    # if ALL logLiks are -Inf)
    finiteLogLiks <- logLiks[is.finite(logLiks)]
    if (length(finiteLogLiks) == 0L) {
      # Degenerate stone — flag as unreliable; set weights to 0
      logRatios[stone] <- NaN
      stoneWeights[[stone]] <- numeric(0)
      stoneMeans[stone] <- NaN
    } else {
      shifted <- deltaBeta * logLiks
      maxShifted <- max(shifted, na.rm = TRUE)
      Lk <- exp(shifted - maxShifted) # centered importance weights
      Lk[!is.finite(Lk)] <- 0         # -Inf logLiks → weight 0
      rk <- mean(Lk)
      logRatios[stone] <- maxShifted + log(rk)
      stoneWeights[[stone]] <- Lk
      stoneMeans[stone] <- rk
    }
  }

  logMarginal <- sum(logRatios)

  # --- SE via delta method (Xie et al. 2011; mcmc3r implementation) ---
  # Var(log r_k) ≈ Var(L_k) / (N_eff * r_k^2)
  varLogRatios <- numeric(nStones)
  unreliable <- character(0)

  for (stone in seq_len(nStones)) {
    Lk <- stoneWeights[[stone]]
    rk <- stoneMeans[stone]

    # Degenerate stone: no finite logLiks were available
    if (length(Lk) == 0L || !is.finite(rk)) {
      varLogRatios[stone] <- NA_real_
      unreliable <- c(unreliable, paste0(stone, " (beta=",
                                          round(betas[stone], 4), ", degenerate)"))
      next
    }

    essK <- .EssVector(Lk)
    vzrK <- var(Lk) / essK
    ratio <- vzrK / rk^2
    varLogRatios[stone] <- ratio

    if (is.finite(ratio) && ratio > 0.1) {
      unreliable <- c(unreliable, paste0(stone, " (beta=",
                                          round(betas[stone], 4), ")"))
    } else if (!is.finite(ratio)) {
      unreliable <- c(unreliable, paste0(stone, " (beta=",
                                          round(betas[stone], 4), ", degenerate)"))
    }
  }

  se <- sqrt(sum(varLogRatios, na.rm = TRUE))

  if (length(unreliable) > 0L && verbose) {
    cli::cli_warn(c(
      "Delta-method SE may be unreliable for {length(unreliable)} stone{?s}.",
      "i" = "Stones with Var(L)/r^2 > 0.1: {unreliable}.",
      "i" = "Consider increasing {.arg nIter} or {.arg nStones}."
    ))
  }

  if (verbose) {
    cli::cli_alert_success(
      "Log marginal likelihood: {round(logMarginal, 2)} (SE: {round(se, 2)})"
    )
  }

  list(
    log_marginal = logMarginal,
    se = se,
    log_ratios = logRatios,
    betas = betas
  )
}


#' Compute effective sample size for a numeric vector
#'
#' Uses \code{coda::effectiveSize} if available, otherwise falls back to an
#' initial positive sequence estimator based on the sample autocorrelation.
#'
#' @param x Numeric vector of MCMC samples.
#' @return Scalar ESS (at least 1).
#' @keywords internal
.EssVector <- function(x) {
  n <- length(x)
  if (n < 2L) return(1)

  s <- sd(x)
  # is.na() does not catch NaN (NaN is not NA in R); use !is.finite() instead
  if (!is.finite(s) || s == 0) return(as.numeric(n))

  if (requireNamespace("coda", quietly = TRUE)) {
    ess <- as.numeric(coda::effectiveSize(coda::mcmc(x)))
    # effectiveSize can return 0 or NA/NaN for degenerate inputs
    if (!is.finite(ess) || ess < 1) return(1)
    return(ess)
  }

  # Fallback: initial positive sequence estimator (Geyer 1992)
  maxLag <- min(n - 1L, floor(10 * log10(n)))
  acfVals <- acf(x, lag.max = maxLag, plot = FALSE)$acf[, , 1]
  # Sum consecutive pairs of autocorrelations; stop when a pair sum
  # is negative (ensures the estimator of tau is monotone)
  tau <- 1
  k <- 2L # acfVals[1] is lag 0 = 1, so start at lag 1
  while (k + 1L <= length(acfVals)) {
    pairSum <- acfVals[k] + acfVals[k + 1L]
    if (pairSum <= 0) break
    tau <- tau + 2 * pairSum
    k <- k + 2L
  }
  max(n / tau, 1)
}
