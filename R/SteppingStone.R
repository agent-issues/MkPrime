#' Marginal likelihood via stepping-stone sampling
#'
#' Estimates the marginal likelihood of the data under the Mk' model using
#' stepping-stone sampling (Xie et al., 2011). This runs MCMC at a series of
#' "power posteriors" interpolating between the prior (beta = 0) and the full
#' posterior (beta = 1), then combines the results into a marginal likelihood
#' estimate.
#'
#' @param data A `phyDat` object or `MkPrimeData` object.
#' @param tree A starting tree (`phylo` object), or `NULL` (default) to use a
#'   neighbour-joining tree built from the data.
#' @param neomorphic Integer vector of neomorphic character indices (if `data`
#'   is `phyDat`). Ignored if `data` is `MkPrimeData`.
#' @param model An `MkPrimeModel` object. Default: `MkPrimeModel()`.
#' @param mcmc `MkPrimeMCMC` object supplying the move schedule and proposal
#'   tuning for each stone; its run-length and stopping-rule settings are
#'   ignored.
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
#' Every stone runs for `warmup + nIter` iterations under the schedule `model`
#' and `mcmc` specify, less the 2D joint proposals (which need a correlation
#' this function does not adapt) and `gibbs_p_marginal` (which is not
#' tempered), and with no weight adaptation. An estimate is therefore
#' comparable only with others computed under the same `model` and `mcmc`.
#'
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
mkp_stepping_stone <- function(data, tree = NULL,
                              neomorphic = integer(0),
                              model = NULL,
                              mcmc = NULL,
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

  if (is.null(tree)) {
    njInput <- if (inherits(data, "phyDat")) data else mkd$phyDat
    tree <- TreeTools::NJTree(njInput, edgeLengths = TRUE)
    if (verbose) {
      cli::cli_alert_info("No starting tree supplied; using neighbour-joining tree.")
    }
  }
  if (is.null(tree$edge.length)) {
    cli::cli_abort(c(
      "{.arg tree} has no branch lengths.",
      "i" = "Supply a tree with edge lengths, e.g. {.code TreeTools::NJTree(data)}."
    ))
  }
  nNeg <- sum(tree$edge.length <= 0)
  if (nNeg > 0L) {
    cli::cli_warn(c(
      "{nNeg} non-positive branch length{?s} clamped to 1e-8.",
      "i" = "Zero or negative lengths arise in NJ trees when taxa are very similar."
    ))
    tree$edge.length[tree$edge.length <= 0] <- 1e-8
  }

  if (is.null(model)) model <- MkPrimeModel()
  if (is.null(mcmc)) mcmc <- MkPrimeMCMC()
  if (isTRUE(mcmc$gibbsPMarginal)) {
    # gibbs_p_marginal accepts on a truncation-normaliser ratio that carries
    # no power, so in a stone it would draw p from the beta = 1 conditional.
    cli::cli_warn(c(
      "{.arg gibbsPMarginal} is ignored by {.fn mkp_stepping_stone}.",
      "i" = "{.code gibbs_p_marginal} is not tempered; {.code mh_logit_p}
             samples {.code p} in each stone instead."
    ))
    mcmc$gibbsPMarginal <- FALSE
  }
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
  # PREORDER INVARIANT: .MkpLogLikelihood() (internal fast-path) requires
  # edges in canonical preorder. Reorder here so .InitState() and all
  # subsequent topology proposals maintain the invariant.
  tree <- TreeTools::Preorder(tree)
  state <- .InitState(tree, mkd, model)

  nEdge <- nrow(tree$edge)
  nTrans <- sum(mkd$type == "transformational")
  hasNeo <- any(mkd$type == "neomorphic")
  # joint2d = FALSE: stepping-stone doesn't adapt rho, so joint moves
  # would always run with rho = 0 (redundant with individual scale moves).
  moves <- .BuildMoves(nEdge, nTrans, hasNeo, mcmc,
                       fixTopology = fixTopology,
                       kPrimePrior = model$kPrimePrior %||% "geometric",
                       qHeterogeneity = isTRUE(model$qHeterogeneity),
                       joint2d = FALSE,
                       likelihoodMode = model$likelihoodMode %||% "sampled_k")
  moveWeights <- vapply(moves, `[[`, numeric(1), "weight")
  names(moveWeights) <- vapply(moves, `[[`, character(1), "name")
  moveWeights <- moveWeights / sum(moveWeights)
  userPins <- .ResolvePinnedWeights(mcmc$moveWeights, names(moveWeights))
  if (!is.null(userPins)) {
    moveWeights <- .NormalizeMoveWeights(moveWeights, userPins)
  }

  mcmcTuning <- mcmc$tuning

  # Build C++ data struct and XPtr state (shared across all stones)
  mcmcData <- .InitMcmcData(mkd, model)
  set_branch_bins(mcmcData, mcmc$nBranchBins)
  transIdx  <- which(mkd$type == "transformational")
  statePtr  <- .InitMcmcChain(state)
  fill_partition_cache(mcmcData, statePtr)
  allocate_cl_workspace(mcmcData, statePtr)  # M-063

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


# .EssVector() is now defined in R/ess.R (native Geyer 1992 + Vehtari 2021).
