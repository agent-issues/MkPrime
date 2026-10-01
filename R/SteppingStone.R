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
#' @param nIter Number of MCMC iterations per stone, each a single move (not
#'   a sweep through every move type). Default 1000.
#' @param warmup Warmup iterations per stone (discarded), counted as for
#'   `nIter`. Default 200.
#' @param alpha Shape parameter for the Beta(alpha, 1) quantile schedule.
#'   Smaller values concentrate more stones near beta = 0 where the integrand
#'   changes fastest. Default 0.3 (Xie et al. recommendation).
#' @param fixTopology If `TRUE`, fix the tree topology and only sample
#'   continuous parameters + branch lengths. Default `FALSE`.
#' @param nRuns Integer specifying the number of independent runs through the
#'   whole ladder; at least two are needed for a standard error. Default 2.
#' @param verbose Logical specifying whether to report per-stone progress;
#'   see [MkPrimeVerbosity()].
#'
#' @return A list with components:
#'   \describe{
#'     \item{log_marginal}{Log marginal likelihood estimate, pooling the
#'       runs.}
#'     \item{se}{Standard error of the estimate (see Details); `NA` when
#'       `nRuns = 1`.}
#'     \item{run_log_marginal}{The estimate from each run alone (length `nRuns`).}
#'     \item{log_ratios}{Per-stone log ratios (length `nStones`).}
#'     \item{betas}{The beta schedule used.}
#'   }
#'   If any stone of any run has no finite log-likelihood, or an infinite one,
#'   `log_marginal` and `se` are `NA`, with a warning.
#'
#' @details
#' Every stone runs for `warmup + nIter` iterations under the schedule `model`
#' and `mcmc` specify, less the 2D joint proposals (which need a correlation
#' this function does not adapt) and `gibbs_p_marginal` (which is not
#' tempered), and with no weight adaptation. An estimate is therefore
#' comparable only with others computed under the same `model` and `mcmc`.
#' One iteration performs one move, so a stone updates a parameter about
#' `nIter` times the weight of the moves that change it; a large tree, whose
#' many branch lengths share those weights, needs a correspondingly larger
#' `nIter`.
#'
#' Under `likelihoodMode = "sampled_k"` and the unconditional `k'` prior, no
#' stone visits `k'_i < kObs_i`, so the prior-only stone is normalised on
#' `k' >= kObs` alone; the estimate adds back the log prior mass of that
#' region, integrated over `p` numerically.
#'
#' Each run is one chain that carries its state from stone to stone, starting
#' from `tree`. The estimate pools the draws of every run at each stone; the
#' standard error is \eqn{s / \sqrt{n}}, where \eqn{s} is the standard
#' deviation of the `n = nRuns` runs' own estimates. Within a run, the
#' autocorrelation of the chain commonly outlasts a stone and carries into
#' the next, so an error computed from each stone's draws alone (the delta
#' method of Xie et al. 2011) can understate the true error several-fold;
#' independent runs capture both. With few runs the standard error is itself
#' uncertain: with the default two it rests on one degree of freedom.
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
                              nRuns = 2L,
                              verbose = MkPrimeVerbosity() > 0L) {

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
  tree <- .PrepareStartTree(tree, mkd)

  classIdx <- vapply(mkd$partitions, function(p) as.integer(p$classIdx %||% 1L),
                     integer(1))
  if (any(classIdx > 1L)) {
    cli::cli_abort(c(
      "{.fn mkp_stepping_stone} does not support partitioned models.",
      "i" = "Supply {.arg data} without a user {.arg partition}."
    ))
  }

  if (is.null(model)) model <- MkPrimeModel()
  if (identical(model$coding, "informative")) {
    mkd <- .DropUninformable(mkd)
  }
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
  model <- .FinalizeModel(model, NULL, mkd)

  nStones <- as.integer(nStones)
  nIter <- as.integer(nIter)
  warmup <- as.integer(warmup)
  nRuns <- as.integer(nRuns)
  if (length(nRuns) != 1L || is.na(nRuns) || nRuns < 1L) {
    cli::cli_abort("{.arg nRuns} must be a positive integer.")
  }

  # --- Beta schedule: quantiles of Beta(alpha, 1) ---
  # beta_k = (k/K)^(1/alpha), k = 0, ..., K
  # This concentrates stones near beta = 0
  betas <- ((seq_len(nStones + 1L) - 1L) / nStones)^(1 / alpha)
  # betas[1] = 0 (prior only), betas[nStones + 1] = 1 (full posterior)

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
                       likelihoodMode = model$likelihoodMode %||% "sampled_k",
                       rateNeoLive = .RateNeoLive(mkd))
  moveWeights <- vapply(moves, `[[`, numeric(1), "weight")
  names(moveWeights) <- vapply(moves, `[[`, character(1), "name")
  moveWeights <- moveWeights / sum(moveWeights)
  pinnedWeights <- .SchedulePins(
    moveWeights, vapply(moves, `[[`, character(1), "type"),
    .ResolvePinnedWeights(mcmc$moveWeights, names(moveWeights))
  )
  if (!is.null(pinnedWeights)) {
    moveWeights <- .NormalizeMoveWeights(moveWeights, pinnedWeights)
  }

  mcmcTuning <- mcmc$tuning
  mcmcData <- .InitMcmcData(mkd, model)
  set_branch_bins(mcmcData, mcmc$nBranchBins)
  transIdx <- which(mkd$type == "transformational")

  # Before the stones, so a failed quadrature costs no sampling.
  logZ0 <- .LogZ0(model, mkd)

  # One chain per run walks the whole ladder, carrying its state from stone
  # to stone; that carry-over is why runs, not stones, are the replicates.
  logLiks <- lapply(seq_len(nRuns), function(run) {
    statePtr <- .InitMcmcChain(state)
    fill_partition_cache(mcmcData, statePtr)
    allocate_cl_workspace(mcmcData, statePtr)  # M-063
    runLogLiks <- matrix(NA_real_, nIter, nStones)
    for (stone in seq_len(nStones)) {
      if (verbose) {
        cli::cli_progress_step(
          "Run {run}/{nRuns}, stone {stone}/{nStones}: beta = \\
           [{round(betas[stone], 4)}, {round(betas[stone + 1L], 4)}]"
        )
      }
      for (iter in seq_len(warmup + nIter)) {
        move <- moves[[sample.int(length(moves), 1L, prob = moveWeights)]]
        .DoMove(move, statePtr, tuning = mcmcTuning,
                beta = betas[stone], transIdx = transIdx, mcmcData = mcmcData)
        if (iter > warmup) {
          runLogLiks[iter - warmup, stone] <- get_state_log_lik(statePtr)
        }
      }
    }
    runLogLiks
  })

  result <- .SteppingStoneEstimate(logLiks, betas, logZ0)
  if (verbose) {
    cli::cli_alert_success(
      "Log marginal likelihood: {round(result$log_marginal, 2)} \\
       (SE: {round(result$se, 2)})"
    )
  }

  # Return:
  c(result, list(betas = betas))
}


# `logLiks` holds one draws x stones matrix per run. The estimate pools every
# run's draws at each stone; its standard error is that of the mean of the
# runs' own estimates, which are independent where a run's stones are not.
.SteppingStoneEstimate <- function(logLiks, betas, logZ0) {
  nStones <- length(betas) - 1L
  nRuns <- length(logLiks)
  runRatios <- matrix(vapply(logLiks, .StoneLogRatios, numeric(nStones),
                             betas), nStones)
  runLogMarginal <- colSums(runRatios) + logZ0
  logRatios <- .StoneLogRatios(do.call(rbind, logLiks), betas)
  degenerate <- which(rowSums(is.na(runRatios)) > 0L)
  if (length(degenerate) > 0L) {
    cli::cli_warn(c(
      "Log marginal likelihood is undefined: {length(degenerate)} stone{?s}
       had no finite log-likelihood, or an infinite one.",
      "i" = "Stone{?s} {degenerate}; check the model and starting values."
    ))
    logMarginal <- NA_real_
    se <- NA_real_
  } else {
    logMarginal <- sum(logRatios) + logZ0
    se <- if (nRuns > 1L) stats::sd(runLogMarginal) / sqrt(nRuns) else NA_real_
  }
  # Return:
  list(log_marginal = logMarginal, se = se, run_log_marginal = runLogMarginal,
       log_ratios = logRatios)
}


# Stepping-stone log ratios, one per stone, from a draws x stones matrix of
# log-likelihoods sampled at betas[-length(betas)]. NA marks a stone with no
# finite draw or one at +Inf / NaN, whose ratio is undefined.
.StoneLogRatios <- function(logLiks, betas) {
  deltaBeta <- diff(betas)
  vapply(seq_along(deltaBeta), function(stone) {
    ll <- logLiks[, stone]
    if (anyNA(ll) || any(ll == Inf) || !any(is.finite(ll))) {
      # Return:
      return(NA_real_)
    }
    shifted <- deltaBeta[stone] * ll
    maxShifted <- max(shifted)
    # Return:
    maxShifted + log(mean(exp(shifted - maxShifted)))
  }, numeric(1))
}


# Under the unconditional (Model A) prior, sampled_k gives k'_i < kObs_i zero
# likelihood and no move proposes it, so the beta = 0 stone samples the prior
# restricted to k' >= kObs and the stones multiply to Z1 / Z0 (#267), with
#   Z0 = E_prior[prod_i P(k'_i >= kObs_i | hyperparameters)].
# Model B already normalises on k' >= kObs, and marginal_k carries that mass in
# its per-character sum, so Z0 = 1 for both.
.LogZ0 <- function(model, mkd) {
  # Resolved as .InitMcmcData resolves it for the sampler.
  priorVariant <- model$priorVariant %||% MkPrimeModel(
    kPrimePrior = model$kPrimePrior %||% "geometric")$priorVariant
  if (!identical(priorVariant, "unconditional") ||
      identical(model$likelihoodMode, "marginal_k")) {
    # Return:
    return(0)
  }
  kObs <- mkd$kObs[mkd$type == "transformational"]
  kObs <- kObs[kObs > 2L]
  if (length(kObs) == 0L) {
    # Return:
    return(0)
  }
  counts <- table(kObs)
  ko <- as.integer(names(counts))
  counts <- as.vector(counts)

  if (identical(model$kPrimePrior, "logseries")) {
    logseriesC <- model$kprimeLogseriesC %||% 0.7
    # Return:
    return(sum(counts * vapply(ko, .LogseriesLogTail, 0, logseriesC)))
  }
  LogTail <- switch(
    model$kPrimePrior,
    geometric = {
      K <- as.integer(model$kprimeTruncK %||% 200L)
      function(k, p, log1mP) {
        (k - 2) * log1mP + .Log1mExp((K - k + 1) * log1mP) -
          .Log1mExp((K - 1) * log1mP)
      }
    },
    empirical_geometric = {
      emp <- model$empiricalNObs %||% .EmpiricalNObs()
      # The two variants differ by exactly the Model B normaliser log Z_i(p).
      function(k, p, log1mP) {
        vapply(p, function(pj) {
          .LogPriorEmpiricalGeometric(k, emp, pj, k, unconditional = TRUE) -
            .LogPriorEmpiricalGeometric(k, emp, pj, k, unconditional = FALSE)
        }, 0)
      }
    },
    cli::cli_abort("No unconditional form for {.val {model$kPrimePrior}}.")
  )

  a <- model$kprimeHyperA
  b <- model$kprimeHyperB
  # Integrate over x = logit(p); p^a (1-p)^b folds in the Jacobian p (1-p).
  LogF <- function(x) {
    p <- stats::plogis(x)
    log1mP <- stats::plogis(-x, log.p = TRUE)
    lf <- a * stats::plogis(x, log.p = TRUE) + b * log1mP - lbeta(a, b)
    for (j in seq_along(ko)) {
      lf <- lf + counts[j] * LogTail(ko[j], p, log1mP)
    }
    # Out at |x| ~ 700, p rounds to 0 or 1 and the tail is 0/0.
    lf[is.nan(lf)] <- -Inf
    # Return:
    lf
  }
  grid <- seq(-40, 40, by = 0.5)
  xMax <- grid[which.max(LogF(grid))]
  xMax <- stats::optimize(LogF, xMax + c(-0.5, 0.5), maximum = TRUE)$maximum
  fMax <- LogF(xMax)
  Integrand <- function(x) exp(LogF(x) - fMax)
  mass <- stats::integrate(Integrand, -Inf, xMax, rel.tol = 1e-10)$value +
    stats::integrate(Integrand, xMax, Inf, rel.tol = 1e-10)$value
  # Return:
  fMax + log(mass)
}


# log P(k' >= ko; logseriesC) under the logseries prior on k' >= 2. The complement is
# summed while it is small; otherwise the tail, whose terms fall by at least
# logseriesC.
.LogseriesLogTail <- function(ko, logseriesC) {
  logC <- log(logseriesC)
  logNorm <- .LogseriesLogNorm(logseriesC)
  below <- seq.int(2L, ko - 1L)
  lower <- sum(exp(below * logC - log(below) - logNorm))
  if (lower < 0.99) {
    # Return:
    return(log1p(-lower))
  }
  kk <- seq.int(ko, ko + ceiling(-60 / logC))
  logTerms <- kk * logC - log(kk)
  mx <- max(logTerms)
  # Return:
  mx + log(sum(exp(logTerms - mx))) - logNorm
}
