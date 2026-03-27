# Model specification for MkPrime MCMC
#
# Defines priors and model options. Used by RunMkPrime() and LogPrior().

#' Specify an MkPrime model
#'
#' @param coding Ascertainment bias correction: `"variable"` (default),
#'   `"informative"`, or `"none"`.
#' @param nCat Number of ACRV rate categories (default 6).
#' @param relabel Apply Mk' relabelling correction for transformational
#'   characters? Default `TRUE`.
#' @param treeLengthShape,treeLengthRate Shape and rate for the Gamma prior
#'   on tree length. Defaults: shape = 2, rate = 2 / `expSteps`.
#' @param expSteps Expected number of character state changes. Used to set
#'   the tree length prior scale. Default `NULL` (computed from parsimony score
#'   of the starting tree).
#' @param rateLossMeanlog,rateLossSdlog Parameters for the LogNormal prior
#'   on `rate_loss` (neomorphic asymmetry). Defaults: meanlog = 0, sdlog = 2.
#' @param rateLogSdShape,rateLogSdRate Shape and rate for the Gamma prior
#'   on `rate_log_sd` (ACRV dispersion). Defaults: shape = 1, rate = 1.
#' @param kPrimePrior Prior distribution for the true number of character states
#'   (`k'`) for transformational characters. One of `"geometric"` (default,
#'   hierarchical geometric with Beta hyperprior on `p`) or `"logseries"`
#'   (logarithmic series with fixed parameter `c`; matches the RevBayes default).
#' @param kprimeHyperA,kprimeHyperB Parameters for the Beta hyperprior on `p`
#'   when `kPrimePrior = "geometric"`. Defaults: a = 1, b = 1 (uniform).
#' @param kprimeLogseriesC The `c` parameter of the log-series prior on `k'`
#'   when `kPrimePrior = "logseries"`. Must be in (0, 1). Default 0.7, matching
#'   the RevBayes `dnMkPrime` default.
#' @param rateNeoMeanlog,rateNeoSdlog Parameters for the LogNormal prior
#'   on the neomorphic partition rate scalar. Defaults: meanlog = 0, sdlog = 2.
#'
#' @return An S3 object of class `MkPrimeModel`.
#' @export
MkPrimeModel <- function(
    coding = "variable",
    nCat = 6L,
    relabel = TRUE,
    treeLengthShape = 2,
    treeLengthRate = NULL,
    expSteps = NULL,
    rateLossMeanlog = 0,
    rateLossSdlog = 2,
    rateLogSdShape = 1,
    rateLogSdRate = 1,
    kPrimePrior = "geometric",
    kprimeHyperA = 1,
    kprimeHyperB = 1,
    kprimeLogseriesC = 0.7,
    rateNeoMeanlog = 0,
    rateNeoSdlog = 2
) {
  coding <- match.arg(coding, c("variable", "informative", "none"))
  kPrimePrior <- match.arg(kPrimePrior, c("geometric", "logseries"))

  # Warn if logseries-specific param is supplied for geometric prior
  if (kPrimePrior == "geometric" && !missing(kprimeLogseriesC)) {
    cli::cli_warn(
      "{.arg kprimeLogseriesC} is ignored when {.arg kPrimePrior = \"geometric\"}."
    )
  }

  # Derive treeLengthRate from expSteps if not provided
  if (is.null(treeLengthRate) && !is.null(expSteps)) {
    treeLengthRate <- 2 / expSteps
  }

  structure(
    list(
      coding = coding,
      nCat = as.integer(nCat),
      relabel = relabel,
      treeLengthShape = treeLengthShape,
      treeLengthRate = treeLengthRate,
      expSteps = expSteps,
      rateLossMeanlog = rateLossMeanlog,
      rateLossSdlog = rateLossSdlog,
      rateLogSdShape = rateLogSdShape,
      rateLogSdRate = rateLogSdRate,
      kPrimePrior = kPrimePrior,
      kprimeHyperA = kprimeHyperA,
      kprimeHyperB = kprimeHyperB,
      kprimeLogseriesC = kprimeLogseriesC,
      rateNeoMeanlog = rateNeoMeanlog,
      rateNeoSdlog = rateNeoSdlog
    ),
    class = "MkPrimeModel"
  )
}


#' Finalize model with data-derived defaults
#'
#' Sets `expSteps` and `treeLengthRate` if not user-specified.
#' Called internally by [RunMkPrime()] before MCMC starts.
#'
#' @param model An `MkPrimeModel` object.
#' @param tree A `phylo` object (starting tree).
#' @param mkd An `MkPrimeData` object.
#' @return Updated `MkPrimeModel` with all defaults resolved.
#' @keywords internal
.FinalizeModel <- function(model, tree, mkd) {
  if (is.null(model$expSteps)) {
    model$expSteps <- max(1, .FitchScore(tree, mkd))
  }
  if (is.null(model$treeLengthRate)) {
    model$treeLengthRate <- 2 / model$expSteps
  }
  model
}


#' Fitch parsimony score on a tree
#'
#' Simple post-order Fitch algorithm. Used to set a data-informed default
#' for `expSteps` (the expected tree length).
#'
#' @param tree A `phylo` object.
#' @param mkd An `MkPrimeData` object.
#' @return Integer parsimony score.
#' @keywords internal
.FitchScore <- function(tree, mkd) {
  tree <- TreeTools::Postorder(tree)
  edge <- tree$edge
  nTip <- length(tree$tip.label)
  nNode <- tree$Nnode

  # Align matrix rows with tree tip order
  tipMat <- mkd$matrix[tree$tip.label, , drop = FALSE]

  total <- 0L
  for (j in seq_len(ncol(tipMat))) {
    # Initialize state sets: list of integer vectors per node
    sets <- vector("list", nTip + nNode)
    for (i in seq_len(nTip)) {
      s <- tipMat[i, j]
      sets[[i]] <- if (is.na(s)) seq.int(0L, mkd$kObs[j] - 1L) else s
    }

    for (i in seq_len(nrow(edge))) {
      p <- edge[i, 1]
      ch <- edge[i, 2]
      if (is.null(sets[[p]])) {
        sets[[p]] <- sets[[ch]]
      } else {
        inter <- intersect(sets[[p]], sets[[ch]])
        if (length(inter) > 0L) {
          sets[[p]] <- inter
        } else {
          sets[[p]] <- union(sets[[p]], sets[[ch]])
          total <- total + 1L
        }
      }
    }
  }

  total
}


#' Compute total log-prior density
#'
#' @param state A list with current parameter values:
#'   `tree_length`, `rel_br_lengths`, `rate_loss`, `rate_log_sd`,
#'   `kPrime` (integer vector). For `kPrimePrior = "geometric"`, also
#'   `p` (hyperprior). For `kPrimePrior = "logseries"`, `p` is absent.
#' @param model An `MkPrimeModel` object (finalized).
#' @param mkd An `MkPrimeData` object (for kObs and character types).
#' @return Scalar log-prior density.
#' @keywords internal
LogPrior <- function(state, model, mkd) {
  # Boundary checks — return -Inf for out-of-support values

  if (state$tree_length <= 0) return(-Inf)
  if (state$rate_log_sd < 0) return(-Inf)
  if (any(state$rel_br_lengths <= 0)) return(-Inf)

  hasNeo <- any(mkd$type == "neomorphic")
  if (hasNeo && state$rate_loss <= 0) return(-Inf)
  if (hasNeo && !is.null(state$rate_neo) && state$rate_neo <= 0) return(-Inf)

  transIdx <- which(mkd$type == "transformational")
  hasTrans <- length(transIdx) > 0L

  if (hasTrans) {
    if (any(state$kPrime[transIdx] < mkd$kObs[transIdx])) return(-Inf)

    if (identical(model$kPrimePrior, "geometric")) {
      if (state$p <= 0 || state$p >= 1) return(-Inf)
    } else {
      # logseries: validate c
      c_ls <- model$kprimeLogseriesC
      if (c_ls <= 0 || c_ls >= 1) return(-Inf)
    }
  }

  lp <- 0.0

  # Tree length: Gamma prior
  lp <- lp + dgamma(state$tree_length,
                     shape = model$treeLengthShape,
                     rate = model$treeLengthRate,
                     log = TRUE)

  # Relative branch lengths: Dirichlet(1, ..., 1) = uniform on simplex
  # log-density is constant: log((n-1)!). Doesn't affect MH ratios but
  # included for correct log-posterior reporting.
  nEdges <- length(state$rel_br_lengths)
  lp <- lp + lfactorial(nEdges - 1L)

  # rate_loss: LogNormal prior (neomorphic characters only)
  if (hasNeo) {
    lp <- lp + dlnorm(state$rate_loss,
                      meanlog = model$rateLossMeanlog,
                      sdlog = model$rateLossSdlog,
                      log = TRUE)

    # rate_neo: partition rate scalar (LogNormal prior)
    if (!is.null(state$rate_neo)) {
      lp <- lp + dlnorm(state$rate_neo,
                         meanlog = model$rateNeoMeanlog,
                         sdlog = model$rateNeoSdlog,
                         log = TRUE)
    }
  }

  # rate_log_sd: Gamma prior (rate_log_sd = 0 is a boundary; dgamma(0) = 0
  # for shape >= 1, but log(0) = -Inf. Treat 0 specially as a valid point.)
  if (state$rate_log_sd > 0) {
    lp <- lp + dgamma(state$rate_log_sd,
                       shape = model$rateLogSdShape,
                       rate = model$rateLogSdRate,
                       log = TRUE)
  }
  # When rate_log_sd == 0 and shape == 1, the density is finite (rate);
  # when shape > 1, density is 0. Handle both:
  if (state$rate_log_sd == 0 && model$rateLogSdShape > 1) return(-Inf)

  if (hasTrans) {
    if (identical(model$kPrimePrior, "geometric")) {
      # k'_i: Geometric(p) shifted by kObs_i
      # P(k'_i = kObs_i + u) = p * (1-p)^u, u = 0, 1, 2, ...
      u <- state$kPrime[transIdx] - mkd$kObs[transIdx]
      lp <- lp + length(transIdx) * log(state$p) + sum(u) * log1p(-state$p)

      # p: Beta hyperprior
      lp <- lp + dbeta(state$p,
                       shape1 = model$kprimeHyperA,
                       shape2 = model$kprimeHyperB,
                       log = TRUE)
    } else {
      # k'_i: Logseries(c)
      # log P(k; c) = k*log(c) - log(k) - log(-log(1-c))
      # Truncation at kObs cancels in MH ratios; constant included here
      # for correct absolute log-posterior reporting.
      c_ls <- model$kprimeLogseriesC
      kp <- state$kPrime[transIdx]
      lp <- lp + sum(kp * log(c_ls) - log(kp)) -
            length(transIdx) * log(-log1p(-c_ls))
    }
  }

  lp
}


#' @export
print.MkPrimeModel <- function(x, ...) {
  cli::cli_h1("MkPrime Model")

  k_prior_str <- if (identical(x$kPrimePrior, "logseries")) {
    "Logseries (c = {x$kprimeLogseriesC})"
  } else {
    "Geometric (Beta hyperprior: a = {x$kprimeHyperA}, b = {x$kprimeHyperB})"
  }

  cli::cli_ul(c(
    "Coding: {x$coding}",
    "ACRV categories: {x$nCat}",
    "Relabelling correction: {x$relabel}",
    "Tree length prior: Gamma({x$treeLengthShape}, {x$treeLengthRate %||% 'auto'})",
    "rate_loss prior: LogNormal({x$rateLossMeanlog}, {x$rateLossSdlog})",
    "rate_log_sd prior: Gamma({x$rateLogSdShape}, {x$rateLogSdRate})",
    paste0("k' prior: ", k_prior_str),
    "rate_neo prior: LogNormal({x$rateNeoMeanlog}, {x$rateNeoSdlog})"
  ))
  invisible(x)
}
