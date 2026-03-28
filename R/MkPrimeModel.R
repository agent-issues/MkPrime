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
#' @param qHeterogeneity Logical. Enable Q-matrix heterogeneity across
#'   characters via discretised Dirichlet-marginal equilibrium frequencies?
#'   Default `FALSE`. When enabled, each character's likelihood is averaged
#'   over a set of F81 rate matrices whose equilibrium frequencies are drawn
#'   from `Beta(beta_scale, (k - 1) * beta_scale)`, where `k` is the number
#'   of states. This is the marginal distribution of one component of a
#'   symmetric `Dirichlet(beta_scale, ..., beta_scale)`. See the
#'   \strong{Q-matrix heterogeneity} section below for details.
#' @param nBetaCat Integer. Number of equal-probability Beta bins for the
#'   heterogeneity discretisation. Default `4L`. Higher values increase
#'   accuracy at the cost of compute time (scales linearly with `nBetaCat`).
#'   Ignored when `qHeterogeneity = FALSE`.
#' @param betaScaleShape,betaScaleRate Shape and rate for the Gamma prior
#'   on `beta_scale` (the symmetric Dirichlet concentration parameter).
#'   Defaults: shape = 1, rate = 1. Ignored when `qHeterogeneity = FALSE`.
#'
#' @section Q-matrix heterogeneity:
#'
#' Standard Mk and Mk' assume all characters share the same (equal-frequency)
#' rate matrix. In reality, some morphological characters may have strongly
#' unequal state frequencies. Q-matrix heterogeneity (`qHeterogeneity = TRUE`)
#' relaxes this by integrating each character's likelihood over a mixture of
#' F81 rate matrices \insertCite{Felsenstein1981}{MkPrime} with varying
#' equilibrium frequencies.
#'
#' The mixture is controlled by a single scalar parameter, `beta_scale`
#' (= \eqn{\alpha}), which acts as the concentration of a symmetric Dirichlet:
#' \itemize{
#'   \item Large \eqn{\alpha}: all characters have nearly equal state
#'     frequencies, recovering the standard Mk/JC model.
#'   \item Small \eqn{\alpha}: characters can have strongly unequal
#'     frequencies (one state dominant, others rare).
#' }
#'
#' The key insight is that the marginal distribution of one component of
#' \eqn{\mathrm{Dirichlet}(\alpha, \ldots, \alpha)} with \eqn{k} components
#' is \eqn{\mathrm{Beta}(\alpha, (k-1)\alpha)}. This allows a unified
#' discretisation for characters of any state count: the scheme adapts
#' automatically to \eqn{k}.
#'
#' For binary characters (\eqn{k = 2}), this simplifies to the symmetric
#' \eqn{\mathrm{Beta}(\alpha, \alpha)}.
#'
#' This feature is intended for **model comparison** (e.g. via
#' [mkp_stepping_stone()]). It multiplies computation time by roughly
#' \eqn{4-5\times}{4-5x} for binary-dominated datasets.
#' See `vignette("het-details", package = "MkPrime")` for the full
#' mathematical derivation.
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
    rateNeoSdlog = 2,
    qHeterogeneity = FALSE,
    nBetaCat = 4L,
    betaScaleShape = 1,
    betaScaleRate = 1
) {
  coding <- match.arg(coding, c("variable", "informative", "none"))
  kPrimePrior <- match.arg(kPrimePrior, c("geometric", "logseries"))

  # Warn if logseries-specific param is supplied for geometric prior
  if (kPrimePrior == "geometric" && !missing(kprimeLogseriesC)) {
    cli::cli_warn(
      "{.arg kprimeLogseriesC} is ignored when {.arg kPrimePrior = \"geometric\"}."
    )
  }

  # M-052: validate Het parameters
  if (isTRUE(qHeterogeneity)) {
    nBetaCat <- as.integer(nBetaCat)
    if (nBetaCat < 1L || nBetaCat > 16L) {
      cli::cli_abort("{.arg nBetaCat} must be between 1 and 16 (got {nBetaCat}).")
    }
    if (betaScaleShape <= 0 || betaScaleRate <= 0) {
      cli::cli_abort(
        "{.arg betaScaleShape} and {.arg betaScaleRate} must be positive."
      )
    }
    # M-100: het_singleton_site_prob() is a stub returning 0, so informative
    # coding produces wrong ascertainment corrections under Het.
    # Remove this guard when F81 singleton correction is implemented.
    if (coding == "informative") {
      cli::cli_abort(c(
        "{.arg qHeterogeneity} cannot be combined with
         {.code coding = \"informative\"} yet.",
        i = "Singleton ascertainment correction under Q-heterogeneity
             is not yet implemented.",
        i = "Use {.code coding = \"variable\"} or disable
             {.code qHeterogeneity}."
      ))
    }
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
      rateNeoSdlog = rateNeoSdlog,
      qHeterogeneity = qHeterogeneity,
      nBetaCat = as.integer(nBetaCat),
      betaScaleShape = betaScaleShape,
      betaScaleRate = betaScaleRate
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
  tree <- TreeTools::Preorder(tree)
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

    # Reverse traversal: preorder edges reversed = bottom-up (Fitch pass)
    for (i in rev(seq_len(nrow(edge)))) {
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
#'   When `qHeterogeneity = TRUE`, also `beta_scale` (positive scalar).
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

  # M-052: beta_scale prior (Q-matrix heterogeneity)
  if (isTRUE(model$qHeterogeneity)) {
    bs <- state$beta_scale
    if (is.null(bs) || bs <= 0) return(-Inf)
    lp <- lp + dgamma(bs,
                       shape = model$betaScaleShape,
                       rate = model$betaScaleRate,
                       log = TRUE)
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

  het_str <- if (isTRUE(x$qHeterogeneity)) {
    "ON ({x$nBetaCat} bins, beta_scale ~ Gamma({x$betaScaleShape}, {x$betaScaleRate}))"
  } else {
    "OFF"
  }

  cli::cli_ul(c(
    "Coding: {x$coding}",
    "ACRV categories: {x$nCat}",
    "Relabelling correction: {x$relabel}",
    "Tree length prior: Gamma({x$treeLengthShape}, {x$treeLengthRate %||% 'auto'})",
    "rate_loss prior: LogNormal({x$rateLossMeanlog}, {x$rateLossSdlog})",
    "rate_log_sd prior: Gamma({x$rateLogSdShape}, {x$rateLogSdRate})",
    paste0("k' prior: ", k_prior_str),
    "rate_neo prior: LogNormal({x$rateNeoMeanlog}, {x$rateNeoSdlog})",
    paste0("Q-matrix heterogeneity: ", het_str)
  ))
  invisible(x)
}
