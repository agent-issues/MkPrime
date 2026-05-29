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
#'   Used as-is in the legacy single-σ path and as the legacy
#'   `priorOnClassRateLogSd = "gamma_independent"` per-class prior.
#'   Not consulted by the default pooled hyperprior (see
#'   `priorOnClassRateLogSd`).
#' @param priorOnClassRateLogSd Prior structure on per-class ACRV
#'   dispersion `σ_c = class_rate_log_sd[c]` when `unlink = "shape"` is
#'   active with two or more user classes. One of:
#'   * `"hyperprior_pooled"` (default): non-centred half-normal hierarchy
#'     `σ_c = τ · z_c`, with `z_c ~ HalfNormal(1)` i.i.d. and
#'     `τ ~ HalfNormal(1)`. Pools information about within-class rate
#'     dispersion across classes; designed for cells where some classes
#'     have few characters and so carry little likelihood signal on σ_c.
#'   * `"gamma_independent"`: each `σ_c` independently
#'     `Gamma(rateLogSdShape, rateLogSdRate)`. Legacy behaviour, kept for
#'     backward compatibility and prior-sensitivity comparisons.
#'
#'   Ignored when `unlink` does not include `"shape"` or when only one
#'   user class is present (the structure is degenerate at K = 1; the
#'   single σ uses the `rateLogSdShape` / `rateLogSdRate` Gamma prior).
#' @param kPrimePrior Prior distribution for the true number of character states
#'   (`k'`) for transformational characters. One of `"empirical_geometric"`
#'   (default; convolution of an empirical pmf on the number of observed
#'   states with a `Geometric(p)` prior on the number of unobserved states,
#'   so `k' = N_obs + N_unobs`),
#'   `"beta_geometric"` (per-character Beta-Geometric with shared
#'   hyperparameters `alpha`, `beta`),
#'   `"geometric"` (hierarchical geometric with Beta hyperprior on `p`),
#'   or `"logseries"` (logarithmic series with fixed parameter `c`; matches
#'   the RevBayes default).
#'   The `"empirical_geometric"` option counters the tendency of inference
#'   to collapse to zero unobserved states by anchoring the observed-state
#'   component to a real-data empirical distribution; the `"beta_geometric"`
#'   option avoids the same over-shrinkage in a different way, by allowing
#'   per-character `p`.
#' @param empiricalNObs An object of class `"MkPrimeEmpiricalPrior"` (see
#'   [MkPrimeEmpiricalPrior()]) giving the empirical pmf on `N_obs` used by
#'   the `"empirical_geometric"` prior. Defaults to the package dataset
#'   [empiricalNObs] tabulated from the `neotrans` corpus. Ignored for other
#'   prior options.
#' @param kprimeHyperA,kprimeHyperB Parameters for the Beta hyperprior on `p`
#'   when `kPrimePrior = "geometric"`. Defaults: a = 1, b = 1 (uniform).
#' @param kprimeAlpha,kprimeBeta Starting values for the shared
#'   hyperparameters of the Beta-Geometric prior
#'   (`kPrimePrior = "beta_geometric"`). Both must be positive.
#'   Defaults: alpha = 1, beta = 1. These are estimated during MCMC with
#'   Exponential(1) hyperpriors.
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
#' @param likelihoodMode How per-character `k'_i` enters the likelihood.
#'   One of:
#'   * `"sampled_k"` (default): `k'_i` is a sampled MCMC state variable;
#'     a per-character integer-walk + Gibbs sweep (case 25) move family
#'     updates it. The trace carries `kPrime_i...` columns and `LogPrior`
#'     evaluates the per-character `P(k'_i | hyperparams)` term.
#'   * `"marginal_k"` (v1 = geometric arm only): `k'_i` is analytically
#'     marginalised out of the likelihood at every evaluation. The
#'     posterior is over `(tree, mu, sigma, p)` only; the slow discrete
#'     coordinate is removed. Required when chain mixing on `k'_i` is
#'     known to be uninformative (per finding EG-003 in
#'     `dev/red-team/findings.md`). v1 supports
#'     `kPrimePrior = "geometric"` only — other arms are §11 follow-ups
#'     in `dev/notes/2026-05-28-marginal-k-plan.md`. Het + marginal-k and
#'     partition-API + marginal-k are deferred (§13 of the plan).
#' @param priorVariant Parameterisation of the geometric `k'` prior under
#'   `likelihoodMode = "marginal_k"`. One of:
#'   * `"conditional"` (default, Model B): `k'_i ~ kObs_i + Geometric(p)`;
#'     the marginal weight for state count `k` is `p (1-p)^(k - kObs_i)`.
#'     This is the shipped marginal-k behaviour.
#'   * `"unconditional"` (Model A): `k'_i ~ 2 + Geometric(p)`, unconditional
#'     on `kObs_i`; the marginal weight for state count `k` is
#'     `p (1-p)^(k - 2)`. The marginal sum still starts at `k = max(2, kObs_i)`
#'     (you cannot have fewer states than observed) but the weight exponent
#'     base is 2 rather than `kObs_i`, with no renormalisation of the truncated
#'     tail. Model A and Model B differ by a per-character factor
#'     `(1-p)^(kObs_i - 2)`. Use `"unconditional"` to match a Model A forward
#'     simulator (e.g. the marginal-k SBC harness).
#'
#'   Only consulted under `likelihoodMode = "marginal_k"` with
#'   `kPrimePrior = "geometric"`; ignored otherwise.
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
    kPrimePrior = "empirical_geometric",
    empiricalNObs = NULL,
    kprimeHyperA = 1,
    kprimeHyperB = 1,
    kprimeAlpha = 1,
    kprimeBeta = 1,
    kprimeLogseriesC = 0.7,
    rateNeoMeanlog = 0,
    rateNeoSdlog = 2,
    qHeterogeneity = FALSE,
    nBetaCat = 4L,
    betaScaleShape = 1,
    betaScaleRate = 1,
    classRateConcentration = 1,
    priorOnClassRateLogSd = c("hyperprior_pooled", "gamma_independent"),
    likelihoodMode = c("sampled_k", "marginal_k"),
    priorVariant = c("conditional", "unconditional")
) {
  coding <- match.arg(coding, c("variable", "informative", "none"))
  kPrimePrior <- match.arg(
    kPrimePrior,
    c("empirical_geometric", "geometric", "beta_geometric", "logseries")
  )
  priorOnClassRateLogSd <- match.arg(priorOnClassRateLogSd)
  likelihoodMode <- match.arg(likelihoodMode)
  priorVariant <- match.arg(priorVariant)

  if (identical(likelihoodMode, "marginal_k")) {
    if (!identical(kPrimePrior, "geometric")) {
      cli::cli_abort(c(
        "{.code likelihoodMode = \"marginal_k\"} requires
         {.code kPrimePrior = \"geometric\"} in v1.",
        i = "Got {.code kPrimePrior = \"{kPrimePrior}\"}.",
        i = "Other arms (empirical_geometric / beta_geometric / logseries)
             are scheduled as §11 follow-ups in
             {.file dev/notes/2026-05-28-marginal-k-plan.md}."
      ))
    }
    if (isTRUE(qHeterogeneity)) {
      cli::cli_abort(c(
        "{.code likelihoodMode = \"marginal_k\"} cannot be combined with
         {.code qHeterogeneity = TRUE}.",
        i = "Het + marginal-k is deferred to v1.x (plan §13).",
        i = "Drop one of the two."
      ))
    }
  }

  # empiricalNObs is only relevant under the empirical_geometric prior; warn
  # if supplied for other priors so the user knows it will be ignored.
  if (kPrimePrior != "empirical_geometric" && !is.null(empiricalNObs)) {
    cli::cli_warn(
      "{.arg empiricalNObs} is ignored when
       {.arg kPrimePrior = \"{kPrimePrior}\"}."
    )
  }
  if (kPrimePrior == "empirical_geometric" && !is.null(empiricalNObs) &&
      !inherits(empiricalNObs, "MkPrimeEmpiricalPrior")) {
    cli::cli_abort(
      "{.arg empiricalNObs} must be an {.cls MkPrimeEmpiricalPrior} object;
       build one with {.fn MkPrimeEmpiricalPrior}."
    )
  }

  # Warn if logseries-specific param is supplied for geometric prior
  if (kPrimePrior != "logseries" && !missing(kprimeLogseriesC)) {
    cli::cli_warn(
      "{.arg kprimeLogseriesC} is ignored when
       {.arg kPrimePrior = \"{kPrimePrior}\"}."
    )
  }

  # Validate beta_geometric hyperparameters
  if (kPrimePrior == "beta_geometric") {
    if (kprimeAlpha <= 0 || kprimeBeta <= 0) {
      cli::cli_abort(
        "{.arg kprimeAlpha} and {.arg kprimeBeta} must be positive."
      )
    }
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
      empiricalNObs = empiricalNObs,
      kprimeHyperA = kprimeHyperA,
      kprimeHyperB = kprimeHyperB,
      kprimeAlpha = kprimeAlpha,
      kprimeBeta = kprimeBeta,
      kprimeLogseriesC = kprimeLogseriesC,
      rateNeoMeanlog = rateNeoMeanlog,
      rateNeoSdlog = rateNeoSdlog,
      qHeterogeneity = qHeterogeneity,
      nBetaCat = as.integer(nBetaCat),
      betaScaleShape = betaScaleShape,
      betaScaleRate = betaScaleRate,
      classRateConcentration = classRateConcentration,
      priorOnClassRateLogSd = priorOnClassRateLogSd,
      likelihoodMode = likelihoodMode,
      priorVariant = priorVariant
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
  if (identical(model$kPrimePrior, "empirical_geometric") &&
      is.null(model$empiricalNObs)) {
    # Lazy-load package data; copy to model so MCMC code can consume it
    # without referencing the package namespace.
    e <- new.env(parent = emptyenv())
    utils::data("empiricalNObs", package = "MkPrime", envir = e)
    model$empiricalNObs <- e$empiricalNObs
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


#' Log P_emp(k) for k = 2, ..., kMax under an `MkPrimeEmpiricalPrior`
#'
#' Returns a numeric vector of length `kMax - 1`; entry `i` is
#' `log P(N_obs = i + 1)`.
#' @param kMax Integer; largest `k` required.
#' @param emp `MkPrimeEmpiricalPrior` object.
#' @return Numeric vector. Entries with zero mass are `-Inf`.
#' @keywords internal
.LogPemp <- function(kMax, emp) {
  bodyLen <- length(emp$body)
  result <- rep_len(-Inf, max(kMax - 1L, 0L))
  bodyEnd <- min(bodyLen, kMax - 1L)
  if (bodyEnd >= 1L) {
    bp <- emp$body[seq_len(bodyEnd)]
    result[seq_len(bodyEnd)] <- ifelse(bp > 0, log(bp), -Inf)
  }
  tailStartK <- emp$tail_start_k
  if (kMax >= tailStartK && emp$tail_start_p > 0 && emp$tail_decay > 0) {
    kk <- seq.int(tailStartK, kMax)
    result[kk - 1L] <- log(emp$tail_start_p) +
                       (kk - tailStartK) * log(emp$tail_decay)
  }
  # Return:
  result
}


#' Log convolution prior on `k'` for the `"empirical_geometric"` option
#'
#' Computes `sum_i log P(k'_i)` where
#' `P(k' = m) = sum_{j=2..m} P_emp(j) * p * (1 - p)^(m - j)`.
#'
#' @param kPrime Integer vector of proposed `k'` values (one per
#'   transformational character).
#' @param emp `MkPrimeEmpiricalPrior` object.
#' @param p Scalar success probability of the geometric on `N_unobs`,
#'   `0 < p < 1`.
#' @return Scalar log density `sum_i log P(k'_i)`. Returns `-Inf` if any
#'   `k'_i < 2`.
#' @keywords internal
.LogPriorEmpiricalGeometric <- function(kPrime, emp, p) {
  if (p <= 0 || p >= 1) {
    # Return:
    return(-Inf)
  }
  if (any(kPrime < 2L)) {
    # Return:
    return(-Inf)
  }
  logP <- log(p)
  log1mP <- log1p(-p)

  kMaxOverall <- max(kPrime)
  logEmp <- .LogPemp(kMaxOverall, emp)

  total <- 0.0
  for (m in kPrime) {
    iVals <- seq.int(2L, m)
    logTerms <- logEmp[iVals - 1L] + logP + (m - iVals) * log1mP
    finite <- is.finite(logTerms)
    if (!any(finite)) {
      # Return:
      return(-Inf)
    }
    mx <- max(logTerms[finite])
    total <- total + mx + log(sum(exp(logTerms[finite] - mx)))
  }
  # Return:
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

    if (identical(model$kPrimePrior, "geometric") ||
        identical(model$kPrimePrior, "empirical_geometric")) {
      if (state$p <= 0 || state$p >= 1) return(-Inf)
    } else if (identical(model$kPrimePrior, "beta_geometric")) {
      ka <- state$kprime_alpha %||% 1.0
      kb <- state$kprime_beta %||% 1.0
      if (ka <= 0 || kb <= 0) return(-Inf)
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

  marginalK <- identical(model$likelihoodMode, "marginal_k")

  if (hasTrans) {
    if (identical(model$kPrimePrior, "geometric")) {
      # k'_i: Geometric(p) shifted by kObs_i
      # P(k'_i = kObs_i + u) = p * (1-p)^u, u = 0, 1, 2, ...
      #
      # Under marginal-k mode, the per-character P(u_i | p) mass is consumed
      # by the marginal-likelihood evaluator (cpp_log_likelihood_marginal),
      # not the prior. The hyperprior on p is kept here unchanged.
      if (!marginalK) {
        u <- state$kPrime[transIdx] - mkd$kObs[transIdx]
        lp <- lp + length(transIdx) * log(state$p) +
              sum(u) * log1p(-state$p)
      }

      # p: Beta hyperprior
      lp <- lp + dbeta(state$p,
                       shape1 = model$kprimeHyperA,
                       shape2 = model$kprimeHyperB,
                       log = TRUE)
    } else if (identical(model$kPrimePrior, "empirical_geometric")) {
      # k'_i = N_obs_i + N_unobs_i, N_obs ~ empirical, N_unobs ~ Geometric(p)
      # Prior on k' is the convolution; truncation k'_i >= kObs_i already
      # enforced above.  N_obs and N_unobs are latent components that sum to
      # k'; we marginalise over their split.
      emp <- model$empiricalNObs
      if (is.null(emp)) {
        # Allow LogPrior to be called on a non-finalised model (e.g. from
        # tests).  Fall back to the package empirical distribution.
        e <- new.env(parent = emptyenv())
        utils::data("empiricalNObs", package = "MkPrime", envir = e)
        emp <- e$empiricalNObs
      }
      lp <- lp + .LogPriorEmpiricalGeometric(
        state$kPrime[transIdx], emp, state$p
      )

      # p: Beta hyperprior
      lp <- lp + dbeta(state$p,
                       shape1 = model$kprimeHyperA,
                       shape2 = model$kprimeHyperB,
                       log = TRUE)
    } else if (identical(model$kPrimePrior, "beta_geometric")) {
      # Per-character p_i marginalized → Beta-Geometric(α, β)
      # log P(k'_i = kObs_i + u | α, β) = lbeta(α+1, β+u) - lbeta(α, β)
      alpha <- state$kprime_alpha %||% 1.0
      beta_ <- state$kprime_beta %||% 1.0
      u <- state$kPrime[transIdx] - mkd$kObs[transIdx]
      lp <- lp + sum(lbeta(alpha + 1, beta_ + u) - lbeta(alpha, beta_))

      # Hyperprior on (α, β): Exponential(1)
      lp <- lp + dexp(alpha, rate = 1, log = TRUE)
      lp <- lp + dexp(beta_, rate = 1, log = TRUE)
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

  # Partition-API (plan v4 §5.1 + §5.4): per-class priors.
  # Only computed when the state carries per-class fields (class_w, etc.).
  # The trivial spec (class_w = 1, length 1) produces zero extra contribution
  # so that the partitioned prior equals the legacy prior (§7b analogue).

  classW <- state$class_w
  if (!is.null(classW)) {
    K <- length(classW)

    # 1. Dirichlet(alpha) on class_w (§5.1)
    #    K == 1: degenerate Dirichlet — contribution is 0.
    if (K > 1) {
      if (any(classW <= 0)) return(-Inf)
      alpha <- model$classRateConcentration %||% 1.0
      logDirConst <- lgamma(K * alpha) - K * lgamma(alpha)
      lp <- lp + logDirConst + (alpha - 1) * sum(log(classW))
    }

    # 2. Per-class ACRV-shape prior (§5.4 / hyperprior extension).
    #    class_rate_log_sd[1] == rate_log_sd: already counted by the scalar
    #    rate_log_sd prior above. Two prior options:
    #
    #      (a) "hyperprior_pooled" (default, K >= 2): non-centred half-normal
    #          hierarchy σ_c = τ · z_c with z_c ~ HN(1) i.i.d. and τ ~ HN(1).
    #          The legacy Gamma(rateLogSdShape, rateLogSdRate) on σ_0 added
    #          above is the wrong prior for the new structure, so we
    #          subtract it and add HN(1) on each z_c (c = 1..K) and on τ.
    #
    #      (b) "gamma_independent" (legacy, K >= 2): each σ_c
    #          independently ~ Gamma. σ_0 prior already counted; loop adds
    #          σ_1..σ_{K-1} contributions.
    #
    #    For K == 1 (linked shape, or unreachable "shape unlinked at K=1"
    #    coerced by .ValidatePartitionArgs): no extra contribution — the
    #    scalar rate_log_sd prior added above covers the single σ.
    classRLS <- state$class_rate_log_sd
    if (!is.null(classRLS) && length(classRLS) > 1L) {
      useHyperprior <-
        identical(model$priorOnClassRateLogSd, "hyperprior_pooled")
      if (useHyperprior) {
        # Subtract legacy Gamma on σ_0 (== rate_log_sd) added above.
        if (state$rate_log_sd > 0) {
          lp <- lp - dgamma(state$rate_log_sd,
                             shape = model$rateLogSdShape,
                             rate = model$rateLogSdRate,
                             log = TRUE)
        }
        # HN(1) on every z_c (c = 1..K): log p(z) = log 2 + dnorm(z;0,1) for z>=0.
        z <- state$class_rate_log_sd_z
        if (is.null(z) || length(z) != length(classRLS)) {
          cli::cli_abort(
            "Hyperprior on {.field class_rate_log_sd} requires
             {.field class_rate_log_sd_z} of matching length."
          )
        }
        if (any(z < 0)) return(-Inf)
        lp <- lp + sum(log(2) + dnorm(z, mean = 0, sd = 1, log = TRUE))
        # HN(1) on τ.
        tau <- state$hyper_tau
        if (is.null(tau)) {
          cli::cli_abort(
            "Hyperprior on {.field class_rate_log_sd} requires
             {.field hyper_tau} scalar."
          )
        }
        if (tau < 0) return(-Inf)
        lp <- lp + log(2) + dnorm(tau, mean = 0, sd = 1, log = TRUE)
      } else {
        # Legacy independent Gamma on each σ_c for c = 2..K.
        extra_sds <- classRLS[-1L]   # classes 2..K
        for (sd_c in extra_sds) {
          if (sd_c < 0) return(-Inf)
          if (sd_c > 0) {
            lp <- lp + dgamma(sd_c,
                               shape = model$rateLogSdShape,
                               rate = model$rateLogSdRate,
                               log = TRUE)
          } else if (model$rateLogSdShape > 1) {
            return(-Inf)
          }
        }
      }
    }
  }

  lp
}


#' @export
print.MkPrimeModel <- function(x, ...) {
  cli::cli_h1("MkPrime Model")

  k_prior_str <- if (identical(x$kPrimePrior, "logseries")) {
    "Logseries (c = {x$kprimeLogseriesC})"
  } else if (identical(x$kPrimePrior, "beta_geometric")) {
    "Beta-Geometric (alpha = {x$kprimeAlpha}, beta = {x$kprimeBeta})"
  } else if (identical(x$kPrimePrior, "empirical_geometric")) {
    paste0(
      "Empirical (N_obs) + Geometric (N_unobs); ",
      "Beta hyperprior on p: a = {x$kprimeHyperA}, b = {x$kprimeHyperB}"
    )
  } else {
    "Geometric (Beta hyperprior: a = {x$kprimeHyperA}, b = {x$kprimeHyperB})"
  }

  het_str <- if (isTRUE(x$qHeterogeneity)) {
    "ON ({x$nBetaCat} bins, beta_scale ~ Gamma({x$betaScaleShape}, {x$betaScaleRate}))"
  } else {
    "OFF"
  }

  lik_mode_str <- x$likelihoodMode %||% "sampled_k"

  cli::cli_ul(c(
    "Coding: {x$coding}",
    "ACRV categories: {x$nCat}",
    "Relabelling correction: {x$relabel}",
    "Tree length prior: Gamma({x$treeLengthShape}, {x$treeLengthRate %||% 'auto'})",
    "rate_loss prior: LogNormal({x$rateLossMeanlog}, {x$rateLossSdlog})",
    "rate_log_sd prior: Gamma({x$rateLogSdShape}, {x$rateLogSdRate})",
    paste0("k' prior: ", k_prior_str),
    "rate_neo prior: LogNormal({x$rateNeoMeanlog}, {x$rateNeoSdlog})",
    paste0("Q-matrix heterogeneity: ", het_str),
    paste0("Likelihood mode: ", lik_mode_str)
  ))
  invisible(x)
}
