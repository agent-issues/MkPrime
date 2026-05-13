
# Compute the total log-likelihood for a MkPrimeData object on a tree.
#
# This is the R-level orchestration function that ties together:
# - JC/MkN pruning (per partition)
# - ACRV rate categories
# - Ascertainment bias correction
# - Mk' relabelling correction
#
# For now this is an internal function used for validation. The MCMC engine
# (Phase 3) will call the C++ functions directly for performance.

#' Compute log-likelihood of morphological data on a tree
#'
#' Evaluates the Mk' likelihood \insertCite{Lewis2001}{MkPrime} via
#' Felsenstein pruning \insertCite{Felsenstein1981}{MkPrime}, with optional
#' ACRV \insertCite{Harrison2014}{MkPrime} and ascertainment bias correction.
#'
#' @references \insertAllCited{}
#' @param tree A `phylo` object (from ape). Must be unrooted or rooted.
#' @param mkd A `MkPrimeData` object.
#' @param kPrime Integer vector of length `mkd$nChar`. The assumed true number
#'   of states for each character. For "known" characters, this is the fixed k.
#'   For "neomorphic", this must be 2. For "transformational", this is the
#'   current MCMC proposal.
#' @param rate_loss Rate asymmetry for neomorphic characters
#'   (1.0 = symmetric). Default 1.0.
#' @param rate_log_sd Standard deviation of the log-normal ACRV distribution.
#'   0 means no rate variation. Default 0.
#' @param nCat Number of ACRV rate categories. Default 6.
#' @param coding Ascertainment bias correction type. `"none"` for no
#'   correction, `"variable"` for conditioning on variable characters,
#'   `"informative"` for conditioning on parsimony-informative characters
#'   (excludes constant + singleton patterns). Default `"variable"`.
#' @param rate_neo Rate scalar for the neomorphic partition relative to
#'   transformational (which is fixed at 1.0). Default 1.0 (equal rates).
#' @param relabel Logical. Apply Mk' relabelling correction for
#'   transformational characters? Set to `FALSE` for standard Mk
#'   likelihood (e.g., for validation against phangorn). Default `TRUE`.
#'
#' @return Scalar log-likelihood.
#' @export
MkpLogLikelihood <- function(tree, mkd,
                               kPrime = NULL,
                               rate_loss = 1.0,
                               rate_log_sd = 0,
                               nCat = 6L,
                               coding = "variable",
                               rate_neo = 1.0,
                               relabel = TRUE) {
  if (!inherits(tree, "phylo")) {
    cli::cli_abort("{.arg tree} must be a {.cls phylo} object.")
  }
  if (!inherits(mkd, "MkPrimeData")) {
    cli::cli_abort("{.arg mkd} must be a {.cls MkPrimeData} object.")
  }
  coding <- match.arg(coding, c("variable", "informative", "none"))

  # Default kPrime: use kObs for transformational, known_k for known, 2 for neo
  if (is.null(kPrime)) {
    kPrime <- mkd$kObs
    knownIdx <- which(mkd$type == "known")
    if (length(knownIdx)) {
      kPrime[knownIdx] <- mkd$known_k[knownIdx]
    }
  }

  # Reorder to canonical preorder for callers without the invariant guarantee.
  # The MCMC hot path calls .MkpLogLikelihood() directly and maintains
  # preorder invariant via .InitState() and all topology proposals.
  tree <- TreeTools::Preorder(tree)

  .MkpLogLikelihood(tree, mkd, kPrime, rate_loss, rate_log_sd,
                    nCat, coding, rate_neo, relabel)
}

# Ecology-aware log-likelihood orchestrator.
#
# Per likelihood call: recompute ecology marginals + edge weights from the
# current tree, then dispatch each partition to the ecology-aware pruning
# variant. The non-ecology fast path remains .MkpLogLikelihood below;
# dispatch between them is the caller's responsibility (typically based on
# model$ecologyAware).
#
# Ascertainment: "variable" coding is supported via per-character mixture-
# aware constant-site probability (one pseudo-character pruning per char).
# "informative" coding (singleton-site exclusion) is a follow-up.
#
# phi: length 1 (magnitudeMode == "global") or kEcology ("per_ecology").
# zMat: nChar x kEcology integer matrix, entries in {0, 1, 2}.
# Per-character constant-site probability under the JC-K ecology mixture.
# zVec is length kEco for one character. Returns scalar P(constant) under
# the mixture induced by zVec, phi, wEdge.
# By JC symmetry the off-diagonal entry of the mixed transition matrix is
# state-independent, so P(constant) = kStates * P(all tips in state 0).
.ConstSiteProbJcEcology <- function(parent, child, edgeLen, nTip, kStates,
                                     rateMultipliers, wEdge,
                                     zVec, phi, mode,
                                     refEcology = 0L, theta = NULL,
                                     pi0 = 0.5) {
  zMat1 <- matrix(as.integer(zVec), nrow = 1, byrow = FALSE)
  states0 <- matrix(0L, nrow = nTip, ncol = 1)
  rootFreqs <- rep(1 / kStates, kStates)
  if (is.null(theta)) theta <- rep(0.5, max(0L, ncol(wEdge) - 1L))
  ll0 <- .PruningJcEcology(parent, child, edgeLen, states0,
                            as.integer(kStates), rootFreqs,
                            rateMultipliers,
                            wEdge, zMat1, phi, as.integer(mode),
                            as.integer(refEcology), as.numeric(theta),
                            as.numeric(pi0))
  kStates * exp(ll0)
}


# Per-character constant-site probability under the MkN ecology mixture.
# The MkN P matrix is asymmetric, so P(all 0) != P(all 1); sum both.
.ConstSiteProbMknEcology <- function(parent, child, edgeLen, nTip,
                                      rateLoss, rateMultipliers, wEdge,
                                      zVec, phi, mode,
                                      refEcology = 0L, theta = NULL,
                                      pi0 = 0.5) {
  zMat1 <- matrix(as.integer(zVec), nrow = 1, byrow = FALSE)
  rootFreqs <- c(rateLoss / (1 + rateLoss), 1 / (1 + rateLoss))
  states0 <- matrix(0L, nrow = nTip, ncol = 1)
  states1 <- matrix(1L, nrow = nTip, ncol = 1)
  if (is.null(theta)) theta <- rep(0.5, max(0L, ncol(wEdge) - 1L))
  refE <- as.integer(refEcology); th <- as.numeric(theta); pi <- as.numeric(pi0)
  ll0 <- .PruningMknEcology(parent, child, edgeLen, states0,
                             rateLoss, rootFreqs, rateMultipliers,
                             wEdge, zMat1, phi, as.integer(mode),
                             refE, th, pi)
  ll1 <- .PruningMknEcology(parent, child, edgeLen, states1,
                             rateLoss, rootFreqs, rateMultipliers,
                             wEdge, zMat1, phi, as.integer(mode),
                             refE, th, pi)
  exp(ll0) + exp(ll1)
}


.MkpEcologyLogLikelihood <- function(tree, mkd, kPrime,
                                       rate_loss, rate_log_sd,
                                       nCat, rate_neo, relabel,
                                       phi, zMat, magnitudeMode = "global",
                                       coding = "none",
                                       refEcology = NULL,
                                       theta = NULL,
                                       pi0 = 0.5) {
  coding <- match.arg(coding, c("none", "variable"))
  # v2 defaults: use mkd$refEcology if not supplied; theta = 0.5 vector.
  if (is.null(refEcology)) refEcology <- mkd$refEcology %||% 0L
  refEcology <- as.integer(refEcology)
  if (is.null(theta)) theta <- rep(0.5, max(0L, mkd$kEcology - 1L))
  theta <- as.numeric(theta)
  pi0 <- as.numeric(pi0)
  # "informative" coding (singleton-site exclusion) is a follow-up.
  if (is.null(mkd$ecology) || is.null(mkd$kEcology)) {
    cli::cli_abort(
      ".MkpEcologyLogLikelihood requires {.arg mkd} with ecology data"
    )
  }
  parent <- tree$edge[, 1]
  child  <- tree$edge[, 2]
  edgeLength <- tree$edge.length
  nTip <- length(tree$tip.label)

  ecologyTip <- as.integer(mkd$ecology)
  ecologyTip[is.na(ecologyTip)] <- -1L
  margMat <- .EcologyNodeMarginals(parent, child, edgeLength,
                                    ecologyTip, mkd$kEcology)
  wEdge <- .EcologyEdgeWeights(margMat, parent, child)

  rates    <- DiscreteLognormalRates(rate_log_sd, nCat)
  modeInt  <- if (identical(magnitudeMode, "global")) 0L else 1L
  totalLoglik <- 0.0

  for (part in mkd$partitions) {
    tipStates <- part$tip_states
    tipStates[is.na(tipStates)] <- -1L
    storage.mode(tipStates) <- "integer"

    zPart <- zMat[part$char_indices, , drop = FALSE]
    storage.mode(zPart) <- "integer"

    if (part$type == "neomorphic") {
      neoEl <- edgeLength * rate_neo
      rootFreqs <- as.numeric(mkn_stationary_freqs(rate_loss))
      ll <- .PruningMknEcology(parent, child, neoEl, tipStates,
                                rate_loss, rootFreqs, rates,
                                wEdge, zPart, phi, modeInt,
                                refEcology, theta, pi0)
      if (coding == "variable") {
        for (c in seq_len(part$nChar)) {
          P_const <- .ConstSiteProbMknEcology(
            parent, child, neoEl, nTip,
            rate_loss, rates, wEdge,
            zPart[c, ], phi, modeInt,
            refEcology, theta, pi0
          )
          ll <- ll - log(1 - P_const)
        }
      }
    } else if (part$type == "known") {
      kStates <- part$k
      rootFreqs <- rep(1.0 / kStates, kStates)
      ll <- .PruningJcEcology(parent, child, edgeLength, tipStates,
                               kStates, rootFreqs, rates,
                               wEdge, zPart, phi, modeInt,
                               refEcology, theta, pi0)
      if (coding == "variable") {
        for (c in seq_len(part$nChar)) {
          P_const <- .ConstSiteProbJcEcology(
            parent, child, edgeLength, nTip, kStates,
            rates, wEdge, zPart[c, ], phi, modeInt,
            refEcology, theta, pi0
          )
          ll <- ll - log(1 - P_const)
        }
      }
    } else {  # transformational: subgroup by kPrime
      kPrimePart <- kPrime[part$char_indices]
      ll <- 0.0
      for (kp in sort(unique(kPrimePart))) {
        cols <- which(kPrimePart == kp)
        subStates <- tipStates[, cols, drop = FALSE]
        rootFreqs <- rep(1.0 / kp, kp)
        zSub <- zPart[cols, , drop = FALSE]
        storage.mode(zSub) <- "integer"
        subLl <- .PruningJcEcology(parent, child, edgeLength, subStates,
                                    kp, rootFreqs, rates,
                                    wEdge, zSub, phi, modeInt,
                                    refEcology, theta, pi0)
        if (coding == "variable") {
          for (c in seq_along(cols)) {
            P_const <- .ConstSiteProbJcEcology(
              parent, child, edgeLength, nTip, kp,
              rates, wEdge, zSub[c, ], phi, modeInt,
              refEcology, theta, pi0
            )
            subLl <- subLl - log(1 - P_const)
          }
        }
        ll <- ll + subLl
      }
      if (relabel) {
        relabelCorr <- mk_prime_relabel_log_batch(
          as.integer(kPrimePart), part$kObs
        )
        ll <- ll + sum(relabelCorr)
      }
    }

    totalLoglik <- totalLoglik + ll
  }

  totalLoglik
}


# Internal fast-path likelihood — no validation, no reorder.
# INVARIANT: tree$edge must already be in canonical preorder. Guaranteed by
# .InitState() and all topology proposals (C++ preorder_weighted_impl).
.MkpLogLikelihood <- function(tree, mkd, kPrime, rate_loss, rate_log_sd,
                               nCat, coding, rate_neo, relabel) {
  parent <- tree$edge[, 1]
  child <- tree$edge[, 2]
  edgeLength <- tree$edge.length
  nTip <- length(tree$tip.label)

  # ACRV rate multipliers
  rates <- DiscreteLognormalRates(rate_log_sd, nCat)

  totalLoglik <- 0.0

  for (part in mkd$partitions) {
    # Prepare tip states: replace NA with -1 for C++
    tipStates <- part$tip_states
    tipStates[is.na(tipStates)] <- -1L
    storage.mode(tipStates) <- "integer"

    nCharPart <- part$nChar

    if (part$type == "neomorphic") {
      # MkN model — apply partition rate scalar
      neoEl <- edgeLength * rate_neo
      rootFreqs <- as.numeric(mkn_stationary_freqs(rate_loss))

      if (rate_log_sd > 0) {
        ll <- pruning_mkn_acrv(parent, child, neoEl,
                               tipStates, rate_loss, rootFreqs, rates)
      } else {
        ll <- pruning_mkn(parent, child, neoEl,
                          tipStates, rate_loss, rootFreqs)
      }

      # Ascertainment correction (uses scaled branches)
      if (coding != "none") {
        puninf <- constant_site_prob_mkn(parent, child, neoEl,
                                         nTip, rate_loss, rootFreqs, rates)
        if (coding == "informative") {
          puninf <- puninf + singleton_site_prob_mkn(
            parent, child, neoEl, nTip, rate_loss, rootFreqs, rates
          )
        }
        ll <- ll - nCharPart * log(1 - puninf)
      }

    } else if (part$type == "known") {
      # Known state space: all chars use the fixed k
      kStates <- part$k
      rootFreqs <- rep(1.0 / kStates, kStates)

      if (rate_log_sd > 0) {
        ll <- pruning_jc_acrv(parent, child, edgeLength,
                              tipStates, kStates, rootFreqs, rates)
      } else {
        ll <- pruning_jc(parent, child, edgeLength,
                         tipStates, kStates, rootFreqs)
      }

      if (coding != "none") {
        puninf <- constant_site_prob_jc(parent, child, edgeLength,
                                        nTip, kStates, rootFreqs, rates)
        if (coding == "informative") {
          puninf <- puninf + singleton_site_prob_jc(
            parent, child, edgeLength, nTip, kStates, rootFreqs, rates
          )
        }
        ll <- ll - nCharPart * log(1 - puninf)
      }

    } else {
      # Transformational: sub-group by current kPrime value.
      # Characters with different k' need different JC(k') matrices.
      kPrimePart <- kPrime[part$char_indices]
      ll <- 0.0

      for (kp in sort(unique(kPrimePart))) {
        cols <- which(kPrimePart == kp)
        subStates <- tipStates[, cols, drop = FALSE]
        nCharSub <- length(cols)

        rootFreqs <- rep(1.0 / kp, kp)

        if (rate_log_sd > 0) {
          subLl <- pruning_jc_acrv(parent, child, edgeLength,
                                    subStates, kp, rootFreqs, rates)
        } else {
          subLl <- pruning_jc(parent, child, edgeLength,
                               subStates, kp, rootFreqs)
        }

        if (coding != "none") {
          puninf <- constant_site_prob_jc(parent, child, edgeLength,
                                          nTip, kp, rootFreqs, rates)
          if (coding == "informative") {
            puninf <- puninf + singleton_site_prob_jc(
              parent, child, edgeLength, nTip, kp, rootFreqs, rates
            )
          }
          subLl <- subLl - nCharSub * log(1 - puninf)
        }

        ll <- ll + subLl
      }

      # Relabelling correction (per character)
      if (relabel) {
        relabelCorr <- mk_prime_relabel_log_batch(
          as.integer(kPrimePart), part$kObs
        )
        ll <- ll + sum(relabelCorr)
      }
    }

    totalLoglik <- totalLoglik + ll
  }

  totalLoglik
}
