
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

  # Reorder to postorder for callers without the invariant guarantee.
  # The MCMC hot path calls .MkpLogLikelihood() directly and maintains
  # the postorder invariant via .InitState() and all topology proposals.
  tree <- ape::reorder.phylo(tree, "postorder")

  .MkpLogLikelihood(tree, mkd, kPrime, rate_loss, rate_log_sd,
                    nCat, coding, rate_neo, relabel)
}

# Internal fast-path likelihood — no validation, no reorder.
# INVARIANT: tree$edge must already be in postorder. This is guaranteed by
# .InitState() and all topology proposals (ProposeNni, ProposeSpr).
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
