
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
#'   (excludes every pattern in which fewer than two states each occur
#'   twice or more). Default `"variable"`.
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

  # Audit Issue 1: RB-style partition-rate normalisation. neoScale and
  # transScale make the nChar-weighted mean partition rate equal 1; this
  # gives tree_length a sampler-independent "expected substitutions per
  # character on the average edge" interpretation that matches RevBayes.
  # Degenerate cases (nNeo == 0 or nTrans == 0): rate_neo has no effect,
  # both scales collapse to 1.0.
  partTypes <- vapply(mkd$partitions, function(p) p$type, character(1))
  nCharByPart <- vapply(mkd$partitions, function(p) as.integer(p$nChar),
                        integer(1))
  nNeo   <- sum(nCharByPart[partTypes == "neomorphic"])
  nTrans <- sum(nCharByPart) - nNeo
  if (nNeo == 0L || nTrans == 0L) {
    neoScale   <- 1.0
    transScale <- 1.0
  } else {
    denom <- 1.0 + rate_neo
    nTotal <- as.numeric(nNeo + nTrans)
    neoScale   <- rate_neo / denom * nTotal / nNeo
    transScale <- 1.0       / denom * nTotal / nTrans
  }
  neoEdge   <- edgeLength * neoScale
  transEdge <- edgeLength * transScale

  totalLoglik <- 0.0

  for (part in mkd$partitions) {
    # Prepare tip states: replace NA with -1 for C++
    tipStates <- part$tip_states
    tipStates[is.na(tipStates)] <- -1L
    storage.mode(tipStates) <- "integer"

    nCharPart <- part$nChar

    if (part$type == "neomorphic") {
      # MkN model — apply RB-style neo partition scale (audit Issue 1)
      neoEl <- neoEdge
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
        ll <- ll - .MaskedAscLog1m(tipStates, puninf, parent, child, neoEl,
                                   nTip, 2L, TRUE, rate_loss, rates, coding)
      }

    } else if (part$type == "known") {
      # Mk: fixed k known per character; collapse all unseen states into a
      # single lumped column.  When the partition has chars of mixed kObs,
      # the SAME collapse with kEff = max_char(kObs) + 1 is exact: JC
      # lumpability holds under any state partition, so chars with smaller
      # kObs simply never populate the upper observed columns. We derive
      # kObsMax from per-char data so future partition refactors that mix
      # kObs (or pack neomorphic + transformational together) keep working.
      kStates <- part$k
      kObsMax <- max(part$kObsPerChar)
      if (kObsMax < kStates) {
        if (rate_log_sd > 0) {
          ll <- pruning_jc_acrv_collapsed(parent, child, transEdge,
                                          tipStates, kStates, kObsMax, rates)
        } else {
          ll <- pruning_jc_collapsed(parent, child, transEdge,
                                     tipStates, kStates, kObsMax)
        }

        if (coding != "none") {
          puninf <- constant_site_prob_jc_collapsed(
            parent, child, transEdge, nTip, kStates, kObsMax, rates
          )
          if (coding == "informative") {
            puninf <- puninf + uninf_nonconst_prob_jc(
              parent, child, transEdge, nTip, kStates, rates
            )
          }
          ll <- ll - .MaskedAscLog1m(tipStates, puninf, parent, child,
                                     transEdge, nTip, kStates, FALSE, 1,
                                     rates, coding)
        }
      } else {
        rootFreqs <- rep(1.0 / kStates, kStates)
        if (rate_log_sd > 0) {
          ll <- pruning_jc_acrv(parent, child, transEdge,
                                tipStates, kStates, rootFreqs, rates)
        } else {
          ll <- pruning_jc(parent, child, transEdge,
                           tipStates, kStates, rootFreqs)
        }

        if (coding != "none") {
          puninf <- constant_site_prob_jc(parent, child, transEdge,
                                          nTip, kStates, rootFreqs, rates)
          if (coding == "informative") {
            puninf <- puninf + uninf_nonconst_prob_jc(
              parent, child, transEdge, nTip, kStates, rates
            )
          }
          ll <- ll - .MaskedAscLog1m(tipStates, puninf, parent, child,
                                     transEdge, nTip, kStates, FALSE, 1,
                                     rates, coding)
        }
      }

    } else {
      # Transformational: sub-group by current kPrime value.
      # Characters with different k' need different JC(k') matrices.
      kPrimePart <- kPrime[part$char_indices]
      ll <- 0.0

      kObsPerCharPart <- part$kObsPerChar
      for (kp in sort(unique(kPrimePart))) {
        cols <- which(kPrimePart == kp)
        subStates <- tipStates[, cols, drop = FALSE]
        nCharSub <- length(cols)

        # Derive kObsMax from per-char kObs of *this sub-batch* so the
        # dispatch survives future partition refactors that mix kObs.
        kObsMaxSub <- max(kObsPerCharPart[cols])
        if (kObsMaxSub < kp) {
          # State-collapse path: kEff = kObsMaxSub + 1.
          if (rate_log_sd > 0) {
            subLl <- pruning_jc_acrv_collapsed(parent, child, transEdge,
                                               subStates, kp, kObsMaxSub, rates)
          } else {
            subLl <- pruning_jc_collapsed(parent, child, transEdge,
                                          subStates, kp, kObsMaxSub)
          }

          if (coding != "none") {
            puninf <- constant_site_prob_jc_collapsed(
              parent, child, transEdge, nTip, kp, kObsMaxSub, rates
            )
            if (coding == "informative") {
              puninf <- puninf + uninf_nonconst_prob_jc(
                parent, child, transEdge, nTip, kp, rates
              )
            }
            subLl <- subLl - .MaskedAscLog1m(subStates, puninf, parent, child,
                                             transEdge, nTip, kp, FALSE, 1,
                                             rates, coding)
          }
        } else {
          rootFreqs <- rep(1.0 / kp, kp)
          if (rate_log_sd > 0) {
            subLl <- pruning_jc_acrv(parent, child, transEdge,
                                      subStates, kp, rootFreqs, rates)
          } else {
            subLl <- pruning_jc(parent, child, transEdge,
                                 subStates, kp, rootFreqs)
          }

          if (coding != "none") {
            puninf <- constant_site_prob_jc(parent, child, transEdge,
                                            nTip, kp, rootFreqs, rates)
            if (coding == "informative") {
              puninf <- puninf + uninf_nonconst_prob_jc(
                parent, child, transEdge, nTip, kp, rates
              )
            }
            subLl <- subLl - .MaskedAscLog1m(subStates, puninf, parent, child,
                                             transEdge, nTip, kp, FALSE, 1,
                                             rates, coding)
          }
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


# Sum over missing-data masks of n_m * log(1 - P_m): the amount the
# ascertainment correction subtracts. A character's ?/- tips are marginalised,
# as it was kept for varying among its observed tips. `puninf` is P for the
# characters with every tip observed.
.MaskedAscLog1m <- function(tipStates, puninf, parent, child, edgeLength,
                            nTip, kStates, neomorphic, rateLoss, rates,
                            coding) {
  missing <- tipStates < 0L
  maskKey <- apply(missing, 2, function(x) paste(which(x), collapse = " "))
  counts <- table(maskKey)
  keys <- names(counts)
  p <- vapply(keys, function(key) {
    if (!nzchar(key)) {
      puninf
    } else {
      asc_site_prob_missing(parent, child, edgeLength, nTip, kStates,
                            neomorphic, rateLoss, rates,
                            missing[, match(key, maskKey)],
                            coding == "informative")
    }
  }, double(1))
  # Return:
  sum(as.integer(counts) * log(1 - p))
}
