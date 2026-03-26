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
#'   correction, `"variable"` for conditioning on variable characters.
#'   Default `"variable"`.
#' @param relabel Logical. Apply Mk' relabelling correction for
#'   transformational characters? Set to `FALSE` for standard Mk
#'   likelihood (e.g., for validation against phangorn). Default `TRUE`.
#'
#' @return Scalar log-likelihood.
#' @export
mkp_loglikelihood <- function(tree, mkd,
                               kPrime = NULL,
                               rate_loss = 1.0,
                               rate_log_sd = 0,
                               nCat = 6L,
                               coding = "variable",
                               relabel = TRUE) {
  if (!inherits(tree, "phylo")) {
    cli::cli_abort("{.arg tree} must be a {.cls phylo} object.")
  }
  if (!inherits(mkd, "MkPrimeData")) {
    cli::cli_abort("{.arg mkd} must be a {.cls MkPrimeData} object.")
  }
  coding <- match.arg(coding, c("variable", "none"))

  # Default kPrime: use kObs for transformational, known_k for known, 2 for neo
  if (is.null(kPrime)) {
    kPrime <- mkd$kObs
    known_idx <- which(mkd$type == "known")
    if (length(known_idx)) {
      kPrime[known_idx] <- mkd$known_k[known_idx]
    }
  }

  # Reorder tree to postorder
  tree <- ape::reorder.phylo(tree, "postorder")
  parent <- tree$edge[, 1]
  child <- tree$edge[, 2]
  edge_length <- tree$edge.length
  nTip <- length(tree$tip.label)

  # ACRV rate multipliers
  rates <- discrete_lognormal_rates(rate_log_sd, nCat)

  total_loglik <- 0.0

  for (part in mkd$partitions) {
    # Prepare tip states: replace NA with -1 for C++
    tip_states <- part$tip_states
    tip_states[is.na(tip_states)] <- -1L
    storage.mode(tip_states) <- "integer"

    nCharPart <- part$nChar

    if (part$type == "neomorphic") {
      # MkN model
      root_freqs <- as.numeric(mkn_stationary_freqs(rate_loss))

      if (rate_log_sd > 0) {
        ll <- pruning_mkn_acrv(parent, child, edge_length,
                               tip_states, rate_loss, root_freqs, rates)
      } else {
        ll <- pruning_mkn(parent, child, edge_length,
                          tip_states, rate_loss, root_freqs)
      }

      # Ascertainment correction
      if (coding == "variable") {
        pconst <- constant_site_prob_mkn(parent, child, edge_length,
                                         nTip, rate_loss, root_freqs, rates)
        ll <- ll - nCharPart * log(1 - pconst)
      }

    } else {
      # JC model (transformational or known)
      kStates <- part$kObs
      if (part$type == "known") {
        kStates <- part$k
      } else {
        # For transformational, use the max kPrime in this partition
        # (all chars share the same kObs but may have different kPrime;
        # however within a partition, we compute with a single k for the
        # rate matrix. For now, use kObs as the matrix size.)
        # TODO: handle per-character kPrime properly in the pruning
        kStates <- part$kObs
      }
      root_freqs <- rep(1.0 / kStates, kStates)

      if (rate_log_sd > 0) {
        ll <- pruning_jc_acrv(parent, child, edge_length,
                              tip_states, kStates, root_freqs, rates)
      } else {
        ll <- pruning_jc(parent, child, edge_length,
                         tip_states, kStates, root_freqs)
      }

      # Ascertainment correction
      if (coding == "variable") {
        pconst <- constant_site_prob_jc(parent, child, edge_length,
                                        nTip, kStates, root_freqs, rates)
        ll <- ll - nCharPart * log(1 - pconst)
      }

      # Relabelling correction for transformational characters (Mk' model)
      if (relabel && part$type == "transformational") {
        kPrime_part <- kPrime[part$char_indices]
        relabel_corr <- mk_prime_relabel_log_batch(
          as.integer(kPrime_part), part$kObs
        )
        ll <- ll + sum(relabel_corr)
      }
    }

    total_loglik <- total_loglik + ll
  }

  total_loglik
}
