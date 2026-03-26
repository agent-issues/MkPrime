# Build partition list from MkPrimeData internals
#
# Groups characters by (type, kObs) for neomorphic/transformational,
# or by (type, known_k) for known-state-space characters.
# Each partition contains tip state data in a format ready for the
# C++ likelihood engine.
#
# @param mkd An MkPrimeData object.
# @return List of partition objects, each with:
#   - type: "neomorphic", "transformational", or "known"
#   - kObs: observed state count for characters in this partition
#   - k: fixed k for "known" partitions; NA for others (inferred during MCMC)
#   - char_indices: integer vector of original character column indices
#   - nChar: number of characters in this partition
#   - tip_states: integer matrix (nTip x nChar), 0-indexed, NA for ambiguous
.build_partitions <- function(mkd) {
  partitions <- list()

  # Neomorphic: all grouped together (all binary, same model)
  neo_idx <- which(mkd$type == "neomorphic")
  if (length(neo_idx)) {
    partitions[[length(partitions) + 1L]] <- list(
      type = "neomorphic",
      kObs = 2L,
      k = NA_integer_,
      char_indices = neo_idx,
      nChar = length(neo_idx),
      tip_states = mkd$matrix[, neo_idx, drop = FALSE]
    )
  }

  # Transformational: group by kObs
  trans_idx <- which(mkd$type == "transformational")
  if (length(trans_idx)) {
    trans_kObs <- mkd$kObs[trans_idx]
    for (ko in sort(unique(trans_kObs))) {
      sel <- trans_idx[trans_kObs == ko]
      partitions[[length(partitions) + 1L]] <- list(
        type = "transformational",
        kObs = ko,
        k = NA_integer_,
        char_indices = sel,
        nChar = length(sel),
        tip_states = mkd$matrix[, sel, drop = FALSE]
      )
    }
  }

  # Known: group by known_k
  known_idx <- which(mkd$type == "known")
  if (length(known_idx)) {
    known_k_vals <- mkd$known_k[known_idx]
    for (kv in sort(unique(known_k_vals))) {
      sel <- known_idx[known_k_vals == kv]
      partitions[[length(partitions) + 1L]] <- list(
        type = "known",
        kObs = max(mkd$kObs[sel]),
        k = kv,
        char_indices = sel,
        nChar = length(sel),
        tip_states = mkd$matrix[, sel, drop = FALSE]
      )
    }
  }

  partitions
}
