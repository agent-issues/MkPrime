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
.BuildPartitions <- function(mkd) {
  partitions <- list()

  # Neomorphic: all grouped together (all binary, same model)
  neoIdx <- which(mkd$type == "neomorphic")
  if (length(neoIdx)) {
    partitions[[length(partitions) + 1L]] <- list(
      type = "neomorphic",
      kObs = 2L,
      k = NA_integer_,
      char_indices = neoIdx,
      nChar = length(neoIdx),
      tip_states = mkd$matrix[, neoIdx, drop = FALSE]
    )
  }

  # Transformational: group by kObs
  transIdx <- which(mkd$type == "transformational")
  if (length(transIdx)) {
    transKObs <- mkd$kObs[transIdx]
    for (ko in sort(unique(transKObs))) {
      sel <- transIdx[transKObs == ko]
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
  knownIdx <- which(mkd$type == "known")
  if (length(knownIdx)) {
    knownKVals <- mkd$known_k[knownIdx]
    for (kv in sort(unique(knownKVals))) {
      sel <- knownIdx[knownKVals == kv]
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
