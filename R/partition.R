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
#   - unique_tip_states: integer matrix (nTip x nUnique), deduplicated columns
#   - pattern_index: integer vector (nChar), 0-based index into unique_tip_states
.BuildPartitions <- function(mkd) {
  partitions <- list()

  # Neomorphic: all grouped together (all binary, same model)
  neoIdx <- which(mkd$type == "neomorphic")
  if (length(neoIdx)) {
    ts <- mkd$matrix[, neoIdx, drop = FALSE]
    pc <- .ComputePatternIndex(ts)
    partitions[[length(partitions) + 1L]] <- list(
      type              = "neomorphic",
      kObs              = 2L,
      k                 = NA_integer_,
      char_indices      = neoIdx,
      nChar             = length(neoIdx),
      tip_states        = ts,
      unique_tip_states = pc$unique_tip_states,
      pattern_index     = pc$pattern_index
    )
  }

  # Transformational: group by kObs
  transIdx <- which(mkd$type == "transformational")
  if (length(transIdx)) {
    transKObs <- mkd$kObs[transIdx]
    for (ko in sort(unique(transKObs))) {
      sel <- transIdx[transKObs == ko]
      ts  <- mkd$matrix[, sel, drop = FALSE]
      pc  <- .ComputePatternIndex(ts)
      partitions[[length(partitions) + 1L]] <- list(
        type              = "transformational",
        kObs              = ko,
        k                 = NA_integer_,
        char_indices      = sel,
        nChar             = length(sel),
        tip_states        = ts,
        unique_tip_states = pc$unique_tip_states,
        pattern_index     = pc$pattern_index
      )
    }
  }

  # Known: group by known_k
  knownIdx <- which(mkd$type == "known")
  if (length(knownIdx)) {
    knownKVals <- mkd$known_k[knownIdx]
    for (kv in sort(unique(knownKVals))) {
      sel <- knownIdx[knownKVals == kv]
      ts  <- mkd$matrix[, sel, drop = FALSE]
      pc  <- .ComputePatternIndex(ts)
      partitions[[length(partitions) + 1L]] <- list(
        type              = "known",
        kObs              = max(mkd$kObs[sel]),
        k                 = kv,
        char_indices      = sel,
        nChar             = length(sel),
        tip_states        = ts,
        unique_tip_states = pc$unique_tip_states,
        pattern_index     = pc$pattern_index
      )
    }
  }

  partitions
}


# M-172: Compute unique tip-state patterns within a character matrix.
#
# For JC/Mk' models, two characters with the same tip-state column (after
# 0-based state remapping) have identical likelihoods under any candidate k'.
# This identifies unique columns so the Gibbs sweep can evaluate each pattern
# once and scatter the result to all characters sharing it.
#
# NA placement is part of the key: c(0,NA,1) != c(0,1,NA).
#
# @param tip_states Integer matrix (nTip x nChar), 0-indexed, NA for ambiguous.
# @return List with:
#   unique_tip_states: integer matrix (nTip x nUnique)
#   pattern_index: integer vector (nChar), 0-based index into unique_tip_states
.ComputePatternIndex <- function(tip_states) {
  nChar <- ncol(tip_states)
  if (nChar <= 1L) {
    return(list(
      unique_tip_states = tip_states,
      pattern_index     = 0L
    ))
  }
  # Serialize each column: NA → -1 so its position is part of the key
  keys <- apply(tip_states, 2L, function(col) {
    paste(ifelse(is.na(col), -1L, as.integer(col)), collapse = "\x1f")
  })
  first_occ     <- !duplicated(keys)
  unique_keys   <- keys[first_occ]
  pattern_index <- match(keys, unique_keys) - 1L  # 0-based
  unique_cols   <- which(first_occ)
  list(
    unique_tip_states = tip_states[, unique_cols, drop = FALSE],
    pattern_index     = pattern_index
  )
}
