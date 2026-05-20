# Build partition list from MkPrimeData internals
#
# Groups characters by (type, kObs) for neomorphic/transformational,
# or by (type, known_k) for known-state-space characters. When a user
# partition is supplied, sub-groups further by (classIdx, type, kObs)
# so each emitted partition belongs to exactly one user class.
# Each partition contains tip state data in a format ready for the
# C++ likelihood engine.
#
# @param mkd An MkPrimeData object.
# @param partition Optional integer vector of length mkd$nChar (post
#   invariant drop) assigning each character to a user class with values
#   forming a contiguous range 1:nClasses. `NULL` (the default) emits a
#   single user class containing every character (classIdx = 1L).
# @return List of partition objects, each with:
#   - type: "neomorphic", "transformational", or "known"
#   - kObs: observed state count for characters in this partition
#   - k: fixed k for "known" partitions; NA for others (inferred during MCMC)
#   - char_indices: integer vector of original character column indices
#   - nChar: number of characters in this partition
#   - tip_states: integer matrix (nTip x nChar), 0-indexed, NA for ambiguous
#   - unique_tip_states: integer matrix (nTip x nUnique), deduplicated columns
#   - pattern_index: integer vector (nChar), 0-based index into unique_tip_states
#   - classIdx: integer scalar, user-class membership (1L when partition = NULL)
.BuildPartitions <- function(mkd, partition = NULL) {
  if (is.null(partition)) {
    partition <- rep(1L, mkd$nChar)
  } else {
    partition <- as.integer(partition)
  }
  nClasses <- max(partition)

  # Sub-group within each (classIdx, type, kObs) bucket.
  # A user class spanning multiple kObs values yields multiple PartInfos
  # carrying the same classIdx; this is necessary because JC pruning needs
  # a fixed k per partition (load-bearing in the C++ loop, plan v4 §2).
  .emit <- function(sel, ptype, kObsVal, kVal, classIdx) {
    if (!length(sel)) return(NULL)
    ts <- mkd$matrix[, sel, drop = FALSE]
    pc <- .ComputePatternIndex(ts)
    list(
      type              = ptype,
      kObs              = kObsVal,
      kObsPerChar       = mkd$kObs[sel],
      k                 = kVal,
      char_indices      = sel,
      nChar             = length(sel),
      tip_states        = ts,
      unique_tip_states = pc$unique_tip_states,
      pattern_index     = pc$pattern_index,
      classIdx          = classIdx
    )
  }

  partitions <- list()
  for (cls in seq_len(nClasses)) {
    classMask <- partition == cls

    # Neomorphic: all chars in this class grouped together (binary, same model)
    neoIdx <- which(mkd$type == "neomorphic" & classMask)
    if (length(neoIdx)) {
      partitions[[length(partitions) + 1L]] <-
        .emit(neoIdx, "neomorphic", 2L, NA_integer_, cls)
    }

    # Transformational: group by kObs within this class
    transIdx <- which(mkd$type == "transformational" & classMask)
    if (length(transIdx)) {
      transKObs <- mkd$kObs[transIdx]
      for (ko in sort(unique(transKObs))) {
        sel <- transIdx[transKObs == ko]
        partitions[[length(partitions) + 1L]] <-
          .emit(sel, "transformational", ko, NA_integer_, cls)
      }
    }

    # Known: group by known_k within this class
    knownIdx <- which(mkd$type == "known" & classMask)
    if (length(knownIdx)) {
      knownKVals <- mkd$known_k[knownIdx]
      for (kv in sort(unique(knownKVals))) {
        sel <- knownIdx[knownKVals == kv]
        partitions[[length(partitions) + 1L]] <-
          .emit(sel, "known", max(mkd$kObs[sel]), kv, cls)
      }
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
