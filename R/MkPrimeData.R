#' Prepare morphological data for Mk' analysis
#'
#' Accepts a `phyDat` object and classifies characters into three types:
#' neomorphic (asymmetric binary), transformational (Mk' with inferred k'),
#' or known state space (standard Mk with user-specified k).
#'
#' @param data A `phyDat` object (from \pkg{ape} or \pkg{TreeTools}).
#' @param neomorphic Integer vector of character indices to treat as neomorphic
#'   (asymmetric binary, 0 = absent, 1 = present). Default: none.
#' @param knownStates Named integer vector specifying the known true number of
#'   states for specific characters. Names are character indices (as strings),
#'   values are the true k. Characters listed here use standard Mk(k) with no
#'   k' inference.
#'
#' @details
#' Characters not listed in `neomorphic` or `knownStates` are classified as
#' **transformational** and will have their true state count k' inferred under
#' the Mk' model.
#'
#' For each character, the number of observed states (kObs) is computed by
#' counting distinct non-ambiguous states across all taxa.
#'
#' @return An object of class `"MkPrimeData"`, a list with components:
#' \describe{
#'   \item{matrix}{Character matrix (nTip x nChar), character states as integers
#'     (0-indexed). `NA` for missing/ambiguous data.}
#'   \item{nTip}{Number of taxa.}
#'   \item{nChar}{Number of characters.}
#'   \item{taxon_names}{Character vector of taxon names.}
#'   \item{type}{Character vector of length nChar: `"neomorphic"`,
#'     `"transformational"`, or `"known"`.}
#'   \item{kObs}{Integer vector: observed state count per character.}
#'   \item{known_k}{Integer vector: known k per character (`NA` for
#'     non-"known" characters).}
#'   \item{levels}{Character state levels from the original phyDat.}
#' }
#'
#' @examples
#' \dontrun{
#' library("TreeTools")
#' data(Lobo.phy)  # example morphological dataset
#' mkd <- MkPrimeData(Lobo.phy)
#' summary(mkd)
#' }
#'
#' @export
MkPrimeData <- function(data,
                         neomorphic = integer(0),
                         knownStates = integer(0)) {
  if (!inherits(data, "phyDat")) {
    cli::cli_abort("{.arg data} must be a {.cls phyDat} object.")
  }

  charMatrix <- .PhyDatToIntMatrix(data)
  nTip <- nrow(charMatrix)
  nChar <- ncol(charMatrix)
  taxonNames <- rownames(charMatrix)
  levels <- attr(data, "levels")

  # Validate neomorphic indices
  neomorphic <- as.integer(neomorphic)
  if (length(neomorphic) && (any(neomorphic < 1) || any(neomorphic > nChar))) {
    cli::cli_abort(
      "{.arg neomorphic} indices must be between 1 and {nChar}."
    )
  }

  # Validate knownStates
  if (length(knownStates)) {
    ksIdx <- as.integer(names(knownStates))
    if (any(is.na(ksIdx)) || any(ksIdx < 1) || any(ksIdx > nChar)) {
      cli::cli_abort(
        "{.arg knownStates} names must be character indices between 1 and
        {nChar}."
      )
    }
    if (any(knownStates < 2L)) {
      cli::cli_abort("{.arg knownStates} values must be >= 2.")
    }
  } else {
    ksIdx <- integer(0)
  }

  # Check for overlap
  overlap <- intersect(neomorphic, ksIdx)
  if (length(overlap)) {
    cli::cli_abort(
      "Character{?s} {overlap} appear{?s/} in both {.arg neomorphic} and
      {.arg knownStates}."
    )
  }

  # Compute kObs per character
  kObs <- apply(charMatrix, 2, function(col) {
    length(unique(col[!is.na(col)]))
  })

  # Drop invariant characters (kObs <= 1) before further validation
  invariant <- which(kObs <= 1L)
  if (length(invariant)) {
    cli::cli_warn(
      "Dropping {length(invariant)} invariant character{?s} (kObs <= 1):
      {?column/columns} {invariant}."
    )
    keep <- setdiff(seq_len(nChar), invariant)
    if (!length(keep)) {
      cli::cli_abort("No variable characters remain after dropping invariants.")
    }
    charMatrix <- charMatrix[, keep, drop = FALSE]
    kObs <- kObs[keep]

    # Remap neomorphic and knownStates indices
    oldToNew <- rep(NA_integer_, nChar)
    oldToNew[keep] <- seq_along(keep)
    neomorphic <- as.integer(na.omit(oldToNew[neomorphic]))
    if (length(ksIdx)) {
      survived <- !is.na(oldToNew[ksIdx])
      newKsIdx <- oldToNew[ksIdx[survived]]
      knownStates <- knownStates[survived]
      names(knownStates) <- as.character(newKsIdx)
      ksIdx <- newKsIdx
    }
    nChar <- length(keep)
  }

  # Validate neomorphic characters are binary (after dropping invariants)
  if (length(neomorphic)) {
    neoNotBinary <- neomorphic[kObs[neomorphic] != 2L]
    if (length(neoNotBinary)) {
      cli::cli_warn(
        "Neomorphic character{?s} {neoNotBinary} {?has/have} kObs != 2.
        Neomorphic model assumes exactly 2 states."
      )
    }
  }

  # Validate knownStates >= kObs (after dropping invariants)
  if (length(knownStates)) {
    tooSmall <- ksIdx[knownStates < kObs[ksIdx]]
    if (length(tooSmall)) {
      cli::cli_abort(
        "Character{?s} {tooSmall}: {.arg knownStates} is less than the
        observed state count (kObs). k must be >= kObs."
      )
    }
  }

  # Classify characters
  type <- rep("transformational", nChar)
  type[neomorphic] <- "neomorphic"
  if (length(ksIdx)) {
    type[ksIdx] <- "known"
  }

  # Build known_k vector (NA for non-known characters)
  knownK <- rep(NA_integer_, nChar)
  if (length(knownStates)) {
    knownK[ksIdx] <- as.integer(knownStates)
  }

  mkd <- structure(
    list(
      matrix = charMatrix,
      nTip = nTip,
      nChar = nChar,
      taxon_names = taxonNames,
      type = type,
      kObs = kObs,
      known_k = knownK,
      levels = levels,
      phyDat = data
    ),
    class = "MkPrimeData"
  )

  mkd$partitions <- .BuildPartitions(mkd)
  mkd
}


#' @export
print.MkPrimeData <- function(x, ...) {
  typeCounts <- table(x$type)
  cli::cli_h2("MkPrimeData: {x$nTip} taxa, {x$nChar} characters")
  for (tp in names(typeCounts)) {
    cli::cli_bullets(c("*" = "{typeCounts[[tp]]} {tp}"))
  }
  nPart <- length(x$partitions)
  cli::cli_bullets(c(
    "i" = "{nPart} partition{?s} (grouped by type and kObs)"
  ))
  kObsRange <- range(x$kObs)
  cli::cli_bullets(c(
    "i" = "kObs range: {kObsRange[1]}\u2013{kObsRange[2]}"
  ))
  invisible(x)
}


#' @export
summary.MkPrimeData <- function(object, ...) {
  cat("MkPrimeData:", object$nTip, "taxa,", object$nChar, "characters\n\n")
  cat("Character types:\n")
  print(table(object$type))
  cat("\nkObs distribution:\n")
  print(table(kObs = object$kObs))
  invisible(object)
}


#' Auto-detect neomorphic character indices
#'
#' Identifies binary characters whose states are `0` (absent) and `1`
#' (present), which are candidates for the asymmetric neomorphic model.
#' Characters that contain state `0` alongside states other than `1` are
#' **not** flagged as neomorphic (they remain transformational by default).
#'
#' @param data A `phyDat` object.
#' @return Integer vector of character indices suitable for the
#'   `neomorphic` argument of [MkPrimeData()].
#' @export
AutoDetectNeomorphic <- function(data) {
  if (!inherits(data, "phyDat")) {
    cli::cli_abort("{.arg data} must be a {.cls phyDat} object.")
  }
  mat <- .PhyDatToIntMatrix(data)
  # Real levels (gap excluded) — matches the 0-based indexing in .PhyDatToIntMatrix
  lvls <- attr(data, "levels")
  lvls <- lvls[lvls != "-"]
  neo <- integer(0)
  for (j in seq_len(ncol(mat))) {
    states <- unique(mat[, j])
    states <- states[!is.na(states)]
    # Map 0-indexed matrix values to level labels
    labels <- lvls[states + 1L]
    # Neomorphic: exactly levels "0" and "1"
    if (length(labels) == 2L && all(sort(labels) == c("0", "1"))) {
      neo <- c(neo, j)
    }
  }
  neo
}


# Convert phyDat to an integer matrix (0-indexed states, NA for ambiguous)
#
# Each row is a taxon, each column is a character. States are integers
# 0, 1, ..., k-1 for the real data levels. Ambiguous states (e.g. "?") and
# the gap character ("-") are NA.
#
# Some phyDat objects (e.g. from TreeTools::ReadAsPhyDat) include "-" as a
# genuine entry in `levels`, where its contrast row maps to a single column
# rather than to all columns. This makes "-" look like a real state (state 0)
# and shifts all actual states up by one, causing out-of-bounds access in the
# C++ pruning code. We handle this by building an explicit mapping from
# contrast column to clean 0-based state, excluding the gap character.
.PhyDatToIntMatrix <- function(data) {
  contrast  <- attr(data, "contrast")
  levels    <- attr(data, "levels")
  weight    <- attr(data, "weight")
  index     <- attr(data, "index")
  nr        <- attr(data, "nr")
  taxa      <- names(data)
  nTip      <- length(taxa)

  # Real data levels: exclude the gap character "-"
  realMask <- levels != "-"
  realCols <- which(realMask)          # 1-indexed positions in `levels`

  # Map: contrast column index (1-based) -> clean 0-based state, NA for gap
  colToState <- rep(NA_integer_, length(levels))
  colToState[realCols] <- seq_along(realCols) - 1L

  # For each allLevel row in the contrast matrix, determine if it resolves to
  # a single real state (exactly one 1 in a real column) or is ambiguous/gap
  singleState <- apply(contrast, 1, function(row) {
    which1 <- which(row == 1)
    if (length(which1) == 1L) colToState[which1] else NA_integer_
  })

  # Build the unique-pattern matrix (nTip x nr)
  patternMat <- matrix(NA_integer_, nrow = nTip, ncol = nr)
  for (i in seq_len(nTip)) {
    patternMat[i, ] <- singleState[data[[i]]]
  }

  # Expand to full character matrix using index
  charMatrix <- patternMat[, index, drop = FALSE]
  rownames(charMatrix) <- taxa

  # Per-character: remap states to contiguous 0-based integers.
  # A character may only use a subset of NEXUS state labels (e.g. {1,2,3}
  # with "0" never observed), leaving state values that exceed kObs - 1 and
  # cause out-of-bounds access in the C++ pruning arrays. For symmetric models
  # (Mk, Mk') the label assignment is arbitrary, so this remapping is valid.
  for (j in seq_len(ncol(charMatrix))) {
    col <- charMatrix[, j]
    notNa <- !is.na(col)
    observed <- sort(unique(col[notNa]))
    if (length(observed) == 0L ||
        identical(observed, seq(0L, length(observed) - 1L))) {
      next  # already contiguous and 0-based
    }
    charMatrix[notNa, j] <- match(col[notNa], observed) - 1L
  }

  charMatrix
}
