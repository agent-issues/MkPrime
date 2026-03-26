#' Prepare morphological data for Mk' analysis
#'
#' Accepts a `phyDat` object and classifies characters into three types:
#' neomorphic (asymmetric binary), transformational (Mk' with inferred k'),
#' or known state space (standard Mk with user-specified k).
#'
#' @param data A `phyDat` object (from ape or TreeTools).
#' @param neomorphic Integer vector of character indices to treat as neomorphic
#'   (asymmetric binary, 0 = absent, 1 = present). Default: none.
#' @param known_states Named integer vector specifying the known true number of
#'   states for specific characters. Names are character indices (as strings),
#'   values are the true k. Characters listed here use standard Mk(k) with no
#'   k' inference.
#'
#' @details
#' Characters not listed in `neomorphic` or `known_states` are classified as
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
#' library(TreeTools)
#' data(Lobo.phy)  # example morphological dataset
#' mkd <- MkPrimeData(Lobo.phy)
#' summary(mkd)
#' }
#'
#' @export
MkPrimeData <- function(data,
                         neomorphic = integer(0),
                         known_states = integer(0)) {
  if (!inherits(data, "phyDat")) {
    cli::cli_abort("{.arg data} must be a {.cls phyDat} object.")
  }

  char_matrix <- .PhyDatToIntMatrix(data)
  nTip <- nrow(char_matrix)
  nChar <- ncol(char_matrix)
  taxon_names <- rownames(char_matrix)
  levels <- attr(data, "levels")

  # Validate neomorphic indices
  neomorphic <- as.integer(neomorphic)
  if (length(neomorphic) && (any(neomorphic < 1) || any(neomorphic > nChar))) {
    cli::cli_abort(
      "{.arg neomorphic} indices must be between 1 and {nChar}."
    )
  }

  # Validate known_states
  if (length(known_states)) {
    ks_idx <- as.integer(names(known_states))
    if (any(is.na(ks_idx)) || any(ks_idx < 1) || any(ks_idx > nChar)) {
      cli::cli_abort(
        "{.arg known_states} names must be character indices between 1 and
        {nChar}."
      )
    }
    if (any(known_states < 2L)) {
      cli::cli_abort("{.arg known_states} values must be >= 2.")
    }
  } else {
    ks_idx <- integer(0)
  }

  # Check for overlap
  overlap <- intersect(neomorphic, ks_idx)
  if (length(overlap)) {
    cli::cli_abort(
      "Character{?s} {overlap} appear{?s/} in both {.arg neomorphic} and
      {.arg known_states}."
    )
  }

  # Compute kObs per character
  kObs <- apply(char_matrix, 2, function(col) {
    length(unique(col[!is.na(col)]))
  })

  # Validate neomorphic characters are binary
  neo_not_binary <- neomorphic[kObs[neomorphic] != 2L]
  if (length(neo_not_binary)) {
    cli::cli_warn(
      "Neomorphic character{?s} {neo_not_binary} {?has/have} kObs != 2.
      Neomorphic model assumes exactly 2 states."
    )
  }

  # Validate known_states >= kObs
  if (length(known_states)) {
    too_small <- ks_idx[known_states < kObs[ks_idx]]
    if (length(too_small)) {
      cli::cli_abort(
        "Character{?s} {too_small}: {.arg known_states} is less than the
        observed state count (kObs). k must be >= kObs."
      )
    }
  }

  # Classify characters
  type <- rep("transformational", nChar)
  type[neomorphic] <- "neomorphic"
  if (length(ks_idx)) {
    type[ks_idx] <- "known"
  }

  # Build known_k vector (NA for non-known characters)
  known_k <- rep(NA_integer_, nChar)
  if (length(known_states)) {
    known_k[ks_idx] <- as.integer(known_states)
  }

  mkd <- structure(
    list(
      matrix = char_matrix,
      nTip = nTip,
      nChar = nChar,
      taxon_names = taxon_names,
      type = type,
      kObs = kObs,
      known_k = known_k,
      levels = levels
    ),
    class = "MkPrimeData"
  )

  mkd$partitions <- .build_partitions(mkd)
  mkd
}


#' @export
print.MkPrimeData <- function(x, ...) {
  type_counts <- table(x$type)
  cli::cli_h2("MkPrimeData: {x$nTip} taxa, {x$nChar} characters")
  for (tp in names(type_counts)) {
    cli::cli_bullets(c("*" = "{type_counts[[tp]]} {tp}"))
  }
  nPart <- length(x$partitions)
  cli::cli_bullets(c(
    "i" = "{nPart} partition{?s} (grouped by type and kObs)"
  ))
  kObs_range <- range(x$kObs)
  cli::cli_bullets(c(
    "i" = "kObs range: {kObs_range[1]}\u2013{kObs_range[2]}"
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


# Convert phyDat to an integer matrix (0-indexed states, NA for ambiguous)
#
# Each row is a taxon, each column is a character. States are integers
# corresponding to positions in `levels`. Ambiguous states (e.g. "?")
# are NA.
.PhyDatToIntMatrix <- function(data) {
  contrast <- attr(data, "contrast")
  levels <- attr(data, "levels")
  weight <- attr(data, "weight")
  index <- attr(data, "index")
  nr <- attr(data, "nr")
  taxa <- names(data)
  nTip <- length(taxa)
  nLevels <- length(levels)

  # For each row of the contrast matrix, determine if it's a single state
  # (exactly one 1) or ambiguous (multiple 1s)
  single_state <- apply(contrast, 1, function(row) {
    which1 <- which(row == 1)
    if (length(which1) == 1L) which1 - 1L else NA_integer_  # 0-indexed
  })

  # Build the unique-pattern matrix (nTip x nr)
  pattern_mat <- matrix(NA_integer_, nrow = nTip, ncol = nr)
  for (i in seq_len(nTip)) {
    pattern_mat[i, ] <- single_state[data[[i]]]
  }

  # Expand to full character matrix using index
  char_matrix <- pattern_mat[, index, drop = FALSE]
  rownames(char_matrix) <- taxa
  char_matrix
}
