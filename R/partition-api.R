# User-class partition API for RunMkPrime — Layer 1 plumbing.
#
# This file owns the new `partition` and `unlink` arguments. It is a sibling
# to (not a replacement of) the existing partition machinery in R/partition.R
# (which groups characters by (type, kObs) for the C++ pruning loop). The
# user-class layer wraps that grouping with an additional `classIdx` axis.
#
# Layer 1 currently implements:
#   - argument validation (partition shape, unlink token resolution)
#   - silent coercion of `unlink` when nothing is partitioned across
# and routes EVERY partition != NULL case to a deferred-implementation
# cli_abort. That keeps the §7a bit-identity guarantee for partition = NULL
# in force while later commits add the user-class state and C++ surface.

# Valid `unlink` tokens for Layer 1. `brlens` is reserved for Layer 2
# (T3 treatment) and rejected here with a "not yet implemented" message —
# it must validate cleanly through the prefix matcher (so users can pass
# uniform `unlink` vectors across treatments) but error before any code
# tries to honour it.
.kPartitionApiUnlinkTokens <- c("shape", "ratemultiplier", "brlens")
.kPartitionApiUnlinkLayer1  <- c("shape", "ratemultiplier")


# Validate and normalise `partition` and `unlink`.
#
# Returns a list with components:
#   partition  - NULL or an integer vector of length mkd$nChar with values
#                that form a contiguous range 1:nClasses.
#   unlink     - character vector of canonical tokens (subset of
#                .kPartitionApiUnlinkLayer1) with no duplicates. Always
#                character(0) when partition is NULL or nClasses == 1.
#   nClasses   - integer scalar (1L when partition is NULL).
#
# Errors hard on: wrong type, length mismatch, NAs, non-contiguous class
# ids, ambiguous unlink prefix, unknown unlink token (with `agrep` "did
# you mean" suggestion). Warns on resolvable partial-prefix matches.
.ValidatePartitionArgs <- function(partition, unlink, mkd) {

  # --- Unlink token resolution -------------------------------------------
  if (!is.null(unlink) && !is.character(unlink)) {
    cli::cli_abort("{.arg unlink} must be a character vector.")
  }
  unlink <- as.character(unlink)
  unlinkIn <- unlink

  resolved <- character(length(unlink))
  for (i in seq_along(unlink)) {
    raw <- unlink[i]
    if (!nzchar(raw)) {
      cli::cli_abort("{.arg unlink} contains an empty string at position {i}.")
    }
    tok <- tolower(raw)
    if (tok %in% .kPartitionApiUnlinkTokens) {
      resolved[i] <- tok
    } else {
      hits <- .kPartitionApiUnlinkTokens[startsWith(.kPartitionApiUnlinkTokens, tok)]
      if (length(hits) == 1L) {
        cli::cli_warn(
          "{.arg unlink} token {.val {unlinkIn[i]}} matched via prefix to {.val {hits}}."
        )
        resolved[i] <- hits
      } else if (length(hits) > 1L) {
        cli::cli_abort(c(
          "{.arg unlink} token {.val {unlinkIn[i]}} is an ambiguous prefix.",
          "i" = "Matches: {.val {hits}}."
        ))
      } else {
        validTokens <- .kPartitionApiUnlinkTokens
        sug <- agrep(tok, validTokens, value = TRUE,
                     max.distance = 0.4, ignore.case = TRUE)
        msgs <- c(
          "{.arg unlink} token {.val {unlinkIn[i]}} is not recognised.",
          "i" = "Valid tokens: {.val {validTokens}}."
        )
        if (length(sug) > 0L) {
          msgs <- c(msgs, "i" = "Did you mean: {.val {sug}}?")
        }
        cli::cli_abort(msgs)
      }
    }
  }
  if (anyDuplicated(resolved)) {
    cli::cli_warn("{.arg unlink} has duplicate tokens; deduplicating.")
    resolved <- unique(resolved)
  }
  unlink <- resolved

  # --- Partition shape validation ---------------------------------------
  nClasses <- 1L
  if (!is.null(partition)) {
    if (!is.numeric(partition)) {
      cli::cli_abort("{.arg partition} must be an integer vector or {.val NULL}.")
    }
    if (anyNA(partition)) {
      cli::cli_abort("{.arg partition} contains {.val NA}; supply an integer vector with no missing values.")
    }
    if (any(partition != floor(partition))) {
      cli::cli_abort("{.arg partition} must contain whole-number class IDs only.")
    }
    partition <- as.integer(partition)
    if (length(partition) != mkd$nChar) {
      cli::cli_abort(c(
        "{.arg partition} has length {length(partition)} but the data has \\
        {mkd$nChar} character{?s} (after invariant-character drop).",
        "i" = "Supply a length-{mkd$nChar} integer vector whose values form a \\
              contiguous range {.val 1}:{.val nClasses}."
      ))
    }
    if (any(partition < 1L)) {
      cli::cli_abort("{.arg partition} must contain values {.val >= 1}; got values < 1.")
    }
    nClasses <- max(partition)
    observed <- sort(unique(partition))
    expected <- seq_len(nClasses)
    if (!identical(observed, expected)) {
      missingCls <- setdiff(expected, observed)
      cli::cli_abort(c(
        "{.arg partition} skips class id{?s} {.val {missingCls}}.",
        "i" = "Class labels must form a contiguous range {.val 1}:{.val nClasses}."
      ))
    }
  }

  # --- Silent coercion when partition is trivial ------------------------
  if ((is.null(partition) || nClasses == 1L) && length(unlink) > 0L) {
    why <- if (is.null(partition)) {
      "no {.arg partition} supplied"
    } else {
      "only one user class in {.arg partition}"
    }
    cli::cli_alert_info(
      "{.arg unlink} ignored ({why}; nothing to unlink across)."
    )
    unlink <- character(0)
  }

  # Return:
  list(partition = partition, unlink = unlink, nClasses = nClasses)
}


# Layer-1 implementation gate. Called from RunMkPrime once validation has
# produced a (possibly non-trivial) partition spec. Until the user-class
# state, eta_neo reparameterisation, and sibling C++ likelihood land in
# follow-up commits, this aborts cleanly so the §7a bit-identity guarantee
# for partition = NULL stays under test.
.RequirePartitionImplemented <- function(spec) {
  if (is.null(spec$partition) && length(spec$unlink) == 0L) {
    # Return:
    return(invisible(NULL))
  }
  cli::cli_abort(c(
    "Partition-aware code path is not yet implemented at this commit.",
    "i" = "Layer 1 of feature/partition-api is in progress; pass \\
          {.code partition = NULL} for now."
  ))
}
