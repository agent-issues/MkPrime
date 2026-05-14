#' Relabel ecology-aware MCMC samples to the phi >= 1 convention
#'
#' @description
#' The ecology-aware MkPrime likelihood is invariant under the joint
#' transformation (phi -> 1/phi, theta_e -> 1 - theta_e, swap encouraged
#' and discouraged z-labels for ecology e), creating a discrete
#' identifiability issue: the posterior has equal mass at the two
#' reflections when the prior on log(phi) and theta_e is symmetric.
#'
#' This function applies the reflection per sample, per ecology, choosing
#' the phi >= 1 representative in each case. When a chain settles in the
#' phi < 1 mode, all samples are relabelled so that phi -> 1/phi,
#' theta_e -> 1 - theta_e, and encouraged/discouraged z-labels are
#' swapped. When a chain straddles both modes (some samples phi > 1, some
#' phi < 1), each sample is corrected independently.
#'
#' For \code{magnitudeMode = "global"}: if \code{log(phi) < 0} for a
#' sample, the entire reflection is applied: phi -> 1/phi, all theta_e ->
#' 1 - theta_e, and all z-matrix columns have codes 1 <-> 2 swapped.
#'
#' For \code{magnitudeMode = "per_ecology"}: each non-reference ecology e
#' is treated independently. If \code{log(phi_e) < 0}, then phi_e ->
#' 1/phi_e, theta_e -> 1 - theta_e, and column e of the z-matrix has
#' codes 1 <-> 2 swapped. The likelihood is invariant cell-by-cell.
#'
#' @param result An \code{MkPosterior} object returned from
#'   \code{\link{RunMkPrime}} with \code{ecologyAware = TRUE}.  Must
#'   contain \code{samples} (numeric matrix with at least a \code{phi} or
#'   \code{phi_1} column, one or more \code{theta_e} columns, and a
#'   \code{pi0} column) and \code{z_samples} (a list of
#'   \code{nChar x (kEco-1)} integer matrices, one per retained sample).
#'   If the run was executed in streaming mode, \code{samples} may be an
#'   empty matrix; in that case \code{result$logFile} must point to a
#'   valid log file, which will be read into memory automatically.
#' @param magnitudeMode Either \code{"global"} (single phi shared across
#'   ecologies) or \code{"per_ecology"} (one phi per ecology).  If
#'   \code{NULL} (default), inferred from column names: \code{"global"} if
#'   a column named \code{phi} exists, \code{"per_ecology"} otherwise.
#'
#' @return The \code{result} object with \code{samples}, \code{z_samples},
#'   and the model's \code{magnitudeMode} field updated in place. All other
#'   fields (\code{trees}, \code{logFile}, \code{model}, etc.) are
#'   unchanged. The returned object carries a new attribute
#'   \code{"relabelled"} set to \code{TRUE}.
#'
#' @examples
#' \dontrun{
#'   res <- RunMkPrime(mkd, model = MkPrimeModel(ecology = 220))
#'   resCanonical <- RelabelEcology(res)
#' }
#' @export
RelabelEcology <- function(result, magnitudeMode = NULL) {
  # Guard: not an MkPosterior
  if (!inherits(result, "MkPosterior")) {
    cli::cli_abort(
      "{.arg result} must be an {.cls MkPosterior} object returned by \\
       {.fn RunMkPrime}."
    )
  }

  # Check ecology-aware
  ecoAware <- isTRUE(result$model$ecologyAware)
  if (!ecoAware) {
    message("No ecology model detected; result returned unchanged.")
    return(result)
  }

  # Load samples from disk if in streaming mode
  samples <- result$samples
  if (nrow(samples) == 0L && !is.null(result$logFile)) {
    samples <- ReadMkLog(result$logFile)
  }
  if (nrow(samples) == 0L) {
    cli::cli_abort(
      "No posterior samples found. Run {.fn ReadMkLog} to load samples \\
       from disk first, or check {.code result$logFile}."
    )
  }

  nSamples  <- nrow(samples)
  colNms    <- colnames(samples)

  # Infer magnitudeMode from column names if not supplied
  if (is.null(magnitudeMode)) {
    magnitudeMode <- if ("phi" %in% colNms) "global" else "per_ecology"
  }
  magnitudeMode <- match.arg(magnitudeMode, c("global", "per_ecology"))

  # Locate phi, theta, and z columns
  if (identical(magnitudeMode, "global")) {
    phiCols  <- "phi"
    if (!phiCols %in% colNms) {
      cli::cli_abort(
        "Column {.val phi} not found. \\
         Did you mean {.arg magnitudeMode = \"per_ecology\"}?"
      )
    }
  } else {
    phiCols <- colNms[startsWith(colNms, "phi_")]
    if (length(phiCols) == 0L) {
      cli::cli_abort(
        "No columns matching {.val phi_*} found. \\
         Did you mean {.arg magnitudeMode = \"global\"}?"
      )
    }
  }

  thetaCols <- colNms[startsWith(colNms, "theta_")]
  if (length(thetaCols) == 0L) {
    cli::cli_abort(
      "No {.val theta_*} columns found. Is this an ecology-aware result?"
    )
  }
  nTheta <- length(thetaCols)

  # Validate z_samples
  zSamples <- result$z_samples
  if (is.null(zSamples) || length(zSamples) == 0L) {
    cli::cli_abort(
      "{.code result$z_samples} is missing or empty."
    )
  }
  if (length(zSamples) != nSamples) {
    cli::cli_abort(
      "{.code result$z_samples} has {length(zSamples)} entries but \\
       {.code samples} has {nSamples} rows; they must align."
    )
  }

  # -------------------------------------------------------------------------
  # Core relabelling
  # -------------------------------------------------------------------------
  #
  # For global mode: one flag per sample (flip if phi < 1).
  # For per_ecology mode: one flag per (sample x non-reference ecology).
  #
  # The j-th theta column (j = 1..nTheta, i.e. theta_j) and the j-th
  # z-matrix column (1-indexed) correspond to the j-th non-reference
  # ecology. In per_ecology mode, the phi for the j-th non-reference
  # ecology must be identified via the model's refEcology.

  if (identical(magnitudeMode, "global")) {
    # -----------------------------------------------------------------------
    # Global mode: single phi
    # -----------------------------------------------------------------------
    phiVec <- samples[, "phi"]
    flip   <- phiVec < 1.0

    if (any(flip)) {
      # Flip phi
      samples[flip, "phi"] <- 1.0 / phiVec[flip]

      # Flip all theta_e columns
      for (tc in thetaCols) {
        samples[flip, tc] <- 1.0 - samples[flip, tc]
      }

      # Swap z codes 1 <-> 2 in all columns for flipped samples
      flipIdx <- which(flip)
      for (i in flipIdx) {
        zm <- zSamples[[i]]
        # Swap in-place: 1->2, 2->1, 0 unchanged
        zSamples[[i]] <- .SwapZ12(zm)
      }
    }

  } else {
    # -----------------------------------------------------------------------
    # Per-ecology mode: one phi per ecology; reflect per non-reference ecology
    # -----------------------------------------------------------------------
    # Determine refEcology (0-indexed). Stored in result$data$refEcology.
    refEco <- result$data$refEcology  # 0-indexed integer
    if (is.null(refEco)) {
      cli::cli_abort(
        "{.code result$data$refEcology} is NULL; cannot determine \\
         reference ecology for per-ecology relabelling."
      )
    }
    kEco <- length(phiCols)  # = kEcology in R (number of phi columns)

    # Build mapping: for each j in 1..nTheta, what R-column-index of phi?
    # Non-reference ecologies in order (0-indexed C++):
    #   j-th (0-indexed) non-reference ecology = j if j < refEco, else j+1
    # In R 1-indexed phi columns (phi_1 .. phi_kEco), that is column e_j+1
    # where e_j = j if j < refEco else j+1 (0-indexed).
    phiColForTheta <- vapply(seq_len(nTheta), function(j) {
      eJ0 <- if ((j - 1L) < refEco) j - 1L else j   # 0-indexed ecology idx
      paste0("phi_", eJ0 + 1L)                        # R column name
    }, character(1L))

    for (j in seq_len(nTheta)) {
      pc  <- phiColForTheta[j]
      tc  <- thetaCols[j]
      phiVec <- samples[, pc]
      flip   <- phiVec < 1.0

      if (any(flip)) {
        samples[flip, pc] <- 1.0 / phiVec[flip]
        samples[flip, tc] <- 1.0 - samples[flip, tc]

        flipIdx <- which(flip)
        for (i in flipIdx) {
          zm <- zSamples[[i]]
          zm[, j] <- .SwapZ12(zm[, j, drop = FALSE])
          zSamples[[i]] <- zm
        }
      }
    }
  }

  # Write back
  result$samples   <- samples
  result$z_samples <- zSamples

  attr(result, "relabelled") <- TRUE
  # Return:
  result
}


# Swap z codes 1 <-> 2, leaving 0 unchanged.
# Accepts an integer vector or integer matrix; returns the same type.
#' @keywords internal
.SwapZ12 <- function(z) {
  # Map: 0->0, 1->2, 2->1. Vectorised via lookup.
  lookup <- c(0L, 2L, 1L)  # index 1=code0, 2=code1, 3=code2
  # z entries are in {0, 1, 2}; use z+1 as index into lookup
  result <- lookup[z + 1L]
  dim(result) <- dim(z)
  result
}
