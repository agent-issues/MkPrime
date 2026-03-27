# Batch-streaming output helpers for MkPrime MCMC
#
# These functions implement the flush-buffer design: scalar MCMC samples
# accumulate in a small in-memory buffer; when the buffer is full it is
# appended to a Tracer-compatible tab-separated log file.  A separate
# circular convergence window keeps recent samples available for PSRF/ESS
# checks without re-reading the file.
#
# None of these functions are exported.  See ReadMkLog() for the
# user-facing counterpart.


# Compute convergence window size from MCMC config.
# The window must be large enough for at least one full checkEvery interval.
# @keywords internal
.ComputeConvWindowSize <- function(mcmc) {
  checkEveryThin <- if (!is.null(mcmc$checkEvery) && mcmc$thin > 0L) {
    max(1L, as.integer(mcmc$checkEvery) %/% mcmc$thin)
  } else {
    100L
  }
  max(mcmc$bufferSize, 4L * checkEveryThin)
}


# Resolve log file paths for nRuns runs.
# Single run: path unchanged.  Multiple runs: "base.log" → "base_1.log", ...
# @keywords internal
.LogFilePaths <- function(logFile, nRuns) {
  if (nRuns == 1L) return(logFile)
  ext  <- tools::file_ext(logFile)
  base <- tools::file_path_sans_ext(logFile)
  if (nzchar(ext)) {
    paste0(base, "_", seq_len(nRuns), ".", ext)
  } else {
    paste0(base, "_", seq_len(nRuns))
  }
}


# Create log file(s) and write the tab-separated header line.
# Returns character vector of resolved paths (length nRuns).
# @keywords internal
.OpenLogFiles <- function(logFile, paramNames, nRuns) {
  paths  <- .LogFilePaths(logFile, nRuns)
  header <- paste(c("Sample", paramNames), collapse = "\t")
  for (p in paths) writeLines(header, p)
  paths
}


# Initialize streaming buffers for one run.
# Returns a named list to be merged into runs[[r]].
# @keywords internal
.InitStreamBuffers <- function(nParams, paramNames, bufferSize, convWindowSize) {
  list(
    flush_buf   = matrix(NA_real_, nrow = bufferSize, ncol = nParams,
                         dimnames = list(NULL, paramNames)),
    flush_iter  = integer(bufferSize),  # MCMC iteration for the Sample column
    flush_idx   = 0L,
    conv_window = matrix(NA_real_, nrow = convWindowSize, ncol = nParams,
                         dimnames = list(NULL, paramNames)),
    conv_head   = 0L,   # next write position in the circular window
    conv_filled = FALSE # TRUE once the window has wrapped around once
  )
}


# Flush nRows rows of the buffer to logFile (tab-separated append).
# @keywords internal
.FlushBuffer <- function(buffer, nRows, iterNums, logFile) {
  if (nRows == 0L) return(invisible(NULL))
  rows  <- buffer[seq_len(nRows), , drop = FALSE]
  lines <- vapply(seq_len(nRows), function(i) {
    paste(c(iterNums[i], rows[i, ]), collapse = "\t")
  }, character(1L))
  cat(paste(lines, collapse = "\n"), "\n", file = logFile, append = TRUE, sep = "")
  invisible(NULL)
}


# Add one sample row to the flush buffer and convergence window.
# Flushes the buffer to logFile when full.  Returns modified run list r.
# @keywords internal
.AddToStreamBuffer <- function(r, row, iterNum, logFile, bufferSize,
                                convWindowSize) {
  # --- flush buffer ---
  r$flush_idx <- r$flush_idx + 1L
  r$flush_buf[r$flush_idx, ] <- row
  r$flush_iter[r$flush_idx]  <- iterNum
  if (r$flush_idx == bufferSize) {
    .FlushBuffer(r$flush_buf, r$flush_idx, r$flush_iter, logFile)
    r$flush_idx <- 0L
  }

  # --- convergence window (circular buffer) ---
  r$conv_head <- (r$conv_head %% convWindowSize) + 1L
  r$conv_window[r$conv_head, ] <- row
  if (!r$conv_filled && r$conv_head == convWindowSize) r$conv_filled <- TRUE

  r
}


# Extract the currently valid rows from a circular convergence window.
# Returns NULL if fewer than minRows rows are available.
# @keywords internal
.ConvWindowRows <- function(r, minRows = 10L) {
  nRows <- if (r$conv_filled) nrow(r$conv_window) else r$conv_head
  if (nRows < minRows) return(NULL)
  if (r$conv_filled) r$conv_window else r$conv_window[seq_len(nRows), , drop = FALSE]
}


#' Read a MkPrime streaming log file
#'
#' Loads the tab-separated parameter log written by [RunMkPrime()] when
#' `logFile` is set in [MkPrimeMCMC()].  The returned matrix has the same
#' column layout as `MkPosterior$samples` and can be assigned directly:
#'
#' ```r
#' result$samples <- ReadMkLog(result$logFile)
#' summary(result)
#' ```
#'
#' @param logFile Character.  Path to one or more log files.  If a vector of
#'   paths is supplied (e.g. the `$logFile` field of an `MkPosterior` from a
#'   multi-run streaming analysis), the samples are combined row-wise.
#'
#' @return A named numeric matrix with one row per thinned sample.  Column
#'   names match those of `MkPosterior$samples`.  The `Sample` column (MCMC
#'   iteration numbers) is stored as the row-names attribute and is not
#'   included as a data column.
#'
#' @export
ReadMkLog <- function(logFile) {
  if (length(logFile) > 1L) {
    mats <- lapply(logFile, ReadMkLog)
    return(do.call(rbind, mats))
  }

  if (!file.exists(logFile)) {
    cli::cli_abort("Log file not found: {.file {logFile}}")
  }

  dat <- utils::read.table(logFile, header = TRUE, sep = "\t",
                           check.names = FALSE, comment.char = "")

  if (!"Sample" %in% names(dat)) {
    cli::cli_abort(
      "Log file {.file {logFile}} does not contain a {.col Sample} column."
    )
  }

  samples <- dat[["Sample"]]
  if (any(diff(samples) <= 0)) {
    cli::cli_warn(
      "The {.col Sample} column in {.file {logFile}} is not strictly increasing."
    )
  }

  mat <- as.matrix(dat[, setdiff(names(dat), "Sample"), drop = FALSE])
  rownames(mat) <- as.character(samples)
  mat
}
