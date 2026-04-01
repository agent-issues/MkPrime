# Batch-streaming output helpers for MkPrime MCMC
#
# These functions implement the flush-buffer design: scalar MCMC samples
# accumulate in a small in-memory buffer; when the buffer is full it is
# appended to a Tracer-compatible tab-separated log file.  A separate
# circular convergence window keeps recent samples available for R-hat/ESS
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
    flushed     = FALSE,  # set TRUE when a flush fires; reset after checkpoint
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
    r$flushed   <- TRUE
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


# Truncate a log file to exactly nDataRows data lines (plus header).
# Used on resume to discard samples written after the last checkpoint.
# @keywords internal
.TruncateLogToN <- function(logFile, nDataRows) {
  lines <- readLines(logFile, warn = FALSE)
  # Separate header/comment lines from data rows.
  # Move-weight log entries start with "#"; the header starts with "Sample".
  isData <- !startsWith(lines, "#") & !startsWith(lines, "Sample\t")
  dataIdx <- which(isData)
  nData <- length(dataIdx)

  if (nData > nDataRows) {
    # Keep header + comments + first nDataRows data lines
    keepIdx <- c(which(!isData), dataIdx[seq_len(nDataRows)])
    keepIdx <- sort(keepIdx)
    nDropped <- nData - nDataRows
    cli::cli_alert_info(
      "Rewinding {.file {logFile}}: discarding {nDropped} post-checkpoint sample{?s}."
    )
    writeLines(lines[keepIdx], logFile)
  } else if (nData < nDataRows) {
    cli::cli_warn(c(
      "Log file {.file {logFile}} has fewer rows than expected.",
      "i" = "Found {nData} data row{?s}, checkpoint recorded {nDataRows}."
    ))
  }
}


#' Recover partial results from an interrupted run
#'
#' When [RunMkPrime()] is interrupted (e.g. by pressing Escape or Ctrl-C),
#' samples already flushed to disk are preserved in a log file.
#' Call `MkPrimeRecover()` to load those samples into an `MkPosterior`
#' object.
#'
#' With no arguments, `MkPrimeRecover()` looks for temporary log files
#' from the most recent interrupted run in the current session.  Pass
#' `logFile` to recover from named log files (e.g. after restarting R).
#'
#' @param logFile Character.  Base path of the log file(s) written by
#'   [RunMkPrime()].  For multi-run analyses, pass the base name
#'   (e.g. `"hyoliths.log"`); the per-run files `hyoliths_1.log`,
#'   `hyoliths_2.log`, etc. are discovered automatically.  A single
#'   file path is also accepted.
#' @param checkpointFile Character.  Path to the `.ckp` checkpoint file.
#'   If `NULL` (default), derived from `logFile` by replacing the
#'   extension with `.ckp`.  The checkpoint supplies `model` and `mcmc`
#'   metadata; recovery still works without it, but the returned object
#'   will have less metadata.
#'
#' @return An [MkPosterior] object containing the partial samples, or
#'   `NULL` (with a message) if no data is available.
#'
#' @seealso [RunMkPrime()], [ReadMkLog()]
#' @export
MkPrimeRecover <- function(logFile = NULL, checkpointFile = NULL) {

  # --- Path 1: recover from named log file(s) on disk ---
  if (!is.null(logFile)) {
    logPaths <- .DiscoverLogFiles(logFile)
    if (is.null(logPaths)) {
      cli::cli_alert_danger("No log files found for {.file {logFile}}.")
      return(invisible(NULL))
    }

    samples <- tryCatch(
      ReadMkLog(logPaths),
      error = function(e) {
        cli::cli_alert_danger(
          "Failed to read log file{?s}: {conditionMessage(e)}"
        )
        return(NULL)
      }
    )

    if (is.null(samples) || nrow(samples) == 0L) {
      nLP <- length(logPaths)
      cli::cli_alert_warning(
        "Log {cli::qty(nLP)}file{?s} contain{?s/} no samples \\
         (run may have been interrupted before any were flushed)."
      )
      return(invisible(NULL))
    }

    # Try to load metadata from checkpoint
    ckpFile <- checkpointFile %||%
      sub("(_\\d+)?\\.[^.]+$", ".ckp", logPaths[1])
    ckp <- NULL
    if (file.exists(ckpFile)) {
      ckp <- tryCatch(readRDS(ckpFile), error = function(e) NULL)
    }

    model <- ckp$model
    mcmc  <- ckp$mcmc

    result <- MkPosterior(
      samples    = samples,
      trees      = list(),
      acceptance = numeric(0),
      model      = model,
      data       = NULL,
      mcmc       = mcmc,
      warmup     = mcmc$warmup %||% 0L,
      tuning     = NULL
    )
    result$partial     <- TRUE
    result$nSamples    <- nrow(samples)
    result$stop_reason <- "recovered"
    result$logFile     <- logPaths

    if (length(logPaths) > 1L) {
      result$nRuns   <- length(logPaths)
      result$per_run <- lapply(logPaths, function(f) {
        s <- tryCatch(ReadMkLog(f), error = function(e) {
          matrix(numeric(0), nrow = 0, ncol = ncol(samples))
        })
        list(samples = s, trees = list(), acceptance = numeric(0),
             saved_idx = nrow(s))
      })
    }

    cli::cli_alert_success(
      "Recovered {nrow(samples)} sample{?s} from \\
       {length(logPaths)} log file{?s}."
    )
    return(result)
  }

  # --- Path 2: recover from session-local temp logs (existing behaviour) ---
  rec <- .mkp_env$recovery
  if (is.null(rec)) {
    cli::cli_alert_info(
      "No interrupted run to recover.
       {.emph Tip: pass {.arg logFile} to recover from a named log file.}"
    )
    return(invisible(NULL))
  }

  # Check that the temp log files still exist
  nFiles <- length(rec$logFiles)
  missing <- !file.exists(rec$logFiles)
  if (all(missing)) {
    cli::cli_alert_danger(
      "Temporary log {cli::qty(nFiles)}file{?s} no longer exist{?s/}. \\
       Cannot recover."
    )
    .mkp_env$recovery <- NULL
    return(invisible(NULL))
  }

  if (any(missing)) {
    cli::cli_alert_warning(
      "Some log files are missing; recovering from {sum(!missing)} of \\
       {nFiles} run{?s}."
    )
    rec$logFiles <- rec$logFiles[!missing]
  }

  # Read samples
  samples <- tryCatch(
    ReadMkLog(rec$logFiles),
    error = function(e) {
      cli::cli_alert_danger(
        "Failed to read log file{?s}: {conditionMessage(e)}"
      )
      return(NULL)
    }
  )

  if (is.null(samples) || nrow(samples) == 0L) {
    cli::cli_alert_warning(
      "Log {cli::qty(length(rec$logFiles))}file{?s} contain{?s/} no samples \\
       (run may have been interrupted before any were flushed)."
    )
    .CleanupTempLogs(rec$logFiles)
    return(invisible(NULL))
  }

  # Build a minimal MkPosterior
  result <- MkPosterior(
    samples    = samples,
    trees      = list(),
    acceptance = numeric(0),
    model      = rec$model,
    data       = rec$data,
    mcmc       = rec$mcmc,
    warmup     = rec$mcmc$warmup,
    tuning     = NULL
  )
  result$partial    <- TRUE
  result$nSamples   <- nrow(samples)
  result$stop_reason <- "interrupted"

  # Clean up
  .CleanupTempLogs(rec$logFiles)

  cli::cli_alert_success(
    "Recovered {nrow(samples)} sample{?s} from interrupted run."
  )
  result
}


#' Discover per-run log files from a base log path
#'
#' Tries the `_N.log` multi-run naming convention first, then falls back
#' to the exact path.  Returns `NULL` if no files are found.
#' @keywords internal
.DiscoverLogFiles <- function(logFile) {
  # Try multi-run pattern: base_1.ext, base_2.ext, ...
  ext  <- tools::file_ext(logFile)
  base <- tools::file_path_sans_ext(logFile)
  first <- if (nzchar(ext)) paste0(base, "_1.", ext) else paste0(base, "_1")

  if (file.exists(first)) {
    paths <- first
    n <- 2L
    repeat {
      nxt <- if (nzchar(ext)) {
        paste0(base, "_", n, ".", ext)
      } else {
        paste0(base, "_", n)
      }
      if (!file.exists(nxt)) break
      paths <- c(paths, nxt)
      n <- n + 1L
    }
    return(paths)
  }

  # Fall back to exact file
 if (file.exists(logFile)) return(logFile)

  NULL
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
                           check.names = FALSE, comment.char = "#")

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
