#' Control MkPrime's console output
#'
#' MkPrime narrates a run as it goes: phase transitions, adapted tuning
#' settings, a progress bar, and notes about files it has written or rewound.
#' `verbosity` selects how much of that narration reaches the console.
#'
#' @details
#' The option `MkPrime.verbosity` sets the level globally; functions that take
#' a `verbosity` argument override it for the duration of that call.
#'
#' Levels are cumulative:
#'
#'  - `0`: silent. Nothing is printed. Warnings and errors are unaffected --
#'    they are conditions, not output, and still propagate.
#'  - `1`: progress bar and the milestones a user running interactively wants
#'    to see (chain stabilisation, frozen move weights, files written).
#'  - `2`: tuning and adaptation diagnostics, e.g. per-move weight decay.
#'
#' The package's own test suite sets `0`, so a test that asserts on console
#' output must raise the level for itself.
#'
#' @return
#' `MkPrimeVerbosity()` returns the current level as an integer.
#'
#' @examples
#' MkPrimeVerbosity()
#' @seealso [RunMkPrime()], which takes a `verbosity` argument.
#' @export
MkPrimeVerbosity <- function() {
  v <- getOption("MkPrime.verbosity", 1L)
  if (!is.numeric(v) || length(v) != 1L || is.na(v)) {
    cli::cli_warn(c(
      "Option {.code MkPrime.verbosity} must be a single number.",
      "i" = "Falling back to {.val 1}."
    ))
    # Return:
    1L
  } else {
    # Return:
    as.integer(v)
  }
}

#' Validate a user-supplied verbosity level
#'
#' @param verbosity Integer giving the requested level.
#' @return `verbosity` as an integer.
#' @keywords internal
.CheckVerbosity <- function(verbosity) {
  if (!is.numeric(verbosity) || length(verbosity) != 1L || is.na(verbosity)) {
    cli::cli_abort(c(
      "{.arg verbosity} must be a single number.",
      "i" = "See {.fn MkPrimeVerbosity} for the available levels."
    ))
  }
  # Return:
  as.integer(verbosity)
}

#' Is MkPrime talkative enough to emit this message?
#'
#' @param level Integer giving the verbosity level at which a message
#'   becomes visible.
#' @return Logical; `TRUE` when the message should be printed.
#' @keywords internal
.Loud <- function(level = 1L) {
  MkPrimeVerbosity() >= level
}

#' Verbosity-gated wrappers around cli output
#'
#' Each forwards to its `cli` counterpart when `.Loud()` permits, and is a
#' no-op otherwise. `.envir` defaults to the calling frame so glue
#' interpolation resolves exactly as it would in a direct `cli` call.
#'
#' @param text,message,... Passed to the corresponding `cli` function.
#' @param level Integer giving the verbosity level at which output appears.
#' @param .envir Environment in which to evaluate glue substitutions.
#' @return `NULL`, invisibly.
#' @keywords internal
#' @name mkp-verbosity-wrappers
NULL

#' @rdname mkp-verbosity-wrappers
.AlertInfo <- function(text, ..., level = 1L, .envir = parent.frame()) {
  if (.Loud(level)) cli::cli_alert_info(text, ..., .envir = .envir)
  invisible(NULL)
}

#' @rdname mkp-verbosity-wrappers
.AlertSuccess <- function(text, ..., level = 1L, .envir = parent.frame()) {
  if (.Loud(level)) cli::cli_alert_success(text, ..., .envir = .envir)
  invisible(NULL)
}

#' @rdname mkp-verbosity-wrappers
.AlertWarning <- function(text, ..., level = 1L, .envir = parent.frame()) {
  if (.Loud(level)) cli::cli_alert_warning(text, ..., .envir = .envir)
  invisible(NULL)
}

#' @rdname mkp-verbosity-wrappers
.AlertDanger <- function(text, ..., level = 1L, .envir = parent.frame()) {
  if (.Loud(level)) cli::cli_alert_danger(text, ..., .envir = .envir)
  invisible(NULL)
}

#' @rdname mkp-verbosity-wrappers
.Inform <- function(message, ..., level = 1L, .envir = parent.frame()) {
  if (.Loud(level)) cli::cli_inform(message, ..., .envir = .envir)
  invisible(NULL)
}

#' @rdname mkp-verbosity-wrappers
.Text <- function(text, ..., level = 1L, .envir = parent.frame()) {
  if (.Loud(level)) cli::cli_text(text, ..., .envir = .envir)
  invisible(NULL)
}

#' @rdname mkp-verbosity-wrappers
.ProgressBar <- function(..., level = 1L, .envir = parent.frame()) {
  if (.Loud(level)) cli::cli_progress_bar(..., .envir = .envir)
  invisible(NULL)
}

#' @rdname mkp-verbosity-wrappers
.ProgressUpdate <- function(..., level = 1L, .envir = parent.frame()) {
  if (.Loud(level)) cli::cli_progress_update(..., .envir = .envir)
  invisible(NULL)
}

#' @rdname mkp-verbosity-wrappers
.ProgressDone <- function(..., level = 1L, .envir = parent.frame()) {
  if (.Loud(level)) cli::cli_progress_done(..., .envir = .envir)
  invisible(NULL)
}
