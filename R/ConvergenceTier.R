# Which sampled columns the diagnostics judge.
#
# The stopping rule, the auto-thin ACT estimate and the tuning bandit each
# carried their own regexp, and the three disagreed. A column is `"gate"`
# unless it is deliberately exempted here, so a new parameter family is judged
# until someone says otherwise.

# Not parameters: `swap_cold` is a swap indicator and `topo_hash` a topology
# hash, so neither has an autocorrelation time to read.
.kBookkeepingCols <- c("iteration", "chain", "swap_cold", "topo_hash")

# Sampled and reported, but not judged. `log_likelihood` is dropped because
# `log_posterior` already covers it; `kPrime_i` are discrete per-character
# nuisance parameters (M-098).
.kNuisanceCols <- "log_likelihood"

#' Convergence tier of each sampled column
#'
#' - `"gate"`: judged by the stopping rule, and minimized over by the
#'   auto-thin ACT estimate and the tuning bandit.
#' - `"nuisance"`: reported in the progress table but not judged.
#' - `"fixed"`: logged, but no move updates it, so it is constant by design.
#' - `"branch"`: the per-edge lengths `br_i`, too numerous to report.
#' - `"bookkeeping"`: not parameters.
#'
#' A gate column whose ESS or R-hat is `NA` makes the stopping rule
#' unassessable; only a `"fixed"` column may be constant without blocking it.
#'
#' @param colNames Character vector naming the sampled columns.
#' @param fixed Character vector naming the columns no move updates, as
#'   given by [.FixedCols()].
#' @return Character vector of tiers, named by `colNames`.
#' @keywords internal
.ConvergenceTier <- function(colNames, fixed = NULL) {
  colNames <- as.character(colNames)
  tier <- rep("gate", length(colNames))
  tier[colNames %in% fixed] <- "fixed"
  tier[grepl("^kPrime_[0-9]+$", colNames)] <- "nuisance"
  tier[colNames %in% .kNuisanceCols] <- "nuisance"
  tier[grepl("^br_[0-9]+$", colNames)] <- "branch"
  tier[colNames %in% .kBookkeepingCols] <- "bookkeeping"
  names(tier) <- colNames
  tier
}

#' Columns the convergence criteria are judged on
#'
#' @inheritParams .ConvergenceTier
#' @return Integer column indices; `.GateCols()` gives the judged columns and
#' `.ReportCols()` those the progress table shows.
#' @keywords internal
.GateCols <- function(colNames, fixed = NULL) {
  unname(which(.ConvergenceTier(colNames, fixed) == "gate"))
}

#' @rdname dot-GateCols
#' @keywords internal
.ReportCols <- function(colNames) {
  unname(which(.ConvergenceTier(colNames) %in% c("gate", "nuisance")))
}

# The k' hyperparameters are logged whatever the data, but only
# transformational characters give them a move.
.kHyperMoveTargets <- list(p = "p",
                           kprime_alpha = c("kprime_s", "kprime_r"),
                           kprime_beta = c("kprime_s", "kprime_r"))

#' Logged columns that no move updates
#'
#' @param colNames Character vector naming the sampled columns.
#' @param moves List of moves, each with a `target` naming what it updates.
#' @return Character vector naming the columns in `colNames` that stay at
#' their initial value by design.
#' @keywords internal
.FixedCols <- function(colNames, moves) {
  targets <- unlist(lapply(moves, `[[`, "target"))
  hyper <- intersect(names(.kHyperMoveTargets), colNames)
  hyper[!vapply(.kHyperMoveTargets[hyper],
                function(x) any(x %in% targets), logical(1))]
}
