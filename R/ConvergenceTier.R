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
#' - `"branch"`: the per-edge lengths `br_i`, too numerous to report.
#' - `"bookkeeping"`: not parameters.
#'
#' @param colNames Character vector naming the sampled columns.
#' @return Character vector of tiers, named by `colNames`.
#' @keywords internal
.ConvergenceTier <- function(colNames) {
  colNames <- as.character(colNames)
  tier <- rep("gate", length(colNames))
  tier[grepl("^kPrime_[0-9]+$", colNames)] <- "nuisance"
  tier[colNames %in% .kNuisanceCols] <- "nuisance"
  tier[grepl("^br_[0-9]+$", colNames)] <- "branch"
  tier[colNames %in% .kBookkeepingCols] <- "bookkeeping"
  names(tier) <- colNames
  tier
}

#' Columns the convergence criteria are judged on
#'
#' @param colNames Character vector naming the sampled columns.
#' @return Integer column indices; `.GateCols()` gives the judged columns and
#' `.ReportCols()` those the progress table shows.
#' @keywords internal
.GateCols <- function(colNames) {
  unname(which(.ConvergenceTier(colNames) == "gate"))
}

#' @rdname dot-GateCols
#' @keywords internal
.ReportCols <- function(colNames) {
  unname(which(.ConvergenceTier(colNames) %in% c("gate", "nuisance")))
}
