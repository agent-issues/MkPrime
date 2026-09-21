# When a candidate is better, and when tuning has stopped paying for itself.
#
# Both rules guard the same failure: a single min-ESS/s window is a noisy
# statistic, so a bare argmax over candidates climbs its own noise, and a fixed
# round budget spends on tuning whatever the target costs to sample.

# Standard errors a candidate must clear. 1 makes adoption roughly a one-sided
# 84% call -- enough to stop noise-climbing without stalling real improvements.
.kTuningGateZ <- 1

# Consecutive rounds that must vote to freeze. A single window can price the
# target anywhere, and an early freeze cannot be recovered from.
.kTuningFreezeStreak <- 2L

#' Does a candidate beat the incumbent by more than measurement noise?
#'
#' The relative standard error of an effective-sample-size estimate is about
#' `sqrt(2 / ESS)` (Geyer 1992); the wall clock contributes negligibly beside
#' it, so the rate inherits the ESS's error.
#'
#' @param candRate,bestRate Numeric min-ESS/s for the candidate and incumbent.
#' @param candEss,bestEss Numeric effective sample sizes the rates came from.
#' @param z Numeric giving the standard errors a candidate must clear.
#' @return `TRUE` where the candidate should be adopted.
#' @keywords internal
.BeatsIncumbent <- function(candRate, candEss, bestRate, bestEss,
                            z = .kTuningGateZ) {
  if (!isTRUE(is.finite(candRate)) || candRate <= 0) {
    return(FALSE)
  }
  if (!isTRUE(is.finite(bestRate)) || bestRate <= 0) {
    return(TRUE)
  }
  rse <- sqrt(2 / max(candEss, 1) + 2 / max(bestEss, 1))
  candRate > bestRate * (1 + z * rse)
}

#' Has tuning stopped paying for itself?
#'
#' Tuning is worth a second only while it saves more than a second of
#' sampling, so the spend so far is compared against the sampling it precedes:
#' the seconds this run needs to collect its share of `targetEss` under the
#' kernel as tuned right now. A run with no ESS target, or no rate to price it
#' with, keeps tuning and is bounded by `tuningRounds` and `tuningBudget` as
#' before.
#'
#' @param streak Integer count of consecutive rounds that have voted to freeze.
#' @param spentSec Numeric seconds spent on warmup and tuning so far.
#' @param bestRate Numeric best min-ESS/s measured this round.
#' @param targetEss Numeric ESS the stopping rule asks of the pooled runs, or
#'   `NULL` where none is set.
#' @param nRuns Integer number of runs sharing that target.
#' @return List with the updated `streak` and a logical `freeze`.
#' @keywords internal
.TuningPayback <- function(streak, spentSec, bestRate, targetEss, nRuns = 1L) {
  vote <- isTRUE(is.finite(bestRate)) && bestRate > 0 &&
    !is.null(targetEss) && isTRUE(is.finite(targetEss)) &&
    spentSec >= (targetEss / max(nRuns, 1L)) / bestRate
  streak <- if (vote) streak + 1L else 0L
  # Return:
  list(streak = streak, freeze = streak >= .kTuningFreezeStreak)
}
