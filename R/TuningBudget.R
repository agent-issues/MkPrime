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
#' @param candFrozen,bestFrozen Logical; `TRUE` where that window never changed
#'   topology, so its rate reflects the scalar gate alone; `FALSE` where it
#'   did; `NA` where no trees were scored. A window that moved the topology
#'   outranks one that froze it whatever their rates.
#' @param topologyCut Logical; `TRUE` where the candidate gives topology moves
#'   less weight than the incumbent.
#' @param everMoved Logical; `TRUE` once this run has adopted a schedule seen
#'   to move the topology.
#' @return `TRUE` where the candidate should be adopted.
#' @details A frozen window's rate prices only the scalar moves, so among
#'   frozen windows the cheapest schedule is the one that spends least on
#'   topology. A frozen candidate is therefore never adopted if it cuts
#'   topology weight, nor once a moving schedule has been seen, whose evidence
#'   one quiet window cannot overturn.
#' @keywords internal
.BeatsIncumbent <- function(candRate, candEss, bestRate, bestEss,
                            z = .kTuningGateZ,
                            candFrozen = NA, bestFrozen = NA,
                            topologyCut = FALSE, everMoved = FALSE) {
  if (!isTRUE(is.finite(candRate)) || candRate <= 0) {
    return(FALSE)
  }
  if (isTRUE(candFrozen) && (isTRUE(topologyCut) || isTRUE(everMoved))) {
    return(FALSE)
  }
  if (!isTRUE(is.finite(bestRate)) || bestRate <= 0) {
    return(TRUE)
  }
  if (!is.na(candFrozen) && !is.na(bestFrozen) && candFrozen != bestFrozen) {
    return(bestFrozen)
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
#' @param frozen Logical; `TRUE` where the best window never changed topology.
#'   Its rate prices the scalar moves alone, so it cannot vote to freeze.
#' @return List with the updated `streak` and a logical `freeze`.
#' @keywords internal
.TuningPayback <- function(streak, spentSec, bestRate, targetEss, nRuns = 1L,
                           frozen = FALSE) {
  vote <- !isTRUE(frozen) && isTRUE(is.finite(bestRate)) && bestRate > 0 &&
    !is.null(targetEss) && isTRUE(is.finite(targetEss)) &&
    spentSec >= (targetEss / max(nRuns, 1L)) / bestRate
  streak <- if (vote) streak + 1L else 0L
  # Return:
  list(streak = streak, freeze = streak >= .kTuningFreezeStreak)
}

# Samples a tuning window collects before it is scored. Admitting a window at
# its first batch scored each candidate on ~29 draws (#78).
.kTuningWindowSamples <- 100L

#' Samples each tuning window collects before it is scored
#'
#' A window spans whole batches, so it is sized in batches: enough to reach
#' `target` samples, but never so many that one round of `nCandidates` windows
#' outruns `budget`.
#'
#' @param budget Integer iterations the tuning phase may spend.
#' @param thin Integer iterations per stored sample.
#' @param batch Integer iterations per tuning batch.
#' @param nCandidates Integer windows per round.
#' @param target Integer samples a window should hold.
#' @return Integer samples at which a window is scored.
#' @keywords internal
.TuningWindowSamples <- function(budget, thin, batch, nCandidates = 4L,
                                 target = .kTuningWindowSamples) {
  perBatch <- max(batch %/% thin, 1L)
  nBatch <- min(ceiling(target / perBatch),
                max(budget %/% (nCandidates * batch), 1L))
  # Return:
  as.integer(max(nBatch * perBatch, 10L))
}

#' Which window does the tuning round run next?
#'
#' The incumbent's window (index 0) is the yardstick every candidate is gated
#' against. If it yields no rate, the first usable candidate would be adopted
#' on one window with no gate at all, so the incumbent is measured again while
#' budget remains, and the round is closed on the incumbent once it does not.
#' A round with no candidates has nothing to gate, so it closes at once.
#'
#' @param candIdx Integer index of the window just scored; 0 is the incumbent.
#' @param rate Numeric min-ESS/s that window gave.
#' @param nCandidates Integer candidate windows in the round.
#' @param budgetLeft Logical; `TRUE` while tuning budget remains.
#' @return Integer index of the next window; above `nCandidates` to end the
#'   round.
#' @keywords internal
.NextTuningWindow <- function(candIdx, rate, nCandidates, budgetLeft) {
  if (candIdx == 0L && nCandidates > 0L && !isTRUE(rate > 0)) {
    return(if (budgetLeft) 0L else nCandidates + 1L)
  }
  candIdx + 1L
}
