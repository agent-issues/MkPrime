#!/usr/bin/env Rscript
# M-131: Post-hoc analysis of warmup stabilisation validation runs.
#
# Loads all 32 RDS files from Hamilton, replays the stabilisation detector
# with a grid of (windowSize, nStableRequired) configs, and compares
# against ground-truth stationarity (last 25% of trace as reference).
#
# Usage: Rscript m131-analyse-warmup.R [results_dir]

library(ggplot2)

results_dir <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(results_dir)) results_dir <- "inst/hamilton/results"

rds_files <- list.files(results_dir, pattern = "^m131_.*\.rds$",
                         full.names = TRUE)
if (length(rds_files) == 0L) stop("No RDS files found in: ", results_dir)

cat(sprintf("Loading %d result files...\n", length(rds_files)))
results <- lapply(rds_files, readRDS)

# --- Replay the stabilisation detector ---

replay_detector <- function(trace, windowSize, zThreshold = 1.5,
                             nStableRequired = 3L, batchSize = 500L) {
  # Returns the iteration at which the detector would have fired, or NA.
  n <- length(trace)
  nStable <- 0L
  for (i in seq_len(n)) {
    if (i < 2L * windowSize) {
      nStable <- 0L
      next
    }
    recent <- trace[(i - windowSize + 1L):i]
    prev   <- trace[(i - 2L * windowSize + 1L):(i - windowSize)]
    meanR <- mean(recent)
    meanP <- mean(prev)
    varR  <- var(recent)
    varP  <- var(prev)
    denom <- sqrt(varR / windowSize + varP / windowSize)
    if (denom < .Machine$double.eps) {
      nStable <- nStable + 1L
    } else {
      z <- (meanR - meanP) / denom
      if (abs(z) < zThreshold) {
        nStable <- nStable + 1L
      } else {
        nStable <- 0L
      }
    }
    if (nStable >= nStableRequired) {
      return(i * batchSize)
    }
  }
  NA_real_
}

# --- Ground truth: equilibrium reference from last 25% ---

ground_truth_mean <- function(trace) {
  n <- length(trace)
  start <- as.integer(0.75 * n) + 1L
  mean(trace[start:n])
}

# --- Grid search ---

windowSizes    <- c(3L, 5L, 7L, 10L, 12L, 15L, 20L)
nStableOptions <- c(2L, 3L, 4L, 5L)

cat("Replaying detector with grid:\n")
cat(sprintf("  windowSize: %s\n", paste(windowSizes, collapse = ", ")))
cat(sprintf("  nStableRequired: %s\n", paste(nStableOptions, collapse = ", ")))

rows <- list()
for (res in results) {
  trace <- res$warmup_trace
  if (is.null(trace) || length(trace) == 0L) next
  eq_mean <- ground_truth_mean(trace)

  for (ws in windowSizes) {
    for (ns in nStableOptions) {
      iter_detected <- replay_detector(trace, ws, nStableRequired = ns)

      # Classify: too-early, on-time, too-late, or never
      if (is.na(iter_detected)) {
        classification <- "never"
        post_mean <- NA_real_
      } else {
        # Mean logP from detection point to end of trace
        snap_idx <- iter_detected / 500L
        if (snap_idx < length(trace)) {
          post_mean <- mean(trace[(snap_idx + 1L):length(trace)])
        } else {
          post_mean <- trace[length(trace)]
        }
        deficit <- eq_mean - post_mean
        if (deficit > 2) {
          classification <- "too_early"
        } else {
          classification <- "on_time"
        }
      }

      # Also compute the adaptive formula result
      adaptive_ws <- max(5L, min(20L, res$nEdge %/% 10L))

      rows[[length(rows) + 1L]] <- data.frame(
        dataset        = res$dataset,
        seed           = res$seed,
        nTip           = res$nTip,
        nEdge          = res$nEdge,
        windowSize     = ws,
        nStableReq     = ns,
        adaptive_ws    = adaptive_ws,
        iter_detected  = iter_detected,
        eq_mean        = eq_mean,
        post_mean      = post_mean,
        classification = classification,
        stringsAsFactors = FALSE
      )
    }
  }
}

grid_results <- do.call(rbind, rows)

# --- Summary table ---

cat("\n=== Classification summary (windowSize x nStableRequired) ===\n\n")
for (ns in nStableOptions) {
  cat(sprintf("nStableRequired = %d:\n", ns))
  sub <- grid_results[grid_results$nStableReq == ns, ]
  tab <- table(sub$windowSize, sub$classification)
  print(tab)
  cat("\n")
}

# --- Adaptive formula comparison ---

cat("\n=== Adaptive formula: windowSize = max(5, min(20, nEdge %/% 10)) ===\n\n")
adaptive_rows <- list()
for (res in results) {
  trace <- res$warmup_trace
  if (is.null(trace) || length(trace) == 0L) next
  eq_mean <- ground_truth_mean(trace)
  ws <- max(5L, min(20L, res$nEdge %/% 10L))

  for (ns in nStableOptions) {
    iter_detected <- replay_detector(trace, ws, nStableRequired = ns)
    if (is.na(iter_detected)) {
      classification <- "never"
    } else {
      snap_idx <- iter_detected / 500L
      post_mean <- if (snap_idx < length(trace)) {
        mean(trace[(snap_idx + 1L):length(trace)])
      } else {
        trace[length(trace)]
      }
      deficit <- eq_mean - post_mean
      classification <- if (deficit > 2) "too_early" else "on_time"
    }
    adaptive_rows[[length(adaptive_rows) + 1L]] <- data.frame(
      dataset        = res$dataset,
      seed           = res$seed,
      nTip           = res$nTip,
      nEdge          = res$nEdge,
      windowSize     = ws,
      nStableReq     = ns,
      iter_detected  = iter_detected,
      classification = classification,
      stringsAsFactors = FALSE
    )
  }
}
adaptive_results <- do.call(rbind, adaptive_rows)

for (ns in nStableOptions) {
  cat(sprintf("nStableRequired = %d:\n", ns))
  sub <- adaptive_results[adaptive_results$nStableReq == ns, ]
  tab <- table(sub$dataset, sub$classification)
  print(tab)
  cat(sprintf("  Median detection iter: %s\n",
              median(sub$iter_detected, na.rm = TRUE)))
  cat("\n")
}

# --- Fixed vs adaptive comparison ---

cat("\n=== Fixed (ws=10) vs Adaptive: median detection iteration ===\n\n")
fixed10 <- grid_results[grid_results$windowSize == 10L &
                         grid_results$nStableReq == 3L, ]
adapt3 <- adaptive_results[adaptive_results$nStableReq == 3L, ]

comparison <- merge(
  fixed10[, c("dataset", "seed", "iter_detected")],
  adapt3[, c("dataset", "seed", "iter_detected")],
  by = c("dataset", "seed"),
  suffixes = c("_fixed10", "_adaptive")
)
comparison$savings <- comparison$iter_detected_fixed10 -
  comparison$iter_detected_adaptive
cat("Per-dataset comparison (nStableRequired = 3):\n")
print(comparison[order(comparison$dataset, comparison$seed), ])
cat(sprintf("\nMedian savings: %.0f iterations\n", median(comparison$savings, na.rm = TRUE)))

# --- Save full results ---

out_file <- file.path(results_dir, "m131_analysis.rds")
saveRDS(list(grid = grid_results, adaptive = adaptive_results,
             comparison = comparison), out_file)
cat(sprintf("\nFull results saved to: %s\n", out_file))

# --- Diagnostic plots ---

pdf(file.path(results_dir, "m131_traces.pdf"), width = 10, height = 8)
for (res in results) {
  trace <- res$warmup_trace
  if (is.null(trace) || length(trace) == 0L) next
  iters <- seq_along(trace) * 500L
  plot(iters, trace, type = "l",
       main = sprintf("%s (seed %d, nTip=%d)", res$dataset, res$seed, res$nTip),
       xlab = "Iteration", ylab = "log P")

  # Mark fixed ws=10, ns=3 detection point
  det_fixed <- replay_detector(trace, 10L, nStableRequired = 3L)
  if (!is.na(det_fixed)) abline(v = det_fixed, col = "red", lty = 2)

  # Mark adaptive detection point
  ws_adapt <- max(5L, min(20L, res$nEdge %/% 10L))
  det_adapt <- replay_detector(trace, ws_adapt, nStableRequired = 3L)
  if (!is.na(det_adapt)) abline(v = det_adapt, col = "blue", lty = 2)

  # Equilibrium reference
  eq <- ground_truth_mean(trace)
  abline(h = eq, col = "green3", lty = 3)

  legend("bottomright",
         c(sprintf("Fixed ws=10: %s", ifelse(is.na(det_fixed), "never", det_fixed)),
           sprintf("Adaptive ws=%d: %s", ws_adapt,
                   ifelse(is.na(det_adapt), "never", det_adapt)),
           "Equilibrium mean"),
         col = c("red", "blue", "green3"), lty = c(2, 2, 3), cex = 0.8)
}
dev.off()
cat(sprintf("Trace plots saved to: %s\n",
            file.path(results_dir, "m131_traces.pdf")))
