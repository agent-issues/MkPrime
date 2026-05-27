#!/usr/bin/env Rscript
# Build a single 4-panel summary plot from the sbc-warmup-trace_*.rds files
# saved during the SBC-WARMUP-001 discriminating check (2026-05-27).
# Panels show the key discriminating runs:
#   (A) sim 1 baseline (tl_true=110.78, data uninformative)
#   (B) sim 1 with N_WARM=20000 — TL move crushed (single unique value)
#   (C) sim 15 baseline (tl_true=4.13, data informative)
#   (D) sim 15 long chain + true topology — still stuck at TL~1.5

ROOT <- "dev/red-team/heavy-tests"
runs <- list(
  A = list(rds = "sbc-warmup-trace_warm2000.rds",
           title = "Sim 1 baseline (tl_true=110.78)"),
  B = list(rds = "sbc-warmup-trace_warm20000_hold.rds",
           title = "Sim 1, N_WARM=20000 — move crushed"),
  C = list(rds = "sbc-warmup-trace_warm2000_edge0.1.rds",
           title = "Sim 15 baseline (tl_true=4.13)"),
  D = list(rds = "sbc-warmup-trace_sim15_warm2000_edge0.1_iter60000_truetopo.rds",
           title = "Sim 15, true topology, 60k iters")
)
data <- lapply(runs, function(r) readRDS(file.path(ROOT, r$rds)))

pngPath <- file.path(ROOT, "sbc-warmup-discriminating-traces.png")
grDevices::png(pngPath, width = 1200, height = 900)
graphics::par(mfrow = c(2L, 2L))
for (lbl in names(runs)) {
  d <- data[[lbl]]
  tr <- d$tl_trace
  ylim_top <- max(c(tr, d$tl_true, 50), na.rm = TRUE) * 1.05
  ylim_bot <- 0
  plot(seq_along(tr), tr, type = "l", col = "steelblue",
       xlab = "kept sample index", ylab = "tree_length",
       main = sprintf("(%s) %s", lbl, runs[[lbl]]$title),
       ylim = c(ylim_bot, ylim_top))
  graphics::abline(h = d$tl_true, col = "red", lwd = 2)
  graphics::abline(h = 50, col = "darkgreen", lty = 2)        # prior mean
  graphics::abline(h = d$start_tl, col = "grey", lty = 3)     # start_tl
  graphics::legend("topright",
                   legend = c(sprintf("tl_trace (rank = %d/%d)",
                                      sum(tr < d$tl_true), length(tr)),
                              sprintf("tl_true = %.2f", d$tl_true),
                              "prior mean = 50",
                              sprintf("start_tl = %.2f", d$start_tl)),
                   col = c("steelblue", "red", "darkgreen", "grey"),
                   lty = c(1, 1, 2, 3), lwd = c(1, 2, 1, 1), cex = 0.8)
}
grDevices::dev.off()
cat(sprintf("Wrote %s\n", pngPath))
