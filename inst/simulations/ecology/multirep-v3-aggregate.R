# multirep-v3-aggregate.R — aggregate sim3-multirep-v3 HPC results.
# Reads per-rep summary.rds from multirep-v3-results/rep0X and prints
# headline aware-vs-blind table.
resDir <- "inst/simulations/ecology/multirep-v3-results"
repDirs <- sort(list.files(resDir, pattern = "^rep\\d{2}$", full.names = TRUE))
summaries <- lapply(repDirs, function(d) readRDS(file.path(d, "summary.rds")))

agg <- do.call(rbind, lapply(summaries, function(s) {
  data.frame(rep = s$repId,
             seed = s$seed,
             pTrue_blind  = s$blind$pTrue,
             pTrue_aware  = s$aware$pTrue,
             pWrong_blind = s$blind$pWrong,
             pWrong_aware = s$aware$pWrong,
             cidTrue_blind = s$blind$cidTrue,
             cidTrue_aware = s$aware$cidTrue)
}))
cat("== Per-replicate ==\n")
print(agg, row.names = FALSE)

cat("\n== Mean +/- sd across 8 reps ==\n")
for (col in setdiff(colnames(agg), c("rep", "seed"))) {
  x <- agg[[col]]
  cat(sprintf("  %-18s %.3f +/- %.3f\n", col, mean(x), stats::sd(x)))
}

cat("\n== Headline aware vs blind ==\n")
cat(sprintf("  Delta P(true)  = %+.3f\n",
            mean(agg$pTrue_aware - agg$pTrue_blind)))
cat(sprintf("  Delta P(wrong) = %+.3f\n",
            mean(agg$pWrong_aware - agg$pWrong_blind)))
cat(sprintf("  Delta CIDtrue  = %+.3f  (lower = closer to truth)\n",
            mean(agg$cidTrue_aware - agg$cidTrue_blind)))
