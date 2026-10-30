#!/usr/bin/env Rscript
# Summarizes mcmc.R (or exact.R) output for the #179 study: per config, the
# mean paired difference empirical_geometric - geometric in each score, with
# its standard error over datasets.
#
# Usage, from the repository root:
#   Rscript dev/red-team/heavy-tests/empirical-benefit/aggregate.R [dir|file]

args <- commandArgs(trailingOnly = TRUE)
src <- if (length(args)) args[[1]] else
  "dev/red-team/heavy-tests/empirical-benefit/results"
rows <- if (dir.exists(src)) {
  files <- list.files(src, "^mcmc-[0-9]+\\.rds$", full.names = TRUE)
  do.call(rbind, lapply(files, function(f) readRDS(f)$rows))
} else {
  readRDS(src)
}

scores <- intersect(
  c("logScore", "exactLogScore", "cover", "bias", "rmse", "cid", "rf"),
  names(rows)
)
keys <- intersect(c("truth", "pTrue", "expSteps"), names(rows))
rows$dataset <- rows$seed
emp <- rows[rows$arm == "empirical_geometric", ]
geo <- rows[rows$arm == "geometric", ]
geo <- geo[match(emp$dataset, geo$dataset), ]
stopifnot(identical(emp$dataset, geo$dataset))

config <- interaction(emp[keys], drop = TRUE)
summary <- do.call(rbind, lapply(split(seq_len(nrow(emp)), config), function(i) {
  out <- data.frame(emp[i[[1]], keys, drop = FALSE], n = length(i))
  for (s in scores) {
    d <- emp[[s]][i] - geo[[s]][i]
    out[[paste0(s, ".EG")]] <- mean(emp[[s]][i])
    out[[paste0(s, ".Geo")]] <- mean(geo[[s]][i])
    out[[paste0("d.", s)]] <- mean(d)
    out[[paste0("se.", s)]] <- stats::sd(d) / sqrt(length(i))
  }
  out
}))
rownames(summary) <- NULL
print(format(summary, digits = 3), row.names = FALSE)
if ("wall" %in% names(rows)) {
  cat(sprintf("\nWall per fit: median %.0f s, max %.0f s; draws per fit %d\n",
              stats::median(rows$wall), max(rows$wall),
              stats::median(rows$nDraw)))
}
