#!/usr/bin/env Rscript
files <- Sys.glob("dev/m9-pilot/syab*/cid_summary.rds")
all <- do.call(rbind, lapply(files, readRDS))
all <- all[order(all$matrix, all$model), ]
cat("=== 6-matrix x 3-model CID-to-WCT summary ===\n\n")
print(all, row.names = FALSE, digits = 4)
cat("\n=== Wide view (mean CID) ===\n")
wide <- reshape(all[, c("matrix", "model", "cid_mean")],
                idvar = "matrix", timevar = "model", direction = "wide")
names(wide) <- sub("cid_mean.", "", names(wide))
print(wide, row.names = FALSE, digits = 4)
saveRDS(all, "dev/m9-pilot/cid_summary_all.rds")
cat("\nSaved: dev/m9-pilot/cid_summary_all.rds\n")
