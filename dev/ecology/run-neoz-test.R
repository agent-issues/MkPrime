# Runner for the EBE neomorphic-only-z red-green test.
# Run RED (current build) before the mask edits, GREEN after recompile.
#   Rscript dev/ecology/run-neoz-test.R
suppressMessages(devtools::load_all(".", quiet = TRUE))
library(testthat)
cat("=== test-ebe-neoz-prior.R ===\n")
res <- test_file("tests/testthat/test-ebe-neoz-prior.R", reporter = "summary")
df  <- as.data.frame(res)
cat(sprintf("\n=== SUMMARY: tests=%d  PASS=%d  FAIL=%d  SKIP=%d ===\n",
            nrow(df), sum(df$passed), sum(df$failed), sum(df$skipped)))
if (sum(df$failed) > 0) cat("RESULT: RED (failures present)\n") else cat("RESULT: GREEN (all pass)\n")
