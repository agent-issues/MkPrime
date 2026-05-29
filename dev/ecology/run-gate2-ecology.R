# GATE-2 ground-truth runner: all ecology + EBE test files (post neo-z fix).
#   Rscript dev/ecology/run-gate2-ecology.R
suppressMessages(devtools::load_all(".", quiet = TRUE))
library(testthat)
cat("=== Ecology + EBE test suite ===\n")
res <- test_dir("tests/testthat", filter = "ecology|ebe",
                reporter = "summary", stop_on_failure = FALSE)
df  <- as.data.frame(res)
cat(sprintf("\n=== TOTAL: tests=%d  PASS=%d  FAIL=%d  WARN=%d  SKIP=%d ===\n",
            nrow(df), sum(df$passed), sum(df$failed),
            sum(df$warning > 0), sum(df$skipped)))
fails <- df[df$failed > 0, c("file", "test", "failed"), drop = FALSE]
if (nrow(fails) > 0) { cat("\n--- FAILING tests ---\n"); print(fails, row.names = FALSE) } else cat("ALL PASS\n")
