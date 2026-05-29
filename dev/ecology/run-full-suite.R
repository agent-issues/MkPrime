# Full-suite GATE-2 regression check (slow tests skip unless MKPRIME_SLOW_TESTS).
#   Rscript dev/ecology/run-full-suite.R
suppressMessages(devtools::load_all(".", quiet = TRUE))
library(testthat)
cat("=== FULL test suite ===\n")
res <- test_dir("tests/testthat", reporter = "summary", stop_on_failure = FALSE)
df  <- as.data.frame(res)
cat(sprintf("\n=== TOTAL: tests=%d  PASS=%d  FAIL=%d  WARN=%d  SKIP=%d ===\n",
            nrow(df), sum(df$passed), sum(df$failed),
            sum(df$warning > 0), sum(df$skipped)))
fails <- df[df$failed > 0, c("file", "test", "failed"), drop = FALSE]
if (nrow(fails) > 0) { cat("\n--- FAILING ---\n"); print(fails, row.names = FALSE) } else cat("ALL PASS (0 failures)\n")
