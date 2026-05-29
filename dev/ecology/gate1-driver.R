# GATE 1 driver: compile + load + run the EBE kernel<->oracle cross-check.
# Run via:  Rscript dev/ecology/gate1-driver.R   (from the package root)
options(warn = 1)
t0 <- Sys.time()
cat("=== GATE 1: compile + load ===\n")
suppressMessages({
  Rcpp::compileAttributes()
  pkgbuild::compile_dll()
  devtools::load_all(".", quiet = TRUE)
})
cat(sprintf("compile+load OK (%.0fs)\n", as.numeric(Sys.time() - t0, units = "secs")))

library(testthat)
cat("\n=== GATE 1: test-ebe-likelihood.R (T1/T2/T3/T4) ===\n")
# test_dir sources helper-ebe-oracle.R; filter selects only the EBE file.
res <- test_dir("tests/testthat", filter = "ebe-likelihood",
                reporter = "summary", stop_on_failure = FALSE)
df <- as.data.frame(res)
cat("\n=== PER-TEST RESULT ===\n")
print(df[, c("test", "nb", "failed", "skipped", "error", "warning")],
      row.names = FALSE)
cat(sprintf("\nTOTAL: %d expectations, %d failed, %d errors, %d skipped\n",
            sum(df$nb), sum(df$failed), sum(df$error), sum(df$skipped)))
cat(if (sum(df$failed) == 0 && sum(df$error) == 0)
      "GATE1: PASS\n" else "GATE1: FAIL\n")
