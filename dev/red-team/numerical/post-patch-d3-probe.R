#!/usr/bin/env Rscript
# Probe degenerate-input behaviour of .ComputeRhat after CONV-002 patch.
# Cases: empty perRun, all-zero-row runs, single-row tail, single-run.

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
})

cat("=== Case A: empty perRun (length 0) ===\n")
out <- tryCatch(
  MkPrime:::.ComputeRhat(list(), integer(0)),
  error = function(e) sprintf("ERROR: %s", conditionMessage(e)),
  warning = function(w) sprintf("WARNING: %s", conditionMessage(w))
)
print(out)

cat("\n=== Case B: one run, zero rows ===\n")
mk_one <- function(nrow_v, nm = "x") {
  m <- matrix(numeric(0), ncol = 1, nrow = nrow_v,
              dimnames = list(NULL, nm))
  list(samples = m)
}
out <- tryCatch(
  MkPrime:::.ComputeRhat(list(mk_one(0)), 1L),
  error = function(e) sprintf("ERROR: %s", conditionMessage(e)),
  warning = function(w) sprintf("WARNING: %s", conditionMessage(w))
)
print(out)

cat("\n=== Case C: two runs, all-zero rows ===\n")
out <- tryCatch(
  MkPrime:::.ComputeRhat(list(mk_one(0), mk_one(0)), 1L),
  error = function(e) sprintf("ERROR: %s", conditionMessage(e)),
  warning = function(w) sprintf("WARNING: %s", conditionMessage(w))
)
print(out)

cat("\n=== Case D: two runs, nKeep collapses to 1 row ===\n")
mk_rand <- function(n) {
  list(samples = matrix(rnorm(n), ncol = 1, dimnames = list(NULL, "x")))
}
set.seed(1)
out <- tryCatch(
  MkPrime:::.ComputeRhat(list(mk_rand(1L), mk_rand(50L)), 1L),
  error = function(e) sprintf("ERROR: %s", conditionMessage(e)),
  warning = function(w) sprintf("WARNING: %s", conditionMessage(w))
)
print(out)

cat("\n=== Case E: two runs, nKeep = 2 ===\n")
out <- tryCatch(
  MkPrime:::.ComputeRhat(list(mk_rand(2L), mk_rand(50L)), 1L),
  error = function(e) sprintf("ERROR: %s", conditionMessage(e)),
  warning = function(w) sprintf("WARNING: %s", conditionMessage(w))
)
print(out)

cat("\n=== Case F: two runs equal-length (sanity) ===\n")
set.seed(2)
out <- tryCatch(
  MkPrime:::.ComputeRhat(list(mk_rand(50L), mk_rand(50L)), 1L),
  error = function(e) sprintf("ERROR: %s", conditionMessage(e)),
  warning = function(w) sprintf("WARNING: %s", conditionMessage(w))
)
print(out)

cat("\n=== Case G: two runs unequal (the bug CONV-002 was meant to fix) ===\n")
set.seed(3)
out <- tryCatch(
  MkPrime:::.ComputeRhat(list(mk_rand(37L), mk_rand(50L)), 1L),
  error = function(e) sprintf("ERROR: %s", conditionMessage(e)),
  warning = function(w) sprintf("WARNING: %s", conditionMessage(w))
)
print(out)

cat("\n[done]\n")
