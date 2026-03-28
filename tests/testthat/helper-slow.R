# Skip guard for slow MCMC integration tests.
# Set MKPRIME_SLOW_TESTS=true to run these (e.g. in CI).
skip_slow_tests <- function() {
  testthat::skip_if_not(
    identical(Sys.getenv("MKPRIME_SLOW_TESTS"), "true"),
    "Slow MCMC tests skipped (set MKPRIME_SLOW_TESTS=true to run)"
  )
}
