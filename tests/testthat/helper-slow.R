# Skip guard for slow MCMC integration tests.
# Set MKPRIME_SLOW_TESTS=true to run these (e.g. in CI).
skip_slow_tests <- function() {
  testthat::skip_if_not(
    identical(Sys.getenv("MKPRIME_SLOW_TESTS"), "true"),
    "Slow MCMC tests skipped (set MKPRIME_SLOW_TESTS=true to run)"
  )
}
# Skip guard for tests whose cost is MCMC iterations rather than code paths.
# mem-check sets MKPRIME_MEMCHECK=true; valgrind costs ~15x AddressSanitizer,
# which puts the full suite ~9 h past GitHub's 6 h job ceiling.  The guarded
# blocks are integration runs whose C++ entry points are already reached by
# cheaper tests in the same file, so skipping them costs depth, not breadth.
# gcc-ASAN still runs every one of them.
skip_under_memcheck <- function() {
  testthat::skip_if(
    identical(Sys.getenv("MKPRIME_MEMCHECK"), "true"),
    "Long MCMC run skipped under valgrind (see .github/workflows/memcheck.yml)"
  )
}
