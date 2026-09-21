test_that("MkPrimeVerbosity() reads the option", {
  withr::local_options(MkPrime.verbosity = NULL)
  expect_equal(MkPrimeVerbosity(), 1L)

  withr::local_options(MkPrime.verbosity = 0)
  expect_equal(MkPrimeVerbosity(), 0L)

  withr::local_options(MkPrime.verbosity = 2.7)
  expect_equal(MkPrimeVerbosity(), 2L)
})

test_that("MkPrimeVerbosity(n) sets the option and returns the old level", {
  withr::local_options(MkPrime.verbosity = 1)

  # Returns the level in force BEFORE the call, not the one it installs.
  expect_equal(MkPrimeVerbosity(0), 1L)
  expect_equal(getOption("MkPrime.verbosity"), 0L)
  expect_equal(MkPrimeVerbosity(), 0L)

  expect_equal(MkPrimeVerbosity(2), 0L)
  expect_equal(MkPrimeVerbosity(), 2L)

  # The option is stored as an integer, whatever was passed.
  expect_identical(getOption("MkPrime.verbosity"), 2L)
  MkPrimeVerbosity(1.9)
  expect_identical(getOption("MkPrime.verbosity"), 1L)
})

test_that("MkPrimeVerbosity(n) returns invisibly, and reads visibly", {
  withr::local_options(MkPrime.verbosity = 1)
  expect_invisible(MkPrimeVerbosity(0))
  expect_visible(MkPrimeVerbosity())
})

test_that("MkPrimeVerbosity(old) restores the previous level", {
  withr::local_options(MkPrime.verbosity = 2)

  Quiet <- function() {
    old <- MkPrimeVerbosity(0)
    on.exit(MkPrimeVerbosity(old), add = TRUE)
    MkPrimeVerbosity()
  }

  expect_equal(Quiet(), 0L)
  expect_equal(MkPrimeVerbosity(), 2L)
})

test_that("MkPrimeVerbosity(n) rejects a bad n without touching the option", {
  withr::local_options(MkPrime.verbosity = 2)

  expect_error(MkPrimeVerbosity("loud"), "must be a single number")
  expect_equal(getOption("MkPrime.verbosity"), 2)

  expect_error(MkPrimeVerbosity(NA), "must be a single number")
  expect_error(MkPrimeVerbosity(1:2), "must be a single number")
  expect_equal(getOption("MkPrime.verbosity"), 2)
})

test_that("MkPrimeVerbosity(n) reports the fallback when the option is junk", {
  withr::local_options(MkPrime.verbosity = "loud")

  # The level in force was effectively the 1 fallback, so that is what a set
  # returns -- restoring it repairs the option instead of reinstating junk.
  expect_warning(old <- MkPrimeVerbosity(0), "must be a single number")
  expect_equal(old, 1L)
  expect_equal(getOption("MkPrime.verbosity"), 0L)

  expect_silent(MkPrimeVerbosity(old))
  expect_equal(getOption("MkPrime.verbosity"), 1L)
})

test_that("MkPrimeVerbosity() warns on a nonsense option", {
  withr::local_options(MkPrime.verbosity = "loud")
  expect_warning(v <- MkPrimeVerbosity(), "must be a single number")
  expect_equal(v, 1L)

  withr::local_options(MkPrime.verbosity = c(1, 2))
  expect_warning(MkPrimeVerbosity(), "must be a single number")

  withr::local_options(MkPrime.verbosity = NA_integer_)
  expect_warning(MkPrimeVerbosity(), "must be a single number")
})

test_that(".Loud() tracks the level", {
  withr::local_options(MkPrime.verbosity = 0)
  expect_false(.Loud())
  expect_false(.Loud(2L))

  withr::local_options(MkPrime.verbosity = 1)
  expect_true(.Loud())
  expect_false(.Loud(2L))

  withr::local_options(MkPrime.verbosity = 2)
  expect_true(.Loud(2L))
})

test_that(".CheckVerbosity() rejects non-numbers", {
  expect_equal(.CheckVerbosity(2), 2L)
  expect_error(.CheckVerbosity("0"), "must be a single number")
  expect_error(.CheckVerbosity(NULL), "must be a single number")
  expect_error(.CheckVerbosity(NA), "must be a single number")
  expect_error(.CheckVerbosity(1:2), "must be a single number")
})

test_that("Output wrappers honour the level", {
  local({
    withr::local_options(MkPrime.verbosity = 0)
    expect_silent(.AlertInfo("info"))
    expect_silent(.AlertSuccess("success"))
    expect_silent(.AlertWarning("warning"))
    expect_silent(.AlertDanger("danger"))
    expect_silent(.Inform("inform"))
    expect_silent(.Text("text"))
    expect_silent(.ProgressBar("bar", total = 2))
    expect_silent(.ProgressUpdate(set = 1))
    expect_silent(.ProgressDone())
  })

  withr::local_options(MkPrime.verbosity = 1)
  expect_message(.AlertInfo("info"), "info")
  expect_message(.AlertSuccess("success"), "success")
  expect_message(.AlertWarning("warning"), "warning")
  expect_message(.AlertDanger("danger"), "danger")
  expect_message(.Inform("inform"), "inform")
  expect_message(.Text("text"), "text")
  expect_silent(.AlertInfo("diagnostic", level = 2L))

  withr::local_options(MkPrime.verbosity = 2)
  expect_message(.AlertInfo("diagnostic", level = 2L), "diagnostic")
})

test_that("Wrappers interpolate in the caller's frame", {
  withr::local_options(MkPrime.verbosity = 1)
  Caller <- function() {
    nTip <- 7L
    .AlertInfo("Tree has {nTip} tips.")
  }
  expect_message(Caller(), "Tree has 7 tips")
})

test_that("verbosity = 0 silences a run", {
  skip_if_not_installed("TreeTools")
  withr::local_options(MkPrime.verbosity = 1)
  set.seed(1)
  # A 200-iteration chain legitimately warns that it never stabilised; that
  # is a condition, not console output, and verbosity does not govern it.
  Run <- function() {
    suppressWarnings(
      RunMkPrime(.mkp_test_pd(), .mkp_test_tree(), nIter = 200L,
                 minWarmup = 50L, maxWarmup = 50L, nRuns = 1L, thin = 10L,
                 autoTune = FALSE, verbosity = 0)
    )
  }
  expect_equal(capture.output(expect_no_message(Run())), character(0))
  # The option is restored, so a second run is audible again.
  expect_equal(MkPrimeVerbosity(), 1L)
})

test_that("RunMkPrime() rejects a nonsense verbosity", {
  expect_error(
    RunMkPrime(.mkp_test_pd(), .mkp_test_tree(), verbosity = "silent"),
    "must be a single number"
  )
})
