test_that("MkBayesianServer initialises with idle status", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")

  shiny::testServer(
    MkBayesianServer,
    args = list(dataset = shiny::reactive(NULL)),
    expr = {
      expect_equal(status(), "idle")
      expect_null(jobFile())
      expect_null(trees())
    }
  )
})


test_that("MkBayesianServer: Run without data shows warning, stays idle", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")

  shiny::testServer(
    MkBayesianServer,
    args = list(dataset = shiny::reactive(NULL)),
    expr = {
      session$setInputs(run = 1L)
      expect_equal(status(), "idle")
    }
  )
})


test_that("MkBayesianServer: Reconnect with missing logDir shows error, stays idle", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")

  shiny::testServer(
    MkBayesianServer,
    args = list(dataset = shiny::reactive(NULL)),
    expr = {
      session$setInputs(logDir = file.path(tempdir(), "no_such_dir_xyz"), reconnect = 1L)
      expect_equal(status(), "idle")
    }
  )
})


test_that(".ParseIntList handles edge cases", {
  expect_equal(.ParseIntList(""),         integer(0))
  expect_equal(.ParseIntList("   "),      integer(0))
  expect_equal(.ParseIntList("1,3,5"),    c(1L, 3L, 5L))
  expect_equal(.ParseIntList("1; 3; 5"),  c(1L, 3L, 5L))
  expect_equal(.ParseIntList("1 3 5"),    c(1L, 3L, 5L))
  expect_equal(.ParseIntList("2,NA,4"),   c(2L, 4L))
})


test_that(".MkLaunchScript produces a runnable-looking script", {
  logDir <- file.path(tempdir(), "mkp_script_test")
  dir.create(logDir, showWarnings = FALSE)
  on.exit(unlink(logDir, recursive = TRUE))

  lines <- .MkLaunchScript(logDir)
  expect_type(lines, "character")
  expect_true(any(grepl("library(MkPrime)", lines, fixed = TRUE)))
  expect_true(any(grepl("RunMkPrime", lines, fixed = TRUE)))
  expect_true(any(grepl("mkp_done.signal", lines, fixed = TRUE)))
  expect_true(any(grepl(".libPaths", lines, fixed = TRUE)))
})


test_that(".ReadAllLogs returns NULL when all files absent", {
  expect_null(.ReadAllLogs(c("/no/such/a.log", "/no/such/b.log")))
})
