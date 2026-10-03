test_that(".MkLaunchScript parses when .libPaths() deparses to several lines", {
  longLibs <- file.path("/very/long/library/path", paste0("site-library-", 1:6))
  longDir  <- file.path(tempdir(), strrep("nested-directory-", 8))
  script <- .MkLaunchScript(longDir, libPaths = longLibs)
  expect_gt(length(deparse(longLibs)), 1L)
  expect_no_error(parse(text = script))
  exprs <- parse(text = script)
  expect_equal(eval(exprs[[1]][[2]]), longLibs)
})


test_that(".MkLaunchScript resumes only when the resume flag exists", {
  script <- .MkLaunchScript("/tmp/x")
  expect_true(any(grepl("overwrite = !file.exists(", script, fixed = TRUE)))
  expect_true(any(grepl("mkp_resume.flag", script, fixed = TRUE)))
  expect_identical(basename(.MkResumeFlagPath("/tmp/x")), "mkp_resume.flag")
})


test_that("GUI log paths match the files MkPrimeMCMC/RunMkPrime write", {
  logDir <- tempfile("gui-logs")
  for (nRuns in 1:2) {
    mcmc <- MkPrimeMCMC(nRuns = nRuns, logFile = .MkLogBase(logDir))
    expect_identical(
      MkLogPaths(.MkLogBase(logDir), nRuns),
      .LogFilePaths(mcmc$logFile, nRuns)
    )
  }
  expect_identical(
    basename(MkLogPaths(.MkLogBase(logDir), 2L)),
    c("run_1.log", "run_2.log")
  )
})


test_that(".MkJobStatus tests the cancel file before the done signal", {
  logDir <- withr::local_tempdir()
  cancelFile <- MkCancelPath(logDir)
  expect_identical(.MkJobStatus(logDir, cancelFile, NA), "running")
  expect_identical(.MkJobStatus(logDir, cancelFile, TRUE), "running")

  file.create(file.path(logDir, "mkp_done.signal"))
  expect_identical(.MkJobStatus(logDir, cancelFile, FALSE), "done")

  file.create(cancelFile)
  expect_identical(.MkJobStatus(logDir, cancelFile, FALSE), "cancelled")
  expect_identical(.MkJobStatus(logDir, cancelFile, NA), "cancelled")
  expect_identical(.MkJobStatus(logDir, cancelFile, TRUE), "running")
})


test_that(".MkJobStatus reports error and dead-without-signal jobs", {
  logDir <- withr::local_tempdir()
  cancelFile <- MkCancelPath(logDir)
  expect_identical(.MkJobStatus(logDir, cancelFile, FALSE), "failed")
  expect_identical(.MkJobStatus(logDir, cancelFile, NA), "running")
  writeLines("boom", file.path(logDir, "mkp_error.txt"))
  expect_identical(.MkJobStatus(logDir, cancelFile, FALSE), "error")
})


test_that(".TailFile returns the last lines, or NULL when absent or empty", {
  f <- withr::local_tempfile()
  expect_null(.TailFile(f))
  file.create(f)
  expect_null(.TailFile(f))
  writeLines(as.character(1:20), f)
  expect_identical(.TailFile(f, 3L), c("18", "19", "20"))
})
