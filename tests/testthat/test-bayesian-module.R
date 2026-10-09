test_that("MkBayesianServer initialises with idle status", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")

  shiny::testServer(
    MkBayesianServer,
    args = list(dataset = shiny::reactive(NULL)),
    expr = {
      expect_equal(rv$status, "idle")
      ret <- session$getReturned()
      expect_null(ret$jobFile())
      expect_null(ret$trees())
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
      expect_equal(rv$status, "idle")
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
      expect_equal(rv$status, "idle")
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
  expect_true(any(grepl('library("MkPrime")', lines, fixed = TRUE)))
  expect_true(any(grepl("RunMkPrime", lines, fixed = TRUE)))
  expect_true(any(grepl("mkp_done.signal", lines, fixed = TRUE)))
  expect_true(any(grepl(".libPaths", lines, fixed = TRUE)))
})


test_that(".ReadAllLogs returns NULL when all files absent", {
  expect_null(.ReadAllLogs(c("/no/such/a.log", "/no/such/b.log")))
})


# ---------------------------------------------------------------------------
# .PidIsAlive (M-093)
# ---------------------------------------------------------------------------

test_that(".PidIsAlive: current process is alive", {
  expect_true(.PidIsAlive(Sys.getpid()))
})

test_that(".PidIsAlive: invalid PIDs return FALSE", {
  expect_false(.PidIsAlive(NA_integer_))
  expect_false(.PidIsAlive(NA_real_))
  expect_false(.PidIsAlive(0L))
  expect_false(.PidIsAlive(-1L))
  expect_false(.PidIsAlive(NULL))
})

test_that(".PidIsAlive: implausibly large PID returns FALSE", {
  # .Machine$integer.max is unlikely to be a running PID on any system
  expect_false(.PidIsAlive(.Machine$integer.max))
})


# ---------------------------------------------------------------------------
# Reconnect signal-file paths (M-082 + M-093)
# ---------------------------------------------------------------------------

# Helper: create a minimal job.rds in a temp dir
.make_test_job <- function(logDir, pid = NA_integer_, checkpointFile = NULL) {
  list(
    logDir         = logDir,
    logFiles       = character(0),
    cancelFile     = file.path(logDir, "mkp_cancel.signal"),
    checkpointFile = checkpointFile,
    scriptFile     = file.path(logDir, "mkp_run.R"),
    nRuns          = 1L,
    startTime      = Sys.time(),
    pid            = pid
  )
}

test_that("Reconnect: done signal → status 'done'", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")

  d <- file.path(tempdir(), "mkp_rc_done")
  dir.create(d, showWarnings = FALSE)
  on.exit(unlink(d, recursive = TRUE))
  saveRDS(.make_test_job(d), file.path(d, "job.rds"))
  file.create(file.path(d, "mkp_done.signal"))

  shiny::testServer(MkBayesianServer,
    args = list(dataset = shiny::reactive(NULL)),
    expr = {
      session$setInputs(logDir = d, reconnect = 1L)
      expect_equal(rv$status, "done")
    }
  )
})

test_that("Reconnect: cancel signal + no checkpoint → status 'cancelled'", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")

  d <- file.path(tempdir(), "mkp_rc_cancel")
  dir.create(d, showWarnings = FALSE)
  on.exit(unlink(d, recursive = TRUE))
  saveRDS(.make_test_job(d), file.path(d, "job.rds"))
  file.create(file.path(d, "mkp_cancel.signal"))

  shiny::testServer(MkBayesianServer,
    args = list(dataset = shiny::reactive(NULL)),
    expr = {
      session$setInputs(logDir = d, reconnect = 1L)
      expect_equal(rv$status, "cancelled")
    }
  )
})

test_that("Reconnect: error file + no checkpoint → status 'error'", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")

  d <- file.path(tempdir(), "mkp_rc_error")
  dir.create(d, showWarnings = FALSE)
  on.exit(unlink(d, recursive = TRUE))
  saveRDS(.make_test_job(d), file.path(d, "job.rds"))
  writeLines("object 'x' not found", file.path(d, "mkp_error.txt"))

  shiny::testServer(MkBayesianServer,
    args = list(dataset = shiny::reactive(NULL)),
    expr = {
      session$setInputs(logDir = d, reconnect = 1L)
      expect_equal(rv$status, "error")
    }
  )
})

test_that("Reconnect: no signals + dead PID + no checkpoint → status 'error'", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")

  d <- file.path(tempdir(), "mkp_rc_deadpid")
  dir.create(d, showWarnings = FALSE)
  on.exit(unlink(d, recursive = TRUE))
  # Use an implausibly large PID so .PidIsAlive() returns FALSE
  saveRDS(.make_test_job(d, pid = .Machine$integer.max),
          file.path(d, "job.rds"))

  shiny::testServer(MkBayesianServer,
    args = list(dataset = shiny::reactive(NULL)),
    expr = {
      session$setInputs(logDir = d, reconnect = 1L)
      expect_equal(rv$status, "error")
    }
  )
})


# ---------------------------------------------------------------------------
# Live-job guards and launch arguments (#421)
# ---------------------------------------------------------------------------

test_that("Reconnect: cancel signal + live PID keeps polling, no relaunch", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")

  d <- file.path(tempdir(), "mkp_rc_cancel_live")
  dir.create(d, showWarnings = FALSE)
  on.exit(unlink(d, recursive = TRUE))
  cp <- file.path(d, "checkpoint.rds")
  file.create(cp)
  writeLines("stop", file.path(d, "mkp_run.R"))
  saveRDS(.make_test_job(d, pid = Sys.getpid(), checkpointFile = cp),
          file.path(d, "job.rds"))
  cancelFile <- file.path(d, "mkp_cancel.signal")
  file.create(cancelFile)

  shiny::testServer(MkBayesianServer,
    args = list(dataset = shiny::reactive(NULL)),
    expr = {
      session$setInputs(logDir = d, reconnect = 1L)
      expect_equal(rv$status, "stopping")
      expect_true(file.exists(cancelFile))
      expect_null(rv$proc)
    }
  )
})

test_that("Run refuses to launch over a live job in the log directory", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("processx")

  d <- file.path(tempdir(), "mkp_run_live")
  dir.create(d, showWarnings = FALSE)
  on.exit(unlink(d, recursive = TRUE))
  saveRDS(.make_test_job(d, pid = Sys.getpid()), file.path(d, "job.rds"))
  cancelFile <- file.path(d, "mkp_cancel.signal")
  file.create(cancelFile)

  shiny::testServer(MkBayesianServer,
    args = list(dataset = shiny::reactive(
      list(a = c("0", "1", "0"), b = c("1", "0", "0"), c = c("0", "0", "1"),
           d = c("1", "1", "0")))),
    expr = {
      session$setInputs(logDir = d, nRuns = 1L, nChains = 1L, heat = 0.2,
                        warmup = 100L, minEss = 10L, maxTimeMins = 1,
                        neomorphic = "", run = 1L)
      expect_equal(rv$status, "idle")
      expect_true(file.exists(cancelFile))
      expect_false(file.exists(file.path(d, "input_mcmc.rds")))
    }
  )
})

test_that(".JobIsLive reflects the recorded PID", {
  d <- file.path(tempdir(), "mkp_joblive")
  dir.create(d, showWarnings = FALSE)
  on.exit(unlink(d, recursive = TRUE))
  expect_false(.JobIsLive(d))
  saveRDS(.make_test_job(d, pid = Sys.getpid()), file.path(d, "job.rds"))
  expect_true(.JobIsLive(d))
  saveRDS(.make_test_job(d, pid = .Machine$integer.max),
          file.path(d, "job.rds"))
  expect_false(.JobIsLive(d))
})

test_that("GUI launch passes maxWarmup, so autoTune stays on", {
  expect_no_warning(
    mcmc <- MkPrimeMCMC(nRuns = 1L, nChains = 1L, maxWarmup = 3000L,
                        minEss = 10L)
  )
  expect_true(mcmc$autoTune)
  src <- paste(deparse(body(MkBayesianServer)), collapse = "\n")
  expect_match(src, "maxWarmup\\s*=\\s*maxWarmup")
  expect_no_match(src, "[^x]warmup\\s*=\\s*warmup")
})

test_that(".TraceTails keeps the last rows of each run", {
  m <- cbind(Sample = 1:7, x = 1:7)
  attr(m, "runRows") <- c(4L, 3L)
  tails <- .TraceTails(m, n = 2L)
  expect_length(tails, 2L)
  expect_equal(tails[[1]][, "Sample"], 3:4)
  expect_equal(tails[[2]][, "Sample"], 6:7)
  expect_length(.TraceTails(m[1:3, ], n = 2L), 1L)
})
