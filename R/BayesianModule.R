# R/BayesianModule.R
# Reusable Shiny module: launch MkPrime MCMC as a detached Rscript process,
# poll TSV log files for live progress, and support Reconnect after a session
# restart.  Depends on M-075 (logFile streaming) and M-081 (cancelFile).

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# Parse a comma/semicolon/space-separated string of integers.
# @keywords internal
.ParseIntList <- function(text) {
  text <- trimws(text)
  if (!nzchar(text)) return(integer(0))
  vals <- suppressWarnings(
    as.integer(strsplit(text, "[,;[:space:]]+")[[1]])
  )
  vals[!is.na(vals)]
}

# Try to read one or more log files and row-bind the results.
# Silently skips missing or unreadable files.
# @keywords internal
.ReadAllLogs <- function(logFiles) {
  mats <- lapply(logFiles, function(f) {
    if (!file.exists(f)) return(NULL)
    tryCatch(ReadMkLog(f), error = function(e) NULL)
  })
  mats <- Filter(Negate(is.null), mats)
  if (!length(mats)) return(NULL)
  do.call(rbind, mats)
}

# Check whether a PID is still alive.
# Uses the ps package when available; falls back to platform-specific methods.
# Returns FALSE for NA/invalid PIDs so callers can treat them as dead.
# @keywords internal
.PidIsAlive <- function(pid) {
  if (is.null(pid) || is.na(pid) || pid <= 0L) return(FALSE)
  if (requireNamespace("ps", quietly = TRUE)) {
    h <- tryCatch(ps::ps_handle(pid = as.integer(pid)), error = function(e) NULL)
    if (is.null(h)) return(FALSE)
    return(tryCatch(ps::ps_is_running(h), error = function(e) FALSE))
  }
  # Windows fallback: tasklist /FI
  if (.Platform$OS.type == "windows") {
    out <- tryCatch(
      system2("tasklist",
              c("/FI", sprintf("PID eq %d", as.integer(pid)), "/NH"),
              stdout = TRUE, stderr = FALSE),
      error = function(e) character(0)
    )
    return(any(grepl(as.character(as.integer(pid)), out, fixed = TRUE)))
  }
  # POSIX fallback: signal 0 checks existence without killing
  tryCatch(tools::pskill(as.integer(pid), 0L) == 0L, error = function(e) FALSE)
}

# Relaunch a stopped run from its checkpoint.
# Clears stale signal files, relaunches the existing mkp_run.R script
# (which calls RunMkPrime with overwrite=FALSE, auto-resuming from the
# checkpoint stored in input_mcmc.rds), and updates rv + job.rds.
# Returns TRUE on success, FALSE if processx is unavailable or scriptFile
# is missing.
# @keywords internal
.RelaunchFromCheckpoint <- function(job, rv) {
  if (!requireNamespace("processx", quietly = TRUE)) return(FALSE)

  scriptFile <- job$scriptFile
  if (is.null(scriptFile) || !file.exists(scriptFile)) return(FALSE)

  # Clear stale signals so the fresh process starts clean
  for (sig in c(job$cancelFile,
                file.path(job$logDir, "mkp_done.signal"),
                file.path(job$logDir, "mkp_error.txt"))) {
    if (!is.null(sig) && file.exists(sig)) file.remove(sig)
  }

  proc <- processx::process$new(
    "Rscript",
    args      = scriptFile,
    supervise = FALSE,
    stdout    = file.path(job$logDir, "mkp_stdout.txt"),
    stderr    = file.path(job$logDir, "mkp_stderr.txt")
  )

  job$pid       <- proc$get_pid()
  job$startTime <- Sys.time()
  saveRDS(job, file.path(job$logDir, "job.rds"))

  rv$job      <- job
  rv$proc     <- proc
  rv$status   <- "running"
  rv$errorMsg <- NULL
  # Keep rv$logSamples: existing partial data stays visible during resume
  TRUE
}

# Build the R script run by the detached Rscript process.
# All inputs are read from RDS files in logDir.
# @keywords internal
.MkLaunchScript <- function(logDir) {
  c(
    sprintf('.libPaths(%s)', deparse(.libPaths())),
    'suppressPackageStartupMessages(library("MkPrime"))',
    sprintf('.d <- %s', deparse(logDir)),
    'data <- readRDS(file.path(.d, "input_data.rds"))',
    'tree <- readRDS(file.path(.d, "input_tree.rds"))',
    'neo  <- readRDS(file.path(.d, "input_neo.rds"))',
    'mcmc <- readRDS(file.path(.d, "input_mcmc.rds"))',
    'result <- tryCatch(',
    '  RunMkPrime(data, tree, neomorphic = neo, mcmc = mcmc),',
    '  error = function(e) {',
    '    writeLines(conditionMessage(e), file.path(.d, "mkp_error.txt"))',
    '    NULL',
    '  }',
    ')',
    'if (!is.null(result)) {',
    '  # Write to a temp file then rename so result.rds is either complete or',
    '  # absent -- prevents a partial file if the process is killed during saveRDS.',
    '  tmp <- paste0(file.path(.d, "result.rds"), ".tmp")',
    '  saveRDS(result, tmp)',
    '  file.rename(tmp, file.path(.d, "result.rds"))',
    '  file.create(file.path(.d, "mkp_done.signal"))',
    '}'
  )
}


# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------

#' Shiny UI module for Bayesian Mk' analysis
#'
#' Renders an accordion with MCMC configuration, output settings, and
#' character-type overrides, together with Run / Stop / Reconnect controls
#' and a live progress area (trace plot + ESS table polled from disk log
#' files).  Pair with [MkBayesianServer()].
#'
#' @param id Shiny module namespace ID (character scalar).
#' @return A `tagList` suitable for embedding in a sidebar or tab panel.
#' @seealso [MkBayesianServer()], [RunMkPrime()], [MkPrimeMCMC()]
#' @export
MkBayesianUi <- function(id) {
  ns <- shiny::NS(id)

  shiny::tagList(
    bslib::accordion(
      id       = ns("settings"),
      open     = c("MCMC", "Output"),
      multiple = TRUE,

      bslib::accordion_panel(
        "MCMC",
        shiny::numericInput(ns("nRuns"),   "Independent runs",   2L,
                            min = 1L, max = 20L, step = 1L),
        shiny::numericInput(ns("nChains"), "Chains per run",     1L,
                            min = 1L, max = 8L,  step = 1L),
        shiny::conditionalPanel(
          paste0("input['", ns("nChains"), "'] > 1"),
          shiny::sliderInput(ns("heat"), "Heat",
                             min = 0.05, max = 0.95, value = 0.2, step = 0.05)
        ),
        shiny::numericInput(ns("warmup"),  "Warmup iterations", 5000L,
                            min = 0L, step = 500L),
        shiny::numericInput(ns("minEss"),  "Target ESS",         200L,
                            min = 10L, step = 50L),
        shiny::numericInput(ns("maxTimeMins"), "Max time (min)",  NA_real_,
                            min = 1)
      ),

      bslib::accordion_panel(
        "Output",
        shiny::textInput(ns("logDir"), "Log directory",
                         placeholder = "Folder for log files and checkpoint")
      ),

      bslib::accordion_panel(
        "Characters",
        shiny::textInput(ns("neomorphic"), "Neomorphic characters",
                         placeholder = "e.g., 1,3,5"),
        shiny::uiOutput(ns("charInfo"))
      )
    ),

    shiny::br(),

    shiny::fluidRow(
      shiny::column(4, shiny::actionButton(
        ns("run"), "Run", class = "btn-primary btn-sm w-100")),
      shiny::column(4, shiny::actionButton(
        ns("stop"), "Stop", class = "btn-outline-danger btn-sm w-100")),
      shiny::column(4, shiny::actionButton(
        ns("reconnect"), "Reconnect", class = "btn-outline-secondary btn-sm w-100"))
    ),

    shiny::uiOutput(ns("statusBadge")),
    shiny::hr(),
    shiny::plotOutput(ns("tracePlot"), height = "260px"),
    shiny::tableOutput(ns("essTable"))
  )
}


# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------

#' Shiny server module for Bayesian Mk' analysis
#'
#' Launches MCMC as a detached `Rscript` process (survives Shiny session
#' close / browser reload), writes a `job.rds` alongside the log files so
#' the "Reconnect" button can resume monitoring after a session restart.
#' Progress is polled every five seconds by reading the streaming TSV log
#' files produced by [RunMkPrime()].
#'
#' @param id Shiny module namespace ID -- must match the `id` passed to
#'   [MkBayesianUi()].
#' @param dataset A *reactive* that returns a `phyDat` object, or `NULL`
#'   when no dataset is loaded.  Supplied by the host app.
#' @param startTree A *reactive* that returns a starting `phylo`, or
#'   `NULL` to generate a random unrooted tree.  Optional.
#'
#' @return A named list of reactives:
#'   \describe{
#'     \item{`jobFile`}{Path to `job.rds` for the current or most recent
#'       run, or `NULL` if no run has been started.}
#'     \item{`trees`}{A `multiPhylo` of post-burnin trees when the run is
#'       complete, otherwise `NULL`.}
#'     \item{`status`}{Character scalar: `"idle"`, `"running"`,
#'       `"done"`, `"cancelled"`, or `"error"`.}
#'   }
#'
#' @section Requirements:
#' The `processx` package must be installed.  Add it via
#' `install.packages("processx")`.
#'
#' @seealso [MkBayesianUi()], [RunMkPrime()], [MkPrimeMCMC()],
#'   [ReadMkLog()], [MkCancelPath()]
#' @export
MkBayesianServer <- function(id, dataset, startTree = NULL) {
  shiny::moduleServer(id, function(input, output, session) {

    rv <- shiny::reactiveValues(
      status     = "idle",  # idle | running | done | cancelled | error
      job        = NULL,    # list: logDir, logFiles, cancelFile, checkpointFile, nRuns, startTime, pid
      proc       = NULL,    # processx::process handle (NULL after session restart)
      logSamples = NULL,    # matrix from .ReadAllLogs() -- latest poll
      errorMsg   = NULL
    )

    # ---- Auto-populate neomorphic from dataset --------------------------------

    shiny::observeEvent(dataset(), {
      dat <- dataset()
      if (is.null(dat)) return()
      neo <- tryCatch(AutoDetectNeomorphic(dat), error = function(e) integer(0))
      shiny::updateTextInput(session, "neomorphic",
                             value = paste(neo, collapse = ", "))
      if (length(neo)) {
        shiny::showNotification(
          sprintf("Auto-detected %d neomorphic character%s.",
                  length(neo), if (length(neo) == 1L) "" else "s"),
          type = "message", duration = 5
        )
      }
    }, ignoreNULL = TRUE)

    # ---- Character info summary -----------------------------------------------

    output$charInfo <- shiny::renderUI({
      dat <- dataset()
      if (is.null(dat)) return(NULL)
      neo <- .ParseIntList(input$neomorphic)
      mkd <- tryCatch(
        MkPrimeData(dat, neomorphic = neo),
        warning = function(w) invokeRestart("muffleWarning"),
        error   = function(e) NULL
      )
      if (is.null(mkd)) return(NULL)
      types <- table(mkd$type)
      shiny::tags$p(
        class = "text-muted small mt-1",
        sprintf("%d characters: %s", mkd$nChar,
                paste(sprintf("%d %s", types, names(types)), collapse = ", "))
      )
    })

    # ---- Run -----------------------------------------------------------------

    shiny::observeEvent(input$run, {
      if (rv$status == "running") {
        shiny::showNotification(
          "MCMC is already running. Stop it before starting a new run.",
          type = "warning")
        return()
      }

      dat <- dataset()
      if (is.null(dat)) {
        shiny::showNotification("Load a dataset before running.", type = "warning")
        return()
      }

      if (!requireNamespace("processx", quietly = TRUE)) {
        shiny::showNotification(
          "The 'processx' package is required. Install it with install.packages('processx').",
          type = "error", duration = 10)
        return()
      }

      logDir <- trimws(input$logDir)
      if (!nzchar(logDir)) {
        shiny::showNotification("Specify a log directory.", type = "warning")
        return()
      }
      if (!dir.exists(logDir)) {
        dir.create(logDir, recursive = TRUE, showWarnings = FALSE)
      }

      # Starting tree
      tree <- if (is.function(startTree) || shiny::is.reactive(startTree)) {
        startTree()
      } else {
        startTree
      }
      if (is.null(tree)) {
        taxa <- names(dat)
        tree <- ape::unroot(ape::rtree(length(taxa), tip.label = taxa))
      }

      # Collect MCMC settings
      nRuns   <- as.integer(input$nRuns)
      nChains <- as.integer(input$nChains)
      heat    <- if (nChains > 1L) input$heat else 0.2
      warmup  <- as.integer(input$warmup)
      minEss  <- as.integer(input$minEss)
      maxTimeMins <- input$maxTimeMins
      maxTime <- if (is.na(maxTimeMins) || is.null(maxTimeMins)) NULL
                 else maxTimeMins * 60
      neo     <- .ParseIntList(input$neomorphic)

      cancelFile     <- MkCancelPath(logDir)
      checkpointFile <- file.path(logDir, "checkpoint.rds")
      logBase        <- file.path(logDir, "run")
      logFiles       <- MkLogPaths(logBase, nRuns)

      # Clear any stale signals from a previous run
      for (sig in c(cancelFile,
                    file.path(logDir, "mkp_done.signal"),
                    file.path(logDir, "mkp_error.txt"))) {
        if (file.exists(sig)) file.remove(sig)
      }

      # Build MCMC config and serialise inputs for the detached script
      mcmc <- MkPrimeMCMC(
        nRuns          = nRuns,
        nChains        = nChains,
        heat           = heat,
        warmup         = warmup,
        minEss         = minEss,
        maxTime        = maxTime,
        logFile        = logBase,
        cancelFile     = cancelFile,
        checkpointFile = checkpointFile
      )
      saveRDS(dat,  file.path(logDir, "input_data.rds"))
      saveRDS(tree, file.path(logDir, "input_tree.rds"))
      saveRDS(neo,  file.path(logDir, "input_neo.rds"))
      saveRDS(mcmc, file.path(logDir, "input_mcmc.rds"))

      # Write the launcher script
      scriptFile <- file.path(logDir, "mkp_run.R")
      writeLines(.MkLaunchScript(logDir), scriptFile)

      # Write job.rds before launching (status inferred from signal files, not stored)
      job <- list(
        logDir         = logDir,
        logFiles       = logFiles,
        cancelFile     = cancelFile,
        checkpointFile = checkpointFile,
        scriptFile     = scriptFile,
        nRuns          = nRuns,
        startTime      = Sys.time(),
        pid            = NA_integer_
      )
      jobFile <- file.path(logDir, "job.rds")
      saveRDS(job, jobFile)

      # Launch detached Rscript
      proc <- processx::process$new(
        "Rscript",
        args      = scriptFile,
        supervise = FALSE,
        stdout    = file.path(logDir, "mkp_stdout.txt"),
        stderr    = file.path(logDir, "mkp_stderr.txt")
      )
      job$pid <- proc$get_pid()
      saveRDS(job, jobFile)

      rv$job        <- job
      rv$proc       <- proc
      rv$status     <- "running"
      rv$logSamples <- NULL
      rv$errorMsg   <- NULL

      shiny::showNotification(
        sprintf("MCMC started (PID %d).", job$pid), type = "message")
    })

    # ---- Stop ----------------------------------------------------------------

    shiny::observeEvent(input$stop, {
      if (rv$status != "running") return()
      job <- rv$job
      if (!is.null(job)) {
        file.create(job$cancelFile)
        shiny::showNotification(
          "Cancel signal sent \u2014 MCMC will stop cleanly at the next checkpoint.",
          type = "warning")
      }
    })

    # ---- Reconnect -----------------------------------------------------------
    # Decision tree:
    #   done signal         -> "done"  (no relaunch)
    #   error file          -> relaunch from checkpoint if present, else "error"
    #   cancel signal       -> relaunch from checkpoint if present, else "cancelled"
    #   no signals, PID alive  -> "running" (re-attach monitoring)
    #   no signals, PID dead   -> relaunch from checkpoint if present, else "error"

    shiny::observeEvent(input$reconnect, {
      logDir  <- trimws(input$logDir)
      jobFile <- file.path(logDir, "job.rds")

      if (!nzchar(logDir) || !file.exists(jobFile)) {
        shiny::showNotification(
          "No job.rds found in the specified log directory.", type = "error")
        return()
      }

      job <- tryCatch(readRDS(jobFile), error = function(e) {
        shiny::showNotification(
          paste("Cannot read job.rds:", conditionMessage(e)), type = "error")
        NULL
      })
      if (is.null(job)) return()

      rv$job        <- job
      rv$proc       <- NULL  # process object is lost after session restart
      rv$logSamples <- .ReadAllLogs(job$logFiles)
      rv$errorMsg   <- NULL

      cpFile   <- job$checkpointFile  # NULL on old job.rds -> file.exists(NULL) = FALSE
      haveCp   <- !is.null(cpFile) && file.exists(cpFile)

      if (file.exists(file.path(job$logDir, "mkp_done.signal"))) {
        # ---- Complete --------------------------------------------------------
        rv$status <- "done"
        shiny::showNotification("Reconnected: analysis complete.", type = "message")

      } else if (file.exists(file.path(job$logDir, "mkp_error.txt"))) {
        # ---- R-level error ---------------------------------------------------
        rv$errorMsg <- paste(
          readLines(file.path(job$logDir, "mkp_error.txt"), warn = FALSE),
          collapse = "\n"
        )
        if (haveCp) {
          shiny::showNotification(
            paste0("Prior run ended with an error \u2014 resuming from checkpoint.\n",
                   rv$errorMsg),
            type = "warning", duration = 10)
          .RelaunchFromCheckpoint(job, rv)
        } else {
          rv$status <- "error"
          shiny::showNotification(
            paste0("Prior run ended with an error; no checkpoint to resume from.\n",
                   rv$errorMsg),
            type = "error", duration = 15)
        }

      } else if (file.exists(job$cancelFile)) {
        # ---- Cancelled -------------------------------------------------------
        if (haveCp) {
          shiny::showNotification(
            "Prior run was cancelled \u2014 resuming from checkpoint.",
            type = "message")
          .RelaunchFromCheckpoint(job, rv)
        } else {
          rv$status <- "cancelled"
          shiny::showNotification(
            "Prior run was cancelled. No checkpoint available to resume from.",
            type = "warning")
        }

      } else {
        # ---- No signal files: check PID liveness (M-093) --------------------
        if (.PidIsAlive(job$pid)) {
          rv$status <- "running"
          shiny::showNotification("Reconnected: monitoring running analysis.",
                                  type = "message")
        } else if (haveCp) {
          shiny::showNotification(
            sprintf(paste0("Process PID %d is no longer running",
                           " \u2014 resuming from checkpoint."),
                    as.integer(job$pid)),
            type = "warning")
          .RelaunchFromCheckpoint(job, rv)
        } else {
          rv$status <- "error"
          rv$errorMsg <- sprintf(
            paste0("Process PID %d is no longer running and left no",
                   " checkpoint. The run may have been killed by the OS.",
                   " Partial log data (if any) is available in:\n%s"),
            as.integer(job$pid), job$logDir)
          shiny::showNotification(rv$errorMsg, type = "error", duration = 20)
        }
      }
    })

    # ---- Poll background process ---------------------------------------------

    shiny::observe({
      if (rv$status != "running") return()
      shiny::invalidateLater(5000)

      job <- rv$job
      if (is.null(job)) return()

      doneSig <- file.path(job$logDir, "mkp_done.signal")
      errFile <- file.path(job$logDir, "mkp_error.txt")

      if (file.exists(doneSig)) {
        rv$logSamples <- .ReadAllLogs(job$logFiles)
        rv$status     <- "done"
        shiny::showNotification("MCMC complete!", type = "message")
        return()
      }

      if (file.exists(errFile)) {
        rv$errorMsg <- paste(
          readLines(errFile, warn = FALSE), collapse = "\n")
        rv$status <- "error"
        shiny::showNotification(
          paste("MCMC error:", rv$errorMsg), type = "error", duration = 15)
        return()
      }

      # If cancel signal exists, wait until the process has actually exited
      if (file.exists(job$cancelFile)) {
        procDead <- if (!is.null(rv$proc)) !rv$proc$is_alive() else TRUE
        if (procDead) {
          rv$logSamples <- .ReadAllLogs(job$logFiles)
          rv$status     <- "cancelled"
          return()
        }
      }

      # Still running: refresh log samples for trace plot / ESS table
      rv$logSamples <- .ReadAllLogs(job$logFiles)
    })

    # ---- Status badge --------------------------------------------------------

    output$statusBadge <- shiny::renderUI({
      cls <- switch(rv$status,
        idle      = "bg-secondary",
        running   = "bg-warning text-dark",
        done      = "bg-success",
        cancelled = "bg-info text-dark",
        error     = "bg-danger",
        "bg-secondary"
      )
      label <- switch(rv$status,
        idle      = "Idle",
        running   = "Running\u2026",
        done      = "Complete",
        cancelled = "Cancelled",
        error     = "Error",
        "Unknown"
      )
      shiny::tags$div(
        class = "mt-2 mb-1 text-center",
        shiny::tags$span(class = paste("badge", cls), label)
      )
    })

    # ---- Trace plot ----------------------------------------------------------

    output$tracePlot <- shiny::renderPlot({
      samp <- rv$logSamples
      if (is.null(samp) || nrow(samp) < 2L) return(NULL)

      # Key scalar parameters only (exclude br_ and Sample columns)
      keep <- !grepl("^(Sample|br_)", colnames(samp))
      if (!any(keep)) return(NULL)
      mat <- samp[, keep, drop = FALSE]

      # Cap at last 500 rows and 6 panels
      nr  <- nrow(mat)
      if (nr > 500L) mat <- mat[(nr - 499L):nr, , drop = FALSE]
      nc  <- min(ncol(mat), 6L)
      mat <- mat[, seq_len(nc), drop = FALSE]

      oldpar <- graphics::par(mfrow = c(nc, 1L), mar = c(2, 4, 1, 1),
                              oma = c(0, 0, 0, 0))
      on.exit(graphics::par(oldpar), add = TRUE)
      for (i in seq_len(nc)) {
        graphics::plot(mat[, i], type = "l",
                       ylab = colnames(mat)[i], xlab = "",
                       col = "#2C7BB6", lwd = 0.8)
      }
    })

    # ---- ESS table -----------------------------------------------------------

    output$essTable <- shiny::renderTable({
      samp <- rv$logSamples
      if (is.null(samp) || nrow(samp) < 10L) return(NULL)

      keep <- !grepl("^(Sample|br_)", colnames(samp))
      mat  <- samp[, keep, drop = FALSE]

      ess <- round(.EssMatrix(mat))

      data.frame(
        Parameter = colnames(mat),
        Samples   = nrow(mat),
        ESS       = ess,
        stringsAsFactors = FALSE
      )
    }, digits = 0)

    # ---- Return value --------------------------------------------------------

    list(
      jobFile = shiny::reactive({
        job <- rv$job
        if (is.null(job)) return(NULL)
        file.path(job$logDir, "job.rds")
      }),
      trees = shiny::reactive({
        if (rv$status != "done") return(NULL)
        job <- rv$job
        if (is.null(job)) return(NULL)
        resultFile <- file.path(job$logDir, "result.rds")
        if (!file.exists(resultFile)) return(NULL)
        tryCatch(readRDS(resultFile)$trees, error = function(e) NULL)
      }),
      status = shiny::reactive(rv$status)
    )
  })
}
