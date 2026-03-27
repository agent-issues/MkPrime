# MkPrime Shiny App — Bayesian Mk' Phylogenetic Inference
#
# Launched by MkPrime::EasyMkPrime(). Runs MCMC in a background process
# with live progress display via PNG polling.

library(shiny)
library(bslib)
library(MkPrime)

# ---------------------------------------------------------------------------
# Helpers (local to app)
# ---------------------------------------------------------------------------

parse_integer_list <- function(text) {
  text <- trimws(text)
  if (!nzchar(text)) return(integer(0))
  vals <- as.integer(strsplit(text, "[,;\\s]+")[[1]])
  vals[!is.na(vals)]
}

read_tree_safe <- function(path, name) {
  ext <- tolower(tools::file_ext(name))
  tr <- if (ext %in% c("nex", "nxs", "nexus")) {
    ape::read.nexus(path)
  } else {
    ape::read.tree(path)
  }
  if (inherits(tr, "multiPhylo")) tr <- tr[[1]]
  if (ape::is.rooted(tr)) tr <- ape::unroot(tr)
  tr
}

format_elapsed <- function(seconds) {
  if (seconds < 60) sprintf("%.0fs", seconds)
  else if (seconds < 3600) sprintf("%.1f min", seconds / 60)
  else sprintf("%.1f h", seconds / 3600)
}

parse_progress_json <- function(path) {
  jt <- tryCatch(readLines(path, warn = FALSE)[1], error = function(e) "")
  if (!nzchar(jt)) return(NULL)
  list(
    iter    = as.integer(sub('.*"iter":([0-9]+).*', '\\1', jt)),
    nIter   = as.integer(sub('.*"nIter":([0-9]+).*', '\\1', jt)),
    warmup  = as.integer(sub('.*"warmup":([0-9]+).*', '\\1', jt)),
    elapsed = as.numeric(sub('.*"elapsed":([0-9.]+).*', '\\1', jt)),
    acc     = as.numeric(sub('.*"recent_acceptance":([0-9.]+).*', '\\1', jt))
  )
}

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------

ui <- page_sidebar(
  title = "MkPrime",
  theme = bs_theme(version = 5),

  sidebar = sidebar(
    width = 310,
    accordion(
      id = "settings",
      open = "Data",

      # ---- Data ----
      accordion_panel(
        "Data",
        fileInput("dataFile", "Character matrix",
                  accept = c(".nex", ".nxs", ".tnt", ".txt")),
        textInput("neomorphic", "Neomorphic characters",
                  placeholder = "e.g., 1,3,5"),
        uiOutput("dataInfo")
      ),

      # ---- Starting tree ----
      accordion_panel(
        "Starting tree",
        radioButtons("treeSource", NULL,
                     choices = c("Random" = "random", "Upload" = "upload"),
                     inline = TRUE),
        conditionalPanel("input.treeSource == 'upload'",
          fileInput("treeFile", "Tree file",
                    accept = c(".nex", ".nxs", ".nwk", ".tre"))
        )
      ),

      # ---- MCMC ----
      accordion_panel(
        "MCMC",
        numericInput("nIter", "Iterations", 10000, min = 100, step = 1000),
        numericInput("warmup", "Warmup", 5000, min = 0, step = 500),
        numericInput("thin", "Thinning", 10, min = 1),
        numericInput("nRuns", "Independent runs", 2, min = 1, max = 10),
        numericInput("nChains", "Chains per run", 1, min = 1, max = 8),
        conditionalPanel("input.nChains > 1",
          sliderInput("heat", "Heat parameter", 0.05, 0.95, 0.2, step = 0.05)
        ),
        input_switch("fixTopology", "Fix topology", value = FALSE)
      )
    ),
    hr(),
    actionButton("run", "Run MCMC", class = "btn-primary w-100"),
    actionButton("stop", "Stop", class = "btn-outline-danger w-100 mt-2"),
    uiOutput("statusBadge")
  ),

  # ---- Main area ----
  navset_card_underline(
    id = "mainTabs",
    full_screen = TRUE,

    nav_panel("Progress",
      uiOutput("progressBar"),
      imageOutput("progressPlot", height = "auto")
    ),

    nav_panel("Traces",
      plotOutput("tracePlot", height = "600px")
    ),

    nav_panel("Summary",
      tableOutput("summaryTable")
    ),

    nav_panel("Consensus",
      plotOutput("consensusPlot", height = "600px")
    )
  )
)

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------

server <- function(input, output, session) {

  rv <- reactiveValues(
    data         = NULL,   # phyDat
    tree         = NULL,   # starting phylo
    bg           = NULL,   # callr::r_bg process handle
    progress_dir = NULL,
    result_file  = NULL,
    result       = NULL,   # MkPosterior
    status       = "idle"  # idle | running | done | error
  )

  # ---- Data loading ----

  observeEvent(input$dataFile, {
    req(input$dataFile)
    tryCatch({
      rv$data <- TreeTools::ReadAsPhyDat(input$dataFile$datapath)
      rv$result <- NULL
      rv$status <- "idle"
      showNotification(
        sprintf("Loaded: %d taxa, %d characters",
                length(rv$data), attr(rv$data, "nr")),
        type = "message"
      )
      # Auto-detect neomorphic characters
      neo <- AutoDetectNeomorphic(rv$data)
      if (length(neo)) {
        updateTextInput(session, "neomorphic",
                        value = paste(neo, collapse = ", "))
        showNotification(
          sprintf("Auto-detected %d neomorphic character%s (binary {0,1}).",
                  length(neo), if (length(neo) == 1L) "" else "s"),
          type = "message", duration = 6
        )
      } else {
        updateTextInput(session, "neomorphic", value = "")
      }
    }, error = function(e) {
      showNotification(paste("Error loading data:", conditionMessage(e)),
                       type = "error", duration = 10)
    })
  })

  output$dataInfo <- renderUI({
    req(rv$data)
    neo <- parse_integer_list(input$neomorphic)
    mkd <- tryCatch(
      MkPrimeData(rv$data, neomorphic = neo),
      warning = function(w) invokeRestart("muffleWarning"),
      error = function(e) NULL
    )
    if (is.null(mkd)) return(NULL)
    types <- table(mkd$type)
    tags$p(
      class = "text-muted small mt-1",
      sprintf("%d characters: %s", mkd$nChar,
              paste(sprintf("%d %s", types, names(types)), collapse = ", "))
    )
  })

  # ---- Run MCMC ----

  observeEvent(input$run, {
    req(rv$data, rv$status != "running")

    taxa <- names(rv$data)
    nTaxa <- length(taxa)

    # Starting tree
    if (input$treeSource == "upload" && !is.null(input$treeFile)) {
      tryCatch({
        rv$tree <- read_tree_safe(input$treeFile$datapath,
                                  input$treeFile$name)
      }, error = function(e) {
        showNotification(paste("Tree error:", conditionMessage(e)),
                         type = "error")
        return()
      })
    } else {
      rv$tree <- ape::unroot(ape::rtree(nTaxa, tip.label = taxa))
    }

    # Set up progress directory
    rv$progress_dir <- tempfile("mkp_shiny_")
    dir.create(rv$progress_dir)
    rv$result_file <- file.path(rv$progress_dir, "result.rds")
    rv$result <- NULL
    rv$status <- "running"

    # Snapshot inputs for background process
    mcmc_args <- list(
      nIter   = as.integer(input$nIter),
      warmup  = as.integer(input$warmup),
      thin    = as.integer(input$thin),
      nRuns   = as.integer(input$nRuns),
      nChains = as.integer(input$nChains),
      heat    = if (input$nChains > 1) input$heat else 0.2,
      plot_every = max(50L, as.integer(input$nIter / 40))
    )
    data_snap <- rv$data
    tree_snap <- rv$tree
    neo_snap  <- parse_integer_list(input$neomorphic)
    fix_snap  <- isTRUE(input$fixTopology)
    prog_dir  <- rv$progress_dir
    res_file  <- rv$result_file

    # Launch background R process
    rv$bg <- callr::r_bg(
      function(data, tree, neomorphic, mcmc_args, fix_topology,
               progress_dir, result_file) {
        library(MkPrime)
        mcmc_args$progress_fn <- mkp_png_progress(progress_dir)
        mcmc <- do.call(MkPrimeMCMC, mcmc_args)
        result <- RunMkPrime(data, tree, neomorphic = neomorphic,
                             mcmc = mcmc, fix_topology = fix_topology)
        saveRDS(result, result_file)
        "done"
      },
      args = list(
        data = data_snap, tree = tree_snap, neomorphic = neo_snap,
        mcmc_args = mcmc_args, fix_topology = fix_snap,
        progress_dir = prog_dir, result_file = res_file
      ),
      libpath = .libPaths(),
      supervise = TRUE
    )
  })

  # ---- Stop ----

  observeEvent(input$stop, {
    req(rv$status == "running", rv$bg)
    if (rv$bg$is_alive()) {
      rv$bg$kill()
      rv$status <- "idle"
      showNotification("MCMC stopped by user.", type = "warning")
    }
  })

  # ---- Poll background process ----

  observe({
    req(rv$status == "running", rv$bg)
    invalidateLater(1500)

    if (!rv$bg$is_alive()) {
      tryCatch({
        rv$bg$get_result()
        if (file.exists(rv$result_file)) {
          rv$result <- readRDS(rv$result_file)
          rv$status <- "done"
          showNotification("MCMC complete!", type = "message")
          nav_select("mainTabs", "Traces")
        } else {
          rv$status <- "error"
          showNotification("MCMC finished but produced no output.",
                           type = "error")
        }
      }, error = function(e) {
        rv$status <- "error"
        showNotification(
          paste("MCMC error:", conditionMessage(e)),
          type = "error", duration = 15
        )
      })
    }
  })

  # ---- Status badge ----

  output$statusBadge <- renderUI({
    cls <- switch(rv$status,
      idle    = if (is.null(rv$data)) "bg-secondary" else "bg-info",
      running = "bg-warning",
      done    = "bg-success",
      error   = "bg-danger"
    )
    label <- switch(rv$status,
      idle    = if (is.null(rv$data)) "Load data to begin" else "Ready",
      running = "Running\u2026",
      done    = "Complete",
      error   = "Error"
    )
    tags$div(class = "mt-2 text-center",
             tags$span(class = paste("badge", cls), label))
  })

  # ---- Progress display ----

  output$progressBar <- renderUI({
    if (rv$status == "running") invalidateLater(1500)
    req(rv$progress_dir)

    json_path <- file.path(rv$progress_dir, "mkp_progress.json")
    if (!file.exists(json_path)) {
      if (rv$status == "running") {
        return(tags$p(class = "text-muted", "Starting MCMC\u2026"))
      }
      if (rv$status == "done") {
        return(tags$p(class = "text-success fw-bold", "Analysis complete."))
      }
      return(NULL)
    }

    info <- parse_progress_json(json_path)
    if (is.null(info)) return(NULL)

    pct <- min(100, round(100 * info$iter / info$nIter))
    phase <- if (info$iter <= info$warmup) "warmup" else "sampling"
    bar_class <- if (rv$status == "done") "bg-success"
                 else "progress-bar-animated"

    tagList(
      tags$div(class = "progress mb-2", style = "height: 22px;",
        tags$div(
          class = paste("progress-bar progress-bar-striped", bar_class),
          role = "progressbar",
          style = paste0("width: ", pct, "%;"),
          paste0(pct, "%")
        )
      ),
      tags$p(class = "text-muted small",
        sprintf("Iteration %s / %s (%s) \u00b7 Acceptance: %.1f%% \u00b7 %s",
                format(info$iter, big.mark = ","),
                format(info$nIter, big.mark = ","),
                phase, info$acc * 100,
                format_elapsed(info$elapsed))
      )
    )
  })

  output$progressPlot <- renderImage({
    if (rv$status == "running") invalidateLater(2500)
    req(rv$progress_dir)
    png_path <- file.path(rv$progress_dir, "mkp_progress.png")
    req(file.exists(png_path))
    list(src = png_path, width = "100%", alt = "MCMC Progress",
         contentType = "image/png")
  }, deleteFile = FALSE)

  # ---- Results: Traces ----

  output$tracePlot <- renderPlot({
    req(rv$result)
    plot(rv$result)
  })

  # ---- Results: Summary ----

  output$summaryTable <- renderTable({
    req(rv$result)
    summary(rv$result)
  }, digits = 4)

  # ---- Results: Consensus ----

  output$consensusPlot <- renderPlot({
    req(rv$result)
    trees <- rv$result$trees
    req(length(trees) >= 1)
    if (length(trees) < 2) {
      ape::plot.phylo(trees[[1]])
      title("Single sampled tree")
    } else {
      cons <- ape::consensus(trees, p = 0.5)
      ape::plot.phylo(cons)
      title("Majority-rule consensus")
    }
  })

  # ---- Cleanup ----

  onSessionEnded(function() {
    if (!is.null(isolate(rv$bg)) && isolate(rv$bg)$is_alive()) {
      isolate(rv$bg)$kill()
    }
  })
}

shinyApp(ui, server)
