# MkPrime Shiny App — Bayesian Mk' Phylogenetic Inference
#
# Launched by MkPrime::EasyMkPrime(). MCMC runs as a detached Rscript
# process (survives browser close / session timeout) via MkBayesianServer.
# Progress is polled from disk log files every 5 s.

library(shiny)
library(bslib)
library(MkPrime)

# ---------------------------------------------------------------------------
# Helpers (local to app)
# ---------------------------------------------------------------------------

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

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------

ui <- page_sidebar(
  title = "MkPrime",
  theme = bs_theme(version = 5),

  sidebar = sidebar(
    width = 330,

    # ---- Data + Starting tree (app-owned) ------------------------------------
    accordion(
      id   = "data_settings",
      open = "Data",

      accordion_panel(
        "Data",
        fileInput("dataFile", "Character matrix",
                  accept = c(".nex", ".nxs", ".tnt", ".txt")),
        uiOutput("dataInfo")
      ),

      accordion_panel(
        "Starting tree",
        radioButtons("treeSource", NULL,
                     choices = c("Random" = "random", "Upload" = "upload"),
                     inline  = TRUE),
        conditionalPanel("input.treeSource == 'upload'",
          fileInput("treeFile", "Tree file",
                    accept = c(".nex", ".nxs", ".nwk", ".tre"))
        )
      )
    ),

    hr(),

    # ---- MCMC module (config, characters, buttons, progress) -----------------
    MkBayesianUi("bayes")
  ),

  # ---- Main area (results only; live progress lives in the sidebar) ----------
  navset_card_underline(
    id = "mainTabs",
    full_screen = TRUE,

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
    data = NULL,  # phyDat
    tree = NULL   # uploaded starting phylo (NULL → random)
  )

  # ---- Data loading ----------------------------------------------------------

  observeEvent(input$dataFile, {
    req(input$dataFile)
    tryCatch({
      rv$data <- TreeTools::ReadAsPhyDat(input$dataFile$datapath)
      showNotification(
        sprintf("Loaded: %d taxa, %d characters",
                length(rv$data), attr(rv$data, "nr")),
        type = "message"
      )
    }, error = function(e) {
      showNotification(paste("Error loading data:", conditionMessage(e)),
                       type = "error", duration = 10)
    })
  })

  output$dataInfo <- renderUI({
    req(rv$data)
    tags$p(class = "text-muted small mt-1",
      sprintf("%d taxa \u00b7 %d characters",
              length(rv$data), attr(rv$data, "nr"))
    )
  })

  # ---- Starting tree ---------------------------------------------------------

  observeEvent(input$treeFile, {
    req(input$treeFile)
    tryCatch({
      rv$tree <- read_tree_safe(input$treeFile$datapath,
                                input$treeFile$name)
    }, error = function(e) {
      showNotification(paste("Tree error:", conditionMessage(e)),
                       type = "error")
    })
  })

  observeEvent(input$treeSource, {
    if (input$treeSource == "random") rv$tree <- NULL
  })

  # ---- Bayesian module -------------------------------------------------------

  module <- MkBayesianServer(
    "bayes",
    dataset   = reactive(rv$data),
    startTree = reactive({
      if (isTRUE(input$treeSource == "upload")) rv$tree else NULL
    })
  )

  # ---- Read full posterior from disk when analysis is complete ---------------

  posterior <- reactive({
    req(module$status() == "done")
    jf <- module$jobFile()
    req(!is.null(jf))
    resultFile <- file.path(dirname(jf), "result.rds")
    if (!file.exists(resultFile)) return(NULL)
    result <- tryCatch(readRDS(resultFile), error = function(e) NULL)
    if (is.null(result)) return(NULL)
    # Streaming mode: RunMkPrime writes samples to TSV log files and returns
    # an MkPosterior with an empty $samples matrix.  Reload from disk so that
    # plot() and summary() work normally.
    if (!is.null(result$logFile) && nrow(result$samples) == 0L) {
      samp <- tryCatch(ReadMkLog(result$logFile), error = function(e) NULL)
      if (!is.null(samp)) result$samples <- samp
    }
    result
  })

  # ---- Traces tab ------------------------------------------------------------

  output$tracePlot <- renderPlot({
    req(posterior())
    plot(posterior())
  })

  # ---- Summary tab -----------------------------------------------------------

  output$summaryTable <- renderTable({
    req(posterior())
    summary(posterior())
  }, digits = 4)

  # ---- Consensus tab ---------------------------------------------------------

  output$consensusPlot <- renderPlot({
    trees <- module$trees()
    req(!is.null(trees), length(trees) >= 1)
    if (length(trees) < 2L) {
      ape::plot.phylo(trees[[1L]])
      title("Single sampled tree")
    } else {
      cons <- ape::consensus(trees, p = 0.5)
      ape::plot.phylo(cons)
      title("Majority-rule consensus")
    }
  })
}

shinyApp(ui, server)
