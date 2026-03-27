# Shiny app launcher for MkPrime

#' Launch the MkPrime graphical interface
#'
#' Opens a Shiny web application for interactive Bayesian phylogenetic
#' analysis under the Mk' model. The app provides data loading, MCMC
#' configuration, live progress display, and result visualization
#' (traces, parameter summaries, consensus tree).
#'
#' @details
#' The MCMC runs as a detached `Rscript` process via \pkg{processx},
#' so the analysis survives browser close and session timeout. Progress
#' is polled from TSV log files written to the chosen log directory.
#' Use the "Reconnect" button to re-attach to a running job after
#' restarting the app.
#'
#' ## Required packages
#'
#' The following packages must be installed (listed in `Suggests`):
#' - **shiny** --- Shiny framework
#' - **bslib** --- Bootstrap 5 UI components
#' - **processx** --- detached process management
#'
#' @return Opens a Shiny app in the default browser; does not return a
#'   value.
#' @seealso [RunMkPrime()], [MkPrimeMCMC()]
#' @export
EasyMkPrime <- function() {
  needed <- c("shiny", "bslib", "processx")
  missing <- needed[!vapply(needed, requireNamespace,
                            logical(1L), quietly = TRUE)]
  if (length(missing)) {
    stop("EasyMkPrime() requires additional packages: ",
         paste(missing, collapse = ", "), ".\n",
         "Install with: install.packages(",
         paste0('"', missing, '"', collapse = ", "), ")",
         call. = FALSE)
  }
  shiny::runApp(system.file("MkPrime", package = "MkPrime"))
}
