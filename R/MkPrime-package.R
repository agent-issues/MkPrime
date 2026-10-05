#' MkPrime: Bayesian Phylogenetic Inference Under the Mk' Model
#'
#' @description
#' Bayesian inference of phylogenetic trees from discrete morphological
#' characters under the Mk' (Mk-prime) model, which infers the true number
#' of character states, including states not observed in the data. Three
#' character types can be combined in a single analysis: asymmetric binary
#' (neomorphic), Mk' with inferred state count (transformational), and
#' standard Mk with a known state space. MCMC sampling supports parallel
#' tempering, among-character rate variation, and ascertainment bias
#' correction, with samples streamed to disk so long runs can be resumed
#' or inspected before they finish.
#'
#' @section Data and model setup:
#' - [MkPrimeData()] classifies characters and builds the data object
#'   `RunMkPrime()` consumes.
#' - [AutoDetectNeomorphic()] flags likely neomorphic (asymmetric binary)
#'   characters.
#' - [MkPrimeModel()] specifies priors and among-character rate variation.
#' - [MkPrimeEmpiricalPrior()] builds an empirical-Bayes prior on k' from a
#'   corpus of comparable characters.
#' - [MkPrimeMCMC()] configures chain counts, tempering, stopping rules,
#'   and output file paths.
#'
#' @section Running MCMC:
#' - [RunMkPrime()] is the main entry point.
#' - [ResumeMkPrime()] continues a run from its checkpoint file.
#' - [MkPrimeRecover()] loads whatever samples an interrupted run already
#'   flushed to disk.
#'
#' @section Monitoring a run in progress:
#' - [ReadMkLog()] reads the streamed parameter log, including while the
#'   run that wrote it is still in progress.
#' - [MkLogPaths()] resolves per-run log file paths for a multi-run analysis.
#' - [MkpWatchLog()] and [MkpTracePlot()] display live trace plots.
#' - [MkpPngProgress()] writes progress plots to disk instead of a device.
#' - [MkCancelPath()] builds the path to a job's cancel-signal file.
#' - [MkPrimeVerbosity()] controls MCMC console output.
#'
#' @section Convergence and summarising results:
#' - [ConvergenceDiagnostics()] computes ESS and rank-normalized R-hat.
#' - [TreeESS()] computes effective sample size for the sampled tree
#'   topologies.
#' - [AutoBurnin()] and [SetBurnin()] detect and set the burn-in period.
#' - [mkp_stepping_stone()] estimates the marginal likelihood for model
#'   comparison.
#' - [MkpLogLikelihood()] evaluates the model log-likelihood for a given
#'   tree and data outside of MCMC.
#'
#' @section Shiny GUI:
#' - [EasyMkPrime()] launches the standalone Shiny app.
#' - [MkBayesianUi()] and [MkBayesianServer()] provide a Shiny module for
#'   the TreeSearch GUI.
#'
#' @seealso
#' `vignette("hyoliths", package = "MkPrime")` for a full worked analysis,
#' including saving, resuming, and monitoring long runs.
#'
#' `vignette("het-details", package = "MkPrime")` for among-character rate
#' variation and Q-matrix heterogeneity.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @useDynLib MkPrime, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom Rdpack reprompt
#' @importFrom graphics abline lines mtext par plot.new title
#' @importFrom stats cor dbeta dexp dgamma dlnorm dnorm fft median na.omit
#' @importFrom stats nextn plogis qlogis qnorm quantile rbeta rnorm runif sd
#' @importFrom stats var
#' @importFrom utils modifyList tail
## usethis namespace: end
NULL

# Package-level environment for interrupt recovery.
# Stores paths to temp log files from the most-recently interrupted run,
# so MkPrimeRecover() can retrieve partial results.
.mkp_env <- new.env(parent = emptyenv())
.mkp_env$recovery <- NULL        # NULL or list(logFiles, paramNames, ...)
.mkp_env$active_temp_logs <- NULL # paths of temp logs from current run
