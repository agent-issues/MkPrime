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
