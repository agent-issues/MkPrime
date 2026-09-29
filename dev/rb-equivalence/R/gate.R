# Equivalence gate for the MkPrime vs RevBayes oracle harness
# (agent-issues/MkPrime#216). Depends only on `posterior`, so it can be tested
# and applied to recorded rds without building MkPrime.

# Every scalar both samplers must log. A missing column is an error, never a
# silently smaller comparison.
ExpectedScalars <- c("tree_length", "rate_log_sd", "rate_loss", "rate_neo")

if (!exists("%||%", mode = "function")) {
  `%||%` <- function(x, y) if (is.null(x)) y else x # nolint: object_name_linter.
}

# Strictly positive with heavy right tails under their LogNormal/Gamma
# priors: compared on the log scale, where a location shift is visible.
LogScaleScalars <- ExpectedScalars

#' Drop the first `frac` of a trace (data frame rows, vector or tree list).
DropBurnin <- function(x, frac = 0.25) {
  n <- if (is.data.frame(x)) nrow(x) else length(x)
  if (n < 4L) return(x)
  keep <- seq.int(floor(n * frac) + 1L, n)
  if (is.data.frame(x)) x[keep, , drop = FALSE] else x[keep]
}

#' Evenly spaced subsample of `n` elements (or rows).
Subsample <- function(x, n) {
  m <- if (is.data.frame(x)) nrow(x) else length(x)
  idx <- round(seq.int(1, m, length.out = n))
  if (is.data.frame(x)) x[idx, , drop = FALSE] else x[idx]
}

#' Stop unless every expected scalar is present in every run.
AssertScalarsPresent <- function(perRun, source, expected = ExpectedScalars) {
  if (!length(perRun)) stop(source, ": no runs")
  for (i in seq_along(perRun)) {
    missingCols <- setdiff(expected, names(perRun[[i]]))
    if (length(missingCols)) {
      stop(source, " run ", i, " lacks expected scalar(s): ",
           paste(missingCols, collapse = ", "))
    }
  }
  invisible(TRUE)
}

#' Iterations x chains matrix of one parameter, each chain thinned evenly to
#' the shortest chain's length (posterior needs equal-length chains).
DrawsMatrix <- function(chains) {
  nMin <- min(lengths(chains))
  if (nMin < 4L) stop("Too few draws per chain (", nMin, ")")
  do.call(cbind, lapply(chains, Subsample, n = nMin))
}

#' One sampler's own convergence on one parameter: rank-normalised split
#' R-hat and bulk ESS across its runs.
#' @param chains List of numeric vectors, one per run.
SamplerConvergence <- function(chains, targetRhat, targetEss) {
  mat <- DrawsMatrix(chains)
  rhat <- posterior::rhat(mat)
  ess <- posterior::ess_bulk(mat)
  list(rhat = rhat, ess = ess,
       passed = is.finite(rhat) && is.finite(ess) &&
         rhat < targetRhat && ess > targetEss)
}

#' Standardised mean-difference test between two samplers.
#'
#' delta = (mean_mk - mean_rb) / pooled posterior SD; its standard error comes
#' from each sampler's Monte Carlo SE of the mean. Fails when |z| exceeds the
#' Bonferroni-corrected two-sided critical value. `mdd` is the smallest shift
#' (in posterior SD) the test detects with probability `power`: a PASS says
#' nothing about shifts smaller than it, so `mdd > tolerance` means the cell
#' is underpowered, not equivalent.
#'
#' @param mkChains,rbChains Lists of numeric vectors, one per run.
#' @param logScale Compare on the log scale (strictly positive parameters).
#' @param alpha Family-wise false-fail rate over `nTests` parameters.
MeanShiftTest <- function(mkChains, rbChains, logScale = FALSE,
                          alpha = 0.05, nTests = 1L, power = 0.8,
                          tolerance = 0.25) {
  if (logScale) {
    allDraws <- unlist(c(mkChains, rbChains))
    if (any(!is.finite(allDraws) | allDraws <= 0)) {
      stop("Non-positive draw in a log-scale parameter")
    }
    mkChains <- lapply(mkChains, log)
    rbChains <- lapply(rbChains, log)
  }
  mk <- DrawsMatrix(mkChains)
  rb <- DrawsMatrix(rbChains)
  pooledSd <- sqrt((stats::var(as.vector(mk)) + stats::var(as.vector(rb))) / 2)
  if (!is.finite(pooledSd) || pooledSd <= 0) stop("Zero posterior variance")
  seDiff <- sqrt(posterior::mcse_mean(mk)^2 + posterior::mcse_mean(rb)^2)
  delta <- (mean(mk) - mean(rb)) / pooledSd
  seDelta <- seDiff / pooledSd
  z <- delta / seDelta
  zCrit <- stats::qnorm(1 - alpha / (2 * nTests))
  mdd <- (zCrit + stats::qnorm(power)) * seDelta
  list(delta = delta, seDelta = seDelta, z = z, zCrit = zCrit, mdd = mdd,
       shifted = !is.finite(z) || abs(z) > zCrit,
       underpowered = !is.finite(mdd) || mdd > tolerance)
}

#' Full gate for one parameter.
#'
#' In order: each sampler's own convergence (on all its post-burn-in draws),
#' cross-sampler rank-normalised R-hat (chains thinned to a common length so
#' neither sampler dominates), and the mean-shift test. `verdict` is "FAIL" if
#' any check fails, else "UNDERPOWERED" if the mean-shift test could miss a
#' `tolerance`-SD shift, else "PASS".
#'
#' @param mkChains,rbChains Lists of numeric vectors, one per run.
#' @return One-row data frame.
GateParam <- function(mkChains, rbChains, targetRhat, targetEss,
                      logScale = FALSE, alpha = 0.05, nTests = 1L,
                      tolerance = 0.25) {
  mkConv <- SamplerConvergence(mkChains, targetRhat, targetEss)
  rbConv <- SamplerConvergence(rbChains, targetRhat, targetEss)
  nCommon <- min(lengths(c(mkChains, rbChains)))
  crossRhat <- posterior::rhat(DrawsMatrix(
    lapply(c(mkChains, rbChains), Subsample, n = nCommon)
  ))
  crossOk <- is.finite(crossRhat) && crossRhat < targetRhat
  shift <- MeanShiftTest(mkChains, rbChains, logScale = logScale,
                         alpha = alpha, nTests = nTests, tolerance = tolerance)

  reasons <- c(
    if (!mkConv$passed) "mkprime_unconverged",
    if (!rbConv$passed) "rb_unconverged",
    if (!crossOk) "cross_rhat",
    if (shift$shifted) "mean_shift"
  )
  verdict <- if (length(reasons)) "FAIL" else if (shift$underpowered) {
    "UNDERPOWERED"
  } else "PASS"

  data.frame(
    mkp_rhat = mkConv$rhat, mkp_ess = mkConv$ess,
    rb_rhat = rbConv$rhat, rb_ess = rbConv$ess,
    cross_rhat = crossRhat,
    scale = if (logScale) "log" else "raw",
    delta_sd = shift$delta, z = shift$z, z_crit = shift$zCrit,
    mdd_sd = shift$mdd, tolerance_sd = tolerance,
    verdict = verdict,
    reasons = paste(reasons, collapse = ";"),
    stringsAsFactors = FALSE
  )
}

#' Stop unless the MkPrime run met the targets the gate applies: it must have
#' stopped on convergence, not on time or iterations, and been asked for at
#' least the gate's R-hat and ESS.
AssertMkRunTargets <- function(posterior, targetRhat, targetEss) {
  reason <- posterior$stop_reason %||% NA_character_
  if (!identical(reason, "converged")) {
    stop("MkPrime run stopped on '", reason, "', not on convergence; ",
         "rerun with a longer --max-time.")
  }
  runRhat <- posterior$mcmc$maxRhat %||% Inf
  runEss <- posterior$mcmc$minEss %||% 0
  if (runRhat > targetRhat || runEss < targetEss) {
    stop(sprintf(paste0(
      "MkPrime run targeted maxRhat = %g, minEss = %g, looser than the ",
      "gate's %g / %g; rerun run_mkprime.R with --rhat/--ess at least as strict."
    ), runRhat, runEss, targetRhat, targetEss))
  }
  invisible(TRUE)
}

#' Model settings both samplers must share, read from a MkPrimeModel object.
MkModelSpec <- function(model) {
  if (!inherits(model, "MkPrimeModel")) stop("MkModelSpec needs a MkPrimeModel")
  list(
    coding = model$coding,
    nCat = as.integer(model$nCat),
    treeLengthShape = as.numeric(model$treeLengthShape),
    treeLengthRate = as.numeric(model$treeLengthRate),
    rateLogSdShape = as.numeric(model$rateLogSdShape),
    rateLogSdRate = as.numeric(model$rateLogSdRate),
    rateLossMeanlog = as.numeric(model$rateLossMeanlog),
    rateLossSdlog = as.numeric(model$rateLossSdlog),
    rateNeoMeanlog = as.numeric(model$rateNeoMeanlog),
    rateNeoSdlog = as.numeric(model$rateNeoSdlog)
  )
}

#' The same settings, parsed from the rendered RevBayes model script that RB
#' actually ran. Stops if any expected declaration is absent or has changed
#' form, so a template edit fails loudly instead of drifting.
#' @param rev Character: the model script's lines or text.
RevModelSpec <- function(rev) {
  rev <- paste(rev, collapse = "\n")
  num <- "([-+]?[0-9.]+(?:[eE][-+]?[0-9]+)?)"
  Grab <- function(pattern, what, n = 1L) {
    m <- regmatches(rev, regexec(pattern, rev, perl = TRUE))[[1]]
    if (length(m) < n + 1L) stop("RevModelSpec: cannot find ", what)
    as.numeric(m[seq_len(n) + 1L])
  }
  shape <- Grab(paste0("(?m)^\\s*gamma_shape\\s*<-\\s*", num, "\\s*$"),
                "gamma_shape")
  steps <- Grab(paste0("(?m)^\\s*exp_steps\\s*<-\\s*", num, "\\s*$"),
                "exp_steps")
  if (!grepl(paste0("tree_length\\s*~\\s*dnGamma\\(\\s*shape\\s*=\\s*gamma_shape",
                    "\\s*,\\s*rate\\s*=\\s*gamma_shape\\s*/\\s*exp_steps\\s*\\)"),
             rev, perl = TRUE)) {
    stop("RevModelSpec: tree_length prior is not ",
         "dnGamma(shape = gamma_shape, rate = gamma_shape / exp_steps)")
  }
  logSd <- Grab(paste0("rate_log_sd\\s*~\\s*dnGamma\\(\\s*", num, "\\s*,\\s*",
                       num, "\\s*\\)"), "rate_log_sd prior", 2L)
  LognormalPrior <- function(par) {
    Grab(paste0(par, "\\s*~\\s*dnLognormal\\(\\s*mean\\s*=\\s*", num,
                "\\s*,\\s*sd\\s*=\\s*", num, "\\s*\\)"),
         paste(par, "prior"), 2L)
  }
  loss <- LognormalPrior("rate_loss")
  neo <- LognormalPrior("rate_neo")
  nCat <- Grab("fnDiscretizeDistribution\\(.*,\\s*([0-9]+)\\s*\\)",
               "fnDiscretizeDistribution category count")
  coding <- unique(regmatches(
    rev, gregexpr('coding\\s*=\\s*"[^"]*"', rev, perl = TRUE)
  )[[1]])
  coding <- unique(sub('.*"([^"]*)"', "\\1", coding))
  if (length(coding) != 1L) {
    stop("RevModelSpec: expected one coding across all dnPhyloCTMC, found ",
         if (length(coding)) paste(coding, collapse = ", ") else "none")
  }
  list(
    coding = coding,
    nCat = as.integer(nCat),
    treeLengthShape = shape,
    treeLengthRate = shape / steps,
    rateLogSdShape = logSd[[1]],
    rateLogSdRate = logSd[[2]],
    rateLossMeanlog = loss[[1]],
    rateLossSdlog = loss[[2]],
    rateNeoMeanlog = neo[[1]],
    rateNeoSdlog = neo[[2]]
  )
}

#' Stop unless both rds carry the same model settings and burn-in fraction.
AssertSameSetup <- function(mk, rb) {
  if (is.null(mk$model_spec) || is.null(rb$model_spec)) {
    stop("model_spec missing from mkprime and/or rb rds; rerun ",
         "run_mkprime.R / post_rb.R (agent-issues/MkPrime#216).")
  }
  same <- all.equal(mk$model_spec, rb$model_spec)
  if (!isTRUE(same)) {
    stop("Model settings differ between samplers: ",
         paste(same, collapse = "; "))
  }
  if (is.null(mk$burnin_frac) || is.null(rb$burnin_frac) ||
      !identical(mk$burnin_frac, rb$burnin_frac)) {
    stop("burnin_frac missing or unequal (mkprime: ",
         format(mk$burnin_frac %||% "NULL"), ", rb: ",
         format(rb$burnin_frac %||% "NULL"), ").")
  }
  invisible(TRUE)
}

#' Row-bind onto an existing CSV, replacing rows with the same key and
#' filling columns either side lacks with NA.
AppendCsv <- function(df, path, keyCols = c("pid", "model", "param", "source")) {
  if (file.exists(path)) {
    existing <- utils::read.csv(path, stringsAsFactors = FALSE)
    keyCols <- intersect(keyCols, intersect(names(df), names(existing)))
    existing <- existing[!do.call(paste, existing[, keyCols, drop = FALSE]) %in%
                         do.call(paste, df[, keyCols, drop = FALSE]), ,
                         drop = FALSE]
    for (col in setdiff(names(df), names(existing))) existing[[col]] <- rep(NA, nrow(existing))
    for (col in setdiff(names(existing), names(df))) df[[col]] <- rep(NA, nrow(df))
    df <- rbind(existing, df[, names(existing)])
  }
  utils::write.csv(df, path, row.names = FALSE)
  invisible(path)
}
