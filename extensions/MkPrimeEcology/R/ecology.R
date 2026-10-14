#' Ecology-biased equilibrium model
#'
#' Builds the ecology-biased equilibrium (EBE) model as an MkPrime engine
#' extension, for `MkPrime::MkPrimeModel(extensions = ...)`.
#'
#' The extension takes over MkPrime's neomorphic (absent / present)
#' characters. On an edge in ecology `s`, the equilibrium frequency of
#' "present" in character `c` is tilted by `logit(pi1) += s_z log(phi)`:
#' towards present when `z[c, s] = 1`, towards absent when `z[c, s] = 2`,
#' and not at all when `z[c, s] = 0` or `s` is the reference ecology. The
#' rate of change is unchanged. Edge ecology is a mixture over ecology
#' states, weighted by the JC(K) reconstruction of tip ecology at the
#' edge's ends. With every `z = 0` the model is MkPrime's neomorphic MkN.
#'
#' Parameters: `phi` (one, or one per non-reference ecology), `pi0` (the
#' prior probability that a character is unaffected by an ecology), and
#' `theta` (per non-reference ecology, the probability that an affected
#' character is encouraged rather than discouraged), all logged with
#' `eco_` names; and the influence categories `z`, logged as
#' `eco_z_<character>_<ecology>`. `z` is updated by a Gibbs move on one
#' character at a time.
#'
#' The tilt makes the process non-stationary, so supply and sample a rooted
#' tree. Use the same `neomorphic` characters in `MkPrime::RunMkPrime()`.
#' Variable and no ascertainment coding are supported; informative coding
#' is not.
#'
#' @param data A `phyDat` object or `MkPrime::MkPrimeData()` result.
#' @param ecology Ecology of each tip: a vector or factor named by tip label.
#'   `NA` marks an unknown ecology.
#' @param neomorphic Indices of neomorphic characters, as for
#'   `MkPrime::RunMkPrime()`. Ignored when `data` is an `MkPrimeData`.
#' @param refEcology The reference ecology, in which no character is
#'   tilted. Default: the first level.
#' @param magnitudeMode `"global"` for one `phi`, or `"per_ecology"` for one
#'   per non-reference ecology.
#' @param pi0Prior,thetaPrior Beta shape parameters of the priors on `pi0`
#'   and each `theta`.
#' @param sigmaPhi Standard deviation of the log-normal prior on `phi`.
#' @param zWeight Relative frequency of the Gibbs move on `z`.
#' @param init Named list of initial values: `phi`, `pi0`, `theta`.
#' @return An `MkPrimeExtension` object.
#' @export
EcologyExtension <- function(data, ecology, neomorphic = integer(0),
                             refEcology = NULL,
                             magnitudeMode = c("global", "per_ecology"),
                             pi0Prior = c(360, 120), thetaPrior = c(2, 2),
                             sigmaPhi = 1.5, zWeight = 5, init = list()) {
  mkd <- if (inherits(data, "MkPrimeData")) {
    data
  } else {
    MkPrime::MkPrimeData(data, neomorphic = neomorphic)
  }
  magnitudeMode <- match.arg(magnitudeMode)
  tips <- mkd$taxon_names
  if (is.null(names(ecology))) {
    stop("`ecology` must be named by tip label.")
  }
  absent <- setdiff(tips, names(ecology))
  if (length(absent)) {
    stop("`ecology` has no entry for tip(s) ", paste(absent, collapse = ", "),
         "; use NA for an unknown ecology.")
  }
  eco <- factor(ecology[tips])
  ecoLevels <- levels(eco)
  kEco <- length(ecoLevels)
  if (kEco < 2L) stop("`ecology` needs at least two ecologies.")
  if (is.null(refEcology)) refEcology <- ecoLevels[1]
  refIdx <- match(as.character(refEcology), ecoLevels)
  if (is.na(refIdx)) stop("`refEcology` is not an ecology in `ecology`.")
  for (nm in c("pi0Prior", "thetaPrior")) {
    pr <- get(nm)
    if (length(pr) != 2L || any(!is.finite(pr)) || any(pr <= 0)) {
      stop("`", nm, "` must be two positive shape parameters.")
    }
  }
  if (length(sigmaPhi) != 1L || !(sigmaPhi > 0)) {
    stop("`sigmaPhi` must be positive.")
  }

  neoParts <- Filter(function(p) identical(p$type, "neomorphic"),
                     mkd$partitions)
  neoChars <- as.integer(unlist(lapply(neoParts, `[[`, "char_indices")))
  nNeo <- length(neoChars)
  if (!nNeo) stop("`data` has no neomorphic characters.")
  nonRef <- ecoLevels[-refIdx]
  label <- gsub("[^A-Za-z0-9]+", "_", nonRef)
  phiNames <- if (magnitudeMode == "global") "eco_phi" else
    paste0("eco_phi_", label)
  paramNames <- c(phiNames, "eco_pi0", paste0("eco_theta_", label))
  phi0 <- init$phi %||% 2
  theta0 <- init$theta %||% 0.5
  paramInit <- c(rep_len(phi0, length(phiNames)), init$pi0 %||% 0.75,
                 rep_len(theta0, length(nonRef)))
  charLabel <- mkd$orig_col[neoChars] %||% neoChars
  zNames <- paste0("eco_z_", rep(charLabel, each = length(nonRef)), "_",
                   rep(label, nNeo))

  pointer <- .EcologyPointerFactory(
    as.integer(eco) - 1L, kEco, refIdx - 1L, nNeo,
    magnitudeMode == "per_ecology", paramNames, as.numeric(paramInit),
    zNames, c(pi0Prior, thetaPrior, sigmaPhi), as.numeric(zWeight)
  )
  ext <- MkPrime::MkPrimeExtension(pointer)
  ext$ecology <- list(levels = ecoLevels, reference = ecoLevels[refIdx],
                      non_reference = nonRef, label = label,
                      magnitude_mode = magnitudeMode,
                      characters = charLabel, tips = tips)
  # Return:
  ext
}

# A closure holding only what the pointer needs, so that it travels cheaply
# to the worker processes that call it.
.EcologyPointerFactory <- function(ecology, kEco, refEco, nNeo, perEcology,
                                   paramNames, paramInit, zNames, priors,
                                   zWeight) {
  force(ecology); force(kEco); force(refEco); force(nNeo); force(perEcology)
  force(paramNames); force(paramInit); force(zNames); force(priors)
  force(zWeight)
  function() {
    .EcologyExtensionPtr(ecology, kEco, refEco, nNeo, perEcology, paramNames,
                         paramInit, zNames, priors, zWeight)
  }
}

`%||%` <- function(x, y) if (is.null(x)) y else x

#' Fold the phi label symmetry
#'
#' The EBE likelihood and prior are unchanged by `phi -> 1/phi` together with
#' `theta -> 1 - theta` and swapping `z` = 1 and 2 (per ecology, under
#' `magnitudeMode = "per_ecology"`). `CanonicaliseEcology()` maps each
#' sample to the representative with `phi >= 1`, so that summaries of
#' `theta` and `z` are interpretable.
#'
#' @param samples Matrix or data frame of samples, e.g. `result$samples`.
#' @param extension The `EcologyExtension()` object used for the run.
#' @return `samples`, folded.
#' @export
CanonicaliseEcology <- function(samples, extension) {
  info <- extension$ecology
  if (is.null(info)) stop("`extension` is not an ecology extension.")
  label <- info$label
  global <- identical(info$magnitude_mode, "global")
  for (j in seq_along(label)) {
    phiCol <- if (global) "eco_phi" else paste0("eco_phi_", label[j])
    flip <- samples[, phiCol] < 1
    if (!global) samples[flip, phiCol] <- 1 / samples[flip, phiCol]
    thetaCol <- paste0("eco_theta_", label[j])
    samples[flip, thetaCol] <- 1 - samples[flip, thetaCol]
    zCols <- intersect(paste0("eco_z_", info$characters, "_", label[j]),
                       colnames(samples))
    for (zc in zCols) {
      z <- samples[flip, zc]
      samples[flip, zc] <- ifelse(z == 1, 2, ifelse(z == 2, 1, z))
    }
  }
  if (global) {
    flip <- samples[, "eco_phi"] < 1
    samples[flip, "eco_phi"] <- 1 / samples[flip, "eco_phi"]
  }
  # Return:
  samples
}
