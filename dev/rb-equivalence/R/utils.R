# Shared helpers for the MkPrime vs RevBayes oracle-validation harness.
#
# These are dev/-only utilities; do not promote to the package surface
# without rethinking the API (see the plan file for rationale).

suppressPackageStartupMessages({
  library(TreeTools)
})

#' Read a single split nex file and return the character matrix (rows = tips,
#' cols = chars). Returns NULL with a warning if the file has zero characters
#' after the standard parsimony-uninformative filtering applied by neotrans's
#' PrepareMatrix().
ReadSplit <- function(path) {
  stopifnot(file.exists(path))
  m <- TreeTools::ReadCharacters(path)
  if (is.null(m) || !length(m) || ncol(m) == 0L) {
    return(NULL)
  }
  m
}

#' Combine separate trans + neo character matrices into a single matrix with
#' all transformational characters appearing FIRST.  Both matrices must share
#' the same taxon set (in any order); the combined matrix is in trans-row order.
#' Returns a list with `matrix`, `nTrans`, `nNeo`.
CombineSplits <- function(trans, neo) {
  stopifnot(!is.null(trans))
  if (is.null(neo)) {
    return(list(matrix = trans, nTrans = ncol(trans), nNeo = 0L))
  }
  common <- intersect(rownames(trans), rownames(neo))
  if (!length(common)) {
    stop("trans and neo have no taxa in common")
  }
  if (length(common) != nrow(trans) || length(common) != nrow(neo)) {
    cli::cli_warn(
      "trans ({nrow(trans)} tips) and neo ({nrow(neo)} tips) taxon sets \\
       differ; using {length(common)} taxa in common"
    )
  }
  trans <- trans[common, , drop = FALSE]
  neo <- neo[common, , drop = FALSE]
  combined <- cbind(trans, neo)
  list(matrix = combined, nTrans = ncol(trans), nNeo = ncol(neo))
}

#' Coerce a character matrix (as returned by TreeTools::ReadCharacters) into a
#' phyDat object suitable for MkPrime::MkPrimeData.  TreeTools::MatrixToPhyDat
#' auto-detects tokens from the matrix.
MatrixToCombinedPhyDat <- function(mat) {
  TreeTools::MatrixToPhyDat(mat)
}

#' Build the knownStates vector for a given (matrices entry, model).
#' Returns a named integer with names = char indices (as char strings) and
#' values = k.  trans chars always occupy indices 1:nTrans by our convention.
BuildKnownStates <- function(matrixInfo, model) {
  nTrans <- matrixInfo$nTrans
  transIdx <- seq_len(nTrans)
  k <- switch(model,
    "by_nt_9v" = rep(9L, nTrans),
    "by_nt_kv" = matrixInfo$kObs,
    stop("Unknown model: ", model)
  )
  ks <- setNames(as.integer(k), as.character(transIdx))
  # Assert contract before passing to MkPrimeData
  ksIdx <- as.integer(names(ks))
  if (any(is.na(ksIdx)) || any(ksIdx < 1L)) {
    stop("knownStates names must be positive integer-coercible char indices")
  }
  ks
}

#' Build the neomorphic-indices vector for a (matrixInfo, model) pair.
#' For trans-only models (not currently in scope but supported), returns
#' integer(0).
BuildNeoIdx <- function(matrixInfo, model) {
  if (matrixInfo$nNeo == 0L) return(integer(0))
  matrixInfo$nTrans + seq_len(matrixInfo$nNeo)
}

#' Construct the MkPrimeModel matching RevBayes priors exactly.
RBMatchedModel <- function() {
  MkPrime::MkPrimeModel(
    coding = "variable",
    nCat = 6L,
    treeLengthShape = 2,
    expSteps = 1,  # treeLengthRate = gamma_shape / exp_steps = 2/1 = 2
    rateLogSdShape = 1,
    rateLogSdRate = 1,
    rateLossMeanlog = 0,
    rateLossSdlog = 2,
    rateNeoMeanlog = 0,
    rateNeoSdlog = 2
  )
}
