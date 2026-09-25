# Shared helpers for the MkPrime vs RevBayes oracle-validation harness.
#
# These are dev/-only utilities; do not promote to the package surface
# without rethinking the API (see the plan file for rationale).

suppressPackageStartupMessages({
  library(TreeTools)
})

# The harness scripts assume `%||%` (base R >= 4.4, the version the Hamilton
# smoke.slurm pins via `module load r/4.5.1`); fall back to rlang's for local
# / cloud runs on an older R rather than silently requiring 4.4+.
if (!exists("%||%", mode = "function")) {
  `%||%` <- rlang::`%||%`
}

#' Read a single split nex file and return the character matrix (rows = tips,
#' cols = chars). Returns NULL with a warning if the file has zero characters
#' after the standard parsimony-uninformative filtering applied by neotrans's
#' PrepareMatrix(). The declared NEXUS `SYMBOLS="..."` string is attached as
#' `attr(., "nexusSymbols")`: RevBayes indexes states by *position in that
#' declaration*, not by the digit's own value or by sorted order, so
#' `RbStateK()` must use it rather than assuming `0-9A-Z`
#' (agent-issues/MkPrime#214).
ReadSplit <- function(path) {
  stopifnot(file.exists(path))
  m <- TreeTools::ReadCharacters(path)
  if (is.null(m) || !length(m) || ncol(m) == 0L) {
    return(NULL)
  }
  attr(m, "nexusSymbols") <- .ReadNexusSymbols(path)
  m
}

#' Extract the declared `SYMBOLS="..."` alphabet from a NEXUS FORMAT line.
#' Falls back to the standard `0-9A-Z` NEXUS default if the file doesn't
#' declare one explicitly.
#' @param path NEXUS file path.
#' @return Character vector of single-character symbols, in declared order.
.ReadNexusSymbols <- function(path) {
  txt <- paste(readLines(path, warn = FALSE), collapse = " ")
  m <- regmatches(
    txt, regexpr('SYMBOLS\\s*=\\s*"[^"]*"', txt, ignore.case = TRUE, perl = TRUE)
  )
  if (!length(m) || !nzchar(m)) {
    return(c(as.character(0:9), LETTERS))
  }
  sym <- sub('(?i)symbols\\s*=\\s*"([^"]*)"', "\\1", m, perl = TRUE)
  strsplit(sym, "", fixed = TRUE)[[1]]
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

#' RevBayes' own state count for one transformational character column.
#'
#' RevBayes sets k to the maximum observed state *index* + 1, indexed by
#' **position in the NEXUS file's own declared `SYMBOLS="..."` string** (not
#' by digit value or sorted order -- a file declaring `SYMBOLS="ABC"` indexes
#' A=0, B=1, C=2), counting each state named inside a polymorphic/uncertain
#' token (`{01}`, `(12)`) individually. This differs from `kObs`
#' (MkPrimeData's distinct-non-ambiguous-token count): a column coded `{1,2}`
#' in a file with the default `SYMBOLS="0123456789"` has kObs = 2 but RB's
#' k = 3 (see RlAbstractHomologousDiscreteCharacterData.cpp:571-622,
#' agent-issues/MkPrime#214).
#'
#' @param col Character vector of raw NEXUS tokens for one character, as
#'   returned by a column of `TreeTools::ReadCharacters()`'s matrix.
#' @param alphabet Character vector giving the file's declared `SYMBOLS`
#'   order, e.g. from `attr(mat, "nexusSymbols")`.
#' @return Integer scalar; `0L` if every cell is missing/gap.
RbStateK <- function(col, alphabet = c(as.character(0:9), LETTERS)) {
  toks <- col[!is.na(col) & col != "?" & col != "-"]
  if (!length(toks)) return(0L)
  syms <- unlist(strsplit(gsub("[{}()]", "", toks), "", fixed = TRUE))
  idx <- match(toupper(syms), toupper(alphabet)) - 1L
  if (anyNA(idx)) {
    stop(
      "Token symbol(s) not in the declared SYMBOLS alphabet (",
      paste(alphabet, collapse = ""), "): ",
      paste(unique(syms[is.na(idx)]), collapse = ", ")
    )
  }
  max(idx) + 1L
}

#' Vectorised `RbStateK()` over every column of a raw character matrix.
#' @param mat Character matrix as returned by `TreeTools::ReadCharacters()`.
#' @param alphabet Character vector; see `RbStateK()`.
#' @return Integer vector of length `ncol(mat)`.
ComputeRbK <- function(mat, alphabet = c(as.character(0:9), LETTERS)) {
  vapply(seq_len(ncol(mat)), function(j) RbStateK(mat[, j], alphabet), integer(1))
}

#' Replace polymorphic/uncertain tokens (`{01}`, `(12)`) with `?` (fully
#' missing).  MkPrimeData already treats a polymorphic cell as fully missing
#' (`R/MkPrimeData.R`'s phyDat conversion has no partial-ambiguity state for
#' it), while RevBayes treats the same token as a partial ambiguity --  a
#' genuine data-content mismatch, not just a labelling one
#' (agent-issues/MkPrime#214). Applying this to both samplers' input before
#' either sees the data means they agree on content, at the cost of RB losing
#' the extra information a partial ambiguity would have given it; this is the
#' conservative fix and the one that requires no change to `R/MkPrimeData.R`.
#'
#' @param mat Character matrix as returned by `TreeTools::ReadCharacters()`.
#' @return `mat`, with `{...}`/`(...)` cells replaced by `"?"`.
SanitizePolymorphic <- function(mat) {
  poly <- grepl("^[{(].*[})]$", mat)
  mat[poly] <- "?"
  mat
}

#' Build the canonical inputs both samplers must be run against for one
#' (pid, model) cell: the common taxon set, per-character state counts
#' matching RevBayes' own rule, and the post-invariant-drop partition sizes
#' MkPrimeData will actually use.  Computing this once and sharing it (via
#' `WriteOrCheckCellInfo()`) is what lets `render_rev.R` and `run_mkprime.R`
#' -- separate script invocations -- provably agree on what they fed their
#' respective sampler (agent-issues/MkPrime#214).
#'
#' Not a gap: MkPrimeData drops columns that become invariant after the
#' taxon intersection and polymorphism sanitisation, while RevBayes keeps
#' them in its data object (`nchar()` still counts them, hence the NCHAR
#' substitution below) but excludes them from the likelihood under
#' `coding="variable"` -- confirmed by a smoke test against a real RevBayes
#' build (2026-09-25, Hamilton, PR agent-issues/MkPrime#231): the
#' log-likelihood is identical with and without a constant column present.
#' The two samplers reach the same likelihood by different bookkeeping, not
#' by different models.
#'
#' @param transMat,neoMat Raw character matrices as returned by `ReadSplit()`
#'   (must carry the `"nexusSymbols"` attribute `ReadSplit()` attaches).
#' @param model One of `"by_nt_9v"`, `"by_nt_kv"`.
#' @return A list with `commonTaxa`, `excludeNeoTaxa`, `excludeTransTaxa`
#'   (taxa present in one split but not the other -- RevBayes must exclude
#'   these explicitly, since its template reads `neo.names()` alone),
#'   `transMatCanonical`, `neoMatCanonical` (taxon-intersected,
#'   polymorphism-sanitised matrices both samplers should read), `nTrans`,
#'   `nNeo`, `kRb` (RB-rule state count per trans character), `neoIdx`,
#'   `knownStates`, `nNeoFinal`, `nTransFinal` (partition sizes after
#'   MkPrimeData's invariant-character drop).
PrepareCellInfo <- function(transMat, neoMat, model) {
  transTaxa <- rownames(transMat)
  neoTaxa <- rownames(neoMat)
  commonTaxa <- intersect(transTaxa, neoTaxa)
  if (!length(commonTaxa)) {
    stop("trans and neo have no taxa in common")
  }

  # `[` subsetting drops non-standard attributes, so capture the declared
  # SYMBOLS alphabets before subsetting and reattach them explicitly --
  # WriteFilteredNexus() must see the SOURCE file's declaration, not the
  # default fallback, or state indices it writes silently shift.
  transSymbols <- attr(transMat, "nexusSymbols") %||% c(as.character(0:9), LETTERS)
  neoSymbols <- attr(neoMat, "nexusSymbols") %||% c(as.character(0:9), LETTERS)

  transSan <- SanitizePolymorphic(transMat[commonTaxa, , drop = FALSE])
  neoSan <- SanitizePolymorphic(neoMat[commonTaxa, , drop = FALSE])
  attr(transSan, "nexusSymbols") <- transSymbols
  attr(neoSan, "nexusSymbols") <- neoSymbols

  kRb <- ComputeRbK(transSan, transSymbols)
  if (any(kRb > 10L)) {
    stop(
      "Character(s) ", paste(which(kRb > 10L), collapse = ", "),
      " have RevBayes state count k > 10; RevBayes silently drops these ",
      "(RlAbstractHomologousDiscreteCharacterData.cpp), so the two ",
      "samplers cannot be fed matching data for this cell as-is."
    )
  }

  comb <- CombineSplits(transSan, neoSan)
  matrixInfo <- list(nTrans = comb$nTrans, nNeo = comb$nNeo, kObs = kRb)
  neoIdx <- BuildNeoIdx(matrixInfo, model)
  knownStates <- BuildKnownStates(matrixInfo, model)

  phy <- MatrixToCombinedPhyDat(comb$matrix)
  mkd <- MkPrime::MkPrimeData(phy, neomorphic = neoIdx, knownStates = knownStates)

  list(
    commonTaxa = commonTaxa,
    excludeNeoTaxa = setdiff(neoTaxa, commonTaxa),
    excludeTransTaxa = setdiff(transTaxa, commonTaxa),
    transMatCanonical = transSan,
    neoMatCanonical = neoSan,
    nTrans = comb$nTrans,
    nNeo = comb$nNeo,
    kRb = kRb,
    neoIdx = neoIdx,
    knownStates = knownStates,
    nNeoFinal = sum(mkd$type == "neomorphic"),
    nTransFinal = sum(mkd$type == "known")
  )
}

#' Write a raw character matrix (tokens as returned by
#' `TreeTools::ReadCharacters()`) out as a minimal NEXUS standard-data file,
#' so RevBayes can be pointed at exactly the taxon set and (sanitised) data
#' MkPrime used (agent-issues/MkPrime#214, RB-119: the harness intersects
#' trans/neo taxa but RB's template previously read the unfiltered splits and
#' clamped to `neo.names()` alone). The `SYMBOLS` declaration is preserved
#' verbatim from the source file (via `symbols`/`attr(mat, "nexusSymbols")`)
#' rather than re-derived from the cells present, so state indices RB assigns
#' do not shift. Round-trips through `TreeTools::ReadCharacters()` unchanged.
#'
#' @param mat Character matrix with taxon names as `rownames()`.
#' @param path File path to write to.
#' @param symbols Character vector giving the `SYMBOLS` declaration order;
#'   defaults to `attr(mat, "nexusSymbols")`.
#' @return `path`, invisibly.
WriteFilteredNexus <- function(mat, path,
                               symbols = attr(mat, "nexusSymbols")) {
  taxa <- rownames(mat)
  stopifnot(!is.null(taxa), !anyNA(taxa), !anyDuplicated(taxa))
  if (is.null(symbols)) symbols <- c(as.character(0:9), LETTERS)
  present <- unlist(strsplit(gsub("[{}()]", "", mat), "", fixed = TRUE))
  present <- present[!is.na(present) & !present %in% c("?", "-") & nzchar(present)]
  missingSym <- setdiff(toupper(unique(present)), toupper(symbols))
  if (length(missingSym)) {
    stop("Symbol(s) present in data but not in the SYMBOLS declaration: ",
        paste(missingSym, collapse = ", "))
  }
  cellStr <- apply(mat, 1, paste, collapse = "")
  taxLabel <- function(x) if (grepl("[^A-Za-z0-9_.]", x)) sprintf("'%s'", x) else x
  lines <- c(
    "#NEXUS",
    "BEGIN DATA;",
    sprintf("  DIMENSIONS NTAX=%d NCHAR=%d;", nrow(mat), ncol(mat)),
    sprintf(
      "  FORMAT DATATYPE=STANDARD MISSING=? GAP=- SYMBOLS=\"%s\" INTERLEAVE=NO;",
      paste(symbols, collapse = "")
    ),
    "  MATRIX",
    sprintf("    %s  %s", vapply(taxa, taxLabel, character(1)), cellStr),
    "  ;",
    "END;"
  )
  writeLines(lines, path)
  invisible(path)
}

#' Persist a `PrepareCellInfo()` result the first time it is computed for a
#' cell, or assert that a later computation (from a separate script
#' invocation) reproduces it exactly.  This is the harness's cross-process
#' equality check that RB and MkPrime were fed the same data
#' (agent-issues/MkPrime#214): if `render_rev.R` and `run_mkprime.R` ever
#' disagree -- e.g. one is run against a stale `matrices.rds`/split file --
#' this stops the run rather than silently comparing mismatched cells.
#'
#' @param path File path for the shared cell-info cache.
#' @param info A `PrepareCellInfo()` result.
#' @return `info`, invisibly.
WriteOrCheckCellInfo <- function(path, info) {
  if (file.exists(path)) {
    prev <- readRDS(path)
    if (!identical(prev, info)) {
      stop(
        "Cell info at ", path, " does not match this run's freshly computed ",
        "value -- MkPrime and RevBayes would be fed different data/state ",
        "counts. Delete the file to regenerate if the input data changed ",
        "intentionally."
      )
    }
  } else {
    dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
    saveRDS(info, path)
  }
  invisible(info)
}

#' Provenance stamp for one harness run: enough to tell a current rds from a
#' stale one (agent-issues/MkPrime#215, RB-107).  A pure mtime rule is
#' unsound (a merge commit's date can predate the rds it should invalidate,
#' `dev/rb-equivalence/` edits were never covered by it, and `fixup_*`
#' scripts bump mtimes); this instead records content identity directly, so
#' `compare.R` can refuse a `templateHash` mismatch outright rather than
#' guess from timestamps.
#'
#' @param scriptDir The `dev/rb-equivalence` directory (for locating the
#'   git checkout root and the Rev templates to hash); `NULL` skips both,
#'   leaving `gitSha`/`gitDirty`/`templateHash` `NA`.
#' @param rbBin Optional path to the RevBayes binary used (e.g. `$RB_BIN`),
#'   recorded verbatim -- this harness cannot itself query an `rb --version`.
#' @return A list: `gitSha`/`gitDirty` of the harness checkout,
#'   `templateHash` (of every `templates/*.Rev` file, sorted), `rbBin`,
#'   package versions, and `timestamp`.
HarnessProvenance <- function(scriptDir = NULL, rbBin = NA_character_) {
  gitSha <- NA_character_
  gitDirty <- NA
  templateHash <- NA_character_
  if (!is.null(scriptDir)) {
    gitSha <- tryCatch(
      system2("git", c("-C", scriptDir, "rev-parse", "HEAD"),
             stdout = TRUE, stderr = FALSE),
      error = function(e) NA_character_
    )
    gitSha <- if (length(gitSha) == 1L) gitSha else NA_character_
    gitDirty <- tryCatch(
      length(system2("git", c("-C", scriptDir, "status", "--porcelain", "."),
                     stdout = TRUE, stderr = FALSE)) > 0L,
      error = function(e) NA
    )
    templateDir <- file.path(scriptDir, "templates")
    templateFiles <- sort(list.files(templateDir, pattern = "\\.Rev$",
                                     full.names = TRUE))
    templateHash <- tryCatch(
      rlang::hash(vapply(templateFiles, function(f) {
        paste(readLines(f, warn = FALSE), collapse = "\n")
      }, character(1))),
      error = function(e) NA_character_
    )
  }
  list(
    gitSha = gitSha,
    gitDirty = gitDirty,
    templateHash = templateHash,
    rbBin = rbBin,
    mkprimeVersion = as.character(utils::packageVersion("MkPrime")),
    treeToolsVersion = as.character(utils::packageVersion("TreeTools")),
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  )
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
