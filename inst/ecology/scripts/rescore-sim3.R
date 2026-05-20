# rescore-sim3.R --------------------------------------------------------------
#
# Re-score every saved v4/v4-cross MCMC result RDS using the root-invariant
# scorer (HasBipartSplits / ScoreTreesUnrooted) in
# inst/ecology/simulations/sim3-scoring.R, and emit a tidy data frame with
# both the corrected and the legacy (prop.part-based) bipartition
# probabilities for side-by-side comparison.
#
# Also reports, per RDS, the root-position bitmask distribution and the
# fraction of post-burnin samples whose root sits INSIDE each target
# tip-set — that fraction bounds how much the legacy scorer could
# disagree with the corrected one for that split.
#
# Outputs:
#   inst/ecology/scripts/rescore-sim3-results.rds  (tidy data frame)
#   inst/ecology/scripts/rescore-sim3-results.csv  (same, human-readable)
#   inst/ecology/scripts/rescore-sim3-root-dist.csv (root-position distributions)
#
# Run from the MkPrime root, e.g.:
#   MKP_LIB=/nobackup/pjjg18/mkp-sim3-multirep-v3/lib \
#   Rscript inst/ecology/scripts/rescore-sim3.R

suppressPackageStartupMessages({
  libPath <- Sys.getenv("MKP_LIB", "")
  if (nzchar(libPath)) .libPaths(c(libPath, .libPaths()))
  library("MkPrime")
  library("TreeTools")
  library("TreeDist")
})

mkpRoot <- Sys.getenv("MKP_REPO_ROOT", getwd())
source(file.path(mkpRoot, "inst/ecology/simulations/sim3-scoring.R"))
source(file.path(mkpRoot, "inst/ecology/simulations/sim3v4-helpers.R"))
source(file.path(mkpRoot, "inst/ecology/simulations/sim3v4cross-helpers.R"))
source(file.path(mkpRoot, "inst/ecology/simulations/sim3v4-params.R"))

# Root of saved results on Hamilton; allow override for local testing.
resultsRoot <- Sys.getenv("MKP_RESULTS_ROOT",
                          file.path("/nobackup", Sys.getenv("USER")))

# Catalogue of every result we want to rescore.
# Each entry: list(runName, rdsPath, design ("v4"|"v4cross"), paramSet,
#                  chain ("blind"|"aware"), repId).
buildCatalogue <- function() {
  out <- list()
  add <- function(...) out[[length(out) + 1L]] <<- list(...)

  for (variant in c("v4a", "v4b", "v4c")) {
    base <- file.path(resultsRoot, sprintf("mkp-sim3-%s", variant), "results")
    for (chain in c("blind", "aware")) {
      add(runName = sprintf("v4-%s-%s", variant, chain),
          rdsPath = file.path(base, sprintf("%s-result.rds", chain)),
          design  = "v4",
          paramSet = variant,
          chain   = chain,
          repId   = NA_integer_)
    }
  }
  for (variant in c("v4cross-b", "v4cross-b-pt", "v4cross-c")) {
    base <- file.path(resultsRoot, sprintf("mkp-sim3-%s", variant), "results")
    paramSet <- switch(variant,
                       "v4cross-b"     = "v4b",
                       "v4cross-b-pt"  = "v4b",
                       "v4cross-c"     = "v4c")
    for (chain in c("blind", "aware")) {
      add(runName = sprintf("%s-%s", variant, chain),
          rdsPath = file.path(base, sprintf("%s-result.rds", chain)),
          design  = "v4cross",
          paramSet = paramSet,
          chain   = chain,
          repId   = NA_integer_)
    }
  }
  # v4cross-b-pt5 — 5-rep array
  pt5Base <- file.path(resultsRoot, "mkp-sim3-v4cross-b-pt5", "results")
  for (rep in 1:5) {
    repDir <- file.path(pt5Base, sprintf("rep%02d", rep))
    for (chain in c("blind", "aware")) {
      add(runName = sprintf("v4cross-b-pt5-rep%02d-%s", rep, chain),
          rdsPath = file.path(repDir,
                              sprintf("%s-rep%02d-result.rds", chain, rep)),
          design  = "v4cross",
          paramSet = "v4b",
          chain   = chain,
          repId   = rep)
    }
  }
  out
}

# Build reference tree + biparts for a given design+paramSet.
designContext <- function(design, paramSet) {
  params <- SIM3V4_PARAMS[[paramSet]]
  if (is.null(params)) stop("Unknown paramSet: ", paramSet)
  if (identical(design, "v4")) {
    tree <- .BuildConvergentTreeV4(params$tipBr, params$stemBrEco,
                                    params$stemBrClade, params$rootBr)
    biparts <- .SimBipartitionsV4()
  } else if (identical(design, "v4cross")) {
    tree <- .BuildConvergentTreeV4Cross(params$tipBr, params$stemBrEco,
                                         params$stemBrClade, params$rootBr)
    biparts <- .SimBipartitionsV4Cross()
  } else {
    stop("Unknown design: ", design)
  }
  list(tree = tree, biparts = biparts)
}

# Process a single catalogue entry. Returns a data frame of scores, plus
# an attached attribute "rootDist" with the root-position counts.
processEntry <- function(entry) {
  if (!file.exists(entry$rdsPath)) return(NULL)
  res <- readRDS(entry$rdsPath)
  trees <- res$trees
  if (is.null(trees) || length(trees) == 0L) return(NULL)
  ctx <- designContext(entry$design, entry$paramSet)
  # Apply 25 % burn-in (same rule as run_v4*.R).
  trBurned <- .DiscardBurnin(trees)
  class(trBurned) <- "multiPhylo"
  nTrees <- length(trBurned)

  scores <- ScoreTreesUnrooted(trees, biparts = ctx$biparts,
                               refTree = ctx$tree,
                               discardBurnin = TRUE)
  # Per-target root-inside-set fractions (only for tip-set biparts).
  rootInside <- vapply(ctx$biparts, function(s) {
    ri <- RootInsideSet(trBurned, s)
    mean(ri, na.rm = TRUE)
  }, numeric(1))
  scoresAug <- merge(scores,
                     data.frame(metric = names(rootInside),
                                fracRootInside = unname(rootInside),
                                stringsAsFactors = FALSE),
                     by = "metric", all.x = TRUE)

  scoresAug$runName  <- entry$runName
  scoresAug$design   <- entry$design
  scoresAug$paramSet <- entry$paramSet
  scoresAug$chain    <- entry$chain
  scoresAug$repId    <- entry$repId
  scoresAug$rdsPath  <- entry$rdsPath

  # Root-position distribution (canonicalised, one side).
  rd <- RootConfigCounts(trBurned)
  rdDf <- data.frame(runName   = entry$runName,
                     chain     = entry$chain,
                     repId     = entry$repId,
                     rootConfig = names(rd),
                     count     = as.integer(unname(rd)),
                     nTrees    = nTrees,
                     stringsAsFactors = FALSE)
  attr(scoresAug, "rootDist") <- rdDf
  scoresAug
}

cat("=== rescore-sim3.R ===\n")
cat("resultsRoot:", resultsRoot, "\n")
cat("MKP_LIB:    ", Sys.getenv("MKP_LIB", "(unset)"), "\n\n")

entries <- buildCatalogue()
cat(sprintf("Catalogue: %d entries\n", length(entries)))

allScores <- list()
allRoot <- list()
nDone <- 0L
nMissing <- 0L
missingPaths <- character(0)
for (i in seq_along(entries)) {
  entry <- entries[[i]]
  if (!file.exists(entry$rdsPath)) {
    nMissing <- nMissing + 1L
    missingPaths <- c(missingPaths, entry$rdsPath)
    cat(sprintf("[%2d/%d] SKIP  %s (missing)\n",
                i, length(entries), entry$runName))
    next
  }
  t0 <- Sys.time()
  res <- tryCatch(processEntry(entry), error = function(e) {
    cat(sprintf("[%2d/%d] ERROR %s: %s\n",
                i, length(entries), entry$runName, conditionMessage(e)))
    NULL
  })
  if (is.null(res)) next
  nDone <- nDone + 1L
  allScores[[length(allScores) + 1L]] <- res
  rd <- attr(res, "rootDist"); attr(res, "rootDist") <- NULL
  allRoot[[length(allRoot) + 1L]] <- rd
  dt <- format(Sys.time() - t0)
  cat(sprintf("[%2d/%d] OK    %s  (n=%d trees, %s)\n",
              i, length(entries), entry$runName,
              res$nTrees[1L], dt))
}

cat(sprintf("\nRescored %d / %d entries  (%d missing)\n",
            nDone, length(entries), nMissing))
if (length(missingPaths)) {
  cat("Missing paths:\n")
  for (p in missingPaths) cat("  ", p, "\n")
}

outDir <- file.path(mkpRoot, "inst/ecology/scripts")
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
if (length(allScores)) {
  scoresDf <- do.call(rbind, allScores)
  scoreCols <- c("runName", "design", "paramSet", "chain", "repId",
                 "metric", "value_corrected", "value_legacy",
                 "fracRootInside", "nTrees", "rdsPath")
  scoresDf <- scoresDf[, scoreCols, drop = FALSE]
  saveRDS(scoresDf, file.path(outDir, "rescore-sim3-results.rds"))
  write.csv(scoresDf, file.path(outDir, "rescore-sim3-results.csv"),
            row.names = FALSE)
} else {
  cat("WARNING: no entries successfully rescored; no output written.\n")
}
if (length(allRoot)) {
  rootDf <- do.call(rbind, allRoot)
  write.csv(rootDf, file.path(outDir, "rescore-sim3-root-dist.csv"),
            row.names = FALSE)
}

cat("\nOutputs:\n")
cat("  ", file.path(outDir, "rescore-sim3-results.rds"), "\n")
cat("  ", file.path(outDir, "rescore-sim3-results.csv"), "\n")
cat("  ", file.path(outDir, "rescore-sim3-root-dist.csv"), "\n")

cat("\nDone.\n")
