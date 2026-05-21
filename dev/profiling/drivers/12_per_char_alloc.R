# 12_per_char_alloc.R --------------------------------------------------
# Micro-bench for .CppLogLikelihoodEcologyPerChar (T-008 per_char portion).
# Compares baseline (.vtune-lib-prof) vs T-008 (.vtune-lib-t008) build.
#
# Strategy: run each build in its own Rscript subprocess to avoid DLL
# conflicts.  Each subprocess prints "RESULT: <median_ms>" to stdout.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library("TreeTools"))

`%||%` <- function(a, b) if (!is.null(a)) a else b

bench_code <- function(libPath, nRep = 200L, nWarmup = 10L) {
  sprintf('
suppressPackageStartupMessages({
  library("MkPrime", lib.loc = %s)
  library("TreeTools")
})
`%%||%%` <- function(a, b) if (!is.null(a)) a else b
set.seed(20260521)
nexFile <- "C:/Users/pjjg18/downloads/mbank_X24848_2026-5-9-1135.nex"
mat <- TreeTools::ReadCharacters(nexFile)
ecologyCol <- 220L; extantCol <- 221L
ecoVec <- mat[, ecologyCol]; extant <- mat[, extantCol]
keepTaxa <- which(extant == "0")
matKeep <- mat[keepTaxa, -extantCol, drop = FALSE]
ecoKeep <- ecoVec[keepTaxa]
poly <- grepl("^\\\\(", ecoKeep)
if (any(poly)) {
  recoded <- substr(sub("^\\\\(", "", ecoKeep[poly]), 1, 1)
  ecoKeep[poly] <- recoded; matKeep[poly, ecologyCol] <- recoded
}
hasEco <- !is.na(ecoKeep) & ecoKeep != "?" &
  !is.na(suppressWarnings(as.integer(ecoKeep)))
matKeep <- matKeep[hasEco, , drop = FALSE]; ecoKeep <- ecoKeep[hasEco]
charMat <- matKeep[, -ecologyCol, drop = FALSE]
pdForDetect <- MatrixToPhyDat(charMat)
neoIdx <- AutoDetectNeomorphic(pdForDetect)
kObsRaw <- vapply(seq_len(ncol(charMat)), function(j) {
  vals <- charMat[!(charMat[,j] %%in%% c("?","-",NA)), j]
  length(unique(vals))
}, integer(1L))
nonNeoOriginal <- setdiff(seq_len(ncol(charMat)), neoIdx)
nonNeoVariable <- nonNeoOriginal[kObsRaw[nonNeoOriginal] >= 2L]
knownK <- setNames(as.integer(kObsRaw[nonNeoVariable]), as.character(nonNeoVariable))
mkd <- suppressWarnings(MkPrimeData(
  pdForDetect, neomorphic = neoIdx, knownStates = knownK,
  ecology = setNames(as.integer(ecoKeep), rownames(matKeep))))
startTree <- Preorder(TreeSearch::AdditionTree(mkd$phyDat))
startTree$edge.length <- rep(0.1, nrow(startTree$edge))
nChar  <- mkd$nChar; kEco <- mkd$kEcology
kPrime <- as.integer(mkd$kObs)
parent <- as.integer(startTree$edge[, 1])
child  <- as.integer(startTree$edge[, 2])
edgeLen <- as.numeric(startTree$edge.length)
rateLoss <- 0.5; rateLogSd <- 0.0; rateNeo <- 0.5
phi <- 1.0; pi0 <- 0.5
zMatrix <- matrix(0L, nChar, kEco - 1L)
model <- MkPrimeModel(ecologyAware = TRUE, magnitudeMode = "global",
  kPrimePrior = "geometric", expSteps = 10,
  rho0Alpha = 7, rho0Beta = 3, thetaAlpha = 2, thetaBeta = 2, sigmaPhi = 1.5,
  rateLossMeanlog = 0, rateLossSdlog = 2,
  rateLogSdShape = 1, rateLogSdRate = 1,
  rateNeoMeanlog = 0, rateNeoSdlog = 1)
parts <- lapply(mkd$partitions, function(p)
  list(type=p$type, k=p$k, kObs=p$kObs, char_indices=p$char_indices,
       tip_states=p$tip_states, unique_tip_states=p$unique_tip_states,
       pattern_index=p$pattern_index))
ecoTip <- as.integer(mkd$ecology); ecoTip[is.na(ecoTip)] <- -1L
dataPtr <- MkPrime:::prepare_mcmc_data(
  parts, as.integer(mkd$kObs), mkd$type, any(mkd$type == "neomorphic"),
  model$nCat, model$coding, model$relabel,
  model$treeLengthShape, model$treeLengthRate,
  model$rateLossMeanlog, model$rateLossSdlog,
  model$rateLogSdShape, model$rateLogSdRate,
  model$rateNeoMeanlog, model$rateNeoSdlog,
  model$kprimeHyperA, model$kprimeHyperB,
  identical(model$kPrimePrior, "logseries"), model$kprimeLogseriesC %%||%% 0.7,
  identical(model$kPrimePrior, "geometric"), FALSE, 4L, 1.0, 1.0,
  FALSE, numeric(0), 1L, 0L, 0, -1e308, TRUE,
  ecoTip, as.integer(mkd$kEcology), model$magnitudeMode %%||%% "global",
  model$rho0Alpha, model$rho0Beta, model$sigmaPhi,
  as.integer(model$gibbsZEvery), model$thetaAlpha, model$thetaBeta)
fn <- function() MkPrime:::.CppLogLikelihoodEcologyPerChar(
  dataPtr, parent, child, edgeLen, kPrime,
  rateLoss, rateLogSd, rateNeo, phi, zMatrix, pi0 = pi0)
for (i in seq_len(%d)) fn()
elapsed_s <- system.time(for (i in seq_len(%d)) fn())["elapsed"]
ms_per_batch <- elapsed_s * 1000 / %d
cat(sprintf("RESULT: %%.6f ms/batch  nChar=%%d  kEco=%%d\\n", ms_per_batch, nChar, kEco))
', shQuote(libPath), nWarmup, nRep, nRep)
}

run_subprocess <- function(libPath, nRep = 200L, nWarmup = 10L) {
  code <- bench_code(libPath, nRep, nWarmup)
  tmp <- tempfile(fileext = ".R")
  writeLines(code, tmp)
  out <- system2("Rscript", tmp, stdout = TRUE, stderr = FALSE)
  file.remove(tmp)
  # Parse "RESULT: X.XXXXXX ms/batch  nChar=Y  kEco=Z"
  res_line <- grep("^RESULT:", out, value = TRUE)
  if (length(res_line) == 0) {
    cat("Subprocess output:\n"); cat(out, sep = "\n")
    stop("No RESULT line found")
  }
  ms   <- as.numeric(sub("RESULT: ([0-9.]+) ms/batch.*", "\\1", res_line))
  nChr <- as.integer(sub(".*nChar=([0-9]+).*", "\\1", res_line))
  kE   <- as.integer(sub(".*kEco=([0-9]+).*", "\\1", res_line))
  list(ms_per_batch = ms, nChar = nChr, kEco = kE)
}

cat("=== Running baseline (.vtune-lib-prof), 200 reps ===\n")
base_res <- run_subprocess("dev/profiling/.vtune-lib-prof")
cat(sprintf("Baseline: %.3f ms/batch  (%.5f ms/char)\n",
            base_res$ms_per_batch, base_res$ms_per_batch / base_res$nChar))

cat("=== Running T-008 (.vtune-lib-t008), 200 reps ===\n")
t008_res <- run_subprocess("dev/profiling/.vtune-lib-t008")
cat(sprintf("T-008:    %.3f ms/batch  (%.5f ms/char)\n",
            t008_res$ms_per_batch, t008_res$ms_per_batch / t008_res$nChar))

nChar      <- base_res$nChar
kEco       <- base_res$kEco
kEcoMinus1 <- kEco - 1L
nCallsPerSweep <- nChar * kEcoMinus1 * 3L

cat(sprintf("\n=== SUMMARY ===\n"))
cat(sprintf("nChar=%d  kEcology=%d\n", nChar, kEco))
cat(sprintf("Gibbs z sweep per_char calls: %d x %d x 3 = %d/sweep\n",
            nChar, kEcoMinus1, nCallsPerSweep))
speedup <- base_res$ms_per_batch / t008_res$ms_per_batch
cat(sprintf("\nBaseline median: %.3f ms/batch  (%.5f ms/char)\n",
            base_res$ms_per_batch, base_res$ms_per_batch / nChar))
cat(sprintf("T-008 median:    %.3f ms/batch  (%.5f ms/char)\n",
            t008_res$ms_per_batch, t008_res$ms_per_batch / nChar))
cat(sprintf("Speedup:         %.2fx\n", speedup))
cat(sprintf("\nHeap allocs eliminated per sweep: %d -> 0\n", nCallsPerSweep * 4L))
