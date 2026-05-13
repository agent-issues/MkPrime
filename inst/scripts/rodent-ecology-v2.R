# rodent-ecology-v2.R --------------------------------------------------------
# Empirical analysis: which characters are associated with which rodent
# ecologies, and how strongly, under the v2 ecology-aware MkPrime model
# (gamma-normalised, reference ecology, asymmetric slab).
#
# Inputs
#   - Rodent character matrix at ~/downloads/mbank_X24848_2026-5-9-1135.nex
#     (102 taxa, 220 morphological characters, col 220 = ecology, col 221 =
#     extant/extinct).  Drop col 221, drop extinct taxa.
#
# Outputs
#   - inst/scripts/rodent-ecology-v2-result.rds : RunMkPrime output + summaries
#   - inst/scripts/rodent-ecology-v2-zPost.csv  : per-(character, ecology)
#       posterior P(z = enc), P(z = disc), with character labels
#   - inst/scripts/rodent-ecology-v2-heatmap.pdf
#
# Headline question: "Per-character, per-ecology effect direction and
# strength."  Posterior P(z != none) × |log phi| gives a unitless effect size
# for ranking characters; sign(P(enc) - P(disc)) gives direction.
suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
  library("TreeTools")
})

args <- commandArgs(trailingOnly = TRUE)
nexFile <- if (length(args) >= 1) args[1] else "~/downloads/mbank_X24848_2026-5-9-1135.nex"
nIter   <- if (length(args) >= 2) as.integer(args[2]) else 500000L

stopifnot(file.exists(nexFile))
raw <- ape::read.nexus.data(nexFile)
taxa <- names(raw)
mat <- do.call(rbind, lapply(raw, function(x) as.character(unlist(x))))
rownames(mat) <- taxa
cat("Raw matrix:", nrow(mat), "tips x", ncol(mat), "chars\n")

ecologyCol <- 220L
extantCol  <- 221L
ecoVec  <- mat[, ecologyCol]
extant  <- mat[, extantCol]
keepTaxa <- which(extant == "0")  # 0 = extant, 1 = extinct
cat("Extant taxa:", length(keepTaxa), "of", nrow(mat), "\n")

matKeep <- mat[keepTaxa, -extantCol, drop = FALSE]  # drop col 221 only
ecoKeep <- ecoVec[keepTaxa]
ecoStates <- sort(unique(ecoKeep[!is.na(ecoKeep) & ecoKeep != "?"]))
cat("Ecology states:", paste(ecoStates, collapse = ", "), "\n")
cat("Ecology tip counts:\n"); print(table(ecoKeep))

# Build MkPrimeData with ecology col (now col 220 in the reduced matrix).
# Auto-detect neomorphic via a temporary phyDat object built from the
# matrix-minus-ecology-column (AutoDetectNeomorphic needs a phyDat).
pdForDetect <- MatrixToPhyDat(matKeep[, -ecologyCol, drop = FALSE])
neoIdx <- AutoDetectNeomorphic(pdForDetect)
cat("Neomorphic chars:", length(neoIdx), " | transformational:",
    ncol(matKeep) - 1L - length(neoIdx), "\n")

mkd <- MkPrimeData(matKeep, ecologyCol = ecologyCol,
                   neomorphic = neoIdx)
cat("MkPrimeData built:",
    " nTip=", mkd$nTip, " nChar=", mkd$nChar,
    " kEcology=", mkd$kEcology,
    " refEcology=", mkd$refEcology, "\n", sep = "")
cat("Ecology tip-fraction:\n"); print(round(mkd$ecoTipFraction, 3))

# Parsimony start tree (per ecology sim pattern; phangorn::optim.parsimony
# crashes on this Windows build).
set.seed(20260512)
randTree <- Preorder(ape::rtree(mkd$nTip, tip.label = mkd$taxon_names))
psTrees   <- TreeSearch::MaximizeParsimony(mkd$phyDat, tree = randTree,
                                           verbosity = 0)
startTree <- Preorder(psTrees[[1]])
startTree$edge.length <- rep(0.1, nrow(startTree$edge))

modelV2 <- MkPrimeModel(
  ecologyAware  = TRUE,
  magnitudeMode = "global",
  kPrimePrior   = "geometric",
  coding        = "variable",
  expSteps      = 10,
  rho0Alpha     = 7,
  rho0Beta      = 3,
  thetaAlpha    = 2,
  thetaBeta     = 2,
  sigmaPhi      = 0.5
)
mcmc <- MkPrimeMCMC(nIter = nIter, nChains = 1L, nRuns = 1L,
                    thin = max(1L, nIter %/% 1000L),
                    treeThin = max(1L, nIter %/% 500L),
                    minWarmup = 10000L, maxWarmup = 50000L,
                    logFile = "inst/scripts/rodent-ecology-v2.log",
                    checkpointFile = NULL)
for (f in c(mcmc$logFile, paste0(tools::file_path_sans_ext(mcmc$logFile),
                                  ".ckp"))) {
  if (file.exists(f)) file.remove(f)
}

t0 <- Sys.time()
res <- RunMkPrime(mkd, tree = startTree, model = modelV2, mcmc = mcmc)
cat("Elapsed:", format(Sys.time() - t0), "\n")

# Post-process posterior on z.
discard <- function(idx) idx[seq.int(ceiling(length(idx) / 4) + 1L,
                                     length(idx))]
samples <- ReadMkLog(mcmc$logFile)
samples <- samples[discard(seq_len(nrow(samples))), , drop = FALSE]

# Z samples may be stored in res$zSamples or similar -- adapt once v2 lands.
# Expected shape: an array nIterSaved x nChar x (kEco - 1) or a list.
if (!is.null(res$zSamples)) {
  zArr <- res$zSamples
  zArr <- zArr[discard(seq_len(dim(zArr)[1])), , , drop = FALSE]
  cat("Posterior z array:", paste(dim(zArr), collapse = " x "), "\n")
  # Per (character, ecology) marginals.
  pEnc <- apply(zArr, c(2, 3), function(x) mean(x == 1L))
  pDisc <- apply(zArr, c(2, 3), function(x) mean(x == 2L))
  pNone <- 1 - pEnc - pDisc
  # Ecology labels (skipping refEcology).
  ecoNames <- setdiff(seq.int(0L, mkd$kEcology - 1L), mkd$refEcology)
  colnames(pEnc) <- colnames(pDisc) <- colnames(pNone) <-
    paste0("eco", ecoNames)
  zPost <- data.frame(
    char = seq_len(mkd$nChar),
    type = mkd$type,
    pEnc, pDisc, pNone
  )
  write.csv(zPost, "inst/scripts/rodent-ecology-v2-zPost.csv",
            row.names = FALSE)
  cat("Wrote z posterior to inst/scripts/rodent-ecology-v2-zPost.csv\n")

  # Effect size per character: max over ecologies of P(z != none) × |log phi|.
  logPhi <- log(samples[, "phi"])
  meanAbsLogPhi <- mean(abs(logPhi))
  effSize <- (1 - pNone) * meanAbsLogPhi  # nChar x (kEco-1)
  topChar <- order(apply(effSize, 1, max), decreasing = TRUE)[1:20]
  cat("\n== Top 20 characters by ecology-effect strength ==\n")
  for (c in topChar) {
    cat(sprintf("  Char %3d (%s):", c, mkd$type[c]))
    for (j in seq_len(ncol(pEnc))) {
      cat(sprintf("  %s: enc=%.2f disc=%.2f", colnames(pEnc)[j],
                  pEnc[c, j], pDisc[c, j]))
    }
    cat("\n")
  }

  # Heatmap.
  pdf("inst/scripts/rodent-ecology-v2-heatmap.pdf", width = 6,
      height = max(4, mkd$nChar / 25))
  par(mar = c(4, 4, 3, 1))
  effSigned <- pEnc - pDisc  # +ve = predominantly encouraged
  image(t(effSigned[seq_len(min(80L, nrow(effSigned))), , drop = FALSE]),
        col = grDevices::colorRampPalette(c("steelblue", "white",
                                            "firebrick"))(64),
        zlim = c(-1, 1), axes = FALSE,
        xlab = "Ecology", ylab = "Character (first 80)",
        main = "P(z = enc) - P(z = disc)")
  axis(1, at = seq(0, 1, length.out = ncol(effSigned)),
       labels = colnames(effSigned))
  dev.off()
  cat("Wrote heatmap to inst/scripts/rodent-ecology-v2-heatmap.pdf\n")
}

# Save full result.
saveRDS(list(mkd = mkd, model = modelV2, mcmc = mcmc, res = res,
             samples = samples,
             pEnc = if (exists("pEnc")) pEnc else NULL,
             pDisc = if (exists("pDisc")) pDisc else NULL),
        "inst/scripts/rodent-ecology-v2-result.rds")
cat("\nDone.\n")
