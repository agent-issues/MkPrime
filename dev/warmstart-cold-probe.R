# Cold-chain probe: start the AWARE model at the blind consensus tree and
# watch the log-likelihood trajectory + topological drift.
#
#   logLik RISES as it leaves blind  -> aware target assigns HIGHER likelihood
#                                        to non-blind trees => target prefers
#                                        biologically-wrong trees (misspec).
#   logLik flat/falls, topology stays -> blind tree is a good region; failure
#                                        to reach it elsewhere is a sampler issue.
#
# Single cold chain (nChains=1) => no PT confound.
suppressPackageStartupMessages({library(ape); library(TreeTools)
  devtools::load_all(".", quiet = TRUE)})

nexFile <- "inst/ecology/data/rodent-X24848.nex"
od <- "inst/ecology/scripts/rodent-comparison"

mat <- TreeTools::ReadCharacters(nexFile)
ecologyCol <- 220L; extantCol <- 221L
ecoVec <- mat[, ecologyCol]; extant <- mat[, extantCol]
keepTaxa <- which(extant == "0")
matKeep <- mat[keepTaxa, -extantCol, drop = FALSE]
ecoKeep <- ecoVec[keepTaxa]
poly <- grepl("^[(]", ecoKeep)
if (any(poly)) { ecoKeep[poly] <- substr(sub("^[(]", "", ecoKeep[poly]), 1, 1)
  matKeep[poly, ecologyCol] <- ecoKeep[poly] }
hasEco <- !is.na(ecoKeep) & ecoKeep != "?" & !is.na(suppressWarnings(as.integer(ecoKeep)))
matKeep <- matKeep[hasEco, , drop = FALSE]; ecoKeep <- ecoKeep[hasEco]

charMat <- matKeep[, -ecologyCol, drop = FALSE]
pdForDetect <- MatrixToPhyDat(charMat)
neoIdx <- AutoDetectNeomorphic(pdForDetect)
.kObsCol <- function(col) { vals <- col[!(col %in% c("?", "-", NA))]
  vals <- unlist(strsplit(gsub("[()]", "", vals), "")); length(unique(vals)) }
kObsRaw <- vapply(seq_len(ncol(charMat)), function(j) .kObsCol(charMat[, j]), integer(1L))
nonNeo <- setdiff(seq_len(ncol(charMat)), neoIdx)
nonNeoVar <- nonNeo[kObsRaw[nonNeo] >= 2L]
knownK <- setNames(as.integer(kObsRaw[nonNeoVar]), as.character(nonNeoVar))
mkd <- MkPrimeData(pdForDetect, neomorphic = neoIdx, knownStates = knownK,
                   ecology = setNames(as.integer(ecoKeep), rownames(matKeep)))

# blind consensus start tree
bl <- do.call(c, lapply(Sys.glob(file.path(od, "rodent-MkNT-blind_trees_*.nwk")), ape::read.tree))
class(bl) <- "multiPhylo"
bl <- bl[seq.int(ceiling(length(bl) / 4) + 1L, length(bl))]
consB <- ape::multi2di(ape::consensus(bl, p = 0.5, rooted = FALSE), random = FALSE)
consB$edge.length <- rep(0.1, nrow(consB$edge))
consB <- TreeTools::Preorder(consB)
storage.mode(consB$edge) <- "integer"
stopifnot(ape::is.binary(consB))
startTree <- consB
blindSplits <- prop.part(consB)

model <- MkPrimeModel(ecologyAware = TRUE, magnitudeMode = "global",
                      kPrimePrior = "geometric", coding = "variable",
                      rho0Alpha = 7, rho0Beta = 3, thetaAlpha = 2, thetaBeta = 2,
                      sigmaPhi = 0.5)

nIter <- 4000L
logFile <- "dev/warmstart_cold_probe.log"
if (file.exists(logFile)) file.remove(logFile)
mcmc <- MkPrimeMCMC(nIter = nIter, nChains = 1L, nRuns = 1L, nCore = 1L,
                    thin = 20L, treeThin = 20L, minWarmup = 500L, maxWarmup = 1500L,
                    logFile = logFile, checkpointFile = NULL)

cat("=== COLD aware chain from BLIND consensus ===\n")
t0 <- Sys.time()
res <- RunMkPrime(mkd, tree = startTree, model = model, mcmc = mcmc)
cat("Elapsed:", format(Sys.time() - t0), "\n")

s <- ReadMkLog(logFile)
if (is.matrix(s)) s <- as.data.frame(s, check.names = FALSE)
ll <- s[["log_likelihood"]]
cat(sprintf("\nlogLik: first=%.1f  q25=%.1f  median=%.1f  max=%.1f  last=%.1f\n",
            ll[1], quantile(ll, .25), median(ll), max(ll), ll[length(ll)]))
cat(sprintf("tree_length: first=%.2f  last=%.2f\n",
            s[["tree_length"]][1], s[["tree_length"]][nrow(s)]))

tf <- sub("\\.log$", "_trees.nwk", logFile)
if (!file.exists(tf)) tf <- paste0(logFile, "_trees.nwk")
if (file.exists(tf)) {
  tr <- ape::read.tree(tf); if (inherits(tr, "phylo")) tr <- list(tr)
  class(tr) <- "multiPhylo"
  n <- length(tr)
  rf_start <- sapply(tr, function(t) ape::dist.topo(ape::unroot(t), ape::unroot(consB)))
  maxrf <- 2 * (Ntip(consB) - 3)
  cat(sprintf("\nRF(sample, blind start) over chain: first=%.0f%%  last=%.0f%%  max=%.0f%% of %d\n",
              100 * rf_start[1] / maxrf, 100 * rf_start[n] / maxrf,
              100 * max(rf_start) / maxrf, maxrf))
  # shared splits between final-quarter consensus and blind start
  post <- tr[seq.int(ceiling(n / 2) + 1L, n)]
  cFin <- ape::consensus(post, p = 0.5, rooted = FALSE)
  shared <- sum(prop.part(cFin) %in% blindSplits)  # crude; report counts
  cat(sprintf("final-half consensus: %d internal splits\n", length(prop.part(cFin)) - 1L))
}
cat("\nDONE-PROBE\n")
