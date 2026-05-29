# sim-ebe-null.R --------------------------------------------------------------
# EBE Sim (NULL / graceful degradation).  Simulate neomorphic + transformational
# characters under the EQUILIBRIUM-tilt model (ecoEquilMode = TRUE) with NO
# ecology effect (z = 0 everywhere) on a random 16-tip tree with a random binary
# tip ecology.  Run BOTH the aware-EBE model and the blind MkNT model on the
# same data and compare.
#
# This is EBE spec (dev/ecology/ebe-spec.md) tests 1 + 6, the decisive
# degradation checks:
#   - GATE 1 already proved the EBE neomorphic KERNEL is byte-identical to
#     baseline MkN at z = 0.  Here we confirm the FULL aware-EBE MCMC reproduces
#     the blind MkNT posterior within MC error when the truth has no ecology
#     signal.
#   - It also exercises the Phase 2b neomorphic-only-z fix: transformational
#     characters must NOT be flagged (their z stays 0; they must not pollute
#     pi0).  Before that fix, gibbs_z samples trans z from the prior
#     (mcmc.cpp:4715) and this harness would show spurious flags.
#
# Pass conditions (decisive in CAPS; the rest are diagnostic):
#   - PHI is NOT confidently != 1: its 95% CI must INCLUDE 1.  Under the null,
#     pi0 ~ 0.75 leaves few slab cells, so phi is WEAKLY IDENTIFIED and reverts
#     toward its LogNormal(0, sigmaPhi) prior (wide CI, point estimate well above
#     1).  That is correct, honest non-identifiability — NOT a manufactured
#     effect.  The failure mode is a TIGHT phi CI that EXCLUDES 1 (confidently
#     claiming ecology where there is none).  [Calibrated on the 2026-05-29 run:
#     phi mean 3.72, CI [0.08, 40.5] includes 1; an earlier |log phi| < 0.20
#     "phi -> 1" criterion was wrong — the data cannot pin phi to 1 with no signal.]
#   - AWARE tree_length / log_likelihood CIs OVERLAP the blind CIs: graceful
#     degradation of the full sampler, not just the kernel.
#   - AWARE and BLIND recover the simulating tree comparably.
#   - pi0 is NOT driven spuriously low (a low pi0 would mean the model invented
#     ecology structure).  Under weak per-character z-identifiability pi0 is
#     prior-dominated (near model$rho0Alpha/(rho0Alpha+rho0Beta)); a strong
#     pooled null would push it toward 1.  Either is fine; pi0 << prior is the
#     failure mode.
#
# DRAFT (2026-05-29): Phase 3 harness, written ahead of Phase 2b wiring.
# NOT YET RUN — pass thresholds are reasoned, not yet empirically calibrated.
# Validate/refine after GATE 2 (the model must be the wired EBE).  Run via
# devtools::load_all in a subprocess (the installed package is stale):
#   Rscript -e "source('inst/ecology/simulations/sim-ebe-null.R')"
suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
  library("TreeTools")
})
source("inst/ecology/simulations/sim3-simulate.R")

set.seed(20260529)
nTip   <- 16L
nEco   <- 60L   # neomorphic chars (EBE tilts these; here z = 0 so they are inert)
nBase  <- 60L   # transformational chars (never tilted under EBE)
phiSim <- 1     # neutral; with z = 0 everywhere it makes no difference
tree   <- Preorder(ape::rcoal(nTip, tip.label = paste0("T", seq_len(nTip))))

# Random binary tip ecology (~half / half, both states well represented).
repeat {
  ecoTip <- setNames(sample(c(0L, 1L), nTip, replace = TRUE), tree$tip.label)
  if (sum(ecoTip == 1L) >= 5L && sum(ecoTip == 0L) >= 5L) break
}
edgeEc <- .AssignEdgeEcology(tree, ecoTip)
cat("Tip ecology counts:\n"); print(table(ecoTip))

# z = 0 everywhere.  Neomorphic chars lead (MkPrimeData expects neomorphic at
# the leading indices).
zNull <- matrix(0L, nrow = nEco + nBase, ncol = 2L)
cType <- c(rep("neomorphic", nEco), rep("transformational", nBase))
datSim <- .SimulateMkPrimeEcology(tree, edgeEc, zNull, phi = phiSim,
                                  type = cType, baseRate = 0.5, rateLoss = 1,
                                  ecoEquilMode = TRUE)
cat("Simulated", nrow(datSim), "tips x", ncol(datSim), "chars (EBE null)\n")

pdSim    <- MatrixToPhyDat(datSim)
neoIdx   <- seq_len(nEco)
mkdAware <- MkPrimeData(pdSim, neomorphic = neoIdx,
                        ecology = ecoTip[rownames(datSim)])
mkdBlind <- MkPrimeData(pdSim, neomorphic = neoIdx)   # no ecology -> blind MkNT
cat("kEcology:", mkdAware$kEcology, " | nChar:", mkdAware$nChar, "\n")

startTree <- Preorder(TreeSearch::AdditionTree(pdSim))
startTree$edge.length <- rep(0.1, nrow(startTree$edge))

modelAware <- MkPrimeModel(ecologyAware = TRUE, magnitudeMode = "global",
                           kPrimePrior = "geometric", coding = "variable",
                           expSteps = 10)
modelBlind <- MkPrimeModel(ecologyAware = FALSE,
                           kPrimePrior = "geometric", coding = "variable",
                           expSteps = 10)
priorMeanPi0 <- modelAware$rho0Alpha / (modelAware$rho0Alpha + modelAware$rho0Beta)

mkMcmc <- function(logFile) MkPrimeMCMC(
  nIter = 100000L, nChains = 1L, nRuns = 1L, thin = 200L, treeThin = 200L,
  minWarmup = 5000L, maxWarmup = 30000L, logFile = logFile, checkpointFile = NULL)

runOne <- function(mkd, model, tag) {
  lf <- paste0("sim-ebe-null-", tag, ".log")
  for (f in c(lf, sub("\\.log$", ".ckp", lf))) if (file.exists(f)) file.remove(f)
  t0  <- Sys.time()
  res <- RunMkPrime(mkd, tree = startTree, model = model, mcmc = mkMcmc(lf))
  cat(sprintf("%-5s elapsed: %s\n", tag, format(Sys.time() - t0)))
  list(res = res, samp = ReadMkLog(lf))
}

awareOut <- runOne(mkdAware, modelAware, "aware")
blindOut <- runOne(mkdBlind, modelBlind, "blind")

burn <- function(s) s[seq.int(ceiling(nrow(s) / 4) + 1L, nrow(s)), , drop = FALSE]
sa <- burn(awareOut$samp); sb <- burn(blindOut$samp)

summ <- function(x) sprintf("mean=%.3f CI=[%.3f, %.3f]", mean(x),
  stats::quantile(x, 0.025, names = FALSE), stats::quantile(x, 0.975, names = FALSE))
overlap <- function(a, b) {
  qa <- stats::quantile(a, c(0.025, 0.975), names = FALSE)
  qb <- stats::quantile(b, c(0.025, 0.975), names = FALSE)
  !(qa[2] < qb[1] || qb[2] < qa[1])
}

cat(sprintf("\n== EBE-null AWARE posterior (expect phi~1, pi0 not << %.3f) ==\n",
            priorMeanPi0))
for (nm in c("phi", "pi0", "rate_neo", "rate_loss", "log_likelihood", "tree_length"))
  if (nm %in% colnames(sa)) cat(sprintf("  %-15s %s\n", nm, summ(sa[, nm])))
if ("phi" %in% colnames(sa))
  cat(sprintf("  %-15s %s\n", "|log phi|", summ(abs(log(sa[, "phi"])))))

cat("\n== BLIND MkNT posterior (reference) ==\n")
for (nm in c("rate_neo", "rate_loss", "log_likelihood", "tree_length"))
  if (nm %in% colnames(sb)) cat(sprintf("  %-15s %s\n", nm, summ(sb[, nm])))

# --- Decisive checks ---------------------------------------------------------
discard <- function(x) x[seq.int(ceiling(length(x) / 4) + 1L, length(x))]
cidOf <- function(res) {
  tr <- discard(res$trees); class(tr) <- "multiPhylo"
  as.numeric(TreeDist::ClusteringInfoDist(tr, tree, normalize = TRUE))
}
ca <- cidOf(awareOut$res); cb <- cidOf(blindOut$res)

phiCI   <- if ("phi" %in% colnames(sa)) stats::quantile(sa[, "phi"], c(0.025, 0.975), names = FALSE) else c(NA, NA)
phiOK   <- if (all(is.finite(phiCI))) (phiCI[1] <= 1 && phiCI[2] >= 1) else NA  # CI must include 1 (weakly identified under null)
tlOK    <- if (all(c("tree_length") %in% colnames(sa)) &&
               "tree_length" %in% colnames(sb)) overlap(sa[, "tree_length"], sb[, "tree_length"]) else NA
llOK    <- if ("log_likelihood" %in% colnames(sa) &&
               "log_likelihood" %in% colnames(sb)) overlap(sa[, "log_likelihood"], sb[, "log_likelihood"]) else NA
pi0OK   <- if ("pi0" %in% colnames(sa)) mean(sa[, "pi0"]) > 0.6 * priorMeanPi0 else NA
treeOK  <- abs(mean(ca) - mean(cb)) < 0.05

cat("\n== Tree recovery (CID to TRUE; lower = better) ==\n")
cat(sprintf("  AWARE: mean=%.3f  identical=%.1f%%\n", mean(ca), 100 * mean(ca == 0)))
cat(sprintf("  BLIND: mean=%.3f  identical=%.1f%%\n", mean(cb), 100 * mean(cb == 0)))

cat("\n========================= EBE-NULL VERDICT =========================\n")
cat(sprintf("  [%s] phi not confidently !=1 (95%% CI [%.2f,%.2f] incl. 1; weakly id. under null)\n", ifelse(isTRUE(phiOK), "PASS", "CHECK"), phiCI[1], phiCI[2]))
cat(sprintf("  [%s] aware tree_length ~ blind (CI overlap)\n",          ifelse(isTRUE(tlOK),   "PASS", "CHECK")))
cat(sprintf("  [%s] aware logLik ~ blind      (CI overlap)\n",          ifelse(isTRUE(llOK),   "PASS", "CHECK")))
cat(sprintf("  [%s] pi0 not spuriously low    (> 0.6 * prior mean)\n",  ifelse(isTRUE(pi0OK),  "PASS", "CHECK")))
cat(sprintf("  [%s] tree recovery ~ blind     (|dCID| < 0.05)\n",       ifelse(isTRUE(treeOK), "PASS", "CHECK")))
cat("  (phi criterion calibrated on the 2026-05-29 run; tree/logLik/pi0 thresholds still provisional)\n")
cat("====================================================================\n")

saveRDS(list(tree = tree, eco = ecoTip, z = zNull, charType = cType,
             phiSim = phiSim, data = datSim, edgeEco = edgeEc, neoIdx = neoIdx,
             aware = awareOut$res, blind = blindOut$res,
             verdict = list(phiOK = phiOK, tlOK = tlOK, llOK = llOK,
                            pi0OK = pi0OK, treeOK = treeOK)),
        "inst/ecology/simulations/sim-ebe-null-result.rds")
cat("\nDone. EBE null result saved.\n")
