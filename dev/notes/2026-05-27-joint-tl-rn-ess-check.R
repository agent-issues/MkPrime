# Empirical ESS/sec check for joint_tl_rn
#
# Note: dev/notes/2026-05-27-rate-neo-ridge-and-joint-moves.md A1 deliverable.
#
# Hypothesis: under the Issue-1 partition-rate fix, (T, r) sit on a tight
# posterior ridge on mixed datasets. The new joint_tl_rn 2D Bactrian should
# beat independent Bactrians on rate_neo ESS/sec.
#
# To run:  Rscript dev/notes/2026-05-27-joint-tl-rn-ess-check.R
#
# Output:  dev/notes/2026-05-27-joint-tl-rn-ess-check.txt

devtools::load_all(quiet = TRUE)
suppressPackageStartupMessages({
  library(coda)
})

# Build a smoke mixed dataset: small TreeSearch dataset with one binary char
# forced neomorphic.
nexFile <- system.file("datasets/Sun2018.nex", package = "TreeSearch")
pd      <- TreeTools::ReadAsPhyDat(nexFile)
mkdBase <- MkPrimeData(pd)
neoIdx  <- which(mkdBase$kObs == 2L)[1]
stopifnot(!is.na(neoIdx))
mkd <- MkPrimeData(pd, neomorphic = neoIdx)

nIter <- 6000L

cfgOn <- MkPrimeMCMC(
  nIter = nIter, nRuns = 1L, nChains = 1L,
  minWarmup = 500L, maxWarmup = 1500L,
  autoTune = TRUE, thin = 5L,
  joint2d = TRUE
)
cfgOff <- MkPrimeMCMC(
  nIter = nIter, nRuns = 1L, nChains = 1L,
  minWarmup = 500L, maxWarmup = 1500L,
  autoTune = TRUE, thin = 5L,
  joint2d = FALSE
)

bench <- function(cfg, seed) {
  set.seed(seed)
  t0  <- Sys.time()
  res <- RunMkPrime(mkd, mcmc = cfg)
  dt  <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  rn  <- res$samples[, "rate_neo"]
  tl  <- res$samples[, "tree_length"]
  list(
    secs        = dt,
    ess_rn      = coda::effectiveSize(rn),
    ess_tl      = coda::effectiveSize(tl),
    ess_per_sec = coda::effectiveSize(rn) / dt,
    acc         = res$acceptance
  )
}

on  <- bench(cfgOn,  9001)
off <- bench(cfgOff, 9001)

out <- c(
  sprintf("joint2d = TRUE:  rate_neo ESS = %.1f, tree_length ESS = %.1f, t = %.1fs, ESS/s = %.2f",
          on$ess_rn,  on$ess_tl,  on$secs,  on$ess_per_sec),
  sprintf("joint2d = FALSE: rate_neo ESS = %.1f, tree_length ESS = %.1f, t = %.1fs, ESS/s = %.2f",
          off$ess_rn, off$ess_tl, off$secs, off$ess_per_sec),
  sprintf("ESS/sec ratio (on/off) on rate_neo: %.2f x",
          on$ess_per_sec / off$ess_per_sec),
  "",
  "Acceptance for joint moves (ON):",
  paste0("  ", names(on$acc), " = ", sprintf("%.3f", on$acc),
         collapse = "\n")
)
cat(out, sep = "\n")
writeLines(out, "dev/notes/2026-05-27-joint-tl-rn-ess-check.txt")
