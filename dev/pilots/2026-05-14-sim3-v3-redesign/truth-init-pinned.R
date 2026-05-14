# Sim 3 v3 truth-init with PINNED ecology move weights
#
# Hypothesis: ecology moves starved at 1% can't counter the rapid drift
# from joint rate moves (joint_tl_rl 7%, rate_neo 6.9%). Pin gibbs_z,
# scale_phi, scale_pi0, scale_theta at higher weights and see if the
# chain can stay at truth.
#
# Run from worktree root:
#   Rscript -e "setwd('C:/Users/pjjg18/GitHub/mkp/.claude/worktrees/ecology-aware'); source('dev/pilots/2026-05-14-sim3-v3-redesign/truth-init-pinned.R')"

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
  library("TreeTools")
})
source("inst/simulations/ecology/sim3-helpers.R")

PILOT_DIR <- "dev/pilots/2026-05-14-sim3-v3-redesign"
saved <- readRDS(file.path(PILOT_DIR, "result.rds"))

tree <- saved$tree
eco  <- saved$eco
datSim <- saved$data
nEco <- saved$config$nEco

pdSim <- MatrixToPhyDat(datSim)
mkdAware <- MkPrimeData(pdSim, neomorphic = seq_len(nEco), ecology = eco)

modelAware <- MkPrimeModel(ecologyAware = TRUE,
                           magnitudeMode = "global",
                           kPrimePrior = "geometric",
                           coding = "variable", expSteps = 10)

# Truth-init patch
origInit <- MkPrime:::.InitState
truthInit <- function(tree, mkd, model) {
  st <- origInit(tree = saved$tree, mkd = mkd, model = model)
  st$phi <- saved$config$phi
  st$pi0 <- 0.75
  st$theta <- 0.999  # avoid Beta(2,2) prior boundary
  st$rate_neo <- 1.0
  st$rate_loss <- 1.0
  st$tree_length <- saved$truthTL
  st$tree <- saved$tree
  st$rel_br_lengths <- saved$tree$edge.length / saved$truthTL
  zTruth <- saved$z[, -1L, drop = FALSE]
  storage.mode(zTruth) <- "integer"
  st$z <- zTruth
  st$log_lik <- MkPrime:::.MkpEcologyLogLikelihood(
    tree = saved$tree, mkd = mkd, kPrime = st$kPrime,
    rate_loss = 1.0, rate_log_sd = st$rate_log_sd,
    nCat = model$nCat, rate_neo = 1.0,
    relabel = model$relabel,
    phi = st$phi, zMat = st$z,
    magnitudeMode = model$magnitudeMode,
    coding = model$coding, refEcology = mkd$refEcology,
    theta = st$theta, pi0 = st$pi0
  )
  st$log_prior <- MkPrime:::LogPrior(st, model, mkd)
  st$log_post <- st$log_lik + st$log_prior
  cat(sprintf("[truth-init] log_lik=%.2f log_prior=%.2f log_post=%.2f\n",
              st$log_lik, st$log_prior, st$log_post))
  st
}
assignInNamespace(".InitState", truthInit, ns = "MkPrime")

# PINNED weights: give ecology moves enough share to counter rate-move drift
pinnedWeights <- c(
  gibbs_z      = 0.08,
  scale_phi    = 0.03,
  scale_pi0    = 0.03,
  scale_theta  = 0.03
)

nIter <- 50000L
thin <- 40L
mcmc <- MkPrimeMCMC(nIter = nIter, nChains = 1L, nRuns = 1L,
                    thin = thin, treeThin = thin,
                    minWarmup = 1000L, maxWarmup = 2000L,
                    moveWeights = pinnedWeights,
                    logFile = NULL, checkpointFile = NULL)
logFile <- file.path(PILOT_DIR, "truth-init-pinned-chain.log")
ckpFile <- file.path(PILOT_DIR, "truth-init-pinned-chain.ckp")
for (f in c(logFile, ckpFile)) if (file.exists(f)) file.remove(f)
mcmc$logFile <- logFile
mcmc$checkpointFile <- ckpFile

cat("== Running truth-init aware chain WITH PINNED ECOLOGY MOVES ==\n")
cat("Pinned weights (total = ", sum(pinnedWeights), "):\n", sep = "")
for (nm in names(pinnedWeights))
  cat(sprintf("  %-12s = %.2f\n", nm, pinnedWeights[[nm]]))

t0 <- Sys.time()
res <- tryCatch(RunMkPrime(mkdAware, tree = saved$tree,
                            model = modelAware, mcmc = mcmc),
                error = function(e) {
                  cat("ERROR:", conditionMessage(e), "\n"); NULL
                })
cat("Elapsed:", format(Sys.time() - t0), "\n")
assignInNamespace(".InitState", origInit, ns = "MkPrime")

if (!is.null(res)) {
  res <- RelabelEcology(res)
  saveRDS(list(res = res, truthTL = saved$truthTL,
               phi_true = saved$config$phi,
               pinnedWeights = pinnedWeights),
          file.path(PILOT_DIR, "truth-init-pinned-result.rds"))
  cat("Saved truth-init-pinned-result.rds\n")
}
