# diagnostic-multirep-mixing.R — apply mixing diagnostic to multirep-v3 reps.
#
# For each rep, compute:
#   (1) chain MAP logPost (chain-reported)
#   (2) fresh-recompute logPost at MAP state (R-side, current code)
#   (3) logPost(truth tree, truth TL) at MAP nuisance (fresh recompute)
#   (4) #samples containing the true AC bipartition
#   (5) unique topology count
#   (6) chain TL trajectory summary
#
# Usage:
#   Rscript inst/simulations/ecology/diagnostic-multirep-mixing.R 01 05

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
  library("TreeTools")
})
source("inst/simulations/ecology/sim3-helpers.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) args <- c("01", "05")

# Truth (matches run_rep.R)
nEco <- 120L; nBase <- 360L; phi <- 4
stemBr <- 0.30; rootBr <- 0.15; tipBr <- 0.5
truthTree <- .BuildConvergentTree(tipBranch = tipBr,
                                   stemBranch = stemBr, rootBranch = rootBr)
truthTL <- sum(truthTree$edge.length)
trueSplit  <- c(paste0("A", 1:4), paste0("C", 1:4))

compute_aware_logpost <- function(tree, mkd, model,
                                  kPrime, rate_loss, rate_log_sd,
                                  rate_neo, p, phi, pi0, theta, zMat) {
  ll <- .MkpEcologyLogLikelihood(
    tree = tree, mkd = mkd, kPrime = kPrime,
    rate_loss = rate_loss, rate_log_sd = rate_log_sd,
    nCat = model$nCat %||% 6L,
    rate_neo = rate_neo, relabel = model$relabel,
    phi = phi, zMat = zMat,
    magnitudeMode = model$magnitudeMode %||% "global",
    coding = model$coding %||% "variable",
    refEcology = mkd$refEcology %||% 0L,
    theta = theta, pi0 = pi0
  )
  tl <- sum(tree$edge.length)
  rel <- if (tl > 0) tree$edge.length / tl
         else rep(1 / nrow(tree$edge), nrow(tree$edge))
  state <- list(tree = tree, tree_length = tl, rel_br_lengths = rel,
                rate_loss = rate_loss, rate_log_sd = rate_log_sd,
                rate_neo = rate_neo, p = p, kPrime = kPrime,
                phi = phi, pi0 = pi0, theta = theta, z = zMat)
  lp <- LogPrior(state = state, model = model, mkd = mkd)
  list(loglik = ll, logprior = lp, logpost = ll + lp)
}

chain_map <- function(res) {
  s <- as.data.frame(res$samples)
  if (nrow(s) == 0) return(NULL)
  idx <- which.max(s$log_posterior)
  kPrimeCols <- grep("^kPrime_", colnames(s), value = TRUE)
  kPrimeIdx <- as.integer(sub("kPrime_", "", kPrimeCols))
  kPrime <- integer(max(kPrimeIdx))
  kPrime[kPrimeIdx] <- as.integer(s[idx, kPrimeCols])
  list(
    idx = idx,
    tree = res$trees[[idx]],
    sample_logpost = s$log_posterior[idx],
    sample_loglik = s$log_likelihood[idx],
    rate_loss = s$rate_loss[idx],
    rate_log_sd = s$rate_log_sd[idx],
    rate_neo = s$rate_neo[idx],
    p = s$p[idx],
    phi = if ("phi" %in% colnames(s)) s$phi[idx] else NA_real_,
    pi0 = if ("pi0" %in% colnames(s)) s$pi0[idx] else NA_real_,
    theta = if ("theta_1" %in% colnames(s)) s$theta_1[idx] else NA_real_,
    kPrime = kPrime,
    tree_length = s$tree_length[idx]
  )
}

hasACBipart <- function(treeList) {
  vapply(treeList, function(tr) {
    cl <- ape::prop.part(tr)
    tips <- attr(cl, "labels")
    splitSet <- which(tips %in% trueSplit)
    any(vapply(cl, function(p) setequal(p, splitSet), logical(1)))
  }, logical(1))
}

diagnose_rep <- function(repId) {
  d <- sprintf("inst/simulations/ecology/multirep-v3-results/rep%s", repId)
  cat(sprintf("\n========== REP %s ==========\n", repId))

  summ <- readRDS(file.path(d, "summary.rds"))
  cat(sprintf("Headline: aware P(true)=%.3f P(wrong)=%.3f CID=%.3f\n",
              summ$aware$pTrue, summ$aware$pWrong, summ$aware$cidTrue))
  cat(sprintf("          blind P(true)=%.3f P(wrong)=%.3f CID=%.3f\n",
              summ$blind$pTrue, summ$blind$pWrong, summ$blind$cidTrue))

  resA <- readRDS(file.path(d, "aware-result.rds"))
  resA <- RelabelEcology(resA)
  mkdA <- resA$data; modelA <- resA$model

  mapA <- chain_map(resA)
  cat(sprintf("\n[AWARE] MAP idx=%d  chain logPost=%.3f  TL=%.3f\n",
              mapA$idx, mapA$sample_logpost, mapA$tree_length))
  cat(sprintf("  MAP params: phi=%.2f pi0=%.3f theta=%.3f rate_loss=%.2f rate_neo=%.2f\n",
              mapA$phi, mapA$pi0, mapA$theta, mapA$rate_loss, mapA$rate_neo))

  # z at MAP
  zMAP <- if (!is.null(resA$z_samples)) resA$z_samples[[mapA$idx]] else NULL
  if (is.null(zMAP)) {
    cat("  (no z_samples; skipping recompute diagnostic for this rep)\n")
  } else {
    # Recompute at MAP state
    rcA <- compute_aware_logpost(
      tree = mapA$tree, mkd = mkdA, model = modelA,
      kPrime = mapA$kPrime,
      rate_loss = mapA$rate_loss, rate_log_sd = mapA$rate_log_sd,
      rate_neo = mapA$rate_neo, p = mapA$p,
      phi = mapA$phi, pi0 = mapA$pi0, theta = mapA$theta, zMat = zMAP
    )
    cat(sprintf("  Recompute at MAP:   logPost=%.3f  (drift vs chain log: %+.3f)\n",
                rcA$logpost, rcA$logpost - mapA$sample_logpost))

    # Truth tree at truth TL, MAP nuisance
    rcT <- compute_aware_logpost(
      tree = truthTree, mkd = mkdA, model = modelA,
      kPrime = mapA$kPrime,
      rate_loss = mapA$rate_loss, rate_log_sd = mapA$rate_log_sd,
      rate_neo = mapA$rate_neo, p = mapA$p,
      phi = mapA$phi, pi0 = mapA$pi0, theta = mapA$theta, zMat = zMAP
    )
    cat(sprintf("  Truth tree (TL=%.2f) + MAP nuisance: logPost=%.3f\n",
                truthTL, rcT$logpost))
    cat(sprintf("  DELTA truth - MAP recompute = %+.3f nats  (positive = chain undersampled)\n",
                rcT$logpost - rcA$logpost))
  }

  hasAC <- hasACBipart(resA$trees)
  cat(sprintf("  Samples containing true AC: %d / %d (%.3f%%)\n",
              sum(hasAC), length(hasAC), 100 * mean(hasAC)))
  tlA <- as.data.frame(resA$samples)$tree_length
  cat(sprintf("  TL trace: min=%.2f med=%.2f max=%.2f last100mean=%.2f  truth=%.2f\n",
              min(tlA), median(tlA), max(tlA), mean(tail(tlA, 100)), truthTL))
  nUniq <- length(unique(resA$samples[, "topo_hash"]))
  cat(sprintf("  Unique topologies: %d / %d (%.1f%%)\n",
              nUniq, nrow(as.matrix(resA$samples)),
              100 * nUniq / nrow(as.matrix(resA$samples))))

  # Acceptance rates
  cat("\n  Aware acceptance (key moves):\n")
  acc <- resA$acceptance
  key <- c("nni", "spr", "tbr", "pspr", "gibbs_spr", "gibbs_subtree_swap",
           "tree_length", "branch_lengths", "scale_phi", "scale_pi0",
           "scale_theta", "gibbs_z", "block_kPrime")
  for (k in key) {
    if (k %in% names(acc))
      cat(sprintf("    %-22s %.3f\n", k, acc[[k]]))
  }

  # Blind
  resB <- readRDS(file.path(d, "blind-result.rds"))
  mapB <- chain_map(resB)
  if (!is.null(mapB)) {
    cat(sprintf("\n[BLIND] MAP idx=%d  chain logPost=%.3f  TL=%.3f\n",
                mapB$idx, mapB$sample_logpost, mapB$tree_length))
    hasACB <- hasACBipart(resB$trees)
    cat(sprintf("  Samples containing true AC: %d / %d (%.3f%%)\n",
                sum(hasACB), length(hasACB), 100 * mean(hasACB)))
    nUniqB <- length(unique(resB$samples[, "topo_hash"]))
    cat(sprintf("  Unique topologies: %d / %d (%.1f%%)\n",
                nUniqB, nrow(as.matrix(resB$samples)),
                100 * nUniqB / nrow(as.matrix(resB$samples))))
  }
}

for (rep in args) diagnose_rep(rep)
cat("\n=== done ===\n")
