# diagnostic-mixing.R — discriminate undersampling from posterior shape.
#
# Reads the post-fix dev pilot result.rds + chain log; answers:
#   (1) Is logPost(truth tree | MAP nuisance) > or < max logPost in chain?
#       > truth: chain *should* have visited but didn't  -> mixing failure
#       < truth: posterior actively disfavours the truth -> not a mixing fix
#   (2) How many unique topologies did each chain visit?
#       low count = stuck; high count = mixing but not converging
#   (3) Per-move-type acceptance rates (already in resXxx$acceptance);
#       flag gibbs_spr / gibbs_subtree_swap which should be 0% under eco gate.
#   (4) Does the true topology hash ever appear in the chain trace?
#
# Run:
#   Rscript dev/pilots/2026-05-14-sim3-v3-redesign/diagnostic-mixing.R

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
  library("TreeTools")
  library("ape")
})

PILOT <- "dev/pilots/2026-05-14-sim3-v3-redesign"
r <- readRDS(file.path(PILOT, "result.rds"))
cfg <- r$config

cat("=== Sim 3 v3 dev pilot diagnostic ===\n")
cat(sprintf("config: nEco=%d nBase=%d phi=%g stem=%g root=%g tip=%g  truthTL=%g\n",
            cfg$nEco, cfg$nBase, cfg$phi, cfg$stemBr, cfg$rootBr, cfg$tipBr,
            r$truthTL))

# ---- Truth: tree + params + z + char type
truthTree <- r$tree
truthZ    <- r$z
truthPhi  <- cfg$phi
truthPi0  <- 0.75
truthTheta <- 1.0
truthEdgeRel <- truthTree$edge.length / sum(truthTree$edge.length)
truthTL   <- r$truthTL

# Helper to compute full log-posterior at an arbitrary state.
compute_logpost <- function(tree, mkd, model, kPrime,
                            rate_loss, rate_log_sd, rate_neo, p,
                            phi, pi0, theta, zMat,
                            ecology_aware) {
  if (ecology_aware) {
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
  } else {
    ll <- .MkpLogLikelihood(
      tree = tree, mkd = mkd, kPrime = kPrime,
      rate_loss = rate_loss, rate_log_sd = rate_log_sd,
      nCat = model$nCat %||% 6L,
      rate_neo = rate_neo, relabel = model$relabel,
      coding = model$coding %||% "variable"
    )
  }
  # Build state for LogPrior
  tl <- sum(tree$edge.length)
  rel <- if (tl > 0) tree$edge.length / tl else rep(1 / nrow(tree$edge), nrow(tree$edge))
  state <- list(
    tree = tree, tree_length = tl, rel_br_lengths = rel,
    rate_loss = rate_loss, rate_log_sd = rate_log_sd, rate_neo = rate_neo,
    p = p, kPrime = kPrime,
    phi = phi, pi0 = pi0, theta = theta, z = zMat
  )
  lp <- LogPrior(state = state, model = model, mkd = mkd)
  list(loglik = ll, logprior = lp, logpost = ll + lp)
}

# Get chain MAP: highest log_posterior sample; assemble its params
chain_map_state <- function(res, mkd) {
  s <- as.data.frame(res$samples)
  idx <- which.max(s$log_posterior)
  tree <- res$trees[[idx]]
  kPrimeCols <- grep("^kPrime_", colnames(s), value = TRUE)
  kPrimeIndices <- as.integer(sub("kPrime_", "", kPrimeCols))
  nCharTotal <- max(kPrimeIndices)
  kPrime <- integer(nCharTotal)
  kPrime[kPrimeIndices] <- as.integer(s[idx, kPrimeCols])
  list(
    idx = idx,
    tree = tree,
    sample_logpost = s$log_posterior[idx],
    sample_loglik  = s$log_likelihood[idx],
    rate_loss   = s$rate_loss[idx],
    rate_log_sd = s$rate_log_sd[idx],
    rate_neo    = s$rate_neo[idx],
    phi         = if ("phi" %in% colnames(s)) s$phi[idx] else NA_real_,
    pi0         = if ("pi0" %in% colnames(s)) s$pi0[idx] else NA_real_,
    theta       = if ("theta_1" %in% colnames(s)) s$theta_1[idx] else NA_real_,
    kPrime      = kPrime,
    topo_hash   = if ("topo_hash" %in% colnames(s)) s$topo_hash[idx] else NA
  )
}

# Construct z from sampled chain at MAP (aware only)
chain_map_z <- function(res, idx, nChar) {
  z_samples <- res$z_samples
  if (is.null(z_samples)) return(NULL)
  # z_samples may be list of matrices per sample, or array
  if (is.list(z_samples)) return(z_samples[[idx]])
  if (length(dim(z_samples)) == 3) return(z_samples[idx, , , drop = TRUE])
  NULL
}

# Hash the truth topology (sorted bipartition signature)
truth_hash <- tryCatch(
  paste(sort(vapply(ape::prop.part(truthTree),
                    function(p) paste(sort(p), collapse = ","),
                    character(1))),
        collapse = "|"),
  error = function(e) NA_character_)

# ---- Aware diagnostic ----
cat("\n--- AWARE chain ---\n")
resA <- r$resAware
mkdA <- resA$data
modelA <- resA$model
mapA <- chain_map_state(resA, mkdA)
zMAP_A <- chain_map_z(resA, mapA$idx, nrow(truthZ))
if (is.null(zMAP_A)) zMAP_A <- truthZ  # fallback

cat(sprintf("MAP idx=%d  chain logPost = %.3f  (loglik %.3f)\n",
            mapA$idx, mapA$sample_logpost, mapA$sample_loglik))
cat(sprintf("MAP params: phi=%.3f pi0=%.3f theta=%.3f rate_loss=%.3f rate_neo=%.3f\n",
            mapA$phi, mapA$pi0, mapA$theta, mapA$rate_loss, mapA$rate_neo))

# Need p (Geometric prior). Pull from sample.
mapA_p <- as.data.frame(resA$samples)$p[mapA$idx]

# (a) logPost at MAP tree with MAP nuisance (sanity recompute)
recompA_MAPtree <- compute_logpost(
  tree = mapA$tree, mkd = mkdA, model = modelA,
  kPrime = mapA$kPrime,
  rate_loss = mapA$rate_loss, rate_log_sd = mapA$rate_log_sd,
  rate_neo = mapA$rate_neo, p = mapA_p,
  phi = mapA$phi, pi0 = mapA$pi0, theta = mapA$theta,
  zMat = zMAP_A, ecology_aware = TRUE
)
cat(sprintf("recompute at MAP state: logPost = %.3f  (loglik %.3f, logprior %.3f)\n",
            recompA_MAPtree$logpost, recompA_MAPtree$loglik, recompA_MAPtree$logprior))
cat(sprintf("  drift vs chain log: %+.3f\n",
            recompA_MAPtree$logpost - mapA$sample_logpost))

# (b) Substitute truth tree, keep MAP nuisance
truthTreeRel <- truthTree
truthTreeRel$edge.length <- truthEdgeRel * sum(mapA$tree$edge.length)  # match TL scale of MAP
# Better: use truth TL too
truthTreeFull <- truthTree

cat("\n--- Counterfactual: truth tree at MAP nuisance (relative-edge match) ---\n")
truthAtMAPp <- compute_logpost(
  tree = truthTreeRel, mkd = mkdA, model = modelA,
  kPrime = mapA$kPrime,
  rate_loss = mapA$rate_loss, rate_log_sd = mapA$rate_log_sd,
  rate_neo = mapA$rate_neo, p = mapA_p,
  phi = mapA$phi, pi0 = mapA$pi0, theta = mapA$theta,
  zMat = zMAP_A, ecology_aware = TRUE
)
cat(sprintf("truth tree (TL=MAP TL=%.3f) + MAP nuisance: logPost = %.3f  (ll %.3f, lp %.3f)\n",
            sum(truthTreeRel$edge.length),
            truthAtMAPp$logpost, truthAtMAPp$loglik, truthAtMAPp$logprior))
cat(sprintf("DELTA (truth - chain MAP) = %+.3f nats\n",
            truthAtMAPp$logpost - mapA$sample_logpost))

cat("\n--- Counterfactual: truth tree FULL (truth TL=13.5) at MAP nuisance ---\n")
truthAtMAP_fullTL <- compute_logpost(
  tree = truthTreeFull, mkd = mkdA, model = modelA,
  kPrime = mapA$kPrime,
  rate_loss = mapA$rate_loss, rate_log_sd = mapA$rate_log_sd,
  rate_neo = mapA$rate_neo, p = mapA_p,
  phi = mapA$phi, pi0 = mapA$pi0, theta = mapA$theta,
  zMat = zMAP_A, ecology_aware = TRUE
)
cat(sprintf("truth tree (TL=%.3f) + MAP nuisance: logPost = %.3f  (ll %.3f, lp %.3f)\n",
            sum(truthTreeFull$edge.length),
            truthAtMAP_fullTL$logpost, truthAtMAP_fullTL$loglik,
            truthAtMAP_fullTL$logprior))
cat(sprintf("DELTA (truth tree truth TL - chain MAP) = %+.3f nats\n",
            truthAtMAP_fullTL$logpost - mapA$sample_logpost))

# Tree length trace: where did the chain converge on TL?
tlA <- as.data.frame(resA$samples)$tree_length
cat(sprintf("\nAware TL trace: min=%.2f  med=%.2f  max=%.2f  last100mean=%.2f  truth=%.2f\n",
            min(tlA), median(tlA), max(tlA),
            mean(tail(tlA, 100)), r$truthTL))

# Does the chain ever visit a tree containing the true bipartition?
trueSplit  <- c(paste0("A", 1:4), paste0("C", 1:4))
hasTrueAC <- vapply(resA$trees, function(tr) {
  cl <- ape::prop.part(tr)
  tips <- attr(cl, "labels")
  splitSet <- which(tips %in% trueSplit)
  any(vapply(cl, function(p) setequal(p, splitSet), logical(1)))
}, logical(1))
cat(sprintf("Aware: %d / %d samples contain the true AC bipartition (%.3f%%)\n",
            sum(hasTrueAC), length(hasTrueAC),
            100 * mean(hasTrueAC)))

# logPost at samples that DO contain AC (if any)
if (any(hasTrueAC)) {
  lpAC <- as.data.frame(resA$samples)$log_posterior[hasTrueAC]
  lpNoAC <- as.data.frame(resA$samples)$log_posterior[!hasTrueAC]
  cat(sprintf("  logPost: AC-containing  mean %.2f / max %.2f  (n=%d)\n",
              mean(lpAC), max(lpAC), sum(hasTrueAC)))
  cat(sprintf("  logPost: non-AC         mean %.2f / max %.2f  (n=%d)\n",
              mean(lpNoAC), max(lpNoAC), sum(!hasTrueAC)))
}

# (c) Unique topology count
nUniqueA <- length(unique(resA$samples[, "topo_hash"]))
cat(sprintf("\nAware unique topologies (out of %d samples): %d\n",
            nrow(as.matrix(resA$samples)), nUniqueA))
cat(sprintf("Most frequent topo_hash freq: %.3f\n",
            max(table(resA$samples[, "topo_hash"])) / nrow(as.matrix(resA$samples))))

# ---- Blind diagnostic ----
cat("\n--- BLIND chain ---\n")
resB <- r$resBlind
mkdB <- resB$data
modelB <- resB$model
if (nrow(as.matrix(resB$samples)) > 0) {
  mapB <- chain_map_state(resB, mkdB)
  cat(sprintf("MAP idx=%d  chain logPost = %.3f  (loglik %.3f)\n",
              mapB$idx, mapB$sample_logpost, mapB$sample_loglik))
  mapB_p <- as.data.frame(resB$samples)$p[mapB$idx]
  truthAtMAPpB <- compute_logpost(
    tree = truthTreeFull, mkd = mkdB, model = modelB,
    kPrime = mapB$kPrime,
    rate_loss = mapB$rate_loss, rate_log_sd = mapB$rate_log_sd,
    rate_neo = mapB$rate_neo, p = mapB_p,
    phi = NA, pi0 = NA, theta = NA, zMat = NULL,
    ecology_aware = FALSE
  )
  cat(sprintf("truth tree (TL=%.3f) + MAP nuisance: logPost = %.3f\n",
              sum(truthTreeFull$edge.length), truthAtMAPpB$logpost))
  cat(sprintf("DELTA (truth - chain MAP) = %+.3f nats\n",
              truthAtMAPpB$logpost - mapB$sample_logpost))
  nUniqueB <- length(unique(resB$samples[, "topo_hash"]))
  cat(sprintf("Blind unique topologies (out of %d): %d  most-freq: %.3f\n",
              nrow(as.matrix(resB$samples)), nUniqueB,
              max(table(resB$samples[, "topo_hash"])) /
                nrow(as.matrix(resB$samples))))
} else {
  cat("** Blind samples table is empty in this pilot file (length-0 matrix);\n")
  cat("** falling back to parsing blind-chain.log directly.\n")
  blog <- read.table(file.path(PILOT, "blind-chain.log"),
                     header = TRUE, sep = "\t", comment.char = "#")
  bIdx <- which.max(blog$log_posterior)
  cat(sprintf("Blind log MAP at sample row %d: log_post = %.3f  TL = %.3f\n",
              bIdx, blog$log_posterior[bIdx], blog$tree_length[bIdx]))
  cat(sprintf("Blind log: %d samples, %d unique topo_hash; final TL = %.3f\n",
              nrow(blog), length(unique(blog$topo_hash)),
              blog$tree_length[nrow(blog)]))
}

# ---- Acceptance rates (already computed) ----
cat("\n--- Acceptance rates (post-warmup) ---\n")
cat("\nAware:\n"); print(round(resA$acceptance, 3))
cat("\nBlind:\n"); print(round(resB$acceptance, 3))

# Flag: gibbs_spr/gibbs_subtree_swap should be 0 in eco mode (gated)
cat("\n*** SANITY: gibbs_spr accept rate in AWARE: ",
    round(resA$acceptance["gibbs_spr"], 4),
    "  -- should be 0 if eco gate is wired ***\n", sep="")
cat("*** SANITY: gibbs_subtree_swap accept rate in AWARE: ",
    round(resA$acceptance["gibbs_subtree_swap"], 4),
    "  -- should be 0 if eco gate is wired ***\n", sep="")

cat("\n=== Diagnostic done ===\n")
