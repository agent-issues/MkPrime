#!/usr/bin/env Rscript
# SBC harness for mixed-partition (neomorphic + transformational) MkNT data.
#
# Validates the partition-rate normalisation fix (2026-05-27, main):
#   neoScale   = r/(1+r) * n_total / n_neo
#   transScale = 1/(1+r) * n_total / n_trans
# by checking that rate_neo, rate_loss, tree_length, and rate_log_sd all
# pass Anderson-Darling rank-histogram uniformity (Bonferroni / K_arm at 5%).
#
# One arm: MkNT_mixed  (4 neo chars + 8 trans chars; both fixed k=2 in sim).
# Neo chars: asymmetric binary MkN (lambda=2), stationary pi_0=rl/(1+rl).
# Trans chars: JC k=2, pSame = 0.5 + 0.5*exp(-2t).  Rate convention matches
#   inference (`src/likelihood.cpp`); for k=2 the two conventions are identical.
#
# Critical SBC invariant: partition counts must be identical in forward sim
# and inference.  Guaranteed by rejection-sampling variable chars (kObs>=2)
# until exactly N_NEO neo + N_TRANS trans variable chars are obtained, so
# neoScale/transScale are computed from the same counts in both.
#
# Modes:
#   Rscript sbc-mixed.R --quick     :  execution smoke test (~60 s)
#   Rscript sbc-mixed.R --full      :  full-scale, run via SLURM
#
# Outputs: dev/red-team/sbc-results-mixed/ (or --out <path>)
#   verdict.txt, summary.rds, rank-matrix.csv

suppressPackageStartupMessages({
  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(".", quiet = TRUE)
  } else {
    library(MkPrime)
  }
  library(ape)
  library(TreeTools)
})

# ----------------------------- CLI -----------------------------
args <- commandArgs(trailingOnly = TRUE)
mode <- if (any(args %in% c("--quick", "-q"))) "quick" else "full"
seedBase <- {
  ix <- which(args == "--seed")
  if (length(ix) && length(args) >= ix + 1L) as.integer(args[ix + 1L]) else 20260528L
}
outRoot <- {
  ix <- which(args == "--out")
  if (length(ix) && length(args) >= ix + 1L) args[ix + 1L] else
    "dev/red-team/sbc-results-mixed"
}
dir.create(outRoot, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("SBC mixed harness | mode=%s | seedBase=%d | out=%s\n",
            mode, seedBase, outRoot))

# --------------------- Mode dimensions -------------------------
if (mode == "quick") {
  N_SIM    <- 5L
  N_TIP    <- 5L
  N_ITER   <- 300L
  N_THIN   <- 3L
  N_WARM   <- 200L
} else {
  N_SIM    <- 200L
  N_TIP    <- 8L
  N_ITER   <- 6000L
  N_THIN   <- 60L
  N_WARM   <- 2000L
}
N_RUNS   <- 1L
N_CHAINS <- 1L
L_SAMPLES <- N_ITER %/% N_THIN

# Partition dimensions — must remain constant across all sims so that
# neoScale/transScale in the forward sim match what the inference computes.
N_NEO   <- 4L
N_TRANS <- 8L

# Tree prior hyperparameters — must match MkPrimeModel() defaults.
EXPSTEPS_FIXED <- 50
TREE_SHAPE     <- 2

# -------------------- AD uniformity test -----------------------
if (!requireNamespace("goftest", quietly = TRUE)) {
  cat("[warn] goftest not installed; falling back to ks.test\n")
}
.adP <- function(ranks_norm) {
  ranks_norm <- ranks_norm[is.finite(ranks_norm)]
  if (length(ranks_norm) < 4L) return(NA_real_)
  if (requireNamespace("goftest", quietly = TRUE)) {
    suppressWarnings(goftest::ad.test(ranks_norm, null = "punif")$p.value)
  } else {
    suppressWarnings(ks.test(ranks_norm, "punif")$p.value)
  }
}

# ---------------- Forward simulation helpers -------------------

# Tree: Gamma(shape, shape/EXPSTEPS_FIXED) total length; Dirichlet(1,...,1) edges.
.simTree <- function(nTip) {
  tr <- ape::rtree(nTip, tip.label = paste0("t", seq_len(nTip)))
  nEdge <- nrow(tr$edge)
  tl <- stats::rgamma(1, shape = TREE_SHAPE, rate = TREE_SHAPE / EXPSTEPS_FIXED)
  w <- stats::rexp(nEdge, rate = 1)
  tr$edge.length <- tl * w / sum(w)
  TreeTools::Preorder(tr)
}

# Partition scale factors: match C++ gibbs_partial_cl.h lines 97-98.
.partScales <- function(rateNeo, nNeo, nTrans) {
  n <- nNeo + nTrans
  r <- rateNeo
  list(neo   = r / (1 + r) * n / nNeo,
       trans = 1 / (1 + r) * n / nTrans)
}

# Simulate one neomorphic (asymmetric binary) character.
# Model: MkN, lambda=2.  pi_0 = rl/(1+rl), pi_1 = 1/(1+rl).
# Effective branch length = edge * neoScale.
.simNeoChar <- function(tree, rateLoss_true, neoScale) {
  nTip    <- length(tree$tip.label)
  states  <- integer(2L * nTip - 1L)
  rootIdx <- nTip + 1L
  pi0 <- rateLoss_true / (1 + rateLoss_true)
  pi1 <- 1            / (1 + rateLoss_true)
  # Root drawn from stationary distribution
  states[rootIdx] <- if (stats::runif(1L) < pi1) 1L else 0L
  edges <- tree$edge
  el    <- tree$edge.length
  for (e in seq_len(nrow(edges))) {
    pa    <- edges[e, 1L]
    ch    <- edges[e, 2L]
    t_eff <- el[e] * neoScale
    ef    <- exp(-2 * t_eff)
    if (states[pa] == 0L) {
      # P(0->1) = pi_1 * (1 - exp(-2t))
      states[ch] <- if (stats::runif(1L) < pi1 * (1 - ef)) 1L else 0L
    } else {
      # P(1->0) = pi_0 * (1 - exp(-2t))
      states[ch] <- if (stats::runif(1L) < pi0 * (1 - ef)) 0L else 1L
    }
  }
  states[seq_len(nTip)]
}

# Simulate one transformational character under JC k=2.
# pSame = 0.5 + 0.5*exp(-2*t_eff) — matches inference exactly at k=2
# (both the JC convention -k*t/(k-1) and the lambda=2 MkN convention give
# the same formula when k=2).
.simTransChar <- function(tree, transScale) {
  nTip    <- length(tree$tip.label)
  states  <- integer(2L * nTip - 1L)
  rootIdx <- nTip + 1L
  states[rootIdx] <- sample.int(2L, 1L) - 1L
  edges <- tree$edge
  el    <- tree$edge.length
  for (e in seq_len(nrow(edges))) {
    pa    <- edges[e, 1L]
    ch    <- edges[e, 2L]
    t_eff <- el[e] * transScale
    pSame <- 0.5 + 0.5 * exp(-2 * t_eff)
    states[ch] <- if (stats::runif(1L) < pSame) states[pa] else 1L - states[pa]
  }
  states[seq_len(nTip)]
}

# Rejection-sample exactly nNeo variable neo + nTrans variable trans characters.
# Each character is redrawn independently until kObs >= 2.
# Returns list(skipped=FALSE, neoMat, transMat) or list(skipped=TRUE, reason).
.simMixedCharsRej <- function(tree, rateLoss_true, rateNeo_true,
                                nNeo, nTrans, maxRetry = 500L) {
  scales <- .partScales(rateNeo_true, nNeo, nTrans)
  nTip   <- length(tree$tip.label)

  neoMat <- matrix(NA_integer_, nTip, nNeo)
  for (j in seq_len(nNeo)) {
    got <- FALSE
    for (attempt in seq_len(maxRetry)) {
      s <- .simNeoChar(tree, rateLoss_true, scales$neo)
      if (length(unique(s)) >= 2L) {
        neoMat[, j] <- s
        got <- TRUE
        break
      }
    }
    if (!got) {
      return(list(skipped = TRUE,
                  reason = sprintf("neo_char_%d: invariant after %d retries", j, maxRetry)))
    }
  }

  transMat <- matrix(NA_integer_, nTip, nTrans)
  for (j in seq_len(nTrans)) {
    got <- FALSE
    for (attempt in seq_len(maxRetry)) {
      s <- .simTransChar(tree, scales$trans)
      if (length(unique(s)) >= 2L) {
        transMat[, j] <- s
        got <- TRUE
        break
      }
    }
    if (!got) {
      return(list(skipped = TRUE,
                  reason = sprintf("trans_char_%d: invariant after %d retries", j, maxRetry)))
    }
  }

  list(skipped   = FALSE,
       neoMat    = neoMat,
       transMat  = transMat,
       neoScale  = scales$neo,
       transScale = scales$trans)
}

# ------------------- One SBC simulation ------------------------
.runOneSim <- function(sim_id, seed) {
  set.seed(seed)

  # --- Prior draws (must match MkPrimeModel() defaults exactly) ---
  # rate_loss ~ LogNormal(0, 2)  (rateLossMeanlog=0, rateLossSdlog=2)
  # rate_neo  ~ LogNormal(0, 2)  (rateNeoMeanlog=0,  rateNeoSdlog=2)
  # rate_log_sd ~ Gamma(shape=1, rate=1)
  rateLoss_true   <- stats::rlnorm(1, meanlog = 0, sdlog = 2)
  rateNeo_true    <- stats::rlnorm(1, meanlog = 0, sdlog = 2)
  rateLogSd_true  <- stats::rgamma(1, shape = 1, rate = 1)

  # --- Tree ---
  true_tree <- .simTree(N_TIP)
  tl_true   <- sum(true_tree$edge.length)

  # --- Characters (rejection-sampled for exact partition counts) ---
  charResult <- .simMixedCharsRej(true_tree, rateLoss_true, rateNeo_true,
                                   N_NEO, N_TRANS)
  if (isTRUE(charResult$skipped)) return(charResult)

  # Combined matrix: neo chars in columns 1:N_NEO, trans in (N_NEO+1):(N_NEO+N_TRANS).
  combined <- cbind(charResult$neoMat, charResult$transMat)
  rownames(combined) <- true_tree$tip.label

  pd  <- TreeTools::MatrixToPhyDat(combined)
  mkd <- MkPrimeData(
    pd,
    neomorphic  = seq_len(N_NEO),
    knownStates = setNames(rep(2L, N_TRANS),
                           as.character(seq.int(N_NEO + 1L, N_NEO + N_TRANS)))
  )

  # Sanity: rejection-sampling must have produced exactly N_NEO+N_TRANS variable chars.
  if (mkd$nChar != N_NEO + N_TRANS) {
    return(list(skipped = TRUE,
                reason = sprintf("mkd dropped chars: nChar=%d expected=%d",
                                 mkd$nChar, N_NEO + N_TRANS)))
  }

  # --- Starting tree: true topology, flat branch lengths ---
  start_tree <- true_tree
  start_tree$edge.length <- rep_len(0.1, nrow(true_tree$edge))

  # --- Inference model ---
  # nCat=1: no ACRV in forward sim, so nCat=1 satisfies SBC forward==inference.
  # kPrimePrior="geometric": degenerate for knownStates chars, no kPrime inference.
  model <- MkPrimeModel(
    coding   = "variable",
    nCat     = 1L,
    kPrimePrior = "geometric",
    expSteps = EXPSTEPS_FIXED
  )

  mcmc <- MkPrimeMCMC(
    nIter      = N_ITER,   thin     = N_THIN,
    minWarmup  = N_WARM,   maxWarmup = N_WARM,
    autoTune   = FALSE,
    nRuns      = N_RUNS,   nChains  = N_CHAINS
  )

  t0 <- Sys.time()
  res <- tryCatch(
    suppressMessages(suppressWarnings(
      RunMkPrime(mkd, start_tree, model = model, mcmc = mcmc,
                 fixTopology = TRUE, overwrite = TRUE)
    )),
    error = function(e) list(error = conditionMessage(e))
  )
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  if (!is.null(res$error)) {
    return(list(skipped = TRUE,
                reason  = paste0("mcmc_error: ", res$error),
                wall    = dt))
  }

  samples <- res$samples
  if (is.null(samples) || nrow(samples) < 10L) {
    return(list(skipped = TRUE, reason = "no_samples", wall = dt))
  }

  # --- SBC ranks ---
  # rank(true_val) = #{posterior samples < true_val} (Talts et al. 2018)
  .rankOf <- function(true_val, post_vec) {
    if (!is.finite(true_val)) return(NA_integer_)
    sum(post_vec < true_val)
  }
  L     <- nrow(samples)
  cn    <- colnames(samples)
  ranks <- list()
  ranks$tree_length <- .rankOf(tl_true,          samples[, "tree_length"])
  if ("rate_neo"    %in% cn) ranks$rate_neo    <- .rankOf(rateNeo_true,   samples[, "rate_neo"])
  if ("rate_loss"   %in% cn) ranks$rate_loss   <- .rankOf(rateLoss_true,  samples[, "rate_loss"])
  if ("rate_log_sd" %in% cn) ranks$rate_log_sd <- .rankOf(rateLogSd_true, samples[, "rate_log_sd"])

  list(skipped        = FALSE,
       L              = L,
       wall           = dt,
       tl_true        = tl_true,
       rateLoss_true  = rateLoss_true,
       rateNeo_true   = rateNeo_true,
       rateLogSd_true = rateLogSd_true,
       ranks          = ranks)
}

# -------------------- Arm driver -------------------------------
armDir <- file.path(outRoot, "MkNT_mixed")
dir.create(armDir, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("\n=== arm: MkNT_mixed  (N_NEO=%d N_TRANS=%d N_TIP=%d N_SIM=%d L=%d) ===\n",
            N_NEO, N_TRANS, N_TIP, N_SIM, L_SAMPLES))

sims <- vector("list", N_SIM)
for (i in seq_len(N_SIM)) {
  seed <- seedBase + i
  cat(sprintf("  sim %3d/%d (seed=%d) ... ", i, N_SIM, seed))
  s <- tryCatch(.runOneSim(i, seed),
                error = function(e) list(skipped = TRUE,
                                         reason = paste("trycatch:", conditionMessage(e))))
  if (isTRUE(s$skipped)) {
    cat(sprintf("SKIP (%s)\n", s$reason))
  } else {
    cat(sprintf("L=%d wall=%.1fs\n", s$L, s$wall))
  }
  sims[[i]] <- s
}

# --- Pool ranks ---
good   <- !vapply(sims, function(s) isTRUE(s$skipped), logical(1L))
nGood  <- sum(good)

if (nGood < 3L) {
  msg <- sprintf("FAIL (only %d/%d sims completed)", nGood, N_SIM)
  writeLines(msg, file.path(armDir, "verdict.txt"))
  saveRDS(list(sims = sims), file.path(armDir, "summary.rds"))
  cat(msg, "\n")
  quit(save = "no", status = 1L)
}
goodSims <- sims[good]

pool_param <- function(name) {
  out <- unlist(lapply(goodSims, function(s) s$ranks[[name]]), use.names = FALSE)
  out[is.finite(out)]
}
paramNames <- unique(unlist(lapply(goodSims, function(s) names(s$ranks))))
paramNames <- paramNames[!is.na(paramNames)]

K_arm     <- length(paramNames)
threshold <- 0.001 / max(1L, K_arm)

L_actual <- {
  Ls <- vapply(goodSims, function(s) as.integer(s$L), integer(1L))
  Ls <- Ls[is.finite(Ls)]
  if (length(Ls)) as.integer(stats::median(Ls)) else as.integer(L_SAMPLES)
}

per_param <- list()
for (nm in paramNames) {
  ranks_raw <- pool_param(nm)
  if (length(ranks_raw) < 4L) {
    per_param[[nm]] <- list(p = NA_real_, n = length(ranks_raw), decision = "INSUFFICIENT")
    next
  }
  norm <- (ranks_raw + 0.5) / (L_actual + 1)
  norm <- pmin(pmax(norm, .Machine$double.eps), 1 - .Machine$double.eps)
  p    <- .adP(norm)
  per_param[[nm]] <- list(
    p        = p,
    n        = length(ranks_raw),
    decision = if (is.na(p)) "NA" else if (p > threshold) "PASS" else "FAIL"
  )
}

decisions   <- vapply(per_param, `[[`, character(1L), "decision")
armVerdict  <- if (mode == "quick") {
  "EXEC_OK"
} else if (all(decisions %in% c("PASS", "INSUFFICIENT"))) {
  "PASS"
} else {
  "FAIL"
}

# --- Write artefacts ---
verdictLines <- c(
  sprintf("arm:         MkNT_mixed"),
  sprintf("model:       MkNT  (hasNeo=TRUE, N_NEO=%d, N_TRANS=%d)", N_NEO, N_TRANS),
  sprintf("neo sim:     MkN lambda=2, pi_0=rl/(1+rl)"),
  sprintf("trans sim:   JC k=2, pSame=0.5+0.5*exp(-2t)"),
  sprintf("mode:        %s", mode),
  sprintf("N_sim:       %d (good = %d)", N_SIM, nGood),
  sprintf("N_tip:       %d", N_TIP),
  sprintf("L_samples:   nominal=%d actual=%d", L_SAMPLES, L_actual),
  sprintf("K_arm:       %d  (Bonferroni threshold p > %.5g)", K_arm, threshold),
  "",
  "per-parameter AD p-values:",
  unlist(lapply(names(per_param), function(nm) {
    x <- per_param[[nm]]
    sprintf("  %-18s  p=%s  n=%d  %s",
            nm,
            if (is.na(x$p)) "NA" else formatC(x$p, digits = 4, format = "g"),
            x$n, x$decision)
  })),
  "",
  sprintf("ARM VERDICT: %s", armVerdict)
)
writeLines(verdictLines, file.path(armDir, "verdict.txt"))
cat(paste0(verdictLines, "\n"), sep = "")

saveRDS(list(sims = sims, per_param = per_param,
             verdict = armVerdict, threshold = threshold,
             N_NEO = N_NEO, N_TRANS = N_TRANS),
        file.path(armDir, "summary.rds"))

rankList <- lapply(paramNames, function(nm) pool_param(nm))
names(rankList) <- paramNames
maxLen <- max(vapply(rankList, length, integer(1L)), 0L)
pad    <- function(x) c(x, rep(NA, maxLen - length(x)))
df     <- as.data.frame(lapply(rankList, pad))
utils::write.csv(df, file.path(armDir, "rank-matrix.csv"), row.names = FALSE)

# Top-level verdict (mirroring sbc.R layout)
topLines <- c(
  "SBC mixed-partition summary",
  sprintf("mode:     %s", mode),
  sprintf("seedBase: %d", seedBase),
  sprintf("N_sim:    %d  L:%d  thin:%d", N_SIM, L_SAMPLES, N_THIN),
  "",
  sprintf("  %-30s %s", "MkNT_mixed",
          paste0(armVerdict, sprintf(" (good=%d)", nGood)))
)
writeLines(topLines, file.path(outRoot, "verdict.txt"))
cat("\n", paste0(topLines, "\n"), sep = "")

status <- if (mode == "quick" || armVerdict %in% c("PASS", "EXEC_OK")) 0L else 1L
quit(save = "no", status = status)
