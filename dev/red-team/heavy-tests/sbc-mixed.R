#!/usr/bin/env Rscript
# SBC harness for mixed-partition (neomorphic + transformational) MkNT data.
#
# Validates the partition-rate normalisation fix (2026-05-27, main):
#   neoScale   = r/(1+r) * n_total / n_neo
#   transScale = 1/(1+r) * n_total / n_trans
# by checking that rate_neo, rate_loss and tree_length all
# pass Anderson-Darling rank-histogram uniformity (Bonferroni: FWER 0.001 / K_arm).
# rate_log_sd is not ranked: nCat = 1, so it never reaches the likelihood.
#
# One arm: MkNT_mixed  (N_NEO neo + N_TRANS trans chars; both fixed k=2 in sim).
# Neo chars: asymmetric binary MkN, mean rate 1 (src/mkn_rates.h): total rate
#   lambda=(1+rl)^2/(2 rl), stationary pi_0=rl/(1+rl).
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
#
# Exit status: 0 PASS / EXEC_OK; 1 FAIL (incl. INSUFFICIENT, a missing
# parameter, or too many skipped sims).

suppressPackageStartupMessages({
  # Race-safe loader: pkgload compiles into src/ in place, so concurrent
  # array tasks corrupt each other's objects unless the tree is pre-built (#15).
  source("dev/red-team/heavy-tests/load-mkprime.R")
  LoadMkPrime(".")
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
# Stamped into the verdict so aggregate-verdicts.R can tell this run's file from
# one a previous submission left in the same outRoot.
runId <- {
  ix <- which(args == "--run-id")
  if (length(ix) && length(args) >= ix + 1L) args[ix + 1L] else
    sprintf("local-%s-%d", format(Sys.time(), "%Y%m%d%H%M%S"), Sys.getpid())
}

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
  N_TIP    <- 16L
  N_ITER   <- 12000L
  N_THIN   <- 60L
  N_WARM   <- 4000L
}
N_RUNS   <- 1L
N_CHAINS <- 1L
L_SAMPLES <- N_ITER %/% N_THIN

# Partition dimensions — must remain constant across all sims so that
# neoScale/transScale in the forward sim match what the inference computes.
# Full mode: 16 tips x 100 characters at mean tree length 1.4, the informative
# regime of marginal-k/T-SBC-marginal-geometric.R. The old 8-tip / expSteps = 50
# regime left tree_length at ~prior sd (posterior sd / prior sd median 0.95).
if (mode == "quick") {
  N_NEO   <- 4L
  N_TRANS <- 8L
} else {
  N_NEO   <- 30L
  N_TRANS <- 70L
}

# Gate thresholds. A sim skipped by the rejection cap or an MCMC error is
# dropped non-randomly in theta, so too many skips invalidate the verdict.
MIN_GOOD_FRAC <- 0.9
MIN_RANKS     <- 50L
EXPECTED_PARAMS <- c("tree_length", "rate_neo", "rate_loss")
FWER <- 0.001

# Tree prior hyperparameters — must match MkPrimeModel() defaults.
EXPSTEPS_FIXED <- if (mode == "quick") 50 else 1.4
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
# Model: MkN with mean rate 1, i.e. total rate lambda = (1+rl)^2/(2 rl)
# (src/mkn_rates.h).  pi_0 = rl/(1+rl), pi_1 = 1/(1+rl).
# Effective branch length = edge * neoScale.
.simNeoChar <- function(tree, rateLoss_true, neoScale) {
  nTip    <- length(tree$tip.label)
  states  <- integer(2L * nTip - 1L)
  rootIdx <- nTip + 1L
  pi0 <- rateLoss_true / (1 + rateLoss_true)
  pi1 <- 1            / (1 + rateLoss_true)
  lambda <- (1 + rateLoss_true)^2 / (2 * rateLoss_true)
  # Root drawn from stationary distribution
  states[rootIdx] <- if (stats::runif(1L) < pi1) 1L else 0L
  edges <- tree$edge
  el    <- tree$edge.length
  for (e in seq_len(nrow(edges))) {
    pa    <- edges[e, 1L]
    ch    <- edges[e, 2L]
    t_eff <- el[e] * neoScale
    ef    <- exp(-lambda * t_eff)
    if (states[pa] == 0L) {
      # P(0->1) = pi_1 * (1 - exp(-lambda t))
      states[ch] <- if (stats::runif(1L) < pi1 * (1 - ef)) 1L else 0L
    } else {
      # P(1->0) = pi_0 * (1 - exp(-lambda t))
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

# Gate: PASS only if enough sims completed and every expected parameter has
# enough ranks and a uniform rank histogram. A missing or under-ranked
# parameter fails; it never passes by omission.
MixedVerdict <- function(sims, expected, L, fwer, minGoodFrac, minRanks, adP) {
  skipped <- vapply(sims, function(s) isTRUE(s$skipped), logical(1L))
  goodSims <- sims[!skipped]
  nGood <- length(goodSims)
  threshold <- fwer / length(expected)
  Ls <- vapply(goodSims, function(s) as.integer(s$L), integer(1L))
  Ls <- Ls[is.finite(Ls)]
  lActual <- if (length(Ls)) as.integer(stats::median(Ls)) else as.integer(L)
  perParam <- lapply(expected, function(nm) {
    ranks <- unlist(lapply(goodSims, function(s) s$ranks[[nm]]), use.names = FALSE)
    ranks <- ranks[is.finite(ranks)]
    if (!length(ranks)) {
      return(list(p = NA_real_, n = 0L, decision = "MISSING"))
    }
    if (length(ranks) < minRanks) {
      return(list(p = NA_real_, n = length(ranks), decision = "INSUFFICIENT"))
    }
    norm <- (ranks + 0.5) / (lActual + 1)
    norm <- pmin(pmax(norm, .Machine$double.eps), 1 - .Machine$double.eps)
    p <- adP(norm)
    decision <- if (is.na(p)) "NA" else if (p > threshold) "PASS" else "FAIL"
    list(p = p, n = length(ranks), decision = decision)
  })
  names(perParam) <- expected
  decisions <- vapply(perParam, `[[`, character(1L), "decision")
  enough <- nGood >= minGoodFrac * length(sims)
  verdict <- if (enough && all(decisions == "PASS")) "PASS" else "FAIL"
  list(verdict = verdict, perParam = perParam, nGood = nGood,
       nSkipped = sum(skipped), enough = enough, threshold = threshold,
       L = lActual)
}

# ------------------- One SBC simulation ------------------------
.runOneSim <- function(sim_id, seed) {
  set.seed(seed)

  # --- Prior draws (must match MkPrimeModel() defaults exactly) ---
  # rate_loss ~ LogNormal(0, 2)  (rateLossMeanlog=0, rateLossSdlog=2)
  # rate_neo  ~ LogNormal(0, 2)  (rateNeoMeanlog=0,  rateNeoSdlog=2)
  rateLoss_true   <- stats::rlnorm(1, meanlog = 0, sdlog = 2)
  rateNeo_true    <- stats::rlnorm(1, meanlog = 0, sdlog = 2)

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
  if ("rate_neo"  %in% cn) ranks$rate_neo  <- .rankOf(rateNeo_true,  samples[, "rate_neo"])
  if ("rate_loss" %in% cn) ranks$rate_loss <- .rankOf(rateLoss_true, samples[, "rate_loss"])

  list(skipped        = FALSE,
       L              = L,
       wall           = dt,
       tl_true        = tl_true,
       rateLoss_true  = rateLoss_true,
       rateNeo_true   = rateNeo_true,
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

# --- Pool ranks and gate ---
V <- MixedVerdict(sims, EXPECTED_PARAMS, L_SAMPLES, FWER, MIN_GOOD_FRAC,
                  MIN_RANKS, .adP)
nGood      <- V$nGood
per_param  <- V$perParam
threshold  <- V$threshold
L_actual   <- V$L
K_arm      <- length(EXPECTED_PARAMS)
paramNames <- names(per_param)
goodSims   <- sims[!vapply(sims, function(s) isTRUE(s$skipped), logical(1L))]
pool_param <- function(name) {
  out <- unlist(lapply(goodSims, function(s) s$ranks[[name]]), use.names = FALSE)
  out[is.finite(out)]
}
armVerdict <- if (mode == "quick") "EXEC_OK" else V$verdict

# --- Write artefacts ---
verdictLines <- c(
  sprintf("arm:         MkNT_mixed"),
  sprintf("model:       MkNT  (hasNeo=TRUE, N_NEO=%d, N_TRANS=%d)", N_NEO, N_TRANS),
  sprintf("neo sim:     MkN mean rate 1, lambda=(1+rl)^2/(2rl), pi_0=rl/(1+rl)"),
  sprintf("trans sim:   JC k=2, pSame=0.5+0.5*exp(-2t)"),
  sprintf("mode:        %s", mode),
  sprintf("N_sim:       %d (good = %d, skipped = %d; need >= %.0f%% good: %s)",
          N_SIM, nGood, V$nSkipped, 100 * MIN_GOOD_FRAC,
          if (V$enough) "ok" else "NOT MET"),
  sprintf("N_tip:       %d", N_TIP),
  sprintf("L_samples:   nominal=%d actual=%d", L_SAMPLES, L_actual),
  sprintf("K_arm:       %d  (Bonferroni FWER %g: p > %.5g)", K_arm, FWER, threshold),
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
  sprintf("run_id:   %s", runId),
  sprintf("mode:     %s", mode),
  sprintf("seedBase: %d", seedBase),
  sprintf("N_sim:    %d  L:%d  thin:%d", N_SIM, L_SAMPLES, N_THIN),
  "",
  sprintf("  %-30s %s", "MkNT_mixed",
          paste0(armVerdict, sprintf(" (good=%d)", nGood)))
)
# Named for the arm, not `verdict.txt`: this script shares `outRoot` with
# sbc.R, so a bare aggregate name collides both with sbc.R's and with any
# concurrent array task's (#16).
writeLines(topLines, file.path(outRoot, "verdict-MkNT_mixed.txt"))
cat("\n", paste0(topLines, "\n"), sep = "")

status <- if (mode == "quick" || armVerdict %in% c("PASS", "EXEC_OK")) 0L else 1L
quit(save = "no", status = status)
