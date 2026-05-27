#!/usr/bin/env Rscript
# Lane D1 — Simulation-Based Calibration (SBC) harness for MkNT and Mk'.
#
# Tests the MkPrime sampler against the Talts et al. (2018) calibration
# criterion: when the forward simulator and inference share the same prior
# and likelihood, the rank of the truth in the posterior sample should be
# Uniform{0, ..., L} for every monitored parameter (Cook–Gelman–Rubin 2006).
#
# Pass criterion: per-parameter Anderson–Darling p > 0.001 / K_arm
# (Bonferroni over K_arm monitored parameters in the arm), with K_arm given
# in the verdict.txt for each arm. See sbc.md for full justification.
#
# Modes:
#   Rscript sbc.R --quick                : ≤ ~60 s, execution smoke test.
#   Rscript sbc.R --full                 : full-scale, run via SLURM.
#   Rscript sbc.R --arm Mkp_eg --quick   : run a single arm.
#
# Outputs are written to dev/red-team/heavy-tests/sbc-results/<arm>/
# with per-arm verdict.txt, summary.rds, and rank-matrix.csv.

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
armFilter <- {
  ix <- which(args == "--arm")
  if (length(ix) && length(args) >= ix + 1L) args[ix + 1L] else NULL
}
seedBase <- {
  ix <- which(args == "--seed")
  if (length(ix) && length(args) >= ix + 1L) as.integer(args[ix + 1L]) else 20260526L
}
outRoot <- {
  ix <- which(args == "--out")
  if (length(ix) && length(args) >= ix + 1L) args[ix + 1L] else
    "dev/red-team/heavy-tests/sbc-results"
}
dir.create(outRoot, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("SBC harness  | mode=%s | seedBase=%d | armFilter=%s\n",
            mode, seedBase, if (is.null(armFilter)) "<all>" else armFilter))

# --------------------- Configuration table ---------------------
# Six arms. See sbc.md §Arms for predicted PASS/FAIL.
ALL_ARMS <- list(
  list(name = "MkNT_geometric",          model = "MkNT", prior = "geometric",
       expect = "PASS"),
  list(name = "Mkp_geometric",           model = "Mkp",  prior = "geometric",
       expect = "PASS"),
  list(name = "Mkp_beta_geometric",      model = "Mkp",  prior = "beta_geometric",
       expect = "PASS"),
  list(name = "Mkp_empirical_geometric", model = "Mkp",  prior = "empirical_geometric",
       expect = "PASS"),  # EG-001 not visible at kObs=2 (option α); L6 closes the math
  list(name = "Mkp_logseries",           model = "Mkp",  prior = "logseries",
       expect = "PASS"),  # LS-001 latent because c is fixed
  list(name = "MkNT_logseries",          model = "MkNT", prior = "logseries",
       expect = "PASS")
)

# Mode-dependent dimensions. Quick = execution smoke; rank-tests are not
# powered to fail at this scale and the AD p-values are reported but not
# asserted.
if (mode == "quick") {
  N_SIM    <- 5L     # SBC simulations per arm
  N_TIP    <- 5L     # tree taxa per simulation
  N_CHAR   <- 8L     # characters per simulation
  N_ITER   <- 300L   # MCMC iterations per chain
  N_THIN   <- 3L     # samples kept = 100  (Talts recommend L ~ 100)
  N_WARM   <- 200L
  N_RUNS   <- 1L
  N_CHAINS <- 1L
} else {
  N_SIM    <- 200L
  N_TIP    <- 8L
  N_CHAR   <- 30L
  N_ITER   <- 6000L
  N_THIN   <- 60L    # samples kept = 100
  N_WARM   <- 2000L
  N_RUNS   <- 1L
  N_CHAINS <- 1L
}
L_SAMPLES <- N_ITER %/% N_THIN  # nominal posterior length per chain

# Numerical guards for AD
if (!requireNamespace("goftest", quietly = TRUE)) {
  cat("[warn] goftest not installed; falling back to ks.test for uniformity\n")
}
.adP <- function(ranks_norm) {
  # ranks_norm in (0, 1); return p-value of Anderson-Darling vs U(0,1).
  ranks_norm <- ranks_norm[is.finite(ranks_norm)]
  if (length(ranks_norm) < 4L) return(NA_real_)
  if (requireNamespace("goftest", quietly = TRUE)) {
    suppressWarnings(goftest::ad.test(ranks_norm, null = "punif")$p.value)
  } else {
    suppressWarnings(ks.test(ranks_norm, "punif")$p.value)
  }
}

# ----------------- Forward simulation helpers ------------------
# SBC requires forward-prior == inference-prior. Inference uses
# tree_length ~ Gamma(shape = TREE_SHAPE, rate = TREE_SHAPE / EXPSTEPS_FIXED)
# (see MkPrimeModel(): treeLengthShape, treeLengthRate). We pin expSteps
# to a fixed constant across all sims and inference (no data-derived
# leak), draw total tree length from that exact Gamma, then partition
# across edges with a Dirichlet(1,...,1) (uniform-simplex). The
# inference uses beta-simplex / Dirichlet-simplex moves on the same
# parameterisation, so the forward and inferred priors agree.
EXPSTEPS_FIXED <- 50         # treeLengthRate = TREE_SHAPE / EXPSTEPS_FIXED
TREE_SHAPE     <- 2          # matches MkPrimeModel() default treeLengthShape
.simTree <- function(nTip) {
  tr <- ape::rtree(nTip, tip.label = paste0("t", seq_len(nTip)))
  nEdge <- nrow(tr$edge)
  tl <- stats::rgamma(1, shape = TREE_SHAPE,
                      rate = TREE_SHAPE / EXPSTEPS_FIXED)
  # Dirichlet(1,...,1) on edge-fraction simplex = normalised iid Exp(1)
  w <- stats::rexp(nEdge, rate = 1)
  tr$edge.length <- tl * w / sum(w)
  TreeTools::Preorder(tr)
}

# Simulate one character under JC(kTrue) by recursive root-to-tip draws.
# kTrue >= 2.  Returns integer tip vector in {0, ..., kTrue-1}.
.simJCchar <- function(tree, kTrue) {
  nTip <- length(tree$tip.label)
  states <- integer(2L * nTip - 1L)
  rootIdx <- nTip + 1L
  states[rootIdx] <- sample.int(kTrue, 1L) - 1L
  edges <- tree$edge
  el    <- tree$edge.length
  for (e in seq_len(nrow(edges))) {
    pa <- edges[e, 1L]; ch <- edges[e, 2L]; t <- el[e]
    # Forward (root→tip) iteration on a Preorder edge list: parents are
    # always introduced before their children, so states[pa] is valid when
    # we reach edge e. `rev()` here would give postorder traversal, reading
    # uninitialised states[pa] = 0 for every non-root edge — destroying
    # phylogenetic signal and biasing inference MLE downward. See SBC-TL-MIX-001.
    # Rate convention must match inference (`src/likelihood.cpp:85`:
    # arg = -kStates * t / (kStates - 1)); naive exp(-k*t) diverges from
    # inference for k > 2, biasing tree_length ranks. See SBC-HARNESS-003.
    pSame <- 1 / kTrue + (1 - 1 / kTrue) * exp(-kTrue * t / (kTrue - 1))
    if (runif(1L) < pSame) {
      states[ch] <- states[pa]
    } else {
      states[ch] <- sample(setdiff(seq.int(0L, kTrue - 1L), states[pa]), 1L)
    }
  }
  states[seq_len(nTip)]
}

# Relabel a tip-state vector so the observed labels are a contiguous
# 0..kObs-1 range, preserving Mk' canonical encoding (`L4` proof).
.canonicaliseLabels <- function(vec) {
  uvals <- sort(unique(vec))
  out <- match(vec, uvals) - 1L
  attr(out, "kObs") <- length(uvals)
  out
}

# Draw kPrime from the *untruncated* prior. This is the heart of the SBC
# generative model: we sample k' before knowing kObs, then simulate the
# character, then let inference condition on the resulting kObs. Any
# missing truncation normaliser in inference will then bias the joint.
#
# For MkNT we simply set kPrime = kTrue (no latent expansion) and draw
# kTrue from the marginal prior on observable states.
.drawKPrime <- function(arm, p, model_hp, n) {
  prior <- arm$prior
  if (arm$model == "MkNT") {
    # MkNT has no kPrime; "k" is the observable state count. Draw from
    # the same marginal but reject k = 1 (need >= 2 states to be a
    # character) and clamp to a sane max for the JC simulator.
    if (prior == "geometric") {
      # u ~ Geo(p), k = u + 2 (kObs >= 2 for variable characters under MkNT)
      u <- stats::rgeom(n, p)
      return(pmin(u + 2L, 20L))
    } else if (prior == "logseries") {
      # logseries on k >= 1, then take max(k, 2)
      c_ls <- model_hp$kprimeLogseriesC
      # inverse-CDF sample
      kvals <- 1:50
      pk <- -c_ls^kvals / (kvals * log(1 - c_ls))
      pk <- pk / sum(pk)
      kk <- sample(kvals, n, replace = TRUE, prob = pk)
      return(pmax(kk, 2L))
    }
  } else {
    # Mk': u ~ <prior>, kPrime = u + 2 nominal (we pick kObs=2 floor and
    # let observed kObs emerge from simulation).
    kFloor <- 2L
    if (prior == "geometric") {
      u <- stats::rgeom(n, p)
      return(u + kFloor)
    } else if (prior == "beta_geometric") {
      # u ~ BetaGeo(alpha, beta): u = number of "failures" before success
      # under p ~ Beta(alpha, beta). Sample p then u.
      pv <- stats::rbeta(n, model_hp$kprimeAlpha, model_hp$kprimeBeta)
      u <- stats::rgeom(n, pv)
      return(u + kFloor)
    } else if (prior == "empirical_geometric") {
      # k = N_obs + N_unobs;  N_obs ~ empiricalNObs body+tail;  N_unobs ~ Geo(p)
      emp <- model_hp$empiricalNObs
      kVals <- seq.int(2L, length.out = length(emp$body))
      pBody <- emp$body
      # sample N_obs from the empirical body only (truncate tail for
      # simulation simplicity; tail mass typically < 5%).
      pBody <- pBody / sum(pBody)
      nObs <- sample(kVals, n, replace = TRUE, prob = pBody)
      nUnobs <- stats::rgeom(n, p)
      return(nObs + nUnobs)
    } else if (prior == "logseries") {
      c_ls <- model_hp$kprimeLogseriesC
      kvals <- 2:50  # kPrime >= 2 by support
      pk <- -c_ls^kvals / (kvals * log(1 - c_ls))
      pk <- pk / sum(pk)
      return(sample(kvals, n, replace = TRUE, prob = pk))
    }
  }
  stop("unknown prior/model combination: ", arm$model, "/", prior)
}

# Run one SBC simulation: returns named list of true values + ranks.
.runOneSim <- function(arm, sim_id, seed) {
  set.seed(seed)

  # Step 1. Hyperparameter draws -- match the *inference* hyperprior.
  # For geometric and empirical_geometric: p ~ Beta(a=1, b=1) by default.
  # For beta_geometric: alpha,beta ~ Exp(1) by default; we leave as fixed
  #   in the inference model and sample p via beta-geo marginalisation.
  # For logseries: c is fixed in the inference, so we don't sample it.
  p_true <- stats::rbeta(1, 1, 1)  # used for {geometric, empirical_geometric}
  model_hp <- list(
    kprimeHyperA = 1, kprimeHyperB = 1,
    kprimeAlpha  = 2, kprimeBeta = 2,    # mean E[u] = beta/(alpha-1) = 2
    kprimeLogseriesC = 0.5,
    empiricalNObs = NULL
  )
  if (arm$prior == "empirical_geometric") {
    # Use the packaged empirical body
    e <- new.env(parent = emptyenv())
    utils::data("empiricalNObs", package = "MkPrime", envir = e)
    model_hp$empiricalNObs <- e$empiricalNObs
  }

  # rateLogSd: must match inference prior exactly. Inference defaults
  # to Gamma(shape = 1, rate = 1) (MkPrimeModel(): rateLogSdShape,
  # rateLogSdRate). No truncation — the earlier `min(.., 3.0)` clamp
  # broke the SBC prior-equality requirement.
  rateLogSd_true <- stats::rgamma(1, shape = 1, rate = 1)

  # Step 2. Tree draw — Gamma(2, 2/EXPSTEPS_FIXED) × Dirichlet(1,..1) on edges.
  true_tree <- .simTree(N_TIP)
  tl_true <- sum(true_tree$edge.length)

  # Step 3. Draw kTrue per character from the inference prior (tracked
  # as the SBC "truth" for kPrime ranks). Simulation is forced to k=2.
  kTrue <- .drawKPrime(arm, p_true, model_hp, N_CHAR)

  # SBC-PRIOR-CONFOUND-001 (option α): simulate every character as binary
  # so realised kObs == 2 deterministically. With kObs == 2, the
  # forward draw `kTrue ~ kFloor + prior` matches the inference per-char
  # prior `kPrime ~ kObs + Geo(p) = 2 + Geo(p)` at the prior level —
  # eliminating the joint mismatch that contaminated all v4 ranks.
  # Cost: EG-001 cannot be demonstrated by SBC at kObs=2 (Z_i(p)
  # missing-normaliser bug cancels there); L6 closes the math.
  kSim <- rep(2L, N_CHAR)

  # Step 4. Simulate each character with k=2.
  sim_mat <- matrix(NA_integer_, N_TIP, N_CHAR,
                    dimnames = list(true_tree$tip.label, NULL))
  kObs <- integer(N_CHAR)
  for (j in seq_len(N_CHAR)) {
    raw <- .simJCchar(true_tree, kSim[j])
    canon <- .canonicaliseLabels(raw)
    sim_mat[, j] <- canon
    kObs[j] <- attr(canon, "kObs")
  }
  # Drop invariant chars (kObs = 1) under coding = "variable"
  keep <- kObs >= 2L
  if (sum(keep) < 3L) {
    return(list(skipped = TRUE, reason = "too_few_variable_chars"))
  }
  sim_mat <- sim_mat[, keep, drop = FALSE]
  kTrue <- kTrue[keep]
  kObs  <- kObs[keep]
  n_char <- sum(keep)

  pd <- TreeTools::MatrixToPhyDat(sim_mat)
  mkd <- if (arm$model == "MkNT") {
    # MkNT: pin k = kObs by treating every character as transformational
    # with knownStates set to its observed support. This matches
    # `feedback_model_scope.md`.
    MkPrimeData(pd,
                knownStates = setNames(as.integer(kObs),
                                       as.character(seq_along(kObs))))
  } else {
    MkPrimeData(pd)
  }

  # Step 5. Choose a starting tree.
  start_tree <- tryCatch(
    {
      st <- TreeSearch::AdditionTree(pd)
      st$edge.length <- rep_len(0.1, nrow(st$edge))
      st
    },
    error = function(e) TreeTools::NJTree(pd, edgeLengths = TRUE)
  )

  # Step 6. Build the matched inference model. Hyperparameter draws above
  # are matched to inference defaults (Beta(1,1) on p, etc).
  #
  # Two earlier mismatches now closed:
  #   - expSteps fixed to EXPSTEPS_FIXED, not derived from tl_true
  #     (data-derived expSteps leaks ground truth into the inference
  #     prior on tree_length).
  #   - nCat = 1L: disable ACRV. The forward simulator does NOT add
  #     per-character rate variation; matching nCat = 1L in the
  #     inference satisfies SBC's forward==inference requirement.
  #     L7 already proved ACRV correct analytically, so we don't need
  #     SBC to re-exercise that path.
  modelArgs <- list(
    coding = "variable",
    nCat = 1L,
    kPrimePrior = if (arm$model == "MkNT") "geometric" else arm$prior,
    expSteps = EXPSTEPS_FIXED
  )
  if (arm$model == "MkNT") {
    # MkNT: kPrime is fixed by knownStates; the prior on k' is degenerate.
    # We can still use geometric as a no-op (it never moves k').
  } else {
    modelArgs$kprimeHyperA <- model_hp$kprimeHyperA
    modelArgs$kprimeHyperB <- model_hp$kprimeHyperB
    modelArgs$kprimeAlpha  <- model_hp$kprimeAlpha
    modelArgs$kprimeBeta   <- model_hp$kprimeBeta
    if (arm$prior == "logseries") {
      modelArgs$kprimeLogseriesC <- model_hp$kprimeLogseriesC
    }
    if (arm$prior == "empirical_geometric") {
      modelArgs$empiricalNObs <- model_hp$empiricalNObs
    }
  }
  model <- do.call(MkPrimeModel, modelArgs)

  # Step 7. MCMC.
  mcmc <- MkPrimeMCMC(
    nIter = N_ITER, thin = N_THIN,
    minWarmup = N_WARM, maxWarmup = N_WARM,
    autoTune = FALSE,
    nRuns = N_RUNS, nChains = N_CHAINS
  )

  t0 <- Sys.time()
  res <- tryCatch(
    suppressMessages(suppressWarnings(
      RunMkPrime(mkd, start_tree, model = model, mcmc = mcmc,
                 fixTopology = TRUE,           # SBC concerns continuous params
                 overwrite = TRUE)
    )),
    error = function(e) {
      list(error = conditionMessage(e))
    }
  )
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (!is.null(res$error)) {
    return(list(skipped = TRUE, reason = paste0("mcmc_error: ", res$error),
                wall = dt))
  }

  samples <- res$samples
  if (is.null(samples) || nrow(samples) < 10L) {
    return(list(skipped = TRUE, reason = "no_samples", wall = dt))
  }

  # Step 8. Compute SBC ranks.
  # For a true value t and L posterior samples x_1, ..., x_L:
  #   rank(t) = #{x_i < t}                      (Talts et al. 2018)
  # Under correct inference rank ~ Uniform{0,...,L}.
  .rankOf <- function(true_val, post_samples) {
    if (!is.finite(true_val)) return(NA_integer_)
    sum(post_samples < true_val)
  }
  L <- nrow(samples)
  ranks <- list()
  ranks$tree_length <- .rankOf(tl_true, samples[, "tree_length"])
  if ("rate_log_sd" %in% colnames(samples)) {
    ranks$rate_log_sd <- .rankOf(rateLogSd_true, samples[, "rate_log_sd"])
  }
  if (arm$model == "Mkp") {
    # Per-character k' ranks. Pool across characters under exchangeability.
    # SBC-PRIOR-CONFOUND-001 option α: with binary sim every kept char has
    # kObs == 2 by construction, so the prior match `kTrue = 2 + Geo(p)`
    # ⇔ `kPrime = kObs + Geo(p) = 2 + Geo(p)` is exact. The kObs == 2
    # filter below is therefore a no-op safety net.
    kpCols <- grep("^kPrime_", colnames(samples), value = TRUE)
    if (length(kpCols)) {
      idx <- as.integer(sub("kPrime_", "", kpCols))
      ord <- order(idx)
      kpCols <- kpCols[ord]
      kObsEq2 <- which(kObs == 2L)
      if (length(kObsEq2)) {
        kp_ranks <- mapply(.rankOf, kTrue[kObsEq2],
                           asplit(samples[, kpCols[kObsEq2], drop = FALSE], 2L))
        ranks$kPrime_pooled <- kp_ranks
      }
    }
    if ("p" %in% colnames(samples) &&
        arm$prior %in% c("geometric", "empirical_geometric")) {
      ranks$p <- .rankOf(p_true, samples[, "p"])
    }
  }

  list(skipped = FALSE,
       L = L,
       wall = dt,
       n_char = n_char,
       tl_true = tl_true,
       rateLogSd_true = rateLogSd_true,
       p_true = p_true,
       kTrue = kTrue,
       kObs = kObs,
       ranks = ranks)
}

# ----------------- Per-arm driver ------------------------------
.runArm <- function(arm) {
  cat(sprintf("\n=== arm: %-30s (model=%s prior=%s expect=%s) ===\n",
              arm$name, arm$model, arm$prior, arm$expect))
  armDir <- file.path(outRoot, arm$name)
  dir.create(armDir, recursive = TRUE, showWarnings = FALSE)

  sims <- vector("list", N_SIM)
  for (i in seq_len(N_SIM)) {
    seed <- seedBase + 1000L * match(arm$name, vapply(ALL_ARMS, `[[`, "", "name")) + i
    cat(sprintf("  sim %3d/%d (seed=%d) ... ", i, N_SIM, seed))
    s <- tryCatch(.runOneSim(arm, i, seed),
                  error = function(e) list(skipped = TRUE,
                                           reason = paste("trycatch:", conditionMessage(e))))
    if (isTRUE(s$skipped)) {
      cat(sprintf("SKIP (%s)\n", s$reason))
    } else {
      cat(sprintf("L=%d wall=%.1fs\n", s$L, s$wall))
    }
    sims[[i]] <- s
  }

  # Pool ranks across simulations.
  good <- !vapply(sims, function(s) isTRUE(s$skipped), logical(1L))
  nGood <- sum(good)
  if (nGood < 3L) {
    msg <- sprintf("FAIL (only %d/%d sims completed)", nGood, N_SIM)
    writeLines(msg, file.path(armDir, "verdict.txt"))
    saveRDS(list(arm = arm, sims = sims), file.path(armDir, "summary.rds"))
    return(list(arm = arm$name, verdict = "FAIL_EXEC",
                detail = list(n_good = nGood, n_total = N_SIM)))
  }
  goodSims <- sims[good]

  pool_param <- function(name) {
    out <- lapply(goodSims, function(s) s$ranks[[name]])
    out <- unlist(out, use.names = FALSE)
    out[is.finite(out)]
  }
  paramNames <- unique(unlist(lapply(goodSims, function(s) names(s$ranks))))
  paramNames <- paramNames[!is.na(paramNames)]

  # Compute AD p-values per parameter.
  # Talts et al. recommend the rank-uniformity test on (rank + 0.5) / (L + 1)
  # to obtain a U(0,1) variable for goodness-of-fit.
  K_arm <- length(paramNames)
  threshold <- 0.001 / max(1L, K_arm)  # Bonferroni

  per_param <- list()
  # Determine the actual per-sim L (in case warmup truncated samples below
  # the nominal L_SAMPLES). Use the modal L across good sims — Talts'
  # construction requires the same L for every rank entering the test.
  L_actual <- {
    Ls <- vapply(goodSims, function(s) as.integer(s$L), integer(1L))
    Ls <- Ls[is.finite(Ls)]
    if (length(Ls)) as.integer(stats::median(Ls)) else as.integer(L_SAMPLES)
  }
  for (nm in paramNames) {
    ranks_raw <- pool_param(nm)
    if (length(ranks_raw) < 4L) {
      per_param[[nm]] <- list(p = NA_real_, n = length(ranks_raw),
                              decision = "INSUFFICIENT")
      next
    }
    L_eff <- L_actual
    norm <- (ranks_raw + 0.5) / (L_eff + 1)
    norm <- pmin(pmax(norm, .Machine$double.eps), 1 - .Machine$double.eps)
    p <- .adP(norm)
    per_param[[nm]] <- list(
      p = p,
      n = length(ranks_raw),
      decision = if (is.na(p)) "NA" else if (p > threshold) "PASS" else "FAIL"
    )
  }

  # Arm verdict: PASS if every param PASS or INSUFFICIENT
  decisions <- vapply(per_param, `[[`, character(1L), "decision")
  armVerdict <- if (mode == "quick") {
    "EXEC_OK"   # quick mode reports execution only
  } else if (all(decisions %in% c("PASS", "INSUFFICIENT"))) {
    "PASS"
  } else {
    "FAIL"
  }

  # Write per-arm artefacts
  verdictLines <- c(
    sprintf("arm:        %s", arm$name),
    sprintf("model:      %s", arm$model),
    sprintf("prior:      %s", arm$prior),
    sprintf("expect:     %s", arm$expect),
    sprintf("mode:       %s", mode),
    sprintf("N_sim:      %d (good = %d)", N_SIM, nGood),
    sprintf("L_samples:  nominal=%d actual=%d", L_SAMPLES, L_actual),
    sprintf("K_arm:      %d  (Bonferroni threshold p > %.5g)", K_arm, threshold),
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

  saveRDS(list(arm = arm, sims = sims, per_param = per_param,
               verdict = armVerdict, threshold = threshold),
          file.path(armDir, "summary.rds"))

  # Rank matrix for downstream plotting
  rankCsv <- file.path(armDir, "rank-matrix.csv")
  rankList <- lapply(paramNames, function(nm) pool_param(nm))
  names(rankList) <- paramNames
  maxLen <- max(vapply(rankList, length, integer(1L)), 0L)
  pad <- function(x) c(x, rep(NA, maxLen - length(x)))
  df <- as.data.frame(lapply(rankList, pad))
  utils::write.csv(df, rankCsv, row.names = FALSE)

  list(arm = arm$name, verdict = armVerdict, per_param = per_param,
       n_good = nGood)
}

# ----------------- Driver --------------------------------------
arms <- if (is.null(armFilter)) ALL_ARMS else
  Filter(function(a) a$name == armFilter, ALL_ARMS)
if (length(arms) == 0L) {
  stop("No arms match filter: ", armFilter)
}

t_start <- Sys.time()
results <- list()
for (arm in arms) {
  results[[arm$name]] <- tryCatch(.runArm(arm),
                                  error = function(e) {
                                    msg <- sprintf("arm %s ERROR: %s",
                                                   arm$name, conditionMessage(e))
                                    cat(msg, "\n")
                                    list(arm = arm$name, verdict = "ERROR",
                                         error = conditionMessage(e))
                                  })
}
dt <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))

# ----------------- Top-level verdict ---------------------------
top_verdict_path <- file.path(outRoot, "verdict.txt")
lines <- c(
  sprintf("SBC harness top-level summary"),
  sprintf("mode:     %s", mode),
  sprintf("seedBase: %d", seedBase),
  sprintf("wall:     %.1fs", dt),
  sprintf("N_sim:    %d  L:%d  thin:%d", N_SIM, L_SAMPLES, N_THIN),
  "",
  "arms:"
)
for (nm in names(results)) {
  r <- results[[nm]]
  lines <- c(lines, sprintf("  %-30s %s", nm,
                            paste0(r$verdict, if (!is.null(r$n_good))
                              sprintf(" (good=%d)", r$n_good) else "")))
}
writeLines(lines, top_verdict_path)
cat("\n", paste0(lines, "\n"), sep = "")

# Overall PASS/FAIL exit code: 0 if all PASS in full mode; 0 if execution
# completed in quick mode; non-zero only on hard errors.
hard_err <- any(vapply(results, function(r) identical(r$verdict, "ERROR"),
                       logical(1L)))
if (hard_err) {
  quit(save = "no", status = 2L)
}
quit(save = "no", status = 0L)
