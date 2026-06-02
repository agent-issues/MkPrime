# =====================================================================
# relabel-marginal-k-driver.R
#
# Decisive numerical audit of the Mk' relabelling term
#   mk_prime_relabel_log(k, kObs) = log(k!/(k-kObs)!)
# added to the per-k log-likelihood at src/mcmc.cpp:4092 (and the R
# reference path R/likelihood.R:262-266).
#
# QUESTION: does the inference's relabel-inclusive per-k likelihood
# P_inf(y | k, tree) track the forward generative model's probability
# of the SAME canonical pattern y as a function of k?  The RATIO across
# k is what drives the k-marginal and hence the geometric-p posterior.
#
# Forward generative model (ground truth, from sbc.R::.simJCchar +
# .canonicaliseLabels):
#   root ~ Uniform{0..k-1}
#   JC(k) transitions, pSame = 1/k + (1-1/k) exp(-k t /(k-1))
#   tip vector CANONICALISED (relabel distinct states 0..kObs-1 in
#   order of appearance).
#
# We estimate P(y | variable, k) = P(y, variable | k) / P(variable | k)
# by brute-force Monte Carlo (simulate, canonicalise, count exact
# matches to y among VARIABLE sims), because the inference call uses
# coding="variable" (constant-site ascertainment), so the comparison
# must condition on variable patterns too.
#
# Compare to exp(per-k LL) from MkpLogLikelihood(..., coding="variable",
# relabel=TRUE), which is the exact relabel-inclusive per-k path used
# by the geometric arm.
#
# OUTPUT: csv table of MC vs inference per-k probabilities and the
# k-to-k ratios, written to relabel-marginal-k-results/.
# =====================================================================

suppressMessages(pkgload::load_all(getwd(), quiet = TRUE))

set.seed(20260529L)

OUTDIR <- file.path(getwd(), "dev", "red-team", "numerical",
                    "relabel-marginal-k-results")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ---- forward simulator (copied verbatim from sbc.R, same conventions) ----
.simJCchar <- function(tree, kTrue) {
  nTip <- length(tree$tip.label)
  states <- integer(2L * nTip - 1L)
  rootIdx <- nTip + 1L
  states[rootIdx] <- sample.int(kTrue, 1L) - 1L
  edges <- tree$edge
  el    <- tree$edge.length
  for (e in seq_len(nrow(edges))) {
    pa <- edges[e, 1L]; ch <- edges[e, 2L]; t <- el[e]
    pSame <- 1 / kTrue + (1 - 1 / kTrue) * exp(-kTrue * t / (kTrue - 1))
    if (runif(1L) < pSame) {
      states[ch] <- states[pa]
    } else {
      states[ch] <- sample(setdiff(seq.int(0L, kTrue - 1L), states[pa]), 1L)
    }
  }
  states[seq_len(nTip)]
}
.canonicalise <- function(vec) {
  uvals <- sort(unique(vec))
  match(vec, uvals) - 1L
}

# ---- fixed tree: preorder, moderate branch lengths so kObs=2,3 occur ----
nTip <- 6L
tr <- TreeTools::PectinateTree(nTip)          # caterpillar; arbitrary fixed topo
tr$edge.length <- rep(0.30, nrow(tr$edge))    # moderate -> multistate patterns
tr <- TreeTools::Preorder(tr)
tipLab <- tr$tip.label

# ---- helper: build single-character MkPrimeData from a canonical vector ----
# kObs is set by the number of distinct symbols actually present.
make_mkd <- function(canonVec) {
  syms <- as.character(canonVec)
  m <- matrix(syms, ncol = 1L, dimnames = list(tipLab, NULL))
  pd <- TreeTools::MatrixToPhyDat(m)
  MkPrime::MkPrimeData(pd)
}

# ---- inference per-k likelihood (relabel-inclusive, variable coding) ----
inf_logp <- function(mkd1, k) {
  tryCatch(
    MkpLogLikelihood(tr, mkd1, kPrime = as.integer(k),
                     rate_log_sd = 0, nCat = 1L,
                     coding = "variable", relabel = TRUE),
    error = function(e) { message("inf err k=", k, ": ", conditionMessage(e)); NA_real_ }
  )
}

# ---- brute-force MC: P(y | variable, k) ----
# Simulate N draws under JC(k); canonicalise; among VARIABLE draws count
# exact matches to the target canonical pattern y.
mc_logp <- function(yCanon, k, N) {
  yKey <- paste(yCanon, collapse = ",")
  nVar <- 0L
  nMatch <- 0L
  # batch for speed
  B <- 50000L
  done <- 0L
  while (done < N) {
    nb <- min(B, N - done)
    for (i in seq_len(nb)) {
      raw <- .simJCchar(tr, k)
      cv  <- .canonicalise(raw)
      nd  <- length(unique(cv))
      if (nd >= 2L) {           # variable
        nVar <- nVar + 1L
        if (paste(cv, collapse = ",") == yKey) nMatch <- nMatch + 1L
      }
    }
    done <- done + nb
  }
  list(p = nMatch / nVar, nMatch = nMatch, nVar = nVar)
}

# =====================================================================
# Run the comparison for a few target patterns.
# =====================================================================
# Target canonical patterns (over 6 tips). Choose patterns reachable
# with moderate frequency so MC has signal.
targets <- list(
  kObs2_split = c(0,0,0,1,1,1),   # kObs = 2, balanced
  kObs2_31    = c(0,0,0,1,0,1),   # kObs = 2, less balanced
  kObs3       = c(0,0,1,1,2,2)    # kObs = 3
)

N_MC <- 4e6L      # MC draws per (pattern, k)
results <- list()

for (tname in names(targets)) {
  y <- targets[[tname]]
  kObs <- length(unique(y))
  mkd1 <- make_mkd(y)
  cat(sprintf("\n=== target %s  kObs=%d  pattern=(%s) ===\n",
              tname, kObs, paste(y, collapse = "")))
  kGrid <- kObs + (0:3)

  # inference per-k LL
  llInf <- vapply(kGrid, function(k) inf_logp(mkd1, k), numeric(1))
  pInf  <- exp(llInf)

  # MC per-k probability (conditioned on variable)
  pMC <- numeric(length(kGrid))
  seMC <- numeric(length(kGrid))
  for (j in seq_along(kGrid)) {
    k <- kGrid[j]
    r <- mc_logp(y, k, N_MC)
    pMC[j] <- r$p
    # binomial SE on conditional probability
    seMC[j] <- sqrt(r$p * (1 - r$p) / r$nVar)
    cat(sprintf("  k=%d  inf=%.6e  MC=%.6e (SE %.2e, nMatch=%d nVar=%d)\n",
                k, pInf[j], r$p, seMC[j], r$nMatch, r$nVar))
  }

  df <- data.frame(
    pattern = tname, kObs = kObs, k = kGrid,
    ll_inf = llInf, p_inf = pInf,
    p_mc = pMC, se_mc = seMC
  )
  # ratios relative to k = kObs (the quantity that drives the k-marginal)
  df$ratio_inf <- df$p_inf / df$p_inf[1]
  df$ratio_mc  <- df$p_mc  / df$p_mc[1]
  df$ratio_se  <- df$ratio_mc * sqrt((df$se_mc / df$p_mc)^2 +
                                     (df$se_mc[1] / df$p_mc[1])^2)
  results[[tname]] <- df

  cat("  --- ratios vs k=kObs (drives the k-marginal / p) ---\n")
  print(df[, c("k", "ratio_inf", "ratio_mc", "ratio_se")], row.names = FALSE)
}

allDf <- do.call(rbind, results)
write.csv(allDf, file.path(OUTDIR, "relabel-mc-vs-inference.csv"),
          row.names = FALSE)
cat("\nWrote", file.path(OUTDIR, "relabel-mc-vs-inference.csv"), "\n")

# ---- summary verdict line ----
allDf$z_ratio <- (allDf$ratio_inf - allDf$ratio_mc) / allDf$ratio_se
cat("\n=== max |z| of (inference ratio - MC ratio) across all (pattern,k>kObs) ===\n")
sub <- allDf[allDf$k > allDf$kObs, ]
cat(sprintf("  max |z| = %.2f  (z>~3 would indicate a real discrepancy)\n",
            max(abs(sub$z_ratio), na.rm = TRUE)))
