# =====================================================================
# marginal-k-pbias-driver.R
#
# DECISIVE localisation of the geometric-p SBC bias (p_true = 0.3 -->
# posterior ~ 0.9, 16/20 rank-0). The per-character marginal likelihood
# has already been exonerated piece-by-piece:
#   * relabel term log(k!/(k-kObs)!)  -- relabel-marginal-k-driver.R
#     showed ratio_inf tracks ratio_mc to within MC SE through k=kObs+3
#     (incl. the collapse-path candidates k >= kObs+2).
#   * root prior 1/k', JC convention, a_k flat across k -- already ruled
#     out in the brief.
#   * truncation: forward caps kTrue at K_MAX_PRIOR=30 with pile-up;
#     inference sums k = kObs..kObs+49 (WIDER) with no pile-up. At p=0.3
#     the dropped/piled tail mass is (1-p)^28 = 0.7^28 ~ 1.6e-5 --
#     numerically negligible. NOT the driver.
#
# This driver asks the one question that a per-k ratio check CANNOT
# answer (ratios cancel any p-independent per-char factor): does the
# PRODUCTION marginal likelihood, summed over k with the geometric
# weights w_k(p), peak at the true p?
#
#   logL(p) = sum_i logSumExp_{k=kObs_i..kObs_i+KCAND-1}
#               [ rawLL_i(k) + ascCorr_i(k) + relabel(k,kObs_i)
#                 + log w_k(p) ]
#
# with, under priorVariant = "unconditional" (Model A, what the harness
# uses to match the forward):
#   log w_k(p) = log p + (k - 2) * log(1 - p)
#
# We build rawLL+asc+relabel via the PRODUCTION per-k path
# (MkpLogLikelihood(..., coding="variable", relabel=TRUE)) -- the exact
# routine the geometric arm calls -- so the test inherits the real
# collapse-kernel dispatch for k >= kObs+2.
#
# THREE EXPERIMENTS:
#   E1. Kernel equivalence: test_persite_collapsed(k) vs
#       test_persite_uncollapsed(k) on the same tree/pattern for
#       k = kObs+2..kObs+4 (the collapse path has no production
#       cross-check). Expect agreement to ~1e-12 by JC lumpability.
#   E2. Marginal logL(p) recovery: forward-simulate N_CHAR chars at
#       p_true = 0.3 on the harness tree (Model A, kTrue=2+Geo, JC,
#       canonicalise), then evaluate logL(p) on a p-grid two ways:
#         (a) Model A weights  log p + (k-2) log(1-p)  -- harness/inference
#         (b) Model B weights  log p + (k-kObs) log(1-p) -- sampled-k R prior
#       Report argmax_p and the curve.
#   E3. Repeat E2 across several p_true to map the bias function
#       p_hat(p_true).
#
# OUTPUT: csv tables + a short console verdict under
#   dev/red-team/numerical/marginal-k-pbias-results/
# =====================================================================

suppressMessages(pkgload::load_all(getwd(), quiet = TRUE))
suppressMessages({ library(ape); library(TreeTools) })

set.seed(20260529L)

OUTDIR <- file.path(getwd(), "dev", "red-team", "numerical",
                    "marginal-k-pbias-results")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

KCAND <- 50L     # inference's kMaxKprimeCand (src/mcmc_state.h:378)

# ---- forward simulator (verbatim from T-SBC-marginal-geometric.R) ----
.simJCchar <- function(tree, kTrue) {
  nTip <- length(tree$tip.label)
  states <- integer(2L * nTip - 1L)
  rootIdx <- nTip + 1L
  states[rootIdx] <- sample.int(kTrue, 1L) - 1L
  edges <- tree$edge; el <- tree$edge.length
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
.canon <- function(v) {
  uvals <- sort(unique(v)); match(v, uvals) - 1L
}

logSumExp <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(-Inf)
  m <- max(x); m + log(sum(exp(x - m)))
}

# ---- fixed tree: the harness "full mode" regime (16 tips, mean TL ~1.4) ----
# Use a moderate caterpillar so kObs in {2,3,4} occur. Keep it fixed and
# deterministic so the experiment is reproducible.
N_TIP <- 16L
tr <- TreeTools::PectinateTree(N_TIP)
set.seed(101L)
# branch lengths: tree_length ~ 1.4 spread over edges via Dirichlet(1)
nE <- nrow(tr$edge)
rel <- as.numeric(rgamma(nE, 1)); rel <- rel / sum(rel)
TL_TRUE <- 1.4
tr$edge.length <- rel * TL_TRUE
tr <- TreeTools::Preorder(tr)
tipLab <- tr$tip.label

# ---- build a 1-char MkPrimeData from a canonical vector ----
make_mkd <- function(canonVec) {
  m <- matrix(as.character(canonVec), ncol = 1L, dimnames = list(tipLab, NULL))
  MkPrime::MkPrimeData(TreeTools::MatrixToPhyDat(m))
}

# ---- production per-k log-lik (relabel + ascertainment, variable coding) ----
inf_logp_k <- function(mkd1, k) {
  tryCatch(
    MkpLogLikelihood(tr, mkd1, kPrime = as.integer(k),
                     rate_log_sd = 0, nCat = 1L,
                     coding = "variable", relabel = TRUE),
    error = function(e) NA_real_)
}

# =====================================================================
# E1. Kernel equivalence (collapse path has no production cross-check)
# =====================================================================
cat("\n===== E1: collapsed vs uncollapsed kernel =====\n")
parent <- tr$edge[, 1]; child <- tr$edge[, 2]; el <- tr$edge.length
# representative canonical patterns
e1_targets <- list(
  kObs2 = .canon(c(0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1)),
  kObs3 = .canon(c(0,0,0,0,0,1,1,1,1,1,2,2,2,2,2,2)),
  kObs4 = .canon(c(0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3))
)
e1_rows <- list()
for (tn in names(e1_targets)) {
  y <- e1_targets[[tn]]; kObs <- length(unique(y))
  ts <- matrix(as.integer(y), ncol = 1L)
  for (k in (kObs + 2L):(kObs + 4L)) {
    ll_un <- test_persite_uncollapsed(parent, child, el, ts, k, c(1.0))
    ll_co <- test_persite_collapsed(parent, child, el, ts, k, kObs, c(1.0))
    e1_rows[[length(e1_rows) + 1L]] <- data.frame(
      pattern = tn, kObs = kObs, k = k,
      ll_uncollapsed = ll_un, ll_collapsed = ll_co,
      abs_diff = abs(ll_un - ll_co))
  }
}
e1 <- do.call(rbind, e1_rows)
print(e1)
write.csv(e1, file.path(OUTDIR, "E1-kernel-equivalence.csv"), row.names = FALSE)
cat(sprintf("E1 max |diff| = %.3e\n", max(e1$abs_diff)))

# =====================================================================
# E2 / E3. Marginal logL(p) recovery.
# =====================================================================
# To make logL(p) cheap and exact, precompute rawLL+asc+relabel ONCE per
# UNIQUE (pattern, k); then logL(p) is pure arithmetic over the p-grid.

P_GRID <- seq(0.02, 0.98, by = 0.01)

# weights: Model A = log p + (k-2) log(1-p); Model B = log p + (k-kObs) log(1-p)
logL_curve <- function(charList, model = c("A", "B")) {
  model <- match.arg(model)
  sapply(P_GRID, function(p) {
    lp <- log(p); l1mp <- log1p(-p)
    tot <- 0.0
    for (ch in charList) {
      kObs <- ch$kObs
      ks <- ch$ks                      # candidate k values
      base <- ch$base                  # rawLL+asc+relabel per k
      w <- if (model == "A") lp + (ks - 2) * l1mp
           else              lp + (ks - kObs) * l1mp
      tot <- tot + logSumExp(base + w)
      if (!is.finite(tot)) break
    }
    tot
  })
}

simulate_chars <- function(p_true, n_char, K_MAX = 30L) {
  u <- rgeom(n_char, p_true)
  kTrue <- pmin(2L + u, K_MAX)
  out <- list()
  for (j in seq_len(n_char)) {
    raw <- .simJCchar(tr, kTrue[j])
    cv  <- .canon(raw)
    if (length(unique(cv)) >= 2L) out[[length(out) + 1L]] <- cv  # variable only
  }
  out
}

# Cache of base[k] keyed by (kObs, pattern-string) to avoid recompute.
base_cache <- new.env(parent = emptyenv())
char_to_eval <- function(cv) {
  kObs <- length(unique(cv))
  key <- paste0(kObs, ":", paste(cv, collapse = ","))
  if (!is.null(base_cache[[key]])) return(base_cache[[key]])
  mkd1 <- make_mkd(cv)
  ks <- kObs + (0:(KCAND - 1L))
  base <- vapply(ks, function(k) inf_logp_k(mkd1, k), numeric(1))
  ok <- is.finite(base)
  obj <- list(kObs = kObs, ks = ks[ok], base = base[ok])
  base_cache[[key]] <- obj
  obj
}

cat("\n===== E2: marginal logL(p) recovery at p_true = 0.3 =====\n")
N_CHAR <- 400L
set.seed(303L)
chars_raw <- simulate_chars(0.30, N_CHAR)
cat(sprintf("simulated %d variable chars; kObs table:\n", length(chars_raw)))
print(table(vapply(chars_raw, function(v) length(unique(v)), integer(1))))
charList <- lapply(chars_raw, char_to_eval)

curveA <- logL_curve(charList, "A")
curveB <- logL_curve(charList, "B")
pHatA <- P_GRID[which.max(curveA)]
pHatB <- P_GRID[which.max(curveB)]
cat(sprintf("p_true = 0.30 | argmax logL  ModelA = %.3f  ModelB = %.3f\n",
            pHatA, pHatB))
e2 <- data.frame(p = P_GRID, logL_modelA = curveA, logL_modelB = curveB)
write.csv(e2, file.path(OUTDIR, "E2-logL-curve-p030.csv"), row.names = FALSE)

# =====================================================================
# E3. p_hat(p_true) bias map.
# =====================================================================
cat("\n===== E3: p_hat(p_true) bias map =====\n")
p_trues <- c(0.10, 0.20, 0.30, 0.50, 0.70)
e3_rows <- list()
for (pt in p_trues) {
  set.seed(1000L + round(pt * 100))
  cr <- simulate_chars(pt, N_CHAR)
  cl <- lapply(cr, char_to_eval)
  cA <- logL_curve(cl, "A"); cB <- logL_curve(cl, "B")
  e3_rows[[length(e3_rows) + 1L]] <- data.frame(
    p_true = pt, n_var = length(cr),
    pHat_modelA = P_GRID[which.max(cA)],
    pHat_modelB = P_GRID[which.max(cB)])
}
e3 <- do.call(rbind, e3_rows)
print(e3)
write.csv(e3, file.path(OUTDIR, "E3-pbias-map.csv"), row.names = FALSE)

cat("\n===== DONE. Results in", OUTDIR, "=====\n")
