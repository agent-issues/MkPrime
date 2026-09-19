# One cell of the gibbsSpr TRUE/FALSE sensitivity experiment.
# Usage: Rscript run-cell.R <tree_idx> <arm:on|off> <seed> <nIter> <outdir>
setTimeLimit(elapsed = 1800, transient = FALSE)
source(file.path(Sys.getenv("GSPR_DIR", "."), "common.R"))

a <- commandArgs(trailingOnly = TRUE)
tree_idx <- as.integer(a[1]); arm <- a[2]; seed <- as.integer(a[3])
nIter <- as.integer(a[4]); outdir <- a[5]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
outfile <- file.path(outdir, sprintf("t%02d_%s_s%03d.rds", tree_idx, arm, seed))
if (file.exists(outfile)) { cat("skip (exists)", outfile, "\n"); quit(save = "no") }

r   <- LoadRep(tree_idx, 1L)
mkd <- MkPrimeData(r$pd)

# Ground truth re-aligned to the sampler's (lexical) character order.
nexf <- sort(list.files(file.path(GT_ROOT, sprintf("tree_%02d/rep_01", tree_idx)),
                        "^chr[0-9]+\\.nex$"))
num  <- as.integer(sub("^chr([0-9]+)\\.nex$", "\\1", nexf))
gtL  <- r$gt[num, ]
stopifnot(all(gtL$kObs == as.integer(mkd$kObs)))

# Start tree: identical across arms and seeds (AdditionTree is deterministic).
st <- TreeSearch::AdditionTree(r$pd)
st$edge.length <- rep(1 / nrow(st$edge), nrow(st$edge))

set.seed(seed)
t0 <- Sys.time()
res <- RunMkPrime(
  mkd, st,
  model = MkPrimeModel(coding = "variable", kPrimePrior = "geometric"),
  mcmc  = MkPrimeMCMC(nIter = nIter, thin = max(1L, nIter %/% 5000L),
                      minWarmup = 2000L, maxWarmup = 5000L,
                      nRuns = 1L, nChains = 1L,
                      gibbsSpr = identical(arm, "on"),
                      maxTime = 1500))
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

s   <- res$samples
kpc <- grep("^kPrime_", colnames(s), value = TRUE)
stopifnot(length(kpc) == nrow(gtL))
kp  <- s[, kpc, drop = FALSE]
k_post <- colMeans(kp)
u_post <- unname(k_post) - gtL$kObs

# Per-character MCMC standard error of k_post via batch means (20 batches).
bm_se <- function(v, nb = 20L) {
  n <- length(v); b <- n %/% nb; if (b < 2L) return(NA_real_)
  m <- colMeans(matrix(v[seq_len(b * nb)], nrow = b))
  sd(m) / sqrt(nb)
}
kp_se <- apply(kp, 2L, bm_se)

saveRDS(list(
  tree_idx = tree_idx, arm = arm, seed = seed, nIter = nIter,
  elapsed = elapsed, stop_reason = res$stop_reason,
  n_samples = nrow(s),
  kObs = gtL$kObs, u_true = gtL$u_true, k_true = gtL$k_true,
  k_post = unname(k_post), u_post = u_post, kp_se = unname(kp_se),
  p_mean = if ("p" %in% colnames(s)) mean(s[, "p"]) else NA_real_,
  tree_length = mean(s[, "tree_length"]),
  log_lik = mean(s[, "log_likelihood"]),
  acceptance = res$acceptance,
  scalar_cols = colnames(s)
), outfile)

cat(sprintf("t%02d %-3s seed=%d  %.0fs  nsamp=%d  mean_u=%+.4f  rho=%+.3f  stop=%s\n",
            tree_idx, arm, seed, elapsed, nrow(s), mean(u_post),
            suppressWarnings(cor(u_post, gtL$u_true, method = "spearman")),
            res$stop_reason))
