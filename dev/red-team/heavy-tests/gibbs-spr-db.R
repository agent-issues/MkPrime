# dev/red-team/heavy-tests/gibbs-spr-db.R
#
# Exact pi-invariance MAGNITUDE measurement and REGRESSION GATE for the
# `gibbs_spr` move (moveType 10).
#
# STATUS OF THE CLAIM
#   That gibbs_spr is not pi-invariant is CONFIRMED analytically and
#   independently of this harness (esjd-allocation.md Q4 + a separate
#   verifier). Under the flat Dirichlet(1,...,1) branch prior
#   (src/mcmc.cpp:323-324) the set E of trees carrying two exactly-equal
#   incident edges is pi-null; every accepted gibbs_spr lands in E; hence
#   (pi K)(E) > 0 = pi(E).
#
#   This script therefore exists for two jobs, NOT to re-demonstrate the bug:
#     (a) MAGNITUDE  -- how far is the sampled distribution from the exact
#                       target, in units someone can act on?
#     (b) GATE       -- fail clearly on current behaviour, pass on a corrected
#                       kernel. The `spr_fixed_surrogate` arm demonstrates the
#                       gate does pass on a correct kernel, so a green gate
#                       after the fix is meaningful rather than vacuous.
#
# WHICH CODE PATH IS DRIVEN
#   The 0.5 * lReg halving is in ALL THREE commit paths:
#     src/mcmc.cpp:1463-1466  main partial-CL  <-- THIS FIXTURE
#     src/mcmc.cpp:1777-1780  Q-heterogeneity
#     src/mcmc.cpp:1921-1924  full-evaluation fallback (the lines Q4 cites)
#   The fixture below uses no Q-heterogeneity and coding = "variable", which
#   routes gibbs_spr_impl (src/mcmc.cpp:1192-1205) to the MAIN PARTIAL-CL
#   path -- the one a default run actually executes. A harness that only
#   exercised the fallback would prove less than it looks like it does.
#
# EXACT TARGET  (see gibbs-spr-db.md)
#   At beta = 0 every Gibbs weight is exp(0) = 1 and every MH log-ratio loses
#   its likelihood term, so the target is the prior. The branch-fraction prior
#   is Dirichlet(1,...,1) -- FLAT (src/mcmc.cpp:323-324) -- and there is no
#   topology prior. No kernel tested here moves treeLength, so conditioning on
#   treeLength = 1:
#       pi = Uniform{labelled unrooted binary topologies}
#            (x) Dirichlet(1,...,1) on the 2n-3 edge fractions
#   pi is exactly i.i.d.-sampleable, so the reference is TRUTH, not a sampler.
#
# TEST  pi K^B = pi from an exact stationary start.
#   Draw x_r ~ pi exactly; apply B sweeps of K; record y_r. The y_r are i.i.d.
#   draws from pi K^B: no burn-in, no thinning, no mixing assumption, valid
#   p-values, and every arm directly comparable.
#
# CONTROLS
#   spr  (code 6)  -> spr_proposal_impl (src/proposals.cpp:110-166):
#                     tau = unif_rand() (:121) and
#                     logHastings = log(lRegraft) - log(lMerge) (:161).
#                     The symmetric fraction density cancels, leaving the
#                     Jacobian as the only required term: exactly correct.
#   tbr  (code 17) -> tbr_proposal_impl (src/tree_moves.cpp:492-494): second
#                     correct reference (corroborating, not gating -- see .md).
#   branch_lengths (code 4): positive control for the branch-fraction stats.
#   spr_fixed_surrogate: R-level MH on spr_proposal accepting at
#                     log U < logHastings. This IS the corrected kernel
#                     (uniform tau + Jacobian + MH), so it doubles as proof
#                     that the gate passes on a fix. GATED.
#   Power ladder: same R-level MH accepting at log U < c * logHastings for
#                     c < 1 -- a SMOOTH, non-atomic bias of tunable size, so
#                     the ladder measures this harness's detection floor.
#
# SECOND ACCUSED (report-only, never enters the gibbs_spr verdict)
#   weighted_spr (code 13) is independently CONFIRMED defective as GSPR-003:
#   it omits the required Jacobian log(lRegraft) - log(lMerge). Included as an
#   independent corroboration arm.
#
# NOTE ON THE MANDATORY TRIPLE-GUARD
#   RunMkPrime() is NEVER called. The driver is do_move_cpp, one move at a
#   time, so maxTime / setTimeLimit() / bash timeout are vacuous here: there is
#   no engine invocation, no streaming output and no convergence criterion that
#   can run away. Cost is the fixed reps x sweeps product printed at startup.
#
# Usage
#   Rscript dev/red-team/heavy-tests/gibbs-spr-db.R --quick     # <= 60 s
#   Rscript dev/red-team/heavy-tests/gibbs-spr-db.R             # full, local
#   ... [--reps N] [--sweeps B] [--n NTIP] [--ref-mult M]
#
# Output
#   dev/red-team/heavy-tests/gibbs-spr-db-results/{verdict.txt,summary.rds}
#   Exit status 1 on FAIL (so it can be used directly as a CI-style gate).
#
# STYLE NOTE
#   R object names are camelCase per the house style. Statistic identifiers
#   ("int_frac", "n_tie", ...) are STRING LITERALS used as matrix column names
#   and are deliberately left snake_case, as are MkPrime's own domain fields
#   (rel_br_lengths, tree_length) reached through the package API.

suppressPackageStartupMessages({
  library(ape)
  library(TreeTools)
})

# ---------------------------------------------------------------------------
# 0. Arguments and configuration
# ---------------------------------------------------------------------------

args  <- commandArgs(trailingOnly = TRUE)
quick <- "--quick" %in% args

.ArgInt <- function(flag, default) {
  ix <- which(args == flag)
  if (length(ix) && length(args) >= ix[1] + 1L)
    as.integer(args[ix[1] + 1L]) else default
}

nTip     <- .ArgInt("--n", 6L)
nReps    <- .ArgInt("--reps",     if (quick) 600L else 20000L)
nSweeps  <- .ArgInt("--sweeps",   if (quick)  40L else   150L)
refMult  <- .ArgInt("--ref-mult", if (quick)   8L else    10L)
nChar    <- 4L

.DfactOdd <- function(m) if (m <= 1L) 1L else prod(seq.int(m, 1L, by = -2L))
nTopos <- .DfactOdd(2L * nTip - 5L)
nEdge  <- 2L * nTip - 3L

# --- pass criterion (fixed BEFORE the code; see the .md) -------------------
gatedArms <- c("gibbs_spr", "gibbs_spr+br_1to1", "gibbs_spr+br_1to10",
               "spr", "branch_lengths", "spr_fixed_surrogate")
nStats    <- 8L                                     # 6 continuous + topo + tie
alphaFam  <- 0.01
alphaTest <- alphaFam / (length(gatedArms) * nStats)
tieMax    <- 1e-3     # floating-point-coincidence guard, not a stat threshold

cat(sprintf(paste0(
  "[gibbs-spr-db] nTip=%d  nEdge=%d  topologies=%d\n",
  "               reps=%d  sweeps=%d  reference=%d  (%s mode)\n",
  "               alphaTest=%.3g  (Bonferroni: %d gated arms x %d stats,",
  " FWER %.3g)\n"),
  nTip, nEdge, nTopos, nReps, nSweeps, refMult * nReps,
  if (quick) "quick" else "full", alphaTest, length(gatedArms), nStats,
  alphaFam))

repo <- tryCatch(system2("git", c("rev-parse", "--show-toplevel"),
                         stdout = TRUE, stderr = FALSE),
                 error = function(e) getwd())
if (!length(repo) || is.na(repo[1]) || !dir.exists(repo[1])) repo <- getwd()
repo <- repo[1]

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(repo, quiet = TRUE)
} else if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(repo, quiet = TRUE)
} else {
  stop("Need {pkgload} or {devtools}")
}

outDir <- file.path(repo, "dev/red-team/heavy-tests/gibbs-spr-db-results")
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# 1. Exact sampler for pi
# ---------------------------------------------------------------------------

# Uniform over labelled unrooted binary topologies, by sequential random-edge
# insertion: attaching taxon k+1 at a uniformly chosen edge of a uniform tree
# on k taxa is uniform over the (2k-3)!! trees on k+1 taxa. The returned
# preorder phylo has root (nTip + 1) of degree 3 always -- exactly MkPrime's
# trifurcating-root representation of an unrooted tree (2n-3 edges).
.RandomUnrootedTopology <- function(n) {
  from <- rep(n + 1L, 3L)
  to   <- 1:3
  if (n > 3L) for (k in 4L:n) {
    w <- n + (k - 2L)
    e <- sample.int(length(from), 1L)
    b <- to[e]
    to[e] <- w
    from  <- c(from, w, w)
    to    <- c(to, b, k)
  }
  root  <- n + 1L
  adj   <- split(c(to, from), c(from, to))
  nE    <- length(from)
  par   <- integer(nE); ch <- integer(nE); np <- 0L
  stack <- root; seen <- root
  while (length(stack)) {
    v  <- stack[length(stack)]; stack <- stack[-length(stack)]
    nb <- adj[[as.character(v)]]
    nb <- nb[!(nb %in% seen)]
    for (x in nb) { np <- np + 1L; par[np] <- v; ch[np] <- x }
    seen  <- c(seen, nb)
    stack <- c(stack, nb[nb > n])
  }
  TreeTools::Preorder(structure(
    list(edge = cbind(par, ch), Nnode = n - 2L,
         tip.label = paste0("t", seq_len(n)), edge.length = rep(1, nE)),
    class = "phylo"))
}

# Dirichlet(1,...,1) on the simplex (exact).
.RandomSimplex <- function(m) { g <- rexp(m); g / sum(g) }

# ---------------------------------------------------------------------------
# 2. Statistics
# ---------------------------------------------------------------------------

# Root-invariant topology key from descendant-tip bitmasks. ape::prop.part is
# NEVER used: it is root-dependent and has caused a real scoring bug in this
# project (project_scoring_bug).
.TopoKey <- function(edge, nTipLocal) {
  mask <- integer(2L * nTipLocal)
  for (i in rev(seq_len(nrow(edge)))) {        # preorder reversed = postorder
    cn <- edge[i, 2L]
    if (cn <= nTipLocal) mask[cn] <- bitwShiftL(1L, cn - 1L)
    pn <- edge[i, 1L]
    mask[pn] <- bitwOr(mask[pn], mask[cn])
  }
  full <- bitwShiftL(1L, nTipLocal) - 1L
  keys <- integer(0)
  for (i in seq_len(nrow(edge))) {
    b  <- mask[edge[i, 2L]]
    nb <- sum(as.integer(intToBits(b))[seq_len(nTipLocal)])
    if (nb < 2L || nb > nTipLocal - 2L) next    # trivial split
    if (bitwAnd(b, 1L) == 1L) b <- bitwAnd(bitwNot(b), full)  # canonical side
    keys <- c(keys, b)
  }
  paste(sort(unique(keys)), collapse = "|")
}

# Adjacency (rooting-invariant): two edge rows are adjacent iff they share a
# node. Two disjoint families, both obtained without combn():
#   sibling pairs      -- rows sharing a parent
#   parent-child pairs -- row j paired with the row whose child is parent[j]
.AdjPairs <- function(edge) {
  pc <- match(edge[, 1L], edge[, 2L])          # NA where parent == root
  ok <- !is.na(pc)
  pcPairs <- cbind(pc[ok], which(ok))
  sibPairs <- do.call(rbind, lapply(split(seq_len(nrow(edge)), edge[, 1L]),
    function(rs) if (length(rs) < 2L) NULL else {
      ii <- rep(seq_along(rs), times = length(rs) - seq_along(rs))
      jj <- unlist(lapply(seq_along(rs), function(a)
              if (a < length(rs)) (a + 1L):length(rs) else integer(0)),
              use.names = FALSE)
      cbind(rs[ii], rs[jj])
    }))
  rbind(pcPairs, sibPairs)
}

# Per-topology invariants are cached: many replicates share an edge matrix, and
# .TopoKey / .AdjPairs are the loopiest part of the harness.
.topoCache <- new.env(hash = TRUE, parent = emptyenv())
.EdgeInfo <- function(edge, nTipLocal) {
  fp  <- paste(edge, collapse = ",")
  got <- .topoCache[[fp]]
  if (!is.null(got)) return(got)
  info <- list(topo = .TopoKey(edge, nTipLocal), ap = .AdjPairs(edge),
               isInt = edge[, 2L] > nTipLocal)
  assign(fp, info, envir = .topoCache)
  info
}

# Column names below are string literals, not R object names: left snake_case.
statCont <- c("int_frac", "ord_min", "ord_med", "ord_max", "simpson",
              "bal_min")

.CollectStats <- function(states, nTipLocal) {
  n <- length(states)
  m <- matrix(NA_real_, n, length(statCont) + 2L,
              dimnames = list(NULL, c(statCont, "n_tie", "n_tie_rel")))
  topo <- character(n)
  for (i in seq_len(n)) {
    edge  <- states[[i]]$edge
    relBr <- states[[i]]$relBr
    eInfo <- .EdgeInfo(edge, nTipLocal)
    a <- relBr[eInfo$ap[, 1L]]
    b <- relBr[eInfo$ap[, 2L]]
    f <- a / (a + b)
    o <- sort(relBr)
    topo[i] <- eInfo$topo
    m[i, ] <- c(sum(relBr[eInfo$isInt]),        # int_frac
                o[1L],                          # ord_min
                o[(length(o) + 1L) %/% 2L],     # ord_med
                o[length(o)],                   # ord_max
                sum(relBr * relBr),             # simpson
                min(abs(f - 0.5)),              # bal_min
                sum(a == b),                    # n_tie      (bit-exact)
                sum(abs(a - b) <=               # n_tie_rel  (<= 1e-12 rel)
                      1e-12 * (abs(a) + abs(b)) / 2))
  }
  list(mat = m, topo = topo)
}

# ---------------------------------------------------------------------------
# 3. Fixture. At beta = 0 the data cannot influence the target; it only has to
#    keep logLik finite, so it is deliberately tiny.
# ---------------------------------------------------------------------------

set.seed(20260812L)
fixTree <- .RandomUnrootedTopology(nTip)
charMat <- matrix(sample(0:1, nTip * nChar, replace = TRUE), nrow = nTip,
                  dimnames = list(fixTree$tip.label, NULL))
mkd     <- suppressWarnings(MkPrimeData(MatrixToPhyDat(charMat)))
model   <- MkPrime:::.FinalizeModel(MkPrimeModel(), fixTree, mkd)
dataPtr <- MkPrime:::.InitMcmcData(mkd, model)

.NewStatePtr <- function(tree, relBr) {
  tree$edge.length <- relBr                    # treeLength == 1
  sp <- MkPrime:::.InitMcmcChain(MkPrime:::.InitState(tree, mkd, model))
  fill_partition_cache(dataPtr, sp)
  allocate_cl_workspace(dataPtr, sp)
  sp
}

# ---------------------------------------------------------------------------
# 4. Kernel drivers
# ---------------------------------------------------------------------------

# One sweep = the move codes in `codes`, applied in order, at beta = 0.
.RunCppArm <- function(codes, nRepsLocal, nSweepsLocal, seed) {
  set.seed(seed)
  states <- vector("list", nRepsLocal); nAcc <- 0L; nProp <- 0L
  for (r in seq_len(nRepsLocal)) {
    sp <- .NewStatePtr(.RandomUnrootedTopology(nTip), .RandomSimplex(nEdge))
    for (s in seq_len(nSweepsLocal)) for (cd in codes) {
      nAcc  <- nAcc + do_move_cpp(dataPtr, sp, cd, 0L, 0.5, 0.5, 1L, 0)
      nProp <- nProp + 1L
    }
    st <- get_mcmc_state(sp)
    states[[r]] <- list(edge = st$edge, relBr = st$relBrLengths)
  }
  list(states = states, accept = nAcc / nProp)
}

# R-level MH on the exported spr_proposal. At beta = 0 the target is the flat
# Dirichlet(1) prior, so the exact acceptance probability is
# min(1, exp(logHastings)). hFrac = 1 is EXACTLY CORRECT (this is the corrected
# kernel: uniform tau, Jacobian, MH). hFrac < 1 is a smooth, tunable,
# non-atomic violation used to calibrate detection power.
.RunRsprArm <- function(hFrac, nRepsLocal, nSweepsLocal, seed) {
  set.seed(seed)
  states <- vector("list", nRepsLocal); nAcc <- 0L; nProp <- 0L
  for (r in seq_len(nRepsLocal)) {
    edge  <- .RandomUnrootedTopology(nTip)$edge
    relBr <- .RandomSimplex(nEdge)
    for (s in seq_len(nSweepsLocal)) {
      pr    <- spr_proposal(edge, nTip, 1.0, relBr)
      nProp <- nProp + 1L
      lh    <- pr$logHastings
      if (is.finite(lh) && all(pr$rel_br_lengths > 0) &&
          log(runif(1L)) < hFrac * lh) {
        edge  <- pr$edge
        relBr <- pr$rel_br_lengths
        nAcc  <- nAcc + 1L
      }
    }
    states[[r]] <- list(edge = edge, relBr = relBr)
  }
  list(states = states, accept = nAcc / nProp)
}

.RunReference <- function(n, seed) {
  set.seed(seed)
  states <- vector("list", n)
  for (r in seq_len(n))
    states[[r]] <- list(edge = .RandomUnrootedTopology(nTip)$edge,
                        relBr = .RandomSimplex(nEdge))
  list(states = states, accept = NA_real_)
}

# ---------------------------------------------------------------------------
# 5. Reference sample + generator self-check
# ---------------------------------------------------------------------------

cat("[gibbs-spr-db] building exact reference sample ...\n")
tStart <- Sys.time()
ref    <- .RunReference(refMult * nReps, seed = 101L)
refS   <- .CollectStats(ref$states, nTip)
cat(sprintf("               %d exact draws in %.1f s\n", length(ref$states),
            as.numeric(difftime(Sys.time(), tStart, units = "secs"))))

# The whole test leans on the reference being exact: make the script prove it
# rather than asserting it.
refTab     <- table(refS$topo)
genBuckets <- length(refTab)
genChi2    <- suppressWarnings(chisq.test(as.numeric(refTab)))
genTies    <- sum(refS$mat[, "n_tie"] > 0)
genOk      <- genBuckets == nTopos && genChi2$p.value > 0.001 && genTies == 0L
cat(sprintf(paste0("               generator self-check: buckets=%d/%d  ",
                   "chi2 p=%.4f  exact-ties=%d  ->  %s\n"),
            genBuckets, nTopos, genChi2$p.value, genTies,
            if (genOk) "OK" else "BROKEN"))

# ---------------------------------------------------------------------------
# 6. Arms
# ---------------------------------------------------------------------------

arms <- list(
  # --- the accused ---------------------------------------------------------
  list(name = "gibbs_spr",          kind = "cpp",  codes = c(10L),
       seed = 2001L, gated = TRUE,  role = "accused"),
  list(name = "gibbs_spr+br_1to1",  kind = "cpp",  codes = c(10L, 4L),
       seed = 2002L, gated = TRUE,  role = "composite"),
  list(name = "gibbs_spr+br_1to10", kind = "cpp",  codes = c(10L, rep(4L, 10L)),
       seed = 2003L, gated = TRUE,  role = "composite"),
  # --- controls -----------------------------------------------------------
  list(name = "spr",                kind = "cpp",  codes = c(6L),
       seed = 2004L, gated = TRUE,  role = "control"),
  list(name = "branch_lengths",     kind = "cpp",  codes = c(4L),
       seed = 2005L, gated = TRUE,  role = "control"),
  list(name = "tbr",                kind = "cpp",  codes = c(17L),
       seed = 2007L, gated = FALSE, role = "control2"),
  list(name = "spr_fixed_surrogate", kind = "rspr", hFrac = 1.00,
       seed = 3001L, gated = TRUE,  role = "fixed"),
  # --- second accused (report only) ---------------------------------------
  list(name = "weighted_spr",       kind = "cpp",  codes = c(13L),
       seed = 2006L, gated = FALSE, role = "gspr003"),
  # --- power ladder (report only) -----------------------------------------
  list(name = "Rspr_c0.95", kind = "rspr", hFrac = 0.95, seed = 3002L,
       gated = FALSE, role = "ladder"),
  list(name = "Rspr_c0.90", kind = "rspr", hFrac = 0.90, seed = 3003L,
       gated = FALSE, role = "ladder"),
  list(name = "Rspr_c0.75", kind = "rspr", hFrac = 0.75, seed = 3004L,
       gated = FALSE, role = "ladder"),
  list(name = "Rspr_c0.50", kind = "rspr", hFrac = 0.50, seed = 3005L,
       gated = FALSE, role = "ladder"),
  list(name = "Rspr_c0.00", kind = "rspr", hFrac = 0.00, seed = 3006L,
       gated = FALSE, role = "ladder")
)

# Branch-mixing dose-response: how much branch-length mixing per gibbs_spr is
# needed before the distortion falls below the detection floor? Full mode only
# (31 moves/sweep is the single most expensive arm).
if (!quick) {
  arms <- c(arms, list(
    list(name = "gibbs_spr+br_1to30", kind = "cpp",
         codes = c(10L, rep(4L, 30L)), seed = 2008L, gated = FALSE,
         role = "composite")))
}

# Arms that never propose a topology change: their topology marginal is the
# exact pi draw by construction, so the chi-squared is uninformative there.
fixedTopoArms <- c("branch_lengths")

# ---------------------------------------------------------------------------
# 7. Assessment: significance AND magnitude
# ---------------------------------------------------------------------------

.Wasserstein1 <- function(x, y) {       # via matched quantiles
  p <- seq_len(999L) / 1000
  mean(abs(as.numeric(quantile(x, p, names = FALSE, type = 7)) -
           as.numeric(quantile(y, p, names = FALSE, type = 7))))
}

.Assess <- function(armName, S, movesTopology) {
  n <- nrow(S$mat)
  nRef <- nrow(refS$mat)
  cont <- lapply(statCont, function(st) {
    x <- S$mat[, st]; y <- refS$mat[, st]
    kk <- suppressWarnings(stats::ks.test(x, y))
    dd <- mean(x) - mean(y)
    se <- sqrt(stats::var(x) / n + stats::var(y) / nRef)
    list(stat = st, p = unname(kk$p.value), dStat = unname(kk$statistic),
         mean = mean(x), refMean = mean(y),
         delta = dd, se = se, ci = dd + c(-1.96, 1.96) * se,
         relBias = dd / mean(y), w1 = .Wasserstein1(x, y),
         pass = unname(kk$p.value) >= alphaTest)
  })
  names(cont) <- statCont

  tab <- table(factor(S$topo, levels = names(refTab)))
  tv  <- 0.5 * sum(abs(as.numeric(tab) / n - 1 / nTopos))
  # Total variation from uniform is heavily noise-dominated at small n/nTopos
  # (E[TV] > 0 even under H0), so it is only interpretable against a matched
  # null: TV of an exact pi subsample of the SAME size. Deterministic seed so
  # the baseline is reproducible and identical across arms.
  tvNull <- {
    savedSeed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    set.seed(917L)
    subTab <- table(factor(sample(refS$topo, n), levels = names(refTab)))
    assign(".Random.seed", savedSeed, envir = .GlobalEnv)
    0.5 * sum(abs(as.numeric(subTab) / n - 1 / nTopos))
  }
  if (movesTopology) {
    chi <- suppressWarnings(chisq.test(as.numeric(tab)))
    # Full occupancy is only evidence of a defect when the expected count per
    # bucket is comfortably large; at n/nTopos < 20 empty buckets are ordinary
    # Poisson noise under H0 and requiring all of them spuriously FAILs
    # correct kernels (cf. SBC-HARNESS-001 in subtree-swap-db.R). chi-squared
    # already encodes the variance, so defer to it below that density.
    occupancyInformative <- (n / nTopos) >= 20
    topo <- list(p = unname(chi$p.value), buckets = sum(tab > 0), tv = tv,
                 tvNull = tvNull, occupancyChecked = occupancyInformative,
                 pass = unname(chi$p.value) >= alphaTest &&
                        (!occupancyInformative || sum(tab > 0) == nTopos))
  } else {
    topo <- list(p = NA_real_, buckets = sum(tab > 0), tv = tv,
                 tvNull = tvNull, occupancyChecked = FALSE, pass = TRUE)
  }

  tie <- list(rate = mean(S$mat[, "n_tie"] > 0),
              rateRel = mean(S$mat[, "n_tie_rel"] > 0),
              meanCount = mean(S$mat[, "n_tie"]))
  tie$pass <- tie$rate <= tieMax

  list(name = armName, n = n, cont = cont, topo = topo, tie = tie,
       pass = all(vapply(cont, `[[`, logical(1), "pass")) && topo$pass &&
              tie$pass,
       minContP = min(vapply(cont, `[[`, numeric(1), "p")),
       maxAbsRelBias = max(abs(vapply(cont, `[[`, numeric(1), "relBias"))))
}

results <- list()
for (arm in arms) {
  cat(sprintf("[gibbs-spr-db]   %-22s ", arm$name))
  tArm <- Sys.time()
  run  <- if (arm$kind == "cpp") .RunCppArm(arm$codes, nReps, nSweeps, arm$seed)
          else .RunRsprArm(arm$hFrac, nReps, nSweeps, arm$seed)
  S    <- .CollectStats(run$states, nTip)
  asr  <- .Assess(arm$name, S, !(arm$name %in% fixedTopoArms))
  dt   <- as.numeric(difftime(Sys.time(), tArm, units = "secs"))
  results[[arm$name]] <- c(asr, list(
    gated = arm$gated, role = arm$role, seed = arm$seed,
    codes = if (arm$kind == "cpp") arm$codes else NA_integer_,
    accept = run$accept, elapsedSec = dt,
    statMat = S$mat, topoTab = table(S$topo)))
  cat(sprintf("%6.1fs acc=%.3f  minP=%-9.3g tie=%.4f |relBias|<=%.4f  %s%s\n",
              dt, run$accept, asr$minContP, asr$tie$rate, asr$maxAbsRelBias,
              if (asr$pass) "indistinguishable" else "DEVIATES",
              if (arm$gated) "" else "  [report-only]"))
}

# ---------------------------------------------------------------------------
# 8. Verdict, magnitude, power
# ---------------------------------------------------------------------------

ctrlOk  <- results[["spr"]]$pass && results[["branch_lengths"]]$pass
fixedOk <- results[["spr_fixed_surrogate"]]$pass
q1Fail  <- !results[["gibbs_spr"]]$pass
q2Fail  <- !results[["gibbs_spr+br_1to1"]]$pass ||
           !results[["gibbs_spr+br_1to10"]]$pass

ladNames <- names(results)[vapply(results, `[[`, character(1), "role") ==
                             "ladder"]
ladder <- data.frame(
  arm    = ladNames,
  cFrac  = as.numeric(sub("^Rspr_c", "", ladNames)),
  minP   = vapply(results[ladNames], `[[`, numeric(1), "minContP"),
  detect = !vapply(results[ladNames], `[[`, logical(1), "pass"),
  dInt   = vapply(results[ladNames], function(r) r$cont$int_frac$delta,
                  numeric(1)),
  relB   = vapply(results[ladNames], `[[`, numeric(1), "maxAbsRelBias"),
  stringsAsFactors = FALSE)
ladder <- ladder[order(-ladder$cFrac), ]
detRows <- ladder[ladder$detect, , drop = FALSE]
powerFloorC <- if (nrow(detRows)) max(detRows$cFrac) else NA_real_
powerFloorB <- if (nrow(detRows))
  detRows$relB[which.max(detRows$cFrac)] else NA_real_

verdict <- if (!genOk || !ctrlOk || !fixedOk) {
  "WARN"
} else if (q1Fail) {
  "FAIL"
} else if (is.na(powerFloorC)) {
  "INCONCLUSIVE"
} else {
  "PASS"
}

.FmtArm <- function(r) {
  hdr <- sprintf("  %-22s n=%d accept=%.3f  %s%s", r$name, r$n, r$accept,
                 if (r$pass) "indistinguishable from pi"
                 else "DEVIATES from pi",
                 if (r$gated) "" else "   [report-only]")
  rows <- vapply(r$cont, function(cc) sprintf(
    paste0("      %-9s KS p=%-10.3g D=%-7.4f mean=%-9.5f exact=%-9.5f ",
           "relBias=%+8.4f W1=%-9.3g %s"),
    cc$stat, cc$p, cc$dStat, cc$mean, cc$refMean, cc$relBias, cc$w1,
    if (cc$pass) "ok" else "FAIL"), character(1))
  topo <- sprintf(
    "      %-9s chi2 p=%-11s buckets=%d/%d  TV=%.5f (H0 baseline %.5f)  %s",
    "topo",
    if (is.na(r$topo$p)) "n/a(fixed)" else
      formatC(r$topo$p, format = "g", digits = 3),
    r$topo$buckets, nTopos, r$topo$tv, r$topo$tvNull,
    if (r$topo$pass) "ok" else "FAIL")
  tie <- sprintf(paste0("      %-9s P(bit-exact tie)=%.5f  P(rel tie)=%.5f  ",
                        "mean ties/state=%.4f  [exact target: 0]  %s"),
                 "n_tie", r$tie$rate, r$tie$rateRel, r$tie$meanCount,
                 if (r$tie$pass) "ok" else "FAIL")
  paste(c(hdr, rows, topo, tie), collapse = "\n")
}

.MagRow <- function(nm) {
  r <- results[[nm]]
  sprintf(paste0("  %-22s pi-null mass=%7.4f  topo TV=%.5f (H0 %.5f)  ",
                 "E[int_frac] %+.5f (%+6.2f%%)  W1(int_frac)=%.3g"),
          nm, r$tie$rate, r$topo$tv, r$topo$tvNull, r$cont$int_frac$delta,
          100 * r$cont$int_frac$relBias, r$cont$int_frac$w1)
}
magArms <- intersect(c("gibbs_spr", "gibbs_spr+br_1to1", "gibbs_spr+br_1to10",
                       "gibbs_spr+br_1to30", "weighted_spr",
                       "spr_fixed_surrogate", "spr", "branch_lengths"),
                     names(results))

lines <- c(
  sprintf("VERDICT: %s", verdict),
  "",
  "Harness: gibbs-spr-db -- exact pi-invariance MAGNITUDE + REGRESSION GATE",
  "Role:    the non-invariance of gibbs_spr is already CONFIRMED analytically",
  "         (esjd-allocation.md Q4). This harness measures HOW MUCH, and gates",
  "         the fix. It is not evidence-of-existence.",
  "Path:    fixture drives the MAIN partial-CL commit path (src/mcmc.cpp:1463-1466)",
  sprintf(paste0("Target:  pi = Uniform{%d unrooted topologies} (x) ",
                 "Dirichlet(1,...,1) on %d fractions"), nTopos, nEdge),
  "Design:  x ~ pi exactly -> B sweeps of K -> compare. i.i.d. replicates; no",
  "         burn-in, no thinning, no mixing assumption.",
  sprintf("Scale:   nTip=%d  reps=%d  sweeps=%d  reference=%d  (%s mode)",
          nTip, nReps, nSweeps, refMult * nReps,
          if (quick) "quick" else "full"),
  sprintf(paste0("Alpha:   %.3g per test (Bonferroni %d gated arms x %d ",
                 "stats, FWER %.3g)"),
          alphaTest, length(gatedArms), nStats, alphaFam),
  sprintf("Tie gate: fail if P(bit-exact adjacent tie) > %.1g  (exact target 0)",
          tieMax),
  "",
  "--- REFERENCE GENERATOR SELF-CHECK ---",
  sprintf("  uniform-topology chi2 p=%.4f over %d/%d buckets; exact ties=%d -> %s",
          genChi2$p.value, genBuckets, nTopos, genTies,
          if (genOk) "OK" else "BROKEN"),
  "",
  "--- CONTROLS (harness validity) ---",
  "  spr = spr_proposal_impl (src/proposals.cpp:121,161): tau ~ U(0,1) and",
  "  logHastings = log(lRegraft) - log(lMerge). Exactly correct.",
  .FmtArm(results[["spr"]]),
  .FmtArm(results[["branch_lengths"]]),
  "  tbr = tbr_proposal_impl (src/tree_moves.cpp:492-494): second correct",
  "  reference; corroborating only (compound move, larger surface).",
  .FmtArm(results[["tbr"]]),
  sprintf("  primary controls %s",
          if (ctrlOk) "PASS -> harness certified"
          else "FAIL -> nothing below is interpretable"),
  "",
  "--- GATE VALIDITY: does the gate pass on a CORRECTED kernel? ---",
  "  spr_fixed_surrogate = uniform tau + Jacobian + MH, i.e. what a fixed",
  "  gibbs_spr must reduce to. If this passes, a green gate after the fix is",
  "  meaningful rather than vacuous.",
  .FmtArm(results[["spr_fixed_surrogate"]]),
  sprintf("  gate-on-correct-kernel: %s",
          if (fixedOk) "PASSES -> gate is meaningful"
          else "FAILS -> gate unusable, harness must be fixed first"),
  "",
  "--- Q1: is the gibbs_spr kernel pi-invariant on its own? ---",
  .FmtArm(results[["gibbs_spr"]]),
  sprintf("  Q1 = %s", if (q1Fail) "FAIL (gibbs_spr is NOT pi-invariant)"
          else "not refuted at this scale"),
  "",
  "--- Q2: do branch-length moves compensate? (composite chains) ---",
  .FmtArm(results[["gibbs_spr+br_1to1"]]),
  .FmtArm(results[["gibbs_spr+br_1to10"]]),
  if ("gibbs_spr+br_1to30" %in% names(results))
    .FmtArm(results[["gibbs_spr+br_1to30"]]) else NULL,
  sprintf("  Q2 = %s", if (q2Fail)
            "FAIL (interleaving branch-length moves does NOT restore pi)"
          else "composite not distinguished from pi at this scale"),
  "",
  "--- MAGNITUDE (the point of this harness) ---",
  "  pi-null mass = fraction of sampled states carrying a bit-exactly equal",
  "  adjacent edge pair; the exact target puts 0 there. topo TV = total",
  "  variation from uniform on the topology marginal, shown against the H0",
  "  baseline (TV of an exact pi subsample of the same size) because TV is",
  "  noise-dominated at small n/nTopos. relBias / W1 quantify distortion of",
  "  the continuous part.",
  paste(vapply(magArms, .MagRow, character(1)), collapse = "\n"),
  "",
  "--- POWER LADDER (R-level MH on spr_proposal, accept at c*logHastings) ---",
  "  c = 1 is exactly correct (that is spr_fixed_surrogate); c < 1 is a smooth,",
  "  non-atomic bias, so this measures the detection floor for biases that do",
  "  NOT leave the atomic tie fingerprint.",
  paste(sprintf("    c=%.2f  min KS p=%-10.3g |relBias|<=%.5f  %s",
                ladder$cFrac, ladder$minP, ladder$relB,
                ifelse(ladder$detect, "DETECTED", "not detected")),
        collapse = "\n"),
  sprintf("  detection floor: smallest detected smooth violation c=%s at |relBias|=%s",
          ifelse(is.na(powerFloorC), "none",
                 formatC(powerFloorC, digits = 3)),
          ifelse(is.na(powerFloorB), "n/a",
                 formatC(powerFloorB, format = "g", digits = 3))),
  "",
  "--- SECOND ACCUSED (report-only): weighted_spr, code 13 = GSPR-003 ---",
  "  Independently CONFIRMED defective: omits the Jacobian",
  "  log(lRegraft) - log(lMerge). Reported here as corroboration only; it never",
  "  enters the gibbs_spr verdict and is not used as a control.",
  .FmtArm(results[["weighted_spr"]]),
  "",
  "--- LIMITATION ---",
  "  beta = 0 makes the target exactly known and exactly sampleable, and the",
  "  branch-length map (merge + 0.5*lReg split) is beta-independent, so a FAIL",
  "  here is a FAIL at any beta. The converse does not hold: the MAGNITUDE",
  "  numbers above are distortions of the beta = 0 target and need not equal",
  "  the beta = 1 posterior distortion. The pi-null mass transfers most",
  "  directly, since it is set by acceptance rate and branch-move dose rather",
  "  than by the likelihood.",
  "",
  "--- INTERPRETATION ---",
  if (verdict == "WARN") {
    paste(c("  Harness suspect; do not read the arms.",
            if (!genOk) "  -> exact-pi generator self-check failed.",
            if (!ctrlOk)
              "  -> a primary positive control failed (spr / branch_lengths).",
            if (!fixedOk)
              "  -> the corrected-kernel surrogate failed, so the gate cannot certify a fix."),
          collapse = "\n")
  } else if (verdict == "FAIL") {
    paste(c(
      "  Controls pass, the corrected-kernel surrogate passes, and gibbs_spr",
      "  deviates from the exact target: the gate is live and currently RED.",
      sprintf("  Q4 confirmed empirically. Mechanism fingerprint: %.2f%% of gibbs_spr",
              100 * results[["gibbs_spr"]]$tie$rate),
      "  states carry a BIT-EXACT adjacent edge tie (exact target 0%), i.e. the",
      "  deterministic 0.5*lReg split at src/mcmc.cpp:1463-1466.",
      sprintf(paste0("  NOTE: the gibbs_spr TOPOLOGY marginal is clean ",
                     "(chi2 p=%.3f, TV=%.5f vs H0 %.5f).
",
                     "  A harness that only checked topology frequencies would ",
                     "have MISSED this entirely;
",
                     "  the defect is confined to the continuous part, exactly ",
                     "as hastings-tree-moves.md
  Sec.5 conceded ",
                     "-- which is why that concession was not a defence."),
              results[["gibbs_spr"]]$topo$p, results[["gibbs_spr"]]$topo$tv,
              results[["gibbs_spr"]]$topo$tvNull),
      sprintf("  Q2: %s", if (q2Fail)
        "the branch-move defence in hastings-tree-moves.md Sec.5 is refuted."
        else "composite arms survive at this scale -- read against the power floor."),
      "  Re-run after the fix: PASS is expected, and is meaningful because the",
      "  corrected-kernel surrogate already passes the same gate."),
      collapse = "\n")
  } else if (verdict == "PASS") {
    sprintf(paste0(
      "  Controls and the corrected-kernel surrogate pass, and gibbs_spr is\n",
      "  indistinguishable from pi at this scale. POWER: smooth violations down\n",
      "  to c=%.2f (|relBias| = %.3g) ARE detected here, and any atomic tie mass\n",
      "  above %.1g would be caught, so this null is informative above that floor."),
      powerFloorC, powerFloorB, tieMax)
  } else {
    paste0("  INCONCLUSIVE: the power ladder detected no smooth violation, so ",
           "a null on gibbs_spr would carry no information. Raise --reps / ",
           "--sweeps.")
  }
)

summaryObj <- list(
  config = list(nTip = nTip, nEdge = nEdge, nTopos = nTopos, reps = nReps,
                sweeps = nSweeps, refN = refMult * nReps, nChar = nChar,
                quick = quick, alphaTest = alphaTest, alphaFamily = alphaFam,
                tieMax = tieMax, gatedArms = gatedArms,
                codePaths = c(mainPartialCl = "src/mcmc.cpp:1463-1466",
                              qHet = "src/mcmc.cpp:1777-1780",
                              fullFallback = "src/mcmc.cpp:1921-1924"),
                drivenPath = "mainPartialCl"),
  generatorCheck = list(buckets = genBuckets, chi2P = genChi2$p.value,
                        exactTies = genTies, ok = genOk),
  reference = list(statMat = refS$mat, topoTab = refTab),
  arms = results,
  powerLadder = ladder,
  powerFloor = list(cFrac = powerFloorC, relBias = powerFloorB),
  q1Fail = q1Fail, q2Fail = q2Fail,
  controlsOk = ctrlOk, gateOnCorrectKernelOk = fixedOk,
  verdict = verdict,
  rVersion = utils::sessionInfo()$R.version$version.string)
saveRDS(summaryObj, file.path(outDir, "summary.rds"))
writeLines(lines, file.path(outDir, "verdict.txt"))

cat("\n", paste(lines, collapse = "\n"), "\n", sep = "")
cat("\nWrote:", file.path(outDir, "verdict.txt"), "\n")
cat("Wrote:", file.path(outDir, "summary.rds"), "\n")

if (verdict == "FAIL") quit(status = 1L)
