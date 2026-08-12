#!/usr/bin/env Rscript
# Topology-multimodal benchmark target for the ESJD proposal-scheduling
# campaign (stage S0 of dev/plans/2026-08-12-esjd-implementation-plan.md).
#
# `BimodalTarget()` constructs a conflicting-signal morphological matrix whose
# topology posterior has two well-separated islands with (by construction)
# comparable mass; `VerifyBimodalTarget()` establishes that claim with
# measurements rather than assertion.  Nothing here touches R/ or src/, and
# nothing here belongs in tests/testthat/.
#
# Usage
#   Rscript dev/benchmarks/bimodal-target.R --quick
#   Rscript dev/benchmarks/bimodal-target.R --stage classifier,valley
#   Rscript dev/benchmarks/bimodal-target.R           # full local verification
#
# `--quick` is an *execution* check (~4 min) and always reports WARN: at
# 40 characters per tree the posterior has no islands to find, so none of the
# thresholds below is powered. The full run is ~2 h on one core and writes
# `summary.rds` / `verdict.txt` (and `summary-merged.rds` across staged
# invocations) to `dev/benchmarks/bimodal-target-results/`.
#
# `source()` it instead to get the constructor without running anything.
#
# ---------------------------------------------------------------------------
# VERDICT, 2026-08-12: this construction does NOT yield a bimodal topology
# posterior, and must not be used as a mixing benchmark until the `valley`
# stage passes. It is committed as a working harness plus the refutation.
#
# Conflicting-signal simulation is the S0 plan's sketch (~60 characters on tree
# A, ~60 on tree B, separated by a distant SPR). Measured, the pooled matrix is
# fitted best by neither generating tree but by the *compromise* topology, with
# the focal clade attached to the central edge between the two contested
# positions:
#
#   * one caterpillar spine, MCMC over 25 placements: valley depth -14.7 nats;
#   * two mirror arms, MCMC at A/B/centre in three regimes: -10.0, -23.9,
#     and +8.9 on max-logL but -0.1 on mean-logL, i.e. no barrier;
#   * two mirror arms, proxy scan of 21 cells (internalLength 0.05-0.70,
#     stateCounts {2,3,4} to {5,6,8}, cladeSize 4-8, arms 5-7, focal stem free
#     or capped): negative in 20, the exception a degenerate cell whose two
#     "islands" are not actually distant.
#
# Decomposed on the default design, committing to the supported island gains
# only +3.0 nats across its own 150 characters while costing +28.1 nats on the
# 150 conflicting ones. The focal clade also goes rogue -- its fitted stem pins
# at the top of the search grid, 1.6-3.2 against a generating 0.35 -- but
# capping the stem makes the barrier *more* negative (-54.7 against -46.4), so
# the rogue stem is a symptom, not the cause. The cause is that a 50/50 mixture
# of two trees is better explained by one intermediate tree than by either
# component. Balance and bimodality pull against each other: the symmetry that
# equalises the two islands is what makes the compromise competitive.
#
# Two design choices are nonetheless settled, and worth keeping if anyone
# revisits this:
#
#  * Two arms, not one caterpillar. On a single spine the candidate placements
#    are nested, so a character uniting the clade with the *large* sister group
#    is equally happy with the clade anywhere deeper, and the conflict is
#    one-sided. With two arms the two positions are mutually exclusive, both
#    character sets discriminate, and every crossing runs through exactly one
#    topology (clade on the central edge) -- which makes the barrier a single
#    measurable number.
#  * A clade, not a rogue tip. Several tips sharing a marker make the focal
#    unit's state reliable per character; a single rogue tip is also trivially
#    relocated by one SPR.
#
# Recommended next step is in the report, not here: prefer an empirical matrix
# already known to be peaky over any simulation under the inference model.
#
# What each stage claims, and what makes it PASS
#
#   classifier  Island labelling is correct and root-independent. Hand-built
#               placements must get their hand-derived labels; a tree with the
#               clade broken up must be "other"; and rooting the same unrooted
#               tree three ways must not change its label. PASS = all of these.
#               (`ape::prop.part` is root-dependent and has caused a real
#               scoring bug here before -- project_scoring_bug -- hence the
#               rooting check.)
#   valley      The islands are separated by a deep log-likelihood valley. A
#               `fixTopology` chain at each of the 25 distinct placements of the
#               clade gives its best fit *with branch lengths adapted*, so the
#               escape above is available and priced in. PASS = min(peak A,
#               peak B) - best non-island placement >= 10 nats, i.e. the
#               crossing topology holds < e^-10 ~ 5e-5 of either island's mass.
#               Also asserts `RunMkPrime(tree = )` is honoured.
#   split       Both islands hold non-trivial mass. Stepping-stone log marginal
#               likelihood under `fixTopology` at A and at B; with a flat
#               topology prior p(T | y) is proportional to p(y | T). PASS =
#               |dlogZ| <= 3, a mode-level split no worse than 95:5. The
#               constructor targets <= 1.5 on a cheaper proxy.
#   hops        The target discriminates samplers. Free-topology chains started
#               in each island, under plain-MH topology moves and under package
#               defaults. PASS = zero island switches under plain MH, so no
#               chain crosses by NNI/SPR/TBR walking alone within the budget.
#               Switches under defaults are reported, not required: they are the
#               signal the ESJD A/B is meant to resolve.
#   reach       The target is a benchmark, not a trap: both islands must be
#               reachable. Metropolis-coupled runs from random starting trees.
#               PASS = both islands occupied somewhere in the set.
#   seeds       The construction generalises. Further data seeds must all reach
#               the balance tolerance with a positive proxy barrier.
#
# `gibbs_spr` is not pi-invariant (it commits a deterministic 0.5 edge split
# with no MH step), so any stage running package defaults -- `hops` in its
# non-plainMh rows, and `reach` -- measures against a slightly wrong target.
# That is fine for relative discrimination, which is all this target needs to
# provide, but every exact claim should come from the plain-MH rows, where
# `nni`/`spr`/`pspr`/`tbr` are correct.
#
# Not covered: neomorphic characters (the target is all-transformational, so
# `rate_loss`/`rate_neo` moves are inert), and crossings that rearrange the
# backbone instead of moving the focal clade.
#
# The A/B driver (`dev/benchmarks/proposal-schedule-ab.R`) is a separate
# deliverable; keep scheduling logic out of this file.

# ---------------------------------------------------------------- environment

.LoadMkPrime <- function() {
  if ("MkPrime" %in% loadedNamespaces()) {
    return(invisible(TRUE))
  }
  # Benchmarks must use an -O2 install, never load_all()'s debug DLL
  # (r-conventions skill, "Building & loading").
  agentLibs <- Sys.glob(".agent-*")
  agentLibs <- agentLibs[dir.exists(file.path(agentLibs, "MkPrime"))]
  if (length(agentLibs)) {
    .libPaths(c(agentLibs[[1]], .libPaths()))
  }
  suppressPackageStartupMessages(library(MkPrime))
  suppressPackageStartupMessages(library(TreeTools))
  invisible(TRUE)
}

# ------------------------------------------------------------ tree scaffolding

.PectinateNewick <- function(labs) {
  n <- length(labs)
  s <- paste0("(", labs[n - 1L], ",", labs[n], ")")
  for (i in seq.int(n - 2L, 1L)) {
    s <- paste0("(", labs[i], ",", s, ")")
  }
  s
}

# Uniform two-value branch lengths: every candidate topology is scored under
# the same rule, so no candidate is handicapped by an arbitrary length draw.
.SetLengths <- function(tr, internal, pendant) {
  tr$edge.length <- ifelse(tr$edge[, 2] <= length(tr$tip.label),
                           pendant, internal)
  TreeTools::Preorder(tr)
}

# Two mirror-image pectinate arms, joined by the central edge. `armRTips` is
# stored in reversed order so that armLTips[j] and armRTips[j] are exchanged by
# the left-right symmetry, making trees A and B isomorphic.
.BackboneTree <- function(design) {
  txt <- paste0("(", .PectinateNewick(design$armLTips), ",",
                .PectinateNewick(design$armRTips), ");")
  .SetLengths(ape::read.tree(text = txt),
              design$internalLength, design$pendantLength)
}

# Graft the focal clade onto the edge above `where` (a node number or tip
# label). Sequential `AddTip` builds the clade pectinately; lengths are then
# reset wholesale so `AddTip`'s own length choices never leak in.
.AttachClade <- function(backbone, where, design) {
  tr <- TreeTools::AddTip(backbone, where = where,
                          label = design$cladeTips[[1]])
  for (i in seq.int(2L, length(design$cladeTips))) {
    tr <- TreeTools::AddTip(tr, where = design$cladeTips[[i - 1L]],
                            label = design$cladeTips[[i]])
  }
  .SetLengths(tr, design$internalLength, design$pendantLength)
}

.ArmNode <- function(backbone, labs) {
  if (length(labs) == 1L) {
    match(labs, backbone$tip.label)
  } else {
    ape::getMRCA(backbone, labs)
  }
}

# Every distinct placement of the focal clade on the backbone: one per edge of
# the *unrooted* backbone. Attaching above the two children of the root gives
# the same unrooted tree, so the list is de-duplicated by split signature --
# an inflated list would double-count island or crossing topologies.
.CandidateTrees <- function(design, backbone) {
  nodes <- setdiff(seq_len(max(backbone$edge)), TreeTools::RootNode(backbone))
  out <- lapply(nodes, .AttachClade, backbone = backbone, design = design)
  sister <- vapply(nodes, function(nd) {
    paste(sort(.DescendantTips(backbone, nd)), collapse = "+")
  }, character(1))
  names(out) <- sister
  tipLabels <- c(design$armLTips, design$armRTips, design$cladeTips)
  sig <- vapply(out, function(tr) {
    paste(sort(.SplitKeys(tr, tipLabels)), collapse = "|")
  }, character(1))
  # Return:
  out[!duplicated(sig)]
}

.DescendantTips <- function(tree, node) {
  nTip <- length(tree$tip.label)
  stack <- node
  out <- integer(0)
  while (length(stack)) {
    nd <- stack[[1]]
    stack <- stack[-1]
    if (nd <= nTip) {
      out <- c(out, nd)
    } else {
      stack <- c(stack, tree$edge[tree$edge[, 1] == nd, 2])
    }
  }
  tree$tip.label[out]
}

# ----------------------------------------------------------- island labelling
#
# `ape::prop.part` is root-dependent and has caused a real scoring bug in this
# project (project_scoring_bug); every split computation below goes through
# `TreeTools::as.Splits` on an explicitly unrooted tree, with polarity
# canonicalized so that a split and its complement share one key.

.SplitKeys <- function(tree, tipLabels) {
  sp <- TreeTools::as.Splits(TreeTools::UnrootTree(tree), tipLabels = tipLabels)
  if (length(sp) == 0L) {
    return(character(0))
  }
  mat <- as.logical(sp)
  mat <- mat[, tipLabels, drop = FALSE]
  apply(mat, 1L, .SplitKeyFromMembership)
}

.SplitKeyFromMembership <- function(inSplit) {
  if (isTRUE(inSplit[[1]])) {
    inSplit <- !inSplit
  }
  paste0(as.integer(inSplit), collapse = "")
}

.MembershipKey <- function(members, tipLabels) {
  .SplitKeyFromMembership(tipLabels %in% members)
}

#' Island label for each sampled tree
#'
#' The focal clade sits in exactly one arm, or on the central edge, or the
#' backbone has been rearranged so neither arm survives.
#'
#' @param trees `phylo` object or list of them.
#' @param target List returned by [BimodalTarget()].
#' @return Character vector of `"A"` (clade inside the left arm), `"B"` (inside
#' the right arm), `"centre"` (both arms intact, so the clade is on the edge
#' between them) or `"other"` (neither arm is monophyletic).
TargetIslands <- function(trees, target) {
  if (inherits(trees, "phylo")) {
    trees <- list(trees)
  }
  tipLabels <- target$tipLabels
  vapply(trees, function(tr) {
    keys <- .SplitKeys(tr, tipLabels)
    armLIntact <- target$keyArmL %in% keys
    armRIntact <- target$keyArmR %in% keys
    if (armLIntact && armRIntact) {
      "centre"
    } else if (armRIntact) {
      "A"
    } else if (armLIntact) {
      "B"
    } else {
      "other"
    }
  }, character(1))
}

# Stem length of the focal clade; NA when the clade is not monophyletic.
# The tree is rooted on a backbone tip first so the clade never spans the root.
.CladeStemLength <- function(tree, target) {
  tr <- tryCatch(TreeTools::RootTree(tree, target$armLTips[[1]]),
                 error = function(e) tree)
  mrca <- tryCatch(ape::getMRCA(tr, target$cladeTips), error = function(e) NULL)
  if (is.null(mrca)) {
    return(NA_real_)
  }
  if (!setequal(.DescendantTips(tr, mrca), target$cladeTips)) {
    return(NA_real_)
  }
  el <- tr$edge.length[tr$edge[, 2] == mrca]
  if (length(el) == 1L) el else NA_real_
}

# -------------------------------------------------------- character simulation

# One character under JC(k). The rate convention must match inference
# (`src/likelihood.cpp`: arg = -k t / (k - 1)); naive exp(-k t) diverges for
# k > 2 and would silently mis-specify the generating branch lengths.
.SimCharacter <- function(tree, k) {
  nTip <- length(tree$tip.label)
  st <- integer(2L * nTip - 1L)
  st[nTip + 1L] <- sample.int(k, 1L) - 1L
  ed <- tree$edge
  el <- tree$edge.length
  # Preorder guarantees each parent is written before it is read.
  for (e in seq_len(nrow(ed))) {
    pa <- ed[e, 1L]
    ch <- ed[e, 2L]
    pSame <- 1 / k + (1 - 1 / k) * exp(-k * el[[e]] / (k - 1))
    st[ch] <- if (stats::runif(1L) < pSame) {
      st[pa]
    } else {
      sample(setdiff(seq.int(0L, k - 1L), st[pa]), 1L)
    }
  }
  st[seq_len(nTip)]
}

.Canonicalize <- function(v) {
  match(v, sort(unique(v))) - 1L
}

.SimMatrix <- function(tree, stateCounts, nChar, tipLabels) {
  ks <- rep(stateCounts, length.out = nChar)
  m <- vapply(ks, function(k) .SimCharacter(tree, k),
              integer(length(tree$tip.label)))
  rownames(m) <- tree$tip.label
  m[tipLabels, , drop = FALSE]
}

# ------------------------------------------------------- mass-balance proxy

# logL maximised over a 3-parameter branch-length family: multipliers on
# internal and on terminal edges, plus the focal clade's stem in absolute
# units. The stem is free because the "uncommitted clade" escape is the whole
# question -- a proxy that held it fixed reported a 20-nat barrier where the
# adapted MCMC found none.
.ProxyMaxLogL <- function(tree, mkd, target,
                          iScale = exp(seq(log(0.6), log(2.0), length.out = 6L)),
                          pScale = exp(seq(log(0.5), log(2.0), length.out = 4L)),
                          stem = c(0.05, 0.1, 0.2, 0.4, 0.8, 1.6, 3.2)) {
  isPendant <- tree$edge[, 2] <= length(tree$tip.label)
  stemEdge <- .CladeStemEdge(tree, target)
  base <- tree$edge.length
  best <- -Inf
  for (i in iScale) {
    for (p in pScale) {
      el <- base * ifelse(isPendant, p, i)
      for (s in stem) {
        if (!is.na(stemEdge)) {
          el[[stemEdge]] <- s
        }
        tree$edge.length <- el
        ll <- MkpLogLikelihood(tree, mkd, nCat = 1L)
        if (is.finite(ll) && ll > best) {
          best <- ll
        }
        if (is.na(stemEdge)) {
          break
        }
      }
    }
  }
  best
}

.CladeStemEdge <- function(tree, target) {
  mrca <- tryCatch(ape::getMRCA(tree, target$cladeTips), error = function(e) NULL)
  if (is.null(mrca)) {
    return(NA_integer_)
  }
  if (!setequal(.DescendantTips(tree, mrca), target$cladeTips)) {
    return(NA_integer_)
  }
  e <- which(tree$edge[, 2] == mrca)
  if (length(e) == 1L) e else NA_integer_
}

.ProxyDelta <- function(mkd, target) {
  .ProxyMaxLogL(target$treeA, mkd, target) -
    .ProxyMaxLogL(target$treeB, mkd, target)
}

# ------------------------------------------------------------- the constructor

#' Construct the bimodal benchmark target
#'
#' Simulates `charsPerTree` characters on each of two topologies that place a
#' focal clade in opposite arms of a symmetric backbone, then keeps the pooled
#' matrix only if the two topologies are near-equally supported.
#'
#' @param seed Integer seeding both simulation and the balance search.
#' @param nArm Integer giving the number of taxa in each backbone arm.
#' @param cladeSize Integer giving the number of taxa in the focal clade.
#' @param armRung Integer giving how deep in its arm the focal clade sits: the
#'   clade is grafted as sister to the arm's outermost `nArm - armRung + 1`
#'   taxa, so larger values sit further from the central edge.
#' @param charsPerTree Integer giving the number of characters simulated on each
#'   of the two topologies.
#' @param stateCounts Integer vector of state-space sizes, cycled over
#'   characters.
#' @param internalLength,pendantLength Numerics giving the generating branch
#'   length of every internal and every terminal edge; the defaults put roughly
#'   6 substitutions per character across the tree, and neither shortening them
#'   towards clean synapomorphies nor lengthening them produced a barrier.
#' @param balanceTol Numeric giving the largest proxy log-likelihood difference
#'   between the two topologies that is accepted.
#' @param maxAttempts Integer capping the number of simulation attempts; the
#'   least unbalanced attempt is returned if none meets `balanceTol`.
#' @param quiet Logical; if `TRUE`, suppress progress messages.
#'
#' @details
#' Both islands must hold comparable mass or the target is useless as a mixing
#' benchmark, and a single simulation replicate misses that badly: the proxy
#' log-likelihood difference between the two generating topologies has a
#' standard deviation of tens of nats, i.e. mass ratios of \eqn{e^{20}}. The
#' constructor therefore *rejection-samples whole datasets*, redrawing the
#' matrix until the two topologies are within `balanceTol` nats. Conditioning on
#' a property of the whole dataset leaves every retained character an untouched
#' draw from the generating process, so no character is selected for agreeing
#' with either tree.
#'
#' Trees A and B are exchanged by the backbone's left-right symmetry, so they
#' are isomorphic, their islands contain the same 12 topologies up to
#' relabelling, and the expected mass split is exactly even; only character
#' noise breaks it, which is what `balanceTol` removes.
#'
#' @return List with the simulated `matrix`, its `phyDat` and `MkPrimeData`
#' representations, ground-truth trees `treeA`/`treeB`, the `backbone`, every
#' distinct placement of the focal clade in `candidates`, island-diagnostic
#' split keys `keyArmL`/`keyArmR`, per-character `charSource`, and the `design`.
BimodalTarget <- function(seed = 20260812L,
                          nArm = 7L,
                          cladeSize = 4L,
                          armRung = 5L,
                          charsPerTree = 150L,
                          stateCounts = c(4L, 5L, 6L),
                          internalLength = 0.35,
                          pendantLength = 0.05,
                          balanceTol = 1.5,
                          maxAttempts = 400L,
                          quiet = FALSE) {
  .LoadMkPrime()
  stopifnot(nArm >= 4L, cladeSize >= 3L, armRung >= 2L, armRung <= nArm)

  design <- list(seed = seed, nArm = nArm, cladeSize = cladeSize,
                 armRung = armRung, charsPerTree = charsPerTree,
                 stateCounts = stateCounts, internalLength = internalLength,
                 pendantLength = pendantLength, balanceTol = balanceTol,
                 armLTips = sprintf("b%02d", seq_len(nArm)),
                 # Reversed so armLTips[j] and armRTips[j] mirror each other.
                 armRTips = rev(sprintf("b%02d", nArm + seq_len(nArm))),
                 cladeTips = sprintf("c%d", seq_len(cladeSize)))
  tipLabels <- c(design$armLTips, rev(design$armRTips), design$cladeTips)

  backbone <- .BackboneTree(design)
  target <- list(design = design, tipLabels = tipLabels,
                 armLTips = design$armLTips, armRTips = design$armRTips,
                 cladeTips = design$cladeTips, backbone = backbone)
  target$treeA <- .AttachClade(
    backbone, .ArmNode(backbone, design$armLTips[armRung:nArm]), design)
  target$treeB <- .AttachClade(
    backbone, .ArmNode(backbone, design$armRTips[armRung:nArm]), design)
  target$candidates <- .CandidateTrees(design, backbone)
  # An arm stays monophyletic exactly while the focal clade is out of it.
  target$keyArmL <- .MembershipKey(design$armLTips, tipLabels)
  target$keyArmR <- .MembershipKey(design$armRTips, tipLabels)

  best <- NULL
  for (attempt in seq_len(maxAttempts)) {
    set.seed(seed + 1000L * attempt)
    mat <- cbind(.SimMatrix(target$treeA, stateCounts, charsPerTree, tipLabels),
                 .SimMatrix(target$treeB, stateCounts, charsPerTree, tipLabels))
    mat <- apply(mat, 2L, .Canonicalize)
    rownames(mat) <- tipLabels
    charSource <- rep(c("A", "B"), each = charsPerTree)
    variable <- apply(mat, 2L, function(v) length(unique(v))) >= 2L
    mat <- mat[, variable, drop = FALSE]
    charSource <- charSource[variable]
    colnames(mat) <- sprintf("%s%03d", charSource, seq_len(ncol(mat)))

    mkd <- MkPrimeData(TreeTools::MatrixToPhyDat(mat))
    delta <- .ProxyDelta(mkd, target)
    if (is.null(best) || abs(delta) < abs(best$proxyDelta)) {
      best <- list(matrix = mat, mkd = mkd, charSource = charSource,
                   proxyDelta = delta, attempt = attempt)
    }
    if (abs(delta) <= balanceTol) {
      break
    }
  }
  if (!quiet) {
    message(sprintf(
      "BimodalTarget: %d chars (%d variable), proxy dlogL = %+.2f after %d attempt(s)",
      2L * charsPerTree, ncol(best$matrix), best$proxyDelta, best$attempt))
  }

  target$matrix <- best$matrix
  target$phyDat <- TreeTools::MatrixToPhyDat(best$matrix)
  target$mkd <- best$mkd
  target$charSource <- best$charSource
  target$proxyDelta <- best$proxyDelta
  target$attempt <- best$attempt
  target$balanced <- abs(best$proxyDelta) <= balanceTol
  # Return:
  target
}

#' Cheap barrier estimate for a constructed target
#'
#' Proxy log-likelihood (see [BimodalTarget()] Details) at the two generating
#' trees and at every non-island placement of the focal clade, maximised over a
#' 3-parameter branch-length family including the clade's own stem. Runs in
#' seconds; use it to screen designs, not to report a barrier.
#'
#' @param target List returned by [BimodalTarget()].
#' @return List with the per-candidate proxy `logL`, island labels, and the
#' implied `barrier`.
ProxyBarrier <- function(target) {
  isl <- vapply(target$candidates, function(tr) TargetIslands(tr, target),
                character(1))
  ll <- vapply(target$candidates, .ProxyMaxLogL, numeric(1),
               mkd = target$mkd, target = target)
  peakA <- max(ll[isl == "A"])
  peakB <- max(ll[isl == "B"])
  # Return:
  list(logL = ll, island = isl, peakA = peakA, peakB = peakB,
       crossing = max(ll[!isl %in% c("A", "B")]),
       barrier = min(peakA, peakB) - max(ll[!isl %in% c("A", "B")]))
}

# ------------------------------------------------------------- MCMC harnesses

.RunFixed <- function(target, tree, nIter, warmup, thin, seed, maxTime) {
  set.seed(seed)
  suppressMessages(RunMkPrime(
    target$mkd, tree = tree, fixTopology = TRUE,
    mcmc = MkPrimeMCMC(nIter = nIter, minWarmup = warmup, maxWarmup = warmup,
                       thin = thin, treeThin = thin, nRuns = 1L, nChains = 1L,
                       autoTune = FALSE, maxTime = maxTime)
  ))
}

.RunFree <- function(target, tree, nIter, warmup, thin, seed, maxTime,
                     plainMh = FALSE, nChains = 1L) {
  set.seed(seed)
  args <- list(nIter = nIter, minWarmup = warmup, maxWarmup = warmup,
               thin = thin, treeThin = thin, nRuns = 1L, nChains = nChains,
               autoTune = FALSE, maxTime = maxTime)
  if (plainMh) {
    args <- c(args, list(gibbsSpr = FALSE, gibbsSubtreeSwap = FALSE,
                         weightedSpr = FALSE, weightedSubtreeSwap = FALSE))
  }
  suppressMessages(RunMkPrime(target$mkd, tree = tree,
                              mcmc = do.call(MkPrimeMCMC, args)))
}

.IslandSeries <- function(post, target) {
  isl <- TargetIslands(post$trees, target)
  runs <- rle(isl[isl %in% c("A", "B")])
  list(islands = isl,
       nTrees = length(isl),
       fracA = mean(isl == "A"),
       fracB = mean(isl == "B"),
       fracCentre = mean(isl == "centre"),
       fracOther = mean(isl == "other"),
       # Switches count A-to-B transitions, ignoring excursions through the
       # crossing topology, so a transient visit is not double-counted.
       nSwitches = max(length(runs$values) - 1L, 0L),
       first = isl[[1]], last = isl[[length(isl)]])
}

.TopologyProposals <- function(post, nSampleIter) {
  w <- post$moveWeights
  if (is.null(w)) {
    return(NA_real_)
  }
  topo <- c("nni", "spr", "pspr", "tbr", "gibbs_spr", "gibbs_subtree_swap",
            "weighted_spr", "weighted_subtree_swap")
  sum(w[intersect(names(w), topo)], na.rm = TRUE) * nSampleIter
}

# ------------------------------------------------------------- verification

.VerifyClassifier <- function(target) {
  d <- target$design
  bb <- target$backbone
  outer <- d$armLTips[seq.int(d$nArm - 1L, d$nArm)]
  outerR <- d$armRTips[seq.int(d$nArm - 1L, d$nArm)]
  cases <- list(
    armL_outer_cherry = list(.ArmNode(bb, outer), "A"),
    armL_deep = list(.ArmNode(bb, d$armLTips[2:d$nArm]), "A"),
    armL_tip = list(d$armLTips[[1]], "A"),
    armR_outer_cherry = list(.ArmNode(bb, outerR), "B"),
    armR_tip = list(d$armRTips[[1]], "B"),
    centre = list(.ArmNode(bb, d$armLTips), "centre")
  )
  got <- vapply(cases, function(cs) {
    TargetIslands(.AttachClade(bb, cs[[1]], d), target)
  }, character(1))
  expected <- vapply(cases, `[[`, character(1), 2L)

  # Break the clade up: one of its tips is moved into the other arm, so
  # neither arm is monophyletic.
  scrambled <- TreeTools::DropTip(target$treeA, d$cladeTips[[1]])
  scrambled <- TreeTools::AddTip(scrambled, d$armRTips[[1]], d$cladeTips[[1]])
  gotScrambled <- TargetIslands(scrambled, target)

  # The known hazard: root-dependent scoring. Rooting the same unrooted tree
  # three different ways must not change its island.
  rootLabels <- c(d$armLTips[[1]], d$armRTips[[1]], d$cladeTips[[1]])
  rootInvariant <- vapply(names(cases), function(nm) {
    tr <- .AttachClade(bb, cases[[nm]][[1]], d)
    lab <- vapply(rootLabels, function(l) {
      TargetIslands(TreeTools::RootTree(tr, l), target)
    }, character(1))
    length(unique(lab)) == 1L && lab[[1]] == got[[nm]]
  }, logical(1))

  islandSizes <- table(vapply(target$candidates,
                              function(tr) TargetIslands(tr, target),
                              character(1)))
  list(pass = all(got == expected) && gotScrambled == "other" &&
         all(rootInvariant) &&
         identical(as.integer(islandSizes[["A"]]),
                   as.integer(islandSizes[["B"]])),
       expected = expected, got = got, scrambled = gotScrambled,
       rootInvariant = rootInvariant, islandSizes = islandSizes,
       rfAB = TreeDist::RobinsonFoulds(target$treeA, target$treeB,
                                       normalize = FALSE))
}

.VerifyValley <- function(target, nIter, warmup, thin, seed, maxTime) {
  cand <- target$candidates
  isl <- vapply(cand, function(tr) TargetIslands(tr, target), character(1))
  rows <- lapply(seq_along(cand), function(i) {
    post <- .RunFixed(target, cand[[i]], nIter, warmup, thin, seed + i, maxTime)
    ll <- post$samples[, "log_likelihood"]
    stem <- vapply(post$trees, .CladeStemLength, numeric(1), target = target)
    # `fixTopology` cannot change the topology, so every sampled tree must
    # carry the island of the tree supplied: a free check that
    # `RunMkPrime(tree = )` is honoured at all.
    data.frame(candidate = names(cand)[[i]], island = isl[[i]],
               startHonoured = all(TargetIslands(post$trees, target) == isl[[i]]),
               maxLogL = max(ll), meanLogL = mean(ll), sdLogL = stats::sd(ll),
               meanTreeLength = mean(post$samples[, "tree_length"]),
               meanCladeStem = mean(stem, na.rm = TRUE),
               nSample = length(ll), stringsAsFactors = FALSE)
  })
  prof <- do.call(rbind, rows)
  peakA <- max(prof$maxLogL[prof$island == "A"])
  peakB <- max(prof$maxLogL[prof$island == "B"])
  cross <- !prof$island %in% c("A", "B")
  depth <- min(peakA, peakB) - max(prof$maxLogL[cross])
  list(pass = depth >= 10 && all(prof$startHonoured),
       profile = prof, valleyDepth = depth, startHonoured = all(prof$startHonoured),
       peakA = peakA, peakB = peakB,
       crossingLogL = max(prof$maxLogL[cross]),
       stemIsland = mean(prof$meanCladeStem[!cross]),
       stemCrossing = mean(prof$meanCladeStem[cross]))
}

.VerifySplit <- function(target, nStones, nIter, warmup, nRep, seed) {
  fit <- function(tree, rep) {
    set.seed(seed + rep)
    suppressMessages(mkp_stepping_stone(
      target$mkd, tree = tree, nStones = nStones, nIter = nIter,
      warmup = warmup, fixTopology = TRUE, verbose = FALSE))$log_marginal
  }
  zA <- vapply(seq_len(nRep), function(r) fit(target$treeA, r), numeric(1))
  zB <- vapply(seq_len(nRep), function(r) fit(target$treeB, r), numeric(1))
  delta <- mean(zA) - mean(zB)
  se <- if (nRep > 1L) {
    sqrt(stats::var(zA) / nRep + stats::var(zB) / nRep)
  } else NA_real_
  list(pass = abs(delta) <= 3, logZA = zA, logZB = zB, delta = delta, se = se,
       massFracA = 1 / (1 + exp(-delta)))
}

.VerifyHops <- function(target, nIter, warmup, thin, seeds, maxTime) {
  grid <- expand.grid(start = c("A", "B"), plainMh = c(TRUE, FALSE),
                      seed = seeds, stringsAsFactors = FALSE)
  rows <- lapply(seq_len(nrow(grid)), function(i) {
    tree <- if (grid$start[[i]] == "A") target$treeA else target$treeB
    post <- .RunFree(target, tree, nIter, warmup, thin, grid$seed[[i]],
                     maxTime, plainMh = grid$plainMh[[i]])
    s <- .IslandSeries(post, target)
    data.frame(start = grid$start[[i]], plainMh = grid$plainMh[[i]],
               seed = grid$seed[[i]], nTrees = s$nTrees, firstIsland = s$first,
               fracA = s$fracA, fracB = s$fracB, fracCentre = s$fracCentre,
               fracOther = s$fracOther, nSwitches = s$nSwitches, last = s$last,
               # Only topology proposals can move a chain between islands, and
               # the default schedule spends most iterations elsewhere: a hop
               # rate per iteration would flatter the barrier.
               topoProposals = .TopologyProposals(post, nIter - warmup),
               stringsAsFactors = FALSE)
  })
  hop <- do.call(rbind, rows)
  plain <- hop[hop$plainMh, ]
  # Fit-for-purpose: plain-MH chains must not leak between islands, or the
  # target cannot discriminate schedules.
  list(pass = all(plain$nSwitches == 0L) &&
         all(plain$firstIsland == plain$start),
       table = hop, plainSwitches = sum(plain$nSwitches),
       defaultSwitches = sum(hop$nSwitches[!hop$plainMh]),
       plainHopsPerTopoProposal =
         sum(plain$nSwitches) / sum(plain$topoProposals),
       defaultHopsPerTopoProposal =
         sum(hop$nSwitches[!hop$plainMh]) / sum(hop$topoProposals[!hop$plainMh]))
}

.VerifyReach <- function(target, nIter, warmup, thin, seeds, maxTime, nChains) {
  rows <- lapply(seeds, function(s) {
    set.seed(s)
    start <- .SetLengths(TreeTools::RandomTree(target$tipLabels, root = FALSE),
                         target$design$internalLength,
                         target$design$pendantLength)
    post <- .RunFree(target, start, nIter, warmup, thin, s, maxTime,
                     nChains = nChains)
    st <- .IslandSeries(post, target)
    data.frame(seed = s, nChains = nChains, fracA = st$fracA,
               fracB = st$fracB, fracCentre = st$fracCentre,
               fracOther = st$fracOther, nSwitches = st$nSwitches,
               last = st$last, stringsAsFactors = FALSE)
  })
  reach <- do.call(rbind, rows)
  list(pass = any(reach$fracA > 0) && any(reach$fracB > 0),
       table = reach,
       islandsReached = unique(c(if (any(reach$fracA > 0)) "A",
                                 if (any(reach$fracB > 0)) "B")))
}

.VerifySeeds <- function(seeds, design, quiet = TRUE) {
  rows <- lapply(seeds, function(s) {
    tg <- BimodalTarget(seed = s, nArm = design$nArm,
                        cladeSize = design$cladeSize, armRung = design$armRung,
                        charsPerTree = design$charsPerTree,
                        stateCounts = design$stateCounts,
                        internalLength = design$internalLength,
                        pendantLength = design$pendantLength,
                        balanceTol = design$balanceTol, quiet = quiet)
    pb <- ProxyBarrier(tg)
    data.frame(seed = s, attempts = tg$attempt, nChar = ncol(tg$matrix),
               proxyDelta = tg$proxyDelta, balanced = tg$balanced,
               proxyBarrier = pb$barrier, stringsAsFactors = FALSE)
  })
  seedTab <- do.call(rbind, rows)
  list(pass = all(seedTab$balanced) && all(seedTab$proxyBarrier > 10),
       table = seedTab)
}

#' Verify the bimodal benchmark target
#'
#' @param target List returned by [BimodalTarget()].
#' @param stages Character vector naming stages to run: `classifier`, `valley`,
#'   `split`, `hops`, `reach`, `seeds`.
#' @param quick Logical; if `TRUE`, run every stage at execution-smoke size, in
#'   which case the verdict is an execution check and not evidence.
#' @param outDir Character path for `summary.rds` and `verdict.txt`.
#' @return List of per-stage results, invisibly.
VerifyBimodalTarget <- function(target = BimodalTarget(),
                                stages = c("classifier", "valley", "split",
                                           "hops", "reach", "seeds"),
                                quick = FALSE,
                                outDir = "dev/benchmarks/bimodal-target-results") {
  .LoadMkPrime()
  dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
  cfg <- if (quick) {
    list(valleyIter = 400L, valleyWarmup = 200L, valleyThin = 20L,
         ssStones = 6L, ssIter = 150L, ssWarmup = 50L, ssRep = 2L,
         hopIter = 1200L, hopWarmup = 400L, hopThin = 10L,
         hopSeeds = 7001L, reachIter = 1200L, reachSeeds = 7101:7102,
         seedList = 20260813L, maxTime = 90)
  } else {
    list(valleyIter = 3000L, valleyWarmup = 1200L, valleyThin = 20L,
         ssStones = 25L, ssIter = 800L, ssWarmup = 300L, ssRep = 3L,
         hopIter = 20000L, hopWarmup = 2000L, hopThin = 10L,
         hopSeeds = c(7001L, 7002L), reachIter = 15000L,
         reachSeeds = 7101:7103, seedList = 20260813:20260816,
         maxTime = 1800)
  }
  res <- list(design = target$design, quick = quick,
              proxyDelta = target$proxyDelta, attempt = target$attempt,
              nChar = ncol(target$matrix), startedAt = Sys.time())

  Report <- function(nm, x) {
    cat(sprintf("[%s] %s\n", if (isTRUE(x$pass)) "PASS" else "FAIL", nm))
    x
  }
  if ("classifier" %in% stages) {
    res$classifier <- Report("classifier", .VerifyClassifier(target))
    cat(sprintf("      RF(A,B) = %d; island sizes %s\n",
                res$classifier$rfAB,
                paste(names(res$classifier$islandSizes),
                      res$classifier$islandSizes, sep = "=", collapse = " ")))
  }
  if ("valley" %in% stages) {
    res$valley <- Report("valley", .VerifyValley(
      target, cfg$valleyIter, cfg$valleyWarmup, cfg$valleyThin, 4100L,
      cfg$maxTime))
    cat(sprintf(
      "      valley depth = %.1f nats (peak A %.1f, B %.1f, crossing %.1f); clade stem %.2f island vs %.2f crossing\n",
      res$valley$valleyDepth, res$valley$peakA, res$valley$peakB,
      res$valley$crossingLogL, res$valley$stemIsland, res$valley$stemCrossing))
  }
  if ("split" %in% stages) {
    res$split <- Report("split", .VerifySplit(
      target, cfg$ssStones, cfg$ssIter, cfg$ssWarmup, cfg$ssRep, 4200L))
    cat(sprintf("      dlogZ = %+.2f (SE %.2f) => P(island A) = %.3f\n",
                res$split$delta, res$split$se, res$split$massFracA))
  }
  if ("hops" %in% stages) {
    res$hops <- Report("hops", .VerifyHops(
      target, cfg$hopIter, cfg$hopWarmup, cfg$hopThin, cfg$hopSeeds,
      cfg$maxTime))
    print(res$hops$table)
  }
  if ("reach" %in% stages) {
    res$reach <- Report("reach", .VerifyReach(
      target, cfg$reachIter, cfg$hopWarmup, cfg$hopThin, cfg$reachSeeds,
      cfg$maxTime, nChains = 4L))
    print(res$reach$table)
  }
  if ("seeds" %in% stages) {
    res$seeds <- Report("seeds", .VerifySeeds(cfg$seedList, target$design))
    print(res$seeds$table)
  }

  res$finishedAt <- Sys.time()
  verdicts <- vapply(res[intersect(names(res), c("classifier", "valley",
                                                 "split", "hops", "reach",
                                                 "seeds"))],
                     function(x) isTRUE(x$pass), logical(1))
  overall <- if (quick) {
    "WARN"  # quick mode is an execution check; thresholds are not powered
  } else if (all(verdicts)) {
    "PASS"
  } else {
    "FAIL"
  }
  res$verdict <- overall
  saveRDS(res, file.path(outDir, "summary.rds"))
  writeLines(c(
    sprintf("VERDICT: %s", overall),
    sprintf("mode: %s",
            if (quick) "quick (execution smoke; NOT powered)" else "full"),
    sprintf("run: %s -- %s", format(res$startedAt), format(res$finishedAt)),
    sprintf("design: %d taxa (2 x %d backbone + %d clade), armRung %d, %d chars",
            length(target$tipLabels), target$design$nArm,
            target$design$cladeSize, target$design$armRung,
            ncol(target$matrix)),
    sprintf("balance: proxy dlogL = %+.2f after %d attempt(s)",
            target$proxyDelta, target$attempt),
    "",
    sprintf("  %-11s %s", names(verdicts), ifelse(verdicts, "PASS", "FAIL"))
  ), file.path(outDir, "verdict.txt"))
  cat(sprintf("\nVERDICT: %s  (%s)\n", overall,
              file.path(outDir, "verdict.txt")))
  invisible(res)
}

# --------------------------------------------------------------------- CLI

.Main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  quick <- any(args %in% c("--quick", "-q"))
  ix <- which(args == "--stage")
  stages <- if (length(ix) && length(args) > ix[[1]]) {
    trimws(strsplit(args[[ix[[1]] + 1L]], ",")[[1]])
  } else {
    c("classifier", "valley", "split", "hops", "reach", "seeds")
  }
  ix <- which(args == "--out")
  outDir <- if (length(ix) && length(args) > ix[[1]]) {
    args[[ix[[1]] + 1L]]
  } else {
    "dev/benchmarks/bimodal-target-results"
  }
  .LoadMkPrime()
  # Cache the constructed target so staged invocations verify one dataset.
  dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
  cache <- file.path(outDir, if (quick) "target-quick.rds" else "target.rds")
  target <- if (file.exists(cache)) {
    readRDS(cache)
  } else {
    tg <- if (quick) {
      BimodalTarget(charsPerTree = 40L, balanceTol = 3, maxAttempts = 40L)
    } else {
      BimodalTarget()
    }
    saveRDS(tg, cache)
    tg
  }
  res <- VerifyBimodalTarget(target, stages = stages, quick = quick,
                             outDir = outDir)
  # Staged runs each write their own summary; merge so a later reader sees
  # every stage that has been run against this target.
  merged <- file.path(outDir, "summary-merged.rds")
  prev <- if (file.exists(merged)) readRDS(merged) else list()
  prev[names(res)] <- res
  saveRDS(prev, merged)
  invisible(res)
}

.thisFile <- sub("^--file=", "",
                 grep("^--file=", commandArgs(), value = TRUE))
if (!interactive() && length(.thisFile) &&
    identical(basename(.thisFile[[1]]), "bimodal-target.R")) {
  .Main()
}
