# parsimony-trap.R --------------------------------------------------------------
#
# Parsimony trap test: do equal-weights (EW) and implied-weights (IW) parsimony
# searches recover the FALSE eco-clade on existing MkPrime simulated datasets?
#
# Hypothesis. If MaximizeParsimony recovers the false eco-clade rather than
# truth on a given dataset, then ML can plausibly be fooled by the same data.
# If parsimony robustly recovers truth even when eco signal is strong,
# no "trap" exists for likelihood either.
#
# Datasets tested
#   v6-realistic rep01-03  (TL=1.16, nChar=~194-200, 8-clade eco)
#   multirep-v3  rep01-08  (TL=13.5,  nChar=480,     8-clade eco)
#   v5break      rep01-03  (TL=1.30, nChar=80,       inner-cherry eco)
#
# Weightings: concavity = Inf (equal weights), 3, 6, 12 (implied weighting).
#
# Per (dataset, rep, weighting) we record:
#   - MPT score (TreeSearch::TreeLength / extended IW)
#   - Whether *any* MPT contains the true sister bipartition
#   - Whether *any* MPT contains the false eco bipartition
#   - Truth-tree score vs. MPT score (Fitch steps, equal weights)
#
# All bipartition checks go through HasBipartSplits (root-invariant). MPT
# searches use 3 random starts; we keep ALL trees within tolerance of the
# best (i.e. the full MP set returned by TreeSearch).
#
# Writes:
#   parsimony-trap-results.csv   (per-rep / per-weighting row)
#   parsimony-trap-summary.csv   (per-dataset / per-weighting aggregate)
#   parsimony-trap-report.md     (narrative)
#

suppressPackageStartupMessages({
  library("TreeTools")
  library("TreeSearch")
  library("ape")
})

scriptDir <- "inst/scripts/parsimony-trap-test"
source("inst/simulations/ecology/sim3-scoring.R")

# Bipartition definitions ------------------------------------------------------
TRUE_SISTER_AC <- c("A1","A2","A3","A4","C1","C2","C3","C4")   # truth = ((A,C),(B,D))
FALSE_ECO_8C   <- c("A1","A2","A3","A4","B1","B2","B3","B4")   # multirep-v3 / v6
FALSE_ECO_4C   <- c("A1","A2","B1","B2")                       # v5break inner cherry

# Datasets ---------------------------------------------------------------------
datasets <- list(
  list(name = "v6-realistic", reps = 1:3,
       path = "inst/simulations/ecology/v6-realistic/rep%02d/blind-result.rds",
       falseEco = FALSE_ECO_8C, TL = 1.16),
  list(name = "multirep-v3",  reps = 1:8,
       path = "inst/simulations/ecology/multirep-v3-results/rep%02d/blind-result.rds",
       falseEco = FALSE_ECO_8C, TL = 13.5),
  list(name = "v5break",      reps = 1:3,
       path = "inst/scripts/parsimony-trap-test/v5break-data/v5break-rep%02d.rds",
       falseEco = FALSE_ECO_4C, TL = 1.30)
)

# Truth tree (same Newick topology for all three datasets — only TL varies,
# which is irrelevant for MP scoring).
.TruthTree <- function(tipLabels) {
  tr <- ape::read.tree(
    text = "((((A1,A2),(A3,A4)),((C1,C2),(C3,C4))),(((B1,B2),(B3,B4)),((D1,D2),(D3,D4))));"
  )
  tr$edge.length <- NULL
  tr
}

weightings <- c(Inf, 3, 6, 12)
nStarts    <- 3L

# Result accumulator -----------------------------------------------------------
allRows <- list()

set.seed(20260520L)

for (ds in datasets) {
  for (rep in ds$reps) {
    f <- sprintf(ds$path, rep)
    if (!file.exists(f)) {
      message("MISSING: ", f); next
    }
    cat(sprintf("\n=== %s rep%02d ===\n", ds$name, rep))
    res <- readRDS(f)
    pd  <- res$data$phyDat
    if (is.null(pd)) {
      message("  no phyDat in result; skipping"); next
    }
    tipLab <- res$data$taxon_names
    nTip   <- res$data$nTip
    nChar  <- res$data$nChar

    truth <- .TruthTree(tipLab)
    # Truth-tree Fitch length on this dataset (equal weights)
    truthLen <- as.numeric(TreeLength(truth, pd, concavity = Inf))

    for (cv in weightings) {
      # Run MaximizeParsimony with nStarts random starting trees, take union
      mpAll <- list()
      bestScore <- Inf
      for (s in seq_len(nStarts)) {
        randTree <- Preorder(ape::rtree(nTip, tip.label = tipLab))
        out <- tryCatch(
          MaximizeParsimony(pd, tree = randTree, concavity = cv,
                            verbosity = 0L,
                            maxReplicates = 24L,
                            targetHits = 10L),
          error = function(e) { message(" MP error: ", conditionMessage(e)); NULL })
        if (is.null(out)) next
        if (inherits(out, "phylo")) out <- list(out)
        for (tr in out) {
          sc <- as.numeric(TreeLength(tr, pd, concavity = cv))
          if (sc < bestScore - 1e-8) {
            mpAll    <- list(tr); bestScore <- sc
          } else if (abs(sc - bestScore) <= 1e-8) {
            mpAll[[length(mpAll) + 1L]] <- tr
          }
        }
      }
      if (length(mpAll) == 0L) {
        cat(sprintf("  cv=%s: NO TREES RETURNED\n", format(cv))); next
      }
      class(mpAll) <- "multiPhylo"
      # Bipartition presence on the MPT set
      hasTrue  <- any(HasBipartSplits(mpAll, TRUE_SISTER_AC))
      hasFalse <- any(HasBipartSplits(mpAll, ds$falseEco))
      # Per-tree counts for richer reporting
      pTrue  <- mean(HasBipartSplits(mpAll, TRUE_SISTER_AC))
      pFalse <- mean(HasBipartSplits(mpAll, ds$falseEco))
      # Equal-weights length of best MPT and of truth, for delta-steps
      bestEW   <- as.numeric(TreeLength(mpAll[[1]], pd, concavity = Inf))

      cat(sprintf(
        "  cv=%-5s nMPT=%-3d  bestScore=%10.3f  truthLen=%10.3f  pTrueSister=%.2f  pFalseEco=%.2f\n",
        format(cv), length(mpAll), bestScore, truthLen, pTrue, pFalse))

      allRows[[length(allRows) + 1L]] <- data.frame(
        dataset     = ds$name,
        rep         = rep,
        TL          = ds$TL,
        nTip        = nTip,
        nChar       = nChar,
        concavity   = ifelse(is.infinite(cv), "Inf", as.character(cv)),
        nMPT        = length(mpAll),
        bestScore   = bestScore,
        bestEW      = bestEW,
        truthLenEW  = truthLen,
        deltaStepsEW = truthLen - bestEW,
        anyMPT_hasTrueSister = hasTrue,
        anyMPT_hasFalseEco   = hasFalse,
        pTrueSister = pTrue,
        pFalseEco   = pFalse,
        stringsAsFactors = FALSE
      )
    }
  }
}

results <- do.call(rbind, allRows)
write.csv(results,
          file = file.path(scriptDir, "parsimony-trap-results.csv"),
          row.names = FALSE)

# Per-dataset / per-weighting summary -----------------------------------------
summ <- aggregate(
  cbind(anyMPT_hasFalseEco, anyMPT_hasTrueSister,
        pFalseEco, pTrueSister, deltaStepsEW) ~ dataset + concavity,
  data = results,
  FUN  = function(x) mean(as.numeric(x))
)
summ$nReps <- aggregate(rep ~ dataset + concavity, data = results,
                        FUN = function(x) length(unique(x)))$rep
summ <- summ[order(summ$dataset, summ$concavity), ]
write.csv(summ,
          file = file.path(scriptDir, "parsimony-trap-summary.csv"),
          row.names = FALSE)

cat("\n\n=== Summary (false-eco-clade hit fraction by dataset x weighting) ===\n")
print(summ, row.names = FALSE)
cat("\nWritten:\n  ",
    file.path(scriptDir, "parsimony-trap-results.csv"), "\n  ",
    file.path(scriptDir, "parsimony-trap-summary.csv"), "\n")
