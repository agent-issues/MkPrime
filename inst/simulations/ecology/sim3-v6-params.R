# sim3-v6-params.R ------------------------------------------------------------
#
# v6-realistic parameter set for the convergent-ecology Sim 3 study.
#
# Motivation
# ----------
# multirep-v3 (truth TL = 13.5, 480 chars at phi = 4) and v5break (truth
# TL = 1.3 with weak eco confound) are now superseded:
#   * multirep-v3 is "junk" — saturated, no real-world interpretability.
#   * v5break truth TL was realistic but blind never failed; eco
#     confound was too weak to break ancestral signal.
#
# v6-realistic uses TL typical for a morphological matrix (~ 1.2) at
# 16 tips x 300 chars, with the v3 8-tip eco-1 topology (full clades A
# and B as ecology = 1). The new auto-expSteps default (commit bc24d52)
# inflates the parsimony score by 1.05x to set the tree-length prior, so
# the prior is on-target for these realistic regimes — not an order of
# magnitude too tight as the old expSteps = 10 default was for v3.
#
# Topology (same builder logic as multirep-v3 but with stemBrEco !=
# stemBrClade): ((A, C), (B, D)). Clades A and B (eco-1) sit on
# *longer* stems (stemBrEco) than clades C and D (stemBrClade). This
# is the only way to give eco-shared parallel substitutions a chance
# to dominate ancestral signal at a realistic global TL: the eco stems
# carry more substitution mass, and (under phi > 1) those subs are
# amplified into parallel synapomorphies linking A and B.
#
# Tree length (16 tips, balanced)
#   Per eco clade : 4 * tipBr (terminals) + 2 * tipBr (cherry-internals)
#                   + stemBrEco
#                = 6 * tipBr + stemBrEco
#   Per non-eco clade : 6 * tipBr + stemBrClade
#   Root        : 2 * rootBr
#   Total       : 12 * tipBr + 2 * stemBrEco + 2 * stemBrClade + 2 * rootBr
#
#   With v6d parameters (tipBr 0.030, stemBrClade 0.035, stemBrEco 0.15,
#   rootBr 0.035):
#     TL = 12 * 0.030 + 2 * 0.15 + 2 * 0.035 + 2 * 0.035
#        = 0.360 + 0.300 + 0.070 + 0.070
#        = 0.80 + 0.36     (terminals + internals/stems)
#        = 1.16
#
# Discriminating local-R Fitch check (5 seeds, run via
# dev/sim-design/v6-discriminate.R):
#   Mean parsimony(true)    : 283
#   Mean parsimony(AB-grouped, false clade) : 292
#   Mean parsimony(MP best) : 283
#   Mean (AB - true) step diff : 9.0
#   MP search recovers the true topology in 100 % of reps.
#   parsimony / (TL * nChar) ratio : 0.81 -- well below saturation.
#
# Interpretation: the AB-grouped tree costs 9 steps over 300 characters
# more than truth. That's parsimony-distinguishable but the *likelihood*
# under the ecology-blind Mk model amplifies eco-driven parallelisms
# beyond what parsimony counts (probabilistic priors on the high-rate
# eco stems concentrate probability mass on synapomorphic shifts), so
# blind MCMC may still wander into the AB-clade region. Aware should
# detect the eco amplification and re-rank — recovering truth.
#
# Predicted auto-expSteps under the new bc24d52 default: ~298 (1.05x
# the parsimony score). This is finally on-target -- old expSteps = 10
# would have put 95 % prior mass on tree lengths 10x below parsimony,
# producing pathological TL collapse.
#
# Use ONLY with the multirep-v3 8-tip-clade ecology assignment
# (.ConvergentEcology returns 1 for tips A1..A4 and B1..B4). See
# .BuildConvergentTreeV6() below for the topology + branch-length
# builder.

#' Build the v6-realistic convergent-clades tree.
#'
#' Eco-1 clades A and B get stem length stemBrEco; eco-0 clades C and D
#' get stem length stemBrClade. All other branches use tipBr (terminals
#' and within-clade cherry-internals). Root edges use rootBr.
#'
#' @param tipBr Numeric. Per-clade terminal & cherry-internal length.
#' @param stemBrClade Numeric. Stem length subtending eco-0 clades C, D.
#' @param stemBrEco Numeric. Stem length subtending eco-1 clades A, B.
#' @param rootBr Numeric. Length of the two root-side internal branches.
#'
#' @return A `phylo` object in Preorder with 16 labelled tips.
.BuildConvergentTreeV6 <- function(tipBr = 0.030,
                                   stemBrClade = 0.035,
                                   stemBrEco = 0.15,
                                   rootBr = 0.035) {
  stopifnot(tipBr > 0, stemBrClade > 0, stemBrEco > 0, rootBr > 0)
  cladeNw <- function(prefix, tb, sb) {
    sprintf(
      "((%s1:%g,%s2:%g):%g,(%s3:%g,%s4:%g):%g):%g",
      prefix, tb, prefix, tb, tb,
      prefix, tb, prefix, tb, tb, sb
    )
  }
  cladeA <- cladeNw("A", tipBr, stemBrEco)    # eco-1
  cladeB <- cladeNw("B", tipBr, stemBrEco)    # eco-1
  cladeC <- cladeNw("C", tipBr, stemBrClade)  # eco-0
  cladeD <- cladeNw("D", tipBr, stemBrClade)  # eco-0
  newick <- sprintf("((%s,%s):%g,(%s,%s):%g);",
                    cladeA, cladeC, rootBr,
                    cladeB, cladeD, rootBr)
  TreeTools::Preorder(ape::read.tree(text = newick))
}

SIM3V6_PARAMS <- list(
  v6 = list(
    name        = "v6_realistic",
    tipBr       = 0.030,
    stemBrClade = 0.035,
    stemBrEco   = 0.150,
    rootBr      = 0.035,
    nNeo        = 80L,
    nTrans      = 220L,
    phi         = 8,
    pi0         = 0.45,
    theta       = 1.0,
    baseRate    = 1.0,
    rateLoss    = 1,
    kStates     = 2L
  )
)
