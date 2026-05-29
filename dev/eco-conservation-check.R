# Quick check: is the ecology covariate phylogenetically conserved on the
# blind consensus tree? If yes, ecology-conditioned rate down-weighting is
# confounded with clade structure.
suppressPackageStartupMessages({library(ape); library(TreeTools); library(phangorn)})
od <- "inst/ecology/scripts/rodent-comparison"
bl <- do.call(c, lapply(Sys.glob(file.path(od, "rodent-MkNT-blind_trees_*.nwk")),
                        ape::read.tree))
class(bl) <- "multiPhylo"
disc <- function(L) L[seq.int(ceiling(length(L) / 4) + 1L, length(L))]
bl <- disc(bl)

mat <- TreeTools::ReadCharacters("inst/ecology/data/rodent-X24848.nex")
ecoVec <- mat[, 220]; extant <- mat[, 221]; keep <- which(extant == "0")
eco <- ecoVec[keep]; nm <- rownames(mat)[keep]
poly <- grepl("^[(]", eco)
eco[poly] <- substr(sub("^[(]", "", eco[poly]), 1, 1)
ok <- !is.na(eco) & eco != "?" & !is.na(suppressWarnings(as.integer(eco)))
eco <- as.integer(eco[ok]); nm <- nm[ok]; names(eco) <- nm

consB <- ape::consensus(bl, p = 0.5, rooted = FALSE)
consB <- multi2di(consB, random = FALSE)
consB$edge.length <- rep(1, nrow(consB$edge))
ed <- eco[consB$tip.label]
lev <- as.character(sort(unique(eco)))
mk_pd <- function(v) phangorn::phyDat(
  matrix(as.character(v), ncol = 1, dimnames = list(consB$tip.label, NULL)),
  type = "USER", levels = lev)

ps <- phangorn::parsimony(consB, mk_pd(ed))
mn <- length(lev) - 1L
cat(sprintf("Ecology states: %s   tip counts: %s\n",
            paste(lev, collapse = ","), paste(table(eco), collapse = "/")))
cat(sprintf("Ecology parsimony on BLIND consensus = %d steps (min possible = %d, CI = %.3f)\n",
            ps, mn, mn / ps))
set.seed(1)
rs <- replicate(500, phangorn::parsimony(consB, mk_pd(sample(ed))))
z <- (ps - mean(rs)) / sd(rs)
cat(sprintf("Random-shuffled ecology on same tree: mean=%.1f sd=%.1f range=[%d,%d]\n",
            mean(rs), sd(rs), min(rs), max(rs)))
cat(sprintf("Observed %d vs random %.1f  (z = %.1f)  => ecology is %s phylogenetically clustered\n",
            ps, mean(rs), z,
            ifelse(ps < mean(rs) - 2 * sd(rs), "STRONGLY", "NOT clearly")))
