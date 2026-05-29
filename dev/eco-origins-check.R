# How many INDEPENDENT origins does each ecology state have on the blind tree?
# Aggregate 17 steps (60 tips, min 3) suggests dispersal; the per-state count
# determines identifiability of an ecology-> trait effect for each regime.
suppressPackageStartupMessages({library(ape); library(TreeTools); library(phangorn)})
od <- "inst/ecology/scripts/rodent-comparison"
bl <- do.call(c, lapply(Sys.glob(file.path(od, "rodent-MkNT-blind_trees_*.nwk")), ape::read.tree))
class(bl) <- "multiPhylo"; bl <- bl[seq.int(ceiling(length(bl)/4)+1L, length(bl))]
mat <- TreeTools::ReadCharacters("inst/ecology/data/rodent-X24848.nex")
eco <- mat[,220]; extant <- mat[,221]; keep <- which(extant=="0")
e <- eco[keep]; nm <- rownames(mat)[keep]
poly <- grepl("^[(]", e); e[poly] <- substr(sub("^[(]","",e[poly]),1,1)
ok <- !is.na(e)&e!="?"&!is.na(suppressWarnings(as.integer(e))); e<-as.integer(e[ok]); nm<-nm[ok]; names(e)<-nm
consB <- ape::multi2di(ape::consensus(bl,p=0.5,rooted=FALSE),random=FALSE)
consB$edge.length <- rep(1,nrow(consB$edge)); consB <- TreeTools::Preorder(consB)
ed <- e[consB$tip.label]; lev <- sort(unique(e))

## ML ancestral reconstruction (ER), MAP state per node, count transitions into each state
fit <- ape::ace(ed, consB, type="discrete", model="ER")
anc <- apply(fit$lik.anc, 1, function(r) lev[which.max(r)])      # internal nodes
tipv <- ed
nTip <- ape::Ntip(consB)
nodeState <- c(tipv[consB$tip.label], setNames(anc, (nTip+1):(nTip+consB$Nnode)))
# index node states by node id
st <- integer(nTip + consB$Nnode)
st[1:nTip] <- tipv[consB$tip.label]
st[(nTip+1):(nTip+consB$Nnode)] <- as.integer(anc)
into <- setNames(integer(length(lev)), lev)
ch <- integer(0)
for (i in seq_len(nrow(consB$edge))) {
  p <- st[consB$edge[i,1]]; c_ <- st[consB$edge[i,2]]
  if (!is.na(p) && !is.na(c_) && p != c_) into[as.character(c_)] <- into[as.character(c_)] + 1L
}
cat("Tip counts per ecology:        "); print(table(ed))
cat("Independent origins per ecology (edges entering state, ML-MAP recon):\n")
print(into)
cat(sprintf("Total transitions (ML-MAP) = %d   (parsimony min steps = %d)\n", sum(into), length(lev)-1L))

## Parsimony cross-check: total steps
pd <- phangorn::phyDat(matrix(as.character(ed),ncol=1,dimnames=list(consB$tip.label,NULL)),type="USER",levels=as.character(lev))
cat(sprintf("Parsimony steps (Fitch) = %d\n", phangorn::parsimony(consB, pd)))
