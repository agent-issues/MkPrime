# Mirror of pspr_proposal_impl steps 1-4: locate the rows the move rewrites and
# enumerate its candidate regraft edges.  Used only to build valid inputs for
# the compiled scorer, so it insists on a bifurcating u rather than resolving a
# polytomy the way the C++ loop does.
.PsprSites <- function(edge, nTip, pruneRow) {
  par <- edge[, 1]
  chi <- edge[, 2]
  u <- par[pruneRow]
  v <- chi[pruneRow]
  parentRow <- which(chi == u)
  sibRow <- which(par == u & chi != v)
  if (length(parentRow) != 1L || length(sibRow) != 1L) {
    return(NULL)
  }
  desc <- v
  repeat {
    extra <- setdiff(chi[par %in% desc], desc)
    if (!length(extra)) break
    desc <- c(desc, extra)
  }
  cand <- which(!(chi %in% desc) & par != u & chi != u)
  if (!length(cand)) {
    return(NULL)
  }
  # Return:
  list(u = u, v = v, pruneRow = pruneRow, parentRow = parentRow,
       sibRow = sibRow, sibNode = chi[sibRow], candidates = cand)
}

# The tree that regrafting u onto candidate row `rr` produces.
.PsprRegraft <- function(edge, sites, rr) {
  out <- edge
  b <- edge[rr, 2]
  out[sites$parentRow, 2] <- sites$sibNode
  out[rr, 2] <- sites$u
  out[sites$sibRow, ] <- c(sites$u, b)
  out[sites$pruneRow, ] <- c(sites$u, sites$v)
  # Return:
  out
}
