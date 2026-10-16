# Q4: does the tuning bandit still reward frozen topology?
# Calls the production .MinEssRate() and .BeatsIncumbent() on synthetic
# 200-sample windows. Scalars are iid in every window, so only the trees
# differ. Run: Rscript dev/red-team/explore-195/q4-bandit.R
source("dev/red-team/explore-195/common.R")
set.seed(195)

n <- 200L
wall <- 10
pool <- lapply(1:50, function(i) RandomTree(8, root = TRUE))
poolKey <- TopoKeys(pool)

Sticky <- function(p) {
  idx <- integer(n); idx[1] <- 1L
  for (i in 2:n) idx[i] <- if (runif(1) < p) idx[i - 1] else sample(50, 1)
  idx
}
windows <- list(
  "iid topologies" = sample(50, n, replace = TRUE),
  "sticky, stay 0.9" = Sticky(0.9),
  "sticky, stay 0.99" = Sticky(0.99),
  "frozen" = rep(1L, n),
  "frozen, one spike at end" = c(rep(1L, n - 1), 2L),
  "frozen, one spike mid" = replace(rep(1L, n), n / 2, 2L)
)

scalars <- matrix(rnorm(n * 2), n, 2,
                  dimnames = list(NULL, c("tree_length", "rate_log_sd")))

rows <- lapply(names(windows), function(nm) {
  idx <- windows[[nm]]
  trees <- pool[idx]
  hash <- match(poolKey[idx], unique(poolKey))
  withHash <- cbind(scalars, topo_hash = hash)
  raw <- mkp$TreeESS(structure(trees, class = "multiPhylo"),
                     dist_fn = TreeDist::RobinsonFoulds,
                     frechet = FALSE)[["medianPseudoESS"]]
  a <- mkp$.MinEssRate(scalars, wall, trees)
  b <- mkp$.MinEssRate(withHash, wall, trees)
  data.frame(window = nm, changes = sum(diff(hash) != 0),
             rawPseudoESS = round(raw, 1),
             rateNoHash = round(a$rate, 2), frozenNoHash = a$frozen,
             rateProd = round(b$rate, 2), frozenProd = b$frozen)
})
tab <- do.call(rbind, rows)
print(tab, row.names = FALSE)

# Head-to-head: does a frozen candidate displace a moving incumbent?
Score <- function(nm) {
  idx <- windows[[nm]]
  hash <- match(poolKey[idx], unique(poolKey))
  mkp$.MinEssRate(cbind(scalars, topo_hash = hash), wall, pool[idx])
}
H2H <- function(cand, best, candScalarSpeedup = 1, everMoved = FALSE,
                topologyCut = FALSE) {
  cs <- Score(cand); bs <- Score(best)
  mkp$.BeatsIncumbent(cs$rate * candScalarSpeedup, cs$ess, bs$rate, bs$ess,
                      candFrozen = cs$frozen, bestFrozen = bs$frozen,
                      topologyCut = topologyCut, everMoved = everMoved)
}
h2h <- data.frame(
  candidate = c("frozen", "frozen (2x faster)", "frozen, one spike at end",
                "frozen, one spike mid", "sticky, stay 0.99",
                "frozen (2x faster)", "frozen (2x faster)",
                "frozen (2x faster)"),
  incumbent = c("sticky, stay 0.9", "sticky, stay 0.99", "sticky, stay 0.9",
                "sticky, stay 0.99", "frozen", "frozen", "frozen",
                "(none yet)"),
  everMoved = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE),
  topologyCut = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE)
)
h2h$adopted <- vapply(seq_len(nrow(h2h)), function(i) {
  cand <- sub(" \\(2x faster\\)", "", h2h$candidate[i])
  speed <- if (grepl("2x", h2h$candidate[i])) 2 else 1
  if (h2h$incumbent[i] == "(none yet)") {
    cs <- Score(cand)
    return(mkp$.BeatsIncumbent(cs$rate * speed, cs$ess, -Inf, NA,
                               candFrozen = cs$frozen,
                               topologyCut = h2h$topologyCut[i]))
  }
  H2H(cand, h2h$incumbent[i], speed, h2h$everMoved[i], h2h$topologyCut[i])
}, logical(1))
print(h2h, row.names = FALSE)
write.csv(tab, "dev/red-team/explore-195/q4-rates.csv", row.names = FALSE)
write.csv(h2h, "dev/red-team/explore-195/q4-h2h.csv", row.names = FALSE)
