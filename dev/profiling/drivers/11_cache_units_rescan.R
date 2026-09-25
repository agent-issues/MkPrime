# #11 (LIKE-003): cache_total_loglik (src/node_cl_cache.h) re-scans
# cache.units once per partition to find each partition's units whenever an
# ascertainment correction applies (coding != "none"). Is that
# O(nParts * nUnits) scan measurable next to the work the same call does?
#
# Isolated micro-benchmark: the scan loop (partIdx compare + skip over a
# CacheUnit-sized struct) against a precomputed per-partition index, at Sun2018
# unit counts and at a pessimistic 10 distinct k' per partition. The floor on
# the rest of one cache_total_loglik call is one constant_site_prob_jc per
# distinct k (jcAscProb memoises per k); that is timed on the same tree.
#
# Usage: Rscript dev/profiling/drivers/11_cache_units_rescan.R [lib]

args <- commandArgs(TRUE)
if (length(args)) .libPaths(c(args[1], .libPaths()))
suppressPackageStartupMessages({
  library(MkPrime)
  library(TreeTools)
})

Rcpp::sourceCpp(code = '
#include <Rcpp.h>
#include <chrono>
using namespace Rcpp;

// sizeof(CacheUnit) is ~200 bytes; pad so each probe touches its own lines.
struct FakeUnit { int partIdx; int kStates; char pad[248]; };

// [[Rcpp::export]]
NumericVector BenchRescan(int nParts, int unitsPerPart, int nRep) {
  std::vector<FakeUnit> units(nParts * unitsPerPart);
  for (size_t i = 0; i < units.size(); ++i) {
    units[i].partIdx = (int)(i % nParts);
    units[i].kStates = 2 + (int)(i / nParts);
  }
  std::vector<std::vector<int>> index(nParts);
  for (int ui = 0; ui < (int)units.size(); ++ui)
    index[units[ui].partIdx].push_back(ui);

  volatile long sinkA = 0, sinkB = 0;
  auto t0 = std::chrono::steady_clock::now();
  for (int r = 0; r < nRep; ++r)
    for (int pi = 0; pi < nParts; ++pi)
      for (int ui = 0; ui < (int)units.size(); ++ui) {
        if (units[ui].partIdx != pi) continue;
        sinkA = sinkA + units[ui].kStates;
      }
  auto t1 = std::chrono::steady_clock::now();
  for (int r = 0; r < nRep; ++r)
    for (int pi = 0; pi < nParts; ++pi)
      for (int ui : index[pi]) sinkB = sinkB + units[ui].kStates;
  auto t2 = std::chrono::steady_clock::now();
  if (sinkA != sinkB) stop("index variant visited different units");
  return NumericVector::create(
    std::chrono::duration<double, std::micro>(t1 - t0).count() / nRep,
    std::chrono::duration<double, std::micro>(t2 - t1).count() / nRep);
}
')

nexFile <- system.file("datasets/Sun2018.nex", package = "TreeSearch")
if (!nzchar(nexFile)) {
  nexFile <- file.path(Sys.getenv("TREESEARCH_SRC"),
                       "inst/datasets/Sun2018.nex")
}
sun <- ReadAsPhyDat(nexFile)
nTip <- length(sun)
mkd <- suppressWarnings(MkPrimeData(sun))
nParts <- length(mkd$partitions)
cat(sprintf("Sun2018: %d partitions, kObs %s\n", nParts,
            paste(sort(unique(mkd$kObs)), collapse = ",")))

nRep <- 1e5L
for (upp in c(1L, 10L)) {
  t <- replicate(7, BenchRescan(nParts, upp, nRep))
  med <- apply(t, 1, median)
  cat(sprintf("%3d units/part: scan %.4f us, index %.4f us, delta %.4f us\n",
              upp, med[1], med[2], med[1] - med[2]))
}

tree <- Preorder(RenumberTips(NJTree(sun, edgeLengths = TRUE), names(sun)))
tree$edge.length[tree$edge.length <= 0] <- 1e-8
rates <- MkPrime:::DiscreteLognormalRates(0.5, 4L)
CspOnce <- function() {
  MkPrime:::constant_site_prob_jc(tree$edge[, 1], tree$edge[, 2],
                                  tree$edge.length, nTip, 2L,
                                  c(0.5, 0.5), rates)
}
cspUs <- median(replicate(7, system.time(for (i in 1:2000) CspOnce())[[3]])) /
  2000 * 1e6
cat(sprintf("one constant_site_prob_jc (k = 2, incl. R call): %.2f us\n",
            cspUs))
