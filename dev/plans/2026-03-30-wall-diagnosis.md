# Diagnosing the topology trap / wall pattern

## Problem statement

In production runs on the hyoliths dataset (54 taxa, 225 chars, `nChains = 4`,
`nRuns = 2`), the cold chain gets stuck around sample 1000. Multiple parameters
(`rate_loss`, `tree_length`, `rate_log_sd`) collapse from proper mixing
distributions to point-mass-like spikes. The stuck regime has ~15 log-posterior
units *better* than the mixing regime, suggesting the chain found a
high-probability basin it cannot escape.

The single-run diagnostic (`nRuns = 1`, `nChains = 4`, 80k iterations) did NOT
show this pattern. This suggests the phenomenon may be seed-dependent, or may
require longer runs to manifest.

## Hypotheses

| # | Hypothesis | Prediction |
|---|-----------|------------|
| H1 | **Tempering swap imports a trap topology.** A heated chain finds a high-posterior tree; swap delivers it to the cold chain; cold-chain topology proposals (SPR/NNI) can't escape the basin. | Transition coincides with a swap event. Topology distance shows a discrete jump at the transition point. |
| H2 | **Rate–topology ridge.** The posterior has a sharp ridge coupling `rate_loss` to topology. Once on the ridge, continuous proposals can't move because likelihood is steep perpendicular to it. | Pairwise scatter of `rate_loss` vs `tree_length` shows an L-shaped or ridge pattern. Topology changes smoothly but parameters freeze. |
| H3 | **kPrime reorganisation.** A collective shift in many `kPrime_i` values creates a new parameter manifold that constrains all scalar parameters. | kPrime distribution (median, spread) differs markedly before vs. after the transition. |
| H4 | **Seed-dependent tuning.** Auto-tuner adapts to the early regime; step sizes become too small to escape once the chain enters the late regime. | The wall only appears with certain seeds. Disabling auto-tuning or resetting tuning at the transition prevents the wall. |

These are not mutually exclusive — H1 could trigger H2 or H3.

## Diagnostic plan

### Phase A: Post-hoc analysis of existing data (no code changes)

Use the production run log + tree file that's currently accumulating
(`hyoliths_prod_1.log`, `hyoliths_prod_trees.nwk`).

**A1. Topology distance time series**

Read the tree file, compute Robinson-Foulds distance between consecutive
samples. Plot as a time series. If H1 is correct, there will be a
discrete jump (large RF distance) at the transition point (~sample 1000).

```r
trees <- ape::read.tree("hyoliths_prod_trees_1.nwk")
rf <- sapply(seq_len(length(trees) - 1), function(i)
  TreeDist::RobinsonFoulds(trees[[i]], trees[[i + 1]]))
plot(rf, type = "l", main = "Consecutive RF distance")
```

Also compute RF distance of each sample to the *first* sample (baseline
topology) and to the *last* sample (trapped topology). If there's a
discrete mode switch, the distance-to-first will jump up and
distance-to-last will drop.

**A2. kPrime shift analysis**

Compare per-character kPrime values before vs. after transition.

```r
kp_cols <- grep("^kPrime_", colnames(r1))
kp_early <- r1[1:800, kp_cols]
kp_late  <- r1[1000:1500, kp_cols]
# Per-character: how many changed their modal value?
mode_early <- apply(kp_early, 2, function(x) as.integer(names(which.max(table(x)))))
mode_late  <- apply(kp_late, 2, function(x) as.integer(names(which.max(table(x)))))
sum(mode_early != mode_late)  # How many characters shifted?
```

**A3. Parameter correlation structure**

Scatter `rate_loss` vs `tree_length`, colored by sample index, to look
for the ridge structure (H2).

### Phase B: Instrumentation (code changes, then re-run)

Add lightweight diagnostic logging to the MCMC engine so we can
distinguish H1 from H2/H3/H4 definitively.

**B1. Swap event log**

Add a per-sample column `swap_cold` to the log file: a binary flag
indicating whether the cold chain participated in a swap that was
*accepted* since the last saved sample. This is cheap (one integer
counter, reset at each sample point) and tells us whether the
transition coincides with a swap.

Implementation: In `run_mcmc_batch_cpp()`:
- Add a counter `coldSwapsSinceSample` initialised to 0.
- After the swap block: if iPair == 0 and swap was accepted,
  increment `coldSwapsSinceSample`.
- At sample extraction: write `coldSwapsSinceSample` to the sample
  row, then reset to 0.

R side: add `"swap_cold"` to `paramNames`. No change to the log format
beyond an extra column.

**B2. Cold-chain log-likelihood per iteration (optional, heavier)**

Track `states[0]->logLik` at every iteration (not just sample points)
into a separate lightweight vector. Return it alongside the batch
result. This gives a high-resolution trace of when the cold chain's
likelihood jumped — a swap signature.

This is heavier (one double per iteration vs. per thinned sample) but
very informative. Store in a raw numeric vector and return to R. Only
enable behind a `diagnosticMode = TRUE` flag to avoid overhead in
production.

**B3. Topology hash per sample**

Compute a cheap hash of the cold chain's topology at each sample point.
Use the sorted edge list as the hash key. Store as an integer column
`topo_hash` in the log. This lets us detect topology transitions without
needing the tree file.

Implementation: hash the parent/child vectors. A simple approach:
sum of `parent[i] * 1000003 + child[i]` modulo a large prime. Not
collision-free but sufficient for detecting *changes*.

### Phase C: Controlled experiments

After instrumentation, run targeted experiments on the hyoliths dataset.

**C1. Tempering on vs. off**

| Config | nChains | nRuns | nIter | Purpose |
|--------|---------|-------|-------|---------|
| A | 4 | 1 | 80k | Tempering ON (reproduce diagnostic run) |
| B | 1 | 1 | 80k | Tempering OFF (no swaps at all) |

Same seed. If the wall appears only with tempering, H1 is strongly
supported. If it appears without tempering too, the problem is in the
posterior geometry (H2/H3).

**C2. Multiple seeds with tempering**

Run config A with 4 different seeds. Check whether the wall appears in
all, some, or none. Seed-dependence supports H4 (tuning) or H1 (rare
swap event).

**C3. Heat parameter sensitivity**

Try `heat = 0.05` (colder hot chains, higher swap rates) and
`heat = 0.5` (hotter, lower swap rates). If H1 is correct, higher
swap rates should make the wall appear sooner/more often.

## Implementation order

1. **Phase A** first (no code changes, immediate results from existing
   data). This narrows the hypotheses.
2. **Phase B1 + B3** (lightweight, add to C++ and R). B2 only if B1
   doesn't settle the question.
3. **Phase C** experiments after instrumentation is in place.

## Deliverables

After the diagnostic runs:

- A clear determination of which hypothesis (or combination) explains
  the wall
- Concrete next steps: e.g., if H1, improve topology proposals or
  tempering; if H2, add parameter-topology joint proposals; if H3,
  investigate kPrime proposal mechanism
