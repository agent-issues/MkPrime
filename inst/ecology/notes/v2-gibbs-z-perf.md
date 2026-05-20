# gibbs_z perf optimization (advisor flag)

The advisor noted during v2 design that γ_e cancels in the Gibbs
categorical ratio across the three z values for a single (c, e) cell.
This is correct — the current `gibbs_z_sweep_impl` over-computes by a
constant factor of ~3.

## Current code (correct but wasteful)

For each (character c, non-ref ecology e):

```cpp
for (int v = 0; v < 3; ++v) {
  zRow[j] = v;
  ll[v] = per_char_log_lik_ecology(  // full per-char pruning
    *data, c, parent, child, edgeLen,
    kp, rateLoss, rateNeo, rates, phi, zRow, wEdge,
    refE, gammaE);
}
```

Each `per_char_log_lik_ecology` call runs a full per-character pruning
pass (nEdge × nCat × kStates²). That's three pruning passes per cell;
with kEco-1 = 1 non-ref ecology and 240 chars × 200,000 iter,
gibbs_z runs many millions of pruning passes per chain.

## Optimization

γ_e doesn't depend on z — only on (π₀, θ_e, φ). So in the Gibbs
ratio

$$
P(z_{c,e} = v \mid \text{rest}) \propto P(\text{data} \mid z_{c,e} = v) \cdot P(z_{c,e} = v),
$$

γ_e enters both the data likelihood (via the rate factor μ(z)/γ_e) and
is the same for all three z values per cell. Pulling 1/γ_e out as a
common multiplicative factor on the per-edge rate (independent of z),
we can precompute the (p_same, p_diff) or (P00, P11, ...) matrices
ONCE per (cat, edge, e) and then sample z without re-running the
pruning three times.

But this requires rewriting `per_char_log_lik_ecology` to expose the
per-z factor independently, which is a non-trivial restructuring of
the C++.

## Simpler partial optimization

Cache the conditional likelihoods at z=0 once per (c, e), then compute
z=1 and z=2 via deltas. For neomorphic the asymmetric matrix changes
significantly between z=1 and z=2, but for transformational (JC-K) the
(p_same, p_diff) shifts predictably with phi → phi/γ_e and phi →
1/(phi·γ_e). For symmetric JC the shift is just rescaling the
effective edge time, so the same pruning buffer can be reused with
re-scaling.

## Recommended priority

LOW. The bigger Sim 3 and rodent runs are 30 min - 3 hr each; cutting
gibbs_z time by ~3× would save 20-60 min per run, but the bigger win is
fixing the MCMC mixing problem (v2 vs v1 logL gap closure) which makes
chains converge in fewer iter overall.

Revisit after v2 is verified end-to-end and the multi-rep paper run is
profiled.
