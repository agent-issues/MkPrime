# Sim 3 v3 dev pilot — mixing diagnostic (2026-05-15)

Investigation prompted by post-fix multirep: blind P(true AC) = 0.001,
aware P(true AC) = 0.002. Question: cherry-picked hard split, or
mixing/posterior shape problem?

Diagnostics on (a) dev pilot result.rds (single rep, pre-`ad25928`
gate) and (b) **post-fix multirep reps 01 + 05** (fully gated
production code). Scripts:
- `dev/pilots/2026-05-14-sim3-v3-redesign/diagnostic-mixing.R`
- `inst/ecology/simulations/diagnostic-multirep-mixing.R`

## Headline

The aware chain **never visits the true topology** (0/2125 samples
contain the true AC bipartition), and a fresh recompute shows the
**truth tree has higher log-posterior than the chain's MAP by +16.8
nats** under the chain's own MAP nuisance parameters. This is a
**topology mixing failure**, not a posterior-shape problem. The chain
explores 1631 distinct topologies out of 2125 samples (77% unique),
so proposals are firing widely — they're just not finding the AC
neighbourhood.

## Key numbers (aware chain, idx 2115 = chain MAP)

| Quantity                                            | Value    |
| --------------------------------------------------- | -------- |
| Chain logged log_posterior at MAP                   | -4774.97 |
| Fresh recompute log_posterior at MAP state          | -4869.21 |
| **Drift (chain log − fresh recompute)**             | **+94.24 nats** |
| Truth tree (TL=13.5) + MAP nuisance, fresh recompute | -4758.18 |
| **Truth − chain MAP (fresh recompute scale)**       | **+16.79 nats** |
| Truth tree (TL=MAP TL=279) + MAP nuisance           | -4882.22 |
| MAP TL vs truth TL                                  | 279 vs 13.5 |
| Median TL across chain                              | 17.4     |
| Last-100 mean TL                                    | 20.4     |
| Unique topologies / total samples                   | 1631 / 2125 (77%) |
| Most-frequent topology share                        | 0.3%     |
| Samples containing true AC bipartition              | **0 / 2125 (0.000%)** |

## Interpretation

### Posterior shape is fine; sampling is broken on AC

When we recompute log-posterior consistently in R (no chain-side
accumulator), the truth tree beats the chain MAP by ~17 nats. The
chain logged 0% AC. Posterior weight on AC is therefore both real and
unsampled. This is **the** definition of a mixing failure.

### Topology proposals fire — they just miss

77% unique topologies and most-frequent ≤ 0.3% rules out a
stuck-chain. The chain wanders broadly across tree space without
finding the AC island. That is consistent with random-MH topology
moves (spr/nni/tbr/pspr) doing local hill-climbing in a posterior
that has a hard-to-cross valley between AB and AC under convergent
data.

### The 94-nat chain-log inflation is separate

The chain's logged log_posterior at the MAP sample is ~94 nats higher
than a fresh recompute of the same state. This is consistent with
accumulator drift / stale state.logLik between resync points. The
periodic eco-resync (every 20 iter) catches drifts > 0.5 nat with
warnings; the post-fix multirep showed zero such warnings.

The 94-nat gap at the MAP sample is likely a snapshot artefact: the
sample row records each parameter's value as of the moment it was
logged, but state.logLik reflects the live accumulator (which moved
ahead of the snapshot via other moves). This decouples
"chain-reported logPost" from "logPost(snapshotted params)" — the two
quantities are not the same thing. **The 17-nat truth-vs-MAP gap is
real and unaffected** because it uses fresh recompute on both sides.

### MAP tree has TL = 279, truth TL = 13.5 — a wandering branch-length problem

The chain MAP sample happens to land at TL=279 (20× truth). Yet the
median TL is 17.4 and last-100 mean is 20.4. The chain reaches the
right neighbourhood of TL on average but occasionally wanders to
TL=200+; the recorded MAP fell at one of those wanderings. High
branch-length acceptance rates (`branch_lengths: 82%`, `tree_length:
77%`) suggest proposal scales are too small — moves are accepted at
high rates because each step is tiny, so the chain takes many small
steps and explores TL slowly.

### Aware chain's nuisance MAP is far from truth on pi0 + rate_loss

- pi0 MAP = 0.027 (truth 0.75)
- theta MAP = 0.989 (truth 1.0) — close
- phi MAP = 3.99 (truth 4.0) — close
- rate_loss MAP = 0.768 (truth 1.0) — close
- rate_neo MAP = 1.048 (truth 1.0) — close

pi0 is the outlier. Truth pi0 = 0.75 means 75% of characters are
"none" (no ecology effect); MAP says only 2.7% are "none". The chain
has settled into a posterior mode that puts almost all characters
into the encouraged/discouraged regime. Even with the Beta(75, 25)
prior (which has mean 0.75 and ESS=100), the chain has been pulled
to 0.03 by the data. This is striking. Either:

1. The data really does prefer pi0 ≈ 0 (over-fitting to noise via
   the slab); or
2. There's a likelihood-prior tension that the moves can't resolve.

This is a separate concern from topology mixing but warrants a
follow-up.

## Acceptance-rate notes

- `gibbs_spr` shows 40.8% accept rate in this pilot. The eco gate at
  `src/mcmc.cpp:833` returns false in eco mode, but the gate landed
  in commit `ad25928` AFTER this pilot finished at 11:06. The pilot
  ran pre-gate code where these moves were silently writing non-eco
  logLik. Current production code (post-fix multirep) is correctly
  gated. **Not a new bug.**
- `block_kPrime: 0%` accept rate in both blind and aware — that move
  never accepts. Could indicate a proposal scale issue but is
  outside the topology-mixing scope.
- `branch_lengths: 82%`, `tree_length: 77%` — very high. Proposal
  scales likely too conservative; chain takes small, almost-always-
  accepted steps.

## Implications for the user's three questions

### 1. Is low P(true AC) cherry-pick, or larger problem?

It's a **mixing problem**, not a hard split or a noisy estimate. The
truth has measurably higher posterior weight (17 nats under MAP
nuisance) but is never sampled. With kEco=2, 8 taxa, and convergent
data placing weight on AB, the chain falls into AB-shaped basins and
the available topology moves (spr/nni/tbr/pspr) can't cross to AC.
Multi-rep variance reflects how often the chain happens to start /
warm into different basins, but no individual rep is doing the right
thing on AC. This is why 0.002 P(true AC) is not "8 reps × random
1/56 baseline" — it's "8 reps × can't escape AB".

### 2. Topology mixing: real, or miscalc?

**Real mixing is wide** (1631 unique topologies) **but
mis-targeted**. The accumulator-drift hypothesis is real (94 nats at
MAP) but doesn't change the diagnostic: the 17-nat truth-vs-chain
gap was computed by fresh R recompute on both sides, so it
sidesteps any accumulator issues. The chain is exploring tree space,
just not in the right direction.

### 3. Red-team

Area #4 completed in parallel: 15 new findings (4 HIGH/MED + 11 lower)
across R-side init and C++ init. Smoking gun for the 113-nat
sample-1 gap was found: `.InitState` doesn't rebuild
`tree$edge.length` after `initOverrides`, so R-side likelihood reads
the old NJ-tree lengths while state stores the override.
`fill_partition_cache` early-returns in eco mode at
`src/mcmc.cpp:497` with no resync, so the mismatch goes silent at
init. See `dev/red-team/findings.md` rows R4-1, R4-7, R4C-1.

## Multirep reps 01 + 05 (post-fix, fully gated)

Pulled aware-result.rds, blind-result.rds, chain logs from Hamilton
for the rep with the lowest aware P(wrong) (rep 01: 0.789) and one
with high aware P(wrong) (rep 05: 0.462) of the 8-rep multirep.
**Two distinct failure modes**:

### Rep 01 — branch-length collapse

| Quantity                                       | Value    |
| ---------------------------------------------- | -------- |
| Aware P(true AC) / P(wrong AB) / CID            | 0.011 / 0.789 / 0.446 |
| Blind P(true AC) / P(wrong AB) / CID            | 0.000 / 0.854 / 0.519 |
| Chain MAP logPost / TL                          | -4819.4 / **2.78** (truth 13.5) |
| Recompute at MAP state                          | -4942.3 |
| Truth tree + MAP nuisance, recompute            | -4972.1 |
| **Delta (truth − MAP recompute)**               | **−29.8 nats** |
| Samples containing true AC                      | 2 / 240 (0.83%) |
| Median chain TL / last-100 mean                 | 3.12 / 3.04 (truth 13.5) |
| Unique topologies                               | 215 / 240 (89.6%) |
| `gibbs_spr` / `gibbs_subtree_swap` accept       | 0.000 / 0.000 (gates fire) |

The chain collapsed branch lengths to ~3 (vs truth 13.5) and found a
short-tree topology with higher logPost than the truth tree at MAP
nuisance. Under collapsed branches, characters become saturated and
the likelihood surface flattens, so multiple wrong topologies look
equally plausible. **Not a topology-mixing failure** — a
branch-length collapse + topology-stuck combination.

### Rep 05 — pure topology mixing failure

| Quantity                                       | Value    |
| ---------------------------------------------- | -------- |
| Aware P(true AC) / P(wrong AB) / CID            | 0.000 / 0.462 / 0.473 |
| Blind P(true AC) / P(wrong AB) / CID            | 0.000 / 0.000 / 0.475 |
| Chain MAP logPost / TL                          | -4770.8 / 15.3 (truth 13.5) |
| Recompute at MAP state                          | -4901.9 |
| Truth tree + MAP nuisance, recompute            | -4795.4 |
| **Delta (truth − MAP recompute)**               | **+106.5 nats** |
| Samples containing true AC                      | 0 / 243 (0.000%) |
| Median chain TL / last-100 mean                 | 18.7 / 20.7 (truth 13.5) |
| Unique topologies                               | 159 / 243 (65.4%) |

Truth tree has logPost 106 nats higher than MAP under MAP nuisance —
a huge undersampling gap. TL is correct (15 vs 13.5). Topology
diversity 65% — chain mixes but doesn't find AC. **Pure topology
mixing failure across a deep posterior valley.**

### Cross-rep observations

- **Gates working in production**: `gibbs_spr: 0.000`,
  `gibbs_subtree_swap: 0.000` in both multirep reps. The dev-pilot
  40% rate was pre-`ad25928`; current code correctly suppresses
  these moves in eco mode.
- **Chain-log drift in both reps**: −123 nats (rep 01), −131 nats
  (rep 05). Chain-logged log_posterior is consistently higher than
  fresh recompute. Snapshot-vs-live decoupling, not catastrophic
  drift (periodic resync logged zero >0.5-nat warnings). Does
  not affect the truth-vs-MAP recompute comparison.
- **pi0 consistently low**: 0.44 (rep 01), 0.38 (rep 05). Truth
  0.75. Beta(75, 25) prior with ESS=100 is being beaten by data
  toward lower pi0. Cross-rep pattern.
- **`block_kPrime` accept = 0 in all reps**: that move never fires.
- **`scale_theta: 0.92` accept rate in both** — proposal too small
  (over-acceptance). Symmetric with `branch_lengths: 55%`,
  `tree_length: 39-57%` (rep-dependent).

### Synthesis: two failure modes

The multirep low P(true) ≈ 0 emerges from at least two distinct
chain pathologies:

1. **Branch-length collapse** (rep 01-type): chain finds a
   high-likelihood region at short TL where saturation washes out
   topology signal. Posterior at that point genuinely beats truth.
   Fix needs to be on the **branch-length sampler / TL prior** so
   the chain doesn't collapse.

2. **Topology mixing failure** (rep 05-type): TL is fine,
   nuisance is fine, but the chain can't cross to the AC island
   that has logPost +106 nats over the current basin. Fix needs to
   be on **topology proposals** with wider span (large-window TBR,
   ecology-aware SPR, simulated tempering on topology).

Both modes have **pi0 ≈ 0.4** (truth 0.75) as a common feature.
Possibly the pi0 sampler over-shrinks under the slab structure when
data signal is moderate — worth a separate audit.

## Next high-leverage moves

1. **Topology mixing proposals** (rep 05-type fix): wider-window TBR
   or simulated tempering on topology specifically. Current
   nni/spr/tbr/pspr at small scales do local hill-climbing only and
   cannot cross 100+ nat valleys.

2. **Branch-length / TL sampling investigation** (rep 01-type fix):
   diagnose why TL collapses to ~3 in some reps despite the TL
   prior favouring mean ≈ expSteps/2. Possible candidates: TL prior
   scale, joint TL-rate proposal, branch_lengths over-acceptance
   driving collapse.

3. **Fix the init bugs identified by red-team R4**: R4-1 (rebuild
   edge.length after initOverrides) and R4C-1 (no fill_partition_cache
   resync in eco mode). Small, safe fixes that close the 113-nat
   sample-1 gap and harden future initOverrides experiments.

4. **Investigate pi0 collapse**: 0.027 / 0.44 / 0.38 across pilot +
   reps 01, 05 vs truth 0.75 under Beta(75, 25) prior is a
   reproducible pattern. Either data legitimately prefers low pi0
   or the pi0/z conditional sampler has a bias. Worth a focused
   audit.

5. **Don't scale to 20 reps** until at least one of (1) or (2)
   lands. More reps of the current chain pathology won't change the
   P(true) ≈ 0 story.
