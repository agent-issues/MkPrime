# Provenance — read before reusing anything in this directory

Audited 2026-09-18 against the open correctness issues on `agent-issues/MkPrime`.

## Produced by a sampler that was not pi-invariant (issue #19)

Driver: `data-raw/hamilton/run_one.R`, which calls `MkPrimeMCMC()` (line 164) with **no
`fixTopology`** and **no `gibbsSpr` override** — free topology with `gibbsSpr = TRUE`, the
configuration issue #19 covers. `gibbs_spr` committed a deterministic `0.5 * lReg` split
with no MH step, so the chain did not target the posterior. Fixed by PR #30.

Distorted: adjacent edge-length fractions at regraft sites (pulled toward equality) and the
topology posterior. `tree_length` second-order only. Magnitude unquantified.

## NOT affected by issues #25 / #26

Checked, not assumed: the streamed logs here carry no `rate_neo` column, so these datasets
contain no neomorphic characters, `compute_partition_scales` returns `(1, 1)`, and the
missing transformational scale is identically zero. The `mk` arms additionally pin `k` via
`knownStates`, so the case-25 Gibbs k-prime sweep never runs for them.

## Also recorded here

`HARNESS-001` (issue #13) applies to these runs: `maxTime` bounded each run separately
rather than the job as a whole, so job 17170442 lost all 260 post-run summary `.rds` files
to the SLURM wall. The streamed `.log` and `_trees.nwk` output survived and is what these
directories mostly contain.

## Before deleting

Roughly 416 MB. It is the streamed output of a multi-day Hamilton campaign and cannot be
cheaply regenerated. Confirm nothing in `data-raw/step1_hamilton_summary.rds`,
`dev/m9-pilot/`, or a manuscript draft still depends on it first.
