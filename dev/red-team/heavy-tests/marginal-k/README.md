# marginal-k heavy tests

Validation drivers for the `likelihoodMode = "marginal_k"` feature
(geometric arm, v1). See:

- `dev/notes/2026-05-28-marginal-k-plan.md` — implementation plan
- `dev/red-team/proofs/marginal-k-geometric.md` — proof note

## Prior variant

Every driver pins `priorVariant` rather than inheriting `MkPrimeModel()`'s
default, which for `geometric` changed from `"conditional"` (Model B) to
`"unconditional"` (Model A) on 2026-06-04. `T-OVL-sampled-vs-marginal.R` and
`gibbs-p-fullchain-check.R` default to `"unconditional"` (what ships) and take
`MARGINAL_K_PRIOR_VARIANT=conditional` to override; the value is written into
their outputs. The two `T-SBC-*-geometric.R` drivers are fixed at
`"unconditional"` because their forward model is Model A.

**The recorded 2026-06 verdicts were Model B:** the T-OVL gated and moves-on
runs and the 200k p-recheck (2026-06-02; `run-phase2-overnight.sh` now pins
Model B to reproduce them) and the "3.2x p-ESS" full-chain check (2026-06-03).
They say nothing about the current default.

## Tests

### `T-OVL-sampled-vs-marginal.R` — plan §7.2 (Rung 2)

Posterior overlap between `sampled_k` and `marginal_k` modes on the same
dataset / seed across a small `(n_tip × n_char)` grid. Equivalence test per
parameter per cell: the difference in mean and in sd must each lie within
`Z_BAR = 3` batch-means standard errors.

**Pass bar.** A cell/parameter below the ESS floor (200 in each mode) is
INCONCLUSIVE. The run PASSes only if no cell/parameter FAILs and at least half
are conclusive; otherwise it is INCONCLUSIVE.

**Quick local run** (no Hamilton submit yet — heavy test):

```bash
Rscript dev/red-team/heavy-tests/marginal-k/T-OVL-sampled-vs-marginal.R
```

### `T-SBC-marginal-geometric.R` — plan §7.3 (Rung 3)

Honest SBC under the geometric Model A forward
(`u_i ~ Geo(p); kTrue = 2 + u`), inference with `likelihoodMode =
"marginal_k"`. Rank tests on continuous parameters only — `kPrime_*` is
no longer a sampled state under marginal-k (proof §2).

**Pass bar (full mode).** Anderson–Darling p > 0.4 on each of
`tree_length`, `rate_log_sd`, `p` at N_SIM = 200, N_ITER = 12000,
N_TIP = 16, N_CHAR = 100. MARGINAL between 0.01 and 0.4. FAIL below.

Both SBC drivers run at `nCat = 1`, so the rate multiplier is 1 for every
character and `rate_log_sd` never enters the likelihood: its rank test checks
only that the chain recovers its prior. It cannot reveal rate-cache bugs.

**Quick local sanity** (env var; ~few minutes):

```bash
MARGINAL_K_SBC_QUICK=1 Rscript dev/red-team/heavy-tests/marginal-k/T-SBC-marginal-geometric.R
```

Quick mode runs N_SIM = 20, N_ITER = 1000 — confirms the harness runs
end-to-end but is **structurally underpowered for AD**. Do not interpret
quick-mode AD p-values as PASS / FAIL.

**Hamilton full run.** See `dev/red-team/heavy-tests/submit-marginal-k-sbc.sh`.
Before sbatch, **pre-build on the login node** per the
`feedback_pkgload_prebuild` memory file:

```bash
cd /nobackup/${USER}/mkp-study/red-team/mkp-source
module load r/4.5.1 && module load gcc/14.2 || true
R_LIBS_USER=/nobackup/${USER}/mkp-study/red-team/lib \
  Rscript -e 'devtools::load_all(".")'
A1=$(sbatch --parsable --export=ALL,BATCH=1 dev/red-team/heavy-tests/submit-marginal-k-sbc.sh)
sbatch --dependency=afterany:${A1} --export=ALL,BATCH=1 dev/red-team/heavy-tests/submit-marginal-k-agg.sh
```

`BATCH=1` (seedBase 20260528) writes to `sbc-results/`, `BATCH=2` (seedBase
20260901) to `sbc-results-b2/`; the driver reads `MARGINAL_K_SBC_OUTDIR` and
`MARGINAL_K_SBC_SEEDBASE`, so concurrent batches never share a shard dir. The
aggregate step mirrors each into `<dir>-hamilton/`. No committed script pools
the two batches for the marginal driver (see `MARGINAL-K-CACHE-002-resume.md`);
`pool-sampled-batches.R` does so for the sampled driver.

## Output layout

```
marginal-k/
├── T-OVL-sampled-vs-marginal.R       # §7.2 driver (PR-B)
├── T-SBC-marginal-geometric.R        # §7.3 driver (PR-C)
├── T-OVL-<tag>-results.rds           # §7.2 results; <tag> = "gated" or the
├── T-OVL-<tag>-verdict.txt           #   MARGINAL_K_OVL_EXTRA moves joined by "+"
├── sbc-results/                      # §7.3 local artefacts
│   ├── sims.rds
│   ├── ranks.rds
│   ├── rank-histograms.png
│   ├── verdict.txt
│   └── verdict-headline.txt          # one-line PASS / MARGINAL / FAIL
├── sbc-results-b2/                   # second seed batch (BATCH=2)
└── sbc-results-hamilton/             # mirror after Hamilton submit
    └── (same shape)
```

## What a failure here would mean

- **T-OVL FAIL.** Marginal evaluator does not produce the same posterior
  as the joint chain on `(tree_length, rate_log_sd, p)`. Most likely
  cause: an arithmetic mismatch in `cpp_log_likelihood_marginal` — the
  placement of the ascertainment / relabelling correction relative to
  the inner `logSumExp` (proof §3). Next: stale per-(char, k) cache
  under `p`-moves (proof §5; the "cache invariance" caveat in proof
  §9 Caveat 2 flags this is implementation-untested by §7.4).

- **T-SBC FAIL on `p`.** The marginal evaluator's pmf weighting is
  off. Compare to the sampled-k Model A SBC on `p` (T6 in
  kprime-viability) to localise: if sampled-k passes and marginal-k
  fails, the bug is in the marginal evaluator alone (not the prior).

- **T-SBC FAIL on `tree_length`.** Either a cache invalidation bug under
  tree moves, or a parametric mismatch
  between the forward simulator and the inference model
  (see SBC-HARNESS-001..006 in `dev/red-team/findings.md` for the
  pattern — the historical SBC harness has had four such mismatches).
  First check the forward sim uses the same JC convention as the
  inference kernel (rate per substitution, not per total tree length).
