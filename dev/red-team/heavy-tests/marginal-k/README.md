# marginal-k heavy tests

Validation drivers for the `likelihoodMode = "marginal_k"` feature
(geometric arm, v1). See:

- `dev/notes/2026-05-28-marginal-k-plan.md` — implementation plan
- `dev/red-team/proofs/marginal-k-geometric.md` — proof note

## Tests

### `T-OVL-sampled-vs-marginal.R` — plan §7.2 (Rung 2)

Posterior overlap between `sampled_k` and `marginal_k` modes on the same
dataset / seed across a small `(n_tip × n_char)` grid. KS test per
parameter per cell.

**Pass bar.** KS p > 0.01 on every parameter (`tree_length`,
`rate_log_sd`, `p`) in every grid cell.

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
N_TIP = 8, N_CHAR = 30. MARGINAL between 0.01 and 0.4. FAIL below.

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
sbatch dev/red-team/heavy-tests/submit-marginal-k-sbc.sh
```

Walltime: 8h initially (per `feedback_resumable_runs`).
Results land in
`dev/red-team/heavy-tests/marginal-k/sbc-results-hamilton/`.

## Output layout

```
marginal-k/
├── T-OVL-sampled-vs-marginal.R       # §7.2 driver (PR-B)
├── T-SBC-marginal-geometric.R        # §7.3 driver (PR-C)
├── T-OVL-results.rds                 # §7.2 results
├── T-OVL-verdict.txt                 # §7.2 verdict
├── sbc-results/                      # §7.3 local artefacts
│   ├── sims.rds
│   ├── ranks.rds
│   ├── rank-histograms.png
│   ├── verdict.txt
│   └── verdict-headline.txt          # one-line PASS / MARGINAL / FAIL
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

- **T-SBC FAIL on `tree_length` or `rate_log_sd`.** Either a cache
  invalidation bug under tree/rate moves, or a parametric mismatch
  between the forward simulator and the inference model
  (see SBC-HARNESS-001..006 in `dev/red-team/findings.md` for the
  pattern — the historical SBC harness has had four such mismatches).
  First check the forward sim uses the same JC convention as the
  inference kernel (rate per substitution, not per total tree length).
