# Red-team rotation log

One entry per round. Bottom of file holds `last_focus:` pointer (single
integer = the area number just reviewed). The skill picks
`(last_focus mod N) + 1` next.

Format per round:

```
## Round R — area #N <name> (YYYY-MM-DD)
- Files reviewed: ...
- Scenarios traced: ...
- Findings filed: <count> (refs to dev/red-team/findings.md rows)
- Trivial fixes applied: ...
- Notes for next reviewer of this area: ...
```

---

## Round 1 — area #1 ecology likelihood + prior math (2026-05-14)

Files reviewed:
- `src/mcmc_ecology.cpp` (1100+ lines)
- `src/mcmc.cpp` ecology block in `cpp_log_prior` (lines 388-433)
- `R/MkPrimeModel.R` `LogPrior` (lines 504-686)
- `inst/simulations/ecology/sim3-simulate.R`
- `vignettes/ecology-details.qmd`

Scenarios traced:
- `gamma_e` formula across simulator / R / C++.
- Per-edge, per-character mixture in pruning loops.
- Reference-ecology skip indexing.
- Boundary cases: theta ∈ {0, 1}, pi0 ∈ {0, 1}.
- wEdge marginal root-edge handling.
- Ecology rate hard-coded at 1.0 in marginals computation.

Findings filed: 5
- B1: wEdge root-edge departs from `1/K` spec, uses `0.5 * (margRoot + margChild)`.
- B2: Simulator scalar theta vs model per-ecology theta_e (kEco > 2 only).
- B3: z column-count contract — simulator has kEco cols, model has kEco - 1.
- B4: R `LogPrior` rejects theta ∈ {0, 1} strict AND has 0 * log(0) = NaN.
- B5: Ecology rate hard-coded at 1.0 in `compute_ecology_node_marginals`.

Trivial fixes applied during round: none (round 1 was research-only).
Fixes since (separate commits): B4 fixed in `b254b98` (theta boundary +
zero-count guard); B1, B2, B3, B5 still open.

Notes for next reviewer of this area:
- After fixing B1 + B5, re-trace `wEdge` on Sim 3's deep stem branch
  vs the simulator's hard-edge labels. Suspect ~50-100 nat impact.
- B2/B3 dormant at kEco = 2; revisit before any kEco > 2 simulation.

## Round 2 — area #2 MCMC proposals + MH ratios (2026-05-14)

Files reviewed:
- `src/mcmc.cpp` move dispatch (cases 30, 34-37), Bactrian kernel (line 78),
  scale_phi/scale_pi0/scale_theta accept blocks
- `R/RunMkPrime.R` `.kMoveTypes` (line 3127), `.BuildMoves` (line 2820+)
- gibbs_z `_impl` at `src/mcmc.cpp:4011-4116`
- log-sum-exp at `src/mcmc.cpp:4023-4045`

Scenarios traced:
- Jacobian correctness for log-phi and logit-pi0/theta moves.
- Boundary reflection (logit always stays in (0,1)).
- gibbs_z behaviour at theta = 1 (zero probability for z = 2).
- scale_phi sampling refE in per_ecology mode.
- `gibbsZEvery` plumbing.
- Move ID renumbering 30 → 34-37 after main merge.

Findings filed: 4
- M1: `scale_phi` perturbs `phi[refEcology]` in per_ecology mode (breaks
  identifiability). Defer until per_ecology mode is actually used.
- M2: `cpp_log_prior` sums LogNormal phi prior over all kEco entries
  (compounds M1). Same defer.
- M3: `gibbsZEvery` plumbed from R → `McmcData` but never read in
  `mcmc.cpp`. Documented behaviour doesn't exist. Move fires purely on
  weighted draw.
- S1 (suspect): `scalarTypes` omits `scale_theta` from 2% weight floor.
  Fixed in `b254b98`.

Trivial fixes applied during round: none.

Notes for next reviewer of this area:
- Adaptive scale tuning during warmup not exhaustively audited — see
  area #7.
- Warmup-end anomaly (chain visible sampling at iter ~7020 with
  `maxWarmup=1000`) explained by Tuning-phase budget; not a moves bug.
  Worth area #7 follow-up.
- M1/M2 dormant until per_ecology mode is used.

## Round 3 — area #3 state synchronisation (2026-05-14)

Files reviewed:
- All `state->logLik =` and `state->logPrior =` assignments in `src/mcmc.cpp`
- `fill_partition_cache` (line 491)
- `gibbs_kprime_sweep_impl`, `block_kprime_shift_impl`, gibbs_z accept,
  scale_phi/pi0/theta accept, slice_scalar_impl
- gibbs_spr / gibbs_subtree_swap and their _het / _full variants
- Pre-proposal drift diagnostic (line 4114)

Scenarios traced:
- Which moves write non-eco partition sums into `state->logLik` in eco mode.
- `wEdgeDirty` lifecycle.
- Per-character Gibbs z cell-by-cell drift.
- Init mismatch between R `.MkpEcologyLogLikelihood` and C++
  `cpp_log_likelihood_ecology`.

Findings filed: 4 (all critical, all fixed)
- S1: `gibbs_kprime_sweep_impl` (case 25, 36% of moves) used
  `cpp_partition_log_likelihood` (non-eco), corrupting `state->logLik`.
  **FIXED `36e9b53`.**
- S2: `block_kprime_shift_impl` (case 26) same bug in `else` branch.
  **FIXED `36e9b53`.**
- S3: Pre-proposal drift diagnostic used `cpp_log_likelihood` in eco
  mode. **FIXED `36e9b53`.**
- S4: `gibbs_spr_impl` + `gibbs_subtree_swap_impl` (cases 10/11, ~34%
  of moves) used partial-CL machinery with no eco awareness. **GATED
  `ad25928`** — these moves now no-op in eco mode.
- S5 (added round 4-ish): `eval_slice_target` + `slice_scalar_impl`
  accept block used non-eco `cpp_log_likelihood`. **FIXED `ad25928`**.

Trivial fixes applied during round: none (large enough to warrant a dedicated commit).

Notes for next reviewer of this area:
- Periodic from-scratch resync now lives in `run_mcmc_batch_cpp` main
  loop (every 20 iter). Drift > 0.5 nats logged to stderr.
- Residual stale-log at sample 1 was still ~113 nats after S1-S5; the
  cause was likely **init mismatch** rather than mid-chain accumulator
  drift — but that diagnostic lives in area #4 (state init).
- If new accept blocks are added in future, force them through
  `compute_full_loglik_at` which already branches on `data.ecologyAware`.

---

last_focus: 3
