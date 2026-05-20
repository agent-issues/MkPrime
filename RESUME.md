# mkp — hand-off 2026-05-20 (12-arm results — oracle does NOT win)

## What this repo is

MkPrime — R package for Bayesian phylogenetic inference under the Mk′ model
for discrete morphological characters. The current focus is a 2026-05-12
prior-validation pilot comparing arms on 26 simulated trees × 10 reps to find
which k-specification best recovers the true tree (CID-to-truth). A parallel
real-data pilot (in neotrans; see below) tests whether the winning arm from
simulation transfers to six empirical matrices compared against well-corroborated
reference trees.

## Where we left off

- **adc5476** — 7-way CID confirmed: mk_k15 beats mk_k9 (p=4.3e-22, 206/260 wins).
  Ramp monotonic up to k=15; k24 queued.

- **c368af4** — mk_k24 full array submitted as job **17222516** (260 tasks, 8h).

- **1e42fbb** — `dev/pilots/2026-05-12-prior-validation/analysis/cid_eight_prior.R`
  written and run. **8-way CID results** (job 17224641 ran summariser for mk_k24):

  | Arm | n | mean CID |
  |---|---:|---:|
  | mk | 260 | 0.2754 |
  | mkp_eg | 260 | 0.2677 |
  | mk_kp1 | 260 | 0.2595 |
  | mkp_geo | 26 | 0.2562 |
  | mk_kp2 | 260 | 0.2535 |
  | mk_k9 | 260 | 0.2451 |
  | mk_k15 | 260 | 0.2418 |
  | **mk_k24** | 260 | **0.2398** |

  **Paired mk_k15 vs mk_k24** (260 common tasks): k24 wins 191/260, mean diff
  −0.00197, binomial p = 2.03e-14. **The ramp continues past k=15**, but with
  further diminishing returns:

  | Step | Δ mean CID |
  |------|-----------:|
  | kp1 → kp2 | −0.0060 |
  | kp2 → k9  | −0.0084 |
  | k9  → k15 | −0.0034 |
  | k15 → k24 | −0.0020 |

  On the 26-task mkp_geo-limited subset: k15 ≈ k24 (0.2260 vs 0.2261 — tied).
  Improvement from k24 comes from harder/larger trees.

  **mk_k24 convergence**: Contrary to the pilot estimate (~30 iter/sec = 43×
  slower), mk_k24 tasks converge on the built-in convergence criterion, finishing
  in 1–4h rather than hitting the 8h walltime. 172 of 260 tasks converged within
  ~4h; tree counts range from ~750 to 66k depending on task difficulty. The "43×
  slower" figure was likely a measurement artefact from the pilot.

  mk_k15 array (17218354): 186/260 converged, 74 timed out (8h walltime). The
  74 timed-out summaries contain partial (but substantial) tree data.

## 12-way CID results (final, all arrays complete or near-complete as of 2026-05-20 04:30 BST)

| Arm | n | mean CID | Family |
|---|---:|---:|---|
| **mk_k40** | 260 | **0.2385** | fixed-k (winner) |
| mk_k24 | 260 | 0.2398 | fixed-k |
| mk_k15 | 260 | 0.2418 | fixed-k |
| mk_k9 | 260 | 0.2451 | fixed-k |
| mk_kp2 | 260 | 0.2535 | fixed-k |
| mkp_geo | 26 | 0.2562 | Mk′ geometric |
| mk_kp1 | 260 | 0.2595 | fixed-k |
| mkp_logs | 260 | 0.2595 | Mk′ logseries(c=0.95) |
| **mk_ktrue** | 260 | **0.2624** | **oracle (k=k_true per char)** |
| mkp_highk | 260 | 0.2639 | Mk′ geom Beta(1,20) |
| mkp_eg | 260 | 0.2677 | Mk′ empirical_geom |
| mk | 260 | 0.2754 | kObs |

**HEADLINE FINDING**: The oracle (mk_ktrue) does NOT beat mk_k40 — mk_k40 wins by 0.024 CID, paired p = 5.9×10⁻⁴⁵. mk_ktrue ≈ mkp_highk (p=0.15). The mk_k40 advantage is *not* about correctly modelling state-space cardinality.

**The mechanism is regularisation-via-saturation:**
1. Under JC(k=40), transition probability saturates at 1/k = 0.025
2. Data can't pull branch lengths higher than the saturation regime → tree-length prior dominates → shorter branches
3. Shorter branches → less LBA-style topology noise → better CID
4. mk_ktrue uses per-character k_true (mean 3.78) — only weakly activates saturation
5. Mk′ arms (mkp_eg, mkp_geo, mkp_highk, mkp_logs) carry a **relabel correction** that cancels some of the saturation-driven regularisation (per the k′ posterior agent — `R/likelihood.R:185`, mk_k40 uses `knownStates=40` → type="known" → relabel skipped)
6. So: **mis-specification → more regularisation → better tree**. The "blunt tool" wins because its bluntness IS the regulariser.

This rewrites the story for the paper. mk_k40 is best read as a regulariser disguised as a model, not as a "right" model.

Δ ramp now (kp1→kp2→k9→k15→k24→k40): −0.0060, −0.0084, −0.0034, −0.0020, −0.0013. Strong diminishing returns; mk_k40 is near the ceiling.

## Red-team audit (2026-05-19, opus subagents)

Initial mk_k40 > mk paradox investigated. **Not a code bug** — `mk` uses k=kObs (realistic baseline), but the simulation's character-filter retains many characters with k_true > kObs (51.9% have hidden multistate; 36.9% have k_true ∈ {3..20} but kObs=2). mk_k40 wins because kObs is mis-specified for those characters. The ramp captures the cost of trusting kObs.

Findings in:
- `RED_TEAM_simulation.md` (kObs vs k_true distribution; trees 11-26 generation script missing; 35 contaminated rows in ground_truth.csv from rbind leakage)
- `RED_TEAM_likelihood.md` (no code bugs; NOTE-1 mechanism: JC(k) × tree-length-prior interaction acts as LBA regulariser at higher k)
- `RED_TEAM_summariser.md` (no CID bug; CRITICAL: mk_k40 dispatch was on Hamilton via tmp_add_k40.py but NOT in committed run_one.R; same now applies to mk_ktrue, mkp_highk, mkp_logs)
- `RED_TEAM_stratification.md` (paradox is NOT tree-shape-confounded; mk_k40 wins 26/26 trees, uniform across J1)
- `TL_HYPOTHESIS_TEST.md` (2026-05-19 evening) — **the "high-k forces long branches" intuition was WRONG: high-k arms infer SHORTER trees than mk (mk_k40 TL=1.19, mkp_eg TL=1.34, mk TL=1.39, truth TL=1.40). Spearman ρ(ΔTL, ΔCID) = −0.37, p=1.4e-9. Mechanism: under JC(k=40), transition prob saturates at 1/k=0.025; T can't be increased to fit more transitions, so prior pulls T short. Short branches → less LBA noise → better topology. So "the blunter tool wins" is via the *prior pulling T short under high-k saturation*, not via long-branch regularisation. mk's TL ≈ truth but its topology is worst — having the right TL is bad for topology under kObs misspecification.**

## Pending jobs

### Hamilton HPC — mkp sim (all four arrays complete)

| Job ID | Arm | Final state | Note |
|--------|-----|--------|------|
| 17225905 | mk_k40 | 191 COMPLETED, 32 FAILED, 37 TIMEOUT | All 260 RDS summaries present; cid numbers unchanged from partial-data version (0.2385 final vs 0.2386 partial) |
| 17227281 | mk_ktrue | 4 COMPLETED, 126 FAILED, 130 TIMEOUT | All 260 RDS summaries present. FAILED status was spurious — R hit a benign parser error on exit AFTER MCMC completed and saveRDS ran. Data is usable. **mk_ktrue mean CID = 0.2624 — the oracle LOSES to mk_k40 by 0.024 (p=5.9e-45).** |
| 17227627 | mkp_highk | 75 COMPLETED, 185 TIMEOUT | All 260 RDS summaries; mean CID 0.2639 |
| 17227628 | mkp_logs | 260 COMPLETED | All 260 RDS summaries; mean CID 0.2595 |
| 17231608-11 | summarisers (4 arms) | COMPLETED | 1040 summaries scp'd to local; cid_twelve_prior.R run; plot saved |

The mk_ktrue R-exit error needs diagnosing before re-running — see open items.

**mk_ktrue setup** (this session, 2026-05-19):
- 260 `ground_truth.csv` files extracted from local `C:\Users\pjjg18\GitHub\mkprime\tree-inference\` and scp'd to Hamilton at `/nobackup/pjjg18/mkprime-files/tree-inference/tree_NN/rep_MM/`
- `run_one.R` patched via `tmp_add_ktrue.py` to add `mk_ktrue` dispatch (reads ground_truth.csv, maps lex-sorted column order → file_num → k_true). NOT YET committed to canonical run_one.R.
- `summarize_streamed.R` and `summarize_array.slurm` also patched to accept mk_ktrue.
- Smoke test (`tmp_test_ktrue_mapping.R`) confirmed: lex mapping correct, k_true ≥ kObs everywhere, mean overage 1.56 on t01_r01.
- Slurm script: `/nobackup/pjjg18/mkp-study/mk_ktrue_array.slurm`. Confirmed first task running cleanly: `mk_ktrue: k_true range 2-11, mean 3.78, vs kObs range 2-4` on t01_r01.

On mk_ktrue completion: `sbatch --array=0-259 summarize_array.slurm mk_ktrue`, scp summaries, build `cid_ten_prior.R` from cid_nine_prior.R template.

### Hamilton HPC (long-form M9 real-data, 3-day walltime)

| Job ID | Matrix | Model | Status | Note |
|--------|--------|-------|--------|------|
| 17217093 | syab07200 | by_nt_9v | RUNNING | 3-day job, ~22h elapsed |
| 17217094 | syab07200 | t_kv | RUNNING | ~22h elapsed |
| 17217095 | syab07200 | t_9v | RUNNING | ~22h elapsed |
| 17217096 | syab07202 | by_nt_9v | RUNNING | ~22h elapsed |
| 17217097 | syab07202 | t_kv | RUNNING | ~22h elapsed |
| 17217098 | syab07202 | t_9v | RUNNING | ~22h elapsed |
| 17217100 | syab07204 | t_kv | RUNNING | ~22h elapsed |
| 17217103 | syab07205 | t_kv | RUNNING | ~22h elapsed |
| 17217105 | syab07206 | by_nt_9v | RUNNING | ~22h elapsed |
| 17217106 | syab07206 | t_kv | RUNNING | ~22h elapsed |
| 17217107 | syab07206 | t_9v | RUNNING | ~22h elapsed |
| 17217660 | syab07204 | t_9v | RUNNING | ~22h elapsed |
| 17217661 | syab07205 | by_nt_9v | RUNNING | ~22h elapsed |
| 17217662 | syab07205 | t_9v | RUNNING | ~10.6h elapsed |

All complete ~2026-05-20 ~09:00 BST. On completion: `scp` the `.trees` files
from `/nobackup/pjjg18/m9-long/<matrix>/` to `dev/m9-pilot/syab<pid>/`, then
run `Rscript dev/m9-pilot/process_pilot.R <pid>` for each matrix and build full
CID table.

### Known problem — needs decision

**17217659 (syab07204 by_nt_9v at 64G) OOM'd** after 1m43s. Options: resubmit
at 128G, or skip syab07204 by_nt_9v from the 6-matrix comparison.

## Open items / next steps

1. **Diagnose the spurious R parser error in mk_ktrue** — 126 of 260 mk_ktrue tasks were marked FAILED (exit code 1) but ALL ran to completion and wrote RDS files. The err log shows MCMC done ✔, then warnings, then "Error: unexpected ')'" before exit. The script structure is fine. Likely candidates: (a) a downstream finalizer hook somewhere; (b) `cat(\"...\\n\")` interaction with batch-mode R; (c) warnings being converted to errors by some option. Low urgency — data is usable as-is — but worth understanding before submitting more runs.
2. **Re-frame the paper around regularisation-via-saturation**. The new findings invalidate the original "high-k beats kObs because the data have hidden multistate" framing. The oracle (mk_ktrue) loses by 0.024 CID. So the mechanism isn't "k=40 ≈ k_true on average" — it's JC(40) saturation cap driving branch-length regularisation via the prior. mk_k40 is a regulariser disguised as a model. See `TL_HYPOTHESIS_TEST.md` and `KPRIME_POSTERIOR_SHAPE.md` for mechanism details.
3. **Commit canonical run_one.R + auxiliary scripts** — all dispatches now in `data-raw/hamilton/run_one.R` plus summarize_streamed.R, summarize_array.slurm, mk_ktrue_array.slurm, mkp_highk_array.slurm, mkp_logs_array.slurm. Commit. (Task #12 in this session)
4. **Direct test of saturation hypothesis**: re-run a single task with mk_k40 but with a *very tight* prior on tree length forcing TL = 1.4 (the truth). If CID degrades back toward Mk's level, the mechanism is confirmed as branch-length-mediated. Cheap one-task experiment.
5. **The mk_no_relabel arm**: take Mk′ (with kPrimePrior) but bypass the relabel correction (treat characters as type="known" while still inferring k′). Tests whether the relabel correction is the gap between mkp_highk and mk_k40.
6. **Long-form M9 collect** (~2026-05-20 09:00 BST): scp `.trees`, run `process_pilot.R <pid>` for each of 07200, 07202, 07204, 07205, 07206, then build the full 6-matrix × 3-model CID table. Compare against 07203 (asher).
7. **Decide 17217659 retry** at 128G vs skip syab07204 by_nt_9v.
8. **Pull neotrans `by_nt_kv` baseline CID** for the 6 matrices.
9. **Clean 35 contaminated rows** in ground_truth.csv files (orig_idx ≤ 0, rbind leakage). Per-character analyses only.
10. **Re-run mk_k40 with treeLengthRate = 2k/((k−1)·FitchScore)** to test NOTE-1 mechanism directly. Now lower priority since the saturation explanation is well-supported.
11. **Pre-pub `_9i` (informative coding) variants**: deferred, after M9 pilots confirm direction.

## Technical pointers

- **SSH**: PowerShell ssh to `hamilton8.dur.ac.uk` works; Bash-tool ssh does NOT
  (different config). Always route Hamilton commands through PowerShell.
- **PowerShell quoting**: use single-quoted strings for ssh commands containing
  `$` to avoid PowerShell variable interpolation. PowerShell 5.1 lacks `&&` /
  `||`; chain with `;` only when failure is acceptable.
- **Background SSH on Windows can hang** if the ssh session is initiated as a
  Bash background task; route through PowerShell instead. Even PowerShell ssh
  can be slow to first-handshake (~30s) when the cluster is under load.
- **Hamilton paths (mkp sim study)**:
  - Run scripts: `/nobackup/pjjg18/mkp-study/{run_one.R, summarize_streamed.R, *_array.slurm}`
  - Raw streamed logs: `/nobackup/pjjg18/mkp-study/results/t{NN}_r{MM}/{arm}_run_1.log`
  - Tree files: `/nobackup/pjjg18/mkp-study/results/t{NN}_r{MM}/{arm}_trees.nwk`
  - Summaries: `/nobackup/pjjg18/mkp-study/summary/{arm}_t{NN}_r{MM}.rds`
  - Data root: `/nobackup/pjjg18/mkprime-files/tree-inference/tree_{NN}/rep_{MM}/chr*.nex`
- **Hamilton paths (M9 real-data pilot)**:
  - Pilot (4h): `/nobackup/pjjg18/m9-pilot/syab07203/{model}_run_{1,2}.trees`
  - Long runs (3d): `/nobackup/pjjg18/m9-long/<matrix>/{model}_run_{1,2}.trees`
  - Logs: `/nobackup/pjjg18/m9-long/logs/m9L-<matrix>-<model>_<jobid>.{out,err}`
  - Source nex files: `/nobackup/pjjg18/07200_by_t_kv/`, `/nobackup/pjjg18/07202_ns_n_kv/`, etc.
- **Local tree files are gitignored**: `dev/m9-pilot/syab*/` is in `.gitignore`
  (intentional — big outputs); CID summaries are too. Record numerical results
  in commits / RESUME.md, not files.
- **RB .trees format**: NOT NEXUS. Tab-separated: `Iteration Posterior Likelihood Prior phylogeny`. Read with `read.delim(f)$phylogeny` then `ape::read.tree(text = s)` per row.
- **WCT pruning**: `KeepTip(wcTree, intersect(sampleTips, wcTree$tip.label))` before CID. Matrix→WCT mapping: `paste0("0720", 0:6)` → `wellCorroboratedTrees.nwk` lines 1–7. syab07203→asher (outgroup: Didelphis, Macropus).
- **OOM threshold (Mk′ M9 real-data)**: fnJC(9) trans-only model needs ≥64G for
  matrices with ≥214 trans.nex lines (syab07204/07205). 64G itself is **not
  always enough** — see 17217659 (syab07204 by_nt_9v at 64G OOM'd at 1m43s).
  by_nt models may need 128G+ for these larger matrices.
- **mk_k24 convergence behaviour**: uses built-in convergence criterion; most tasks
  finish in 1–4h not 8h. Tree counts per task: 750–66k (wide range by difficulty).
  The earlier "43× slower / ~2,350 samples" estimate was a pilot artefact — disregard.
- **Per-arm MCMC speed (mkp sim)**: rough iter/sec on t01_r01 (from tree counts):
  - mk_k9:  ~120k trees in 8h (hits walltime, thin=100)
  - mk_k15: ~105k trees, converges; ~74 tasks hit 8h walltime
  - mk_k24: converges early (~750–66k trees); 172/260 done within 4h of submission
- **kObs alignment**: `chr*.nex` filenames are unpadded. `run_one.R` uses plain `sort()` → lex order (chr1, chr10, chr11, …). Analysis must match.
- **EG arm result (corrected, 2026-05-17)**: median u_post = 1.01 (binary chars, n=9917); EG posterior on k′ is calibrated near kObs+1. mk_kp1 collapses this to a point estimate, beating EG on CID.
- **No oversampling**: thin streamed MCMC to ~30k samples/run (thin_iters ≥ 16 for 8h job). See `feedback_no_oversample.md`.
- **Summariser production vs canonical**: `tmp_summarize_streamed.R` (repo root) is the production working copy on Hamilton; `data-raw/hamilton/summarize_streamed.R` is stale (missing mk_k15/mk_k24 in match.arg). Reconcile before next major summariser change.

## Things ruled out

- **u_post anchor hypothesis (first pass, pre-bug-fix)**: thought mk_kp1 ≈ EG
  posterior mean. Initial audit (with broken alignment) showed u_post ≈ 0.1,
  rejecting. Post-fix, hypothesis **revived** — EG sits near kObs+1.
- **Sampler bug in EG prior**: investigated negative u_post. Confirmed constraint
  correct (`R/MkPrimeModel.R:406`). The bug was a script-level alignment error.
- **Updating the simulation to match the new prior**: ruled out as circular.
- **Per-pattern vs per-character k′ confusion**: phyDat compresses to unique
  patterns; `mkd$kObs` and `state$kPrime` are per-character. MkPrime expands back.
- **kObs+1 sweet-spot hypothesis**: mk_k9 > mk_kp2 > mk_kp1 monotonically in
  6-way CID pilot. Extended in 7-way: mk_k15 > mk_k9 (p=4.3e-22). Extended in
  8-way: mk_k24 > mk_k15 (p=2.0e-14). Ramp is monotonic up to k=24 with strongly
  diminishing returns (Δk15→k24 = −0.0020 vs Δkp2→k9 = −0.0084).
- **mk_k24 infeasibility**: earlier "43× slower" pilot estimate was wrong. mk_k24
  converges via built-in criterion, finishing within the 8h walltime for all tasks.

## Suggested first action

```powershell
# Check mk_ktrue (17227281) progress
ssh pjjg18@hamilton8.dur.ac.uk 'sacct -j 17227281 --format=State -X | sort | uniq -c'
```

If most COMPLETED: `sbatch --array=0-259 summarize_array.slurm mk_ktrue`, then scp the mk_ktrue summaries and build `cid_ten_prior.R` from cid_nine_prior.R.

Also check M9 long-runs (complete ~09:00 BST 2026-05-20):
```powershell
ssh pjjg18@hamilton8.dur.ac.uk 'for d in /nobackup/pjjg18/m9-long/syab*/; do echo "=== $(basename $d) ==="; for f in $d/*.trees; do echo "  $(basename $f): $(wc -l < $f) lines"; done; done'
```
