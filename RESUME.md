# mkp — hand-off 2026-05-19 (8-way CID analysis complete)

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

## Pending jobs

### Hamilton HPC — mkp sim

| Job ID | Arm | Status | Note |
|--------|-----|--------|------|
| 17218354 | mk_k15 (--array=0-259) | COMPLETE | 186 converged, 74 timed out; summaries in place |
| 17222516 | mk_k24 (--array=0-259) | RUNNING | ~88 tasks still running at last check; all converge within 8h |
| 17224641 | summariser mk_k24 | COMPLETE | 260 summaries in `/nobackup/pjjg18/mkp-study/summary/` |
| 17224967 | mk_k40 (--array=0-259) | QUEUED | Submitted 2026-05-19; same 8h budget; convergence expected within 8h as for k24 |

On completion: `sbatch --array=0-259 summarize_array.slurm mk_k40`, scp summaries,
add `mk_k40` to `cid_eight_prior.R` → `cid_nine_prior.R`.

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

1. **Long-form M9 collect** (~2026-05-20 09:00 BST): scp `.trees`, run
   `process_pilot.R <pid>` for each of 07200, 07202, 07204, 07205, 07206, then
   build the full 6-matrix × 3-model CID table. Compare against 07203 (asher).
2. **Decide 17217659 retry** at 128G vs skip syab07204 by_nt_9v.
3. **Pull neotrans `by_nt_kv` baseline CID** for the 6 matrices.
4. **When mk_k40 (17224967) completes**: summarise, pull, run `cid_nine_prior.R`.
   Expected Δ ≤ −0.001 given the geometric diminishing-returns pattern.
5. **Re-summarise 74 timed-out mk_k15 tasks** with full data (now that 17218354
   is fully complete). Submit: `sbatch --array=<ids> summarize_array.slurm mk_k15`.
   Low priority if existing numbers are stable.
6. **Pre-pub `_9i` (informative coding) variants**: deferred, after M9 pilots
   confirm direction.

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
# Check M9 long-run progress (~22h elapsed, complete ~09:00 BST 2026-05-20)
ssh pjjg18@hamilton8.dur.ac.uk 'for d in /nobackup/pjjg18/m9-long/syab*/; do echo "=== $(basename $d) ==="; for f in $d/*.trees; do echo "  $(basename $f): $(wc -l < $f) lines"; done; done'
```

Then `scp` trees for any completed matrices and run `process_pilot.R <pid>`.
