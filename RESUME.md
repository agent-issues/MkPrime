# mkp — hand-off 2026-05-18

Branch: `claude/wonderful-tereshkova-08a8b5` (worktree)

## What this repo is

MkPrime — R package for Bayesian phylogenetic inference under the Mk′ model
for discrete morphological characters. The current focus is a 2026-05-12
prior-validation pilot comparing six arms on 26 simulated trees × 10 reps:

- `mk`         — M(kObs) — naive plug-in (one k per char = observed states)
- `mk_kp1`     — M(kObs+1) — fixed +1 hidden state per char
- `mk_kp2`     — M(kObs+2) — fixed +2 hidden states per char
- `mk_k9`      — M(9) — fixed k=9 across all variable chars
- `mkp_eg`     — Mk′ with empirical_geometric prior on k′ ("EG")
- `mkp_geo`    — Mk′ with plain geometric prior on k′ (pilot variant)

The headline live finding (pre-correction): `mk_kp1` beat all three Bayesian
arms (mk / EG / geo) on CID-to-truth across 260 tasks. Post-correction of an
alignment bug, EG's posterior on k′ turns out to concentrate near kObs+1
(median u_post = 1.01 for binary chars), making `mk_kp1` an oracle-ish point
estimate for the posterior median.

## Where we left off

- **c112167** — committed this session's Hamilton scripts (run_one.R updated
  with mk_kp2 + mk_k9 arms; mkp_geo thin bumped to 100 per no-oversample
  rule), new slurm scripts, dispatched-summariser fix, and
  `dev/pilots/2026-05-12-prior-validation/analysis/u_post_eg.R`.

Key files touched:
- [data-raw/hamilton/run_one.R](data-raw/hamilton/run_one.R) — mk_kp2, mk_k9, mkp_geo thin
- [data-raw/hamilton/summarize_array.slurm](data-raw/hamilton/summarize_array.slurm) — all six arms supported
- [tmp_summarize_streamed.R](tmp_summarize_streamed.R) — fixed lex-sort kObs alignment; now stores `$kObs` in each summary RDS (production copy is at `/nobackup/pjjg18/mkp-study/summarize_streamed.R` on Hamilton)
- [dev/pilots/2026-05-12-prior-validation/analysis/u_post_eg.R](dev/pilots/2026-05-12-prior-validation/analysis/u_post_eg.R) — corrected per-char u_post audit
- [dev/pilots/2026-05-12-prior-validation/analysis/cid_four_prior.R](dev/pilots/2026-05-12-prior-validation/analysis/cid_four_prior.R) — existing 4-way CID figure (unaffected by alignment bug)

## Pending jobs

None. All Hamilton arrays completed during this session:

| Arm | Run job | Outcome | Summarise job | Result |
|---|---|---|---|---|
| `mkp_geo` (rep 01 re-run, thin=100) | 17205133 | 8 COMPLETED, 18 TIMEOUT (logs intact) | 17214743 | 26/26 RDS on Hamilton |
| `mk_kp2` | 17204958 | 106 COMPLETED, 153 TIMEOUT | 17214745 | 260/260 RDS on Hamilton |
| `mk_k9` | 17204959 | 141 COMPLETED, 118 TIMEOUT | 17215030 | 260/260 RDS on Hamilton |

TIMEOUT just means SIGKILL at 8 h wall — streamed logs were complete and the
post-hoc summariser produced valid RDS. Disk usage: 242 GB / 600 GB quota.

## Open items / next steps

1. **Pull new summaries to local** (mkp_geo, mk_kp2, mk_k9 — all in `/nobackup/pjjg18/mkp-study/summary/` on Hamilton, target `dev/pilots/2026-05-12-prior-validation/summary/` locally).
2. **Redo u_post audit on mkp_geo** with `dev/pilots/2026-05-12-prior-validation/analysis/u_post_eg.R` adapted to mkp_geo prefix. This tests whether geo's posterior also concentrates near kObs+1 or differs (would explain why mk_kp1 also beats geo).
3. **6-way CID comparison**: extend `cid_four_prior.R` to include mk_kp2 and mk_k9. The natural discriminator question: does mk_k9 win (flexibility story) or does mk_kp1 still win (kObs+1 is genuinely special)?
4. **Write up corrected u_post + 6-way CID into `dev/red-team/findings.md` and `dev/pilots/2026-05-12-prior-validation/campaign-log.md`.**

## Technical pointers

- **SSH**: PowerShell ssh to `hamilton8.dur.ac.uk` works; Bash-tool ssh does NOT (different config). Always route Hamilton commands through PowerShell.
- **PowerShell quoting**: use single-quoted strings for ssh commands containing `$` to avoid PowerShell variable interpolation. PowerShell 5.1 lacks `&&` / `||`; chain with `;` only when failure is acceptable.
- **Hamilton paths**:
  - Run scripts: `/nobackup/pjjg18/mkp-study/{run_one.R, summarize_streamed.R, *_array.slurm}`
  - Raw streamed logs: `/nobackup/pjjg18/mkp-study/results/t{NN}_r{MM}/{arm}_run_1.log`
  - Tree files: `/nobackup/pjjg18/mkp-study/results/t{NN}_r{MM}/{arm}_trees.nwk`
  - Summaries: `/nobackup/pjjg18/mkp-study/summary/{arm}_t{NN}_r{MM}.rds`
  - Data root: `/nobackup/pjjg18/mkprime-files/tree-inference/tree_{NN}/rep_{MM}/chr*.nex`
- **kObs alignment**: `chr*.nex` filenames are unpadded. `run_one.R` uses plain `sort()` → lex order (chr1, chr10, chr11, …, chr2, chr20, …). The MCMC log column `kPrime_<i>` refers to position i in this lex-sorted character list. Any post-hoc analysis must either (a) load chars in lex order to match, or (b) read `$kObs` directly from the summary RDS (now stored).
- **EG arm result (corrected, 2026-05-17)**: median u_post = 1.01 (binary chars, n=9917); EG posterior on k′ is calibrated near kObs+1 for the dominant binary class. mk_kp1 collapses this to a point estimate and beats EG on CID by reducing nuisance-parameter integration variance.
- **No oversampling**: thin streamed MCMC to ~30k samples/run (thin_iters ≥ 16 for an 8 h job). Default thin=10 in `make_mcmc` is too dense; new arms explicitly pass thin_iters=100. See `feedback_no_oversample.md`.
- **Mk' study HARNESS-001**: nRuns=2L + 8h SLURM wall = TIMEOUT before saveRDS; post-hoc summariser recovers data from streamed logs. Mostly tolerable.

## Things ruled out

- **u_post anchor hypothesis (first pass, pre-bug-fix)**: thought mk_kp1 ≈ EG posterior mean. Initial audit (with broken alignment) showed u_post ≈ 0.1, rejecting the hypothesis. Post-fix, the hypothesis is **revived** — EG actually does sit near kObs+1.
- **Sampler bug in EG prior**: investigated whether negative u_post values were from samples with k′ < kObs. Verified locally that the constraint is correctly enforced (`R/MkPrimeModel.R:406` returns -Inf for k′ < kObs). The "negative u_post" was a script-level alignment bug, not a sampler bug.
- **Updating the simulation to match the new prior**: ruled out as circular (would be tuning the DGP to the inference choice).
- **Per-pattern vs per-character k′ confusion**: phyDat compresses chars to unique patterns (`nr = 37` patterns for 50 chars in t01_r01), but `mkd$kObs` and `state$kPrime` are both per-character (length nChar). MkPrime expands pattern-level structure back to char-level before MCMC.

## Suggested first action

```powershell
# Pull all new summaries to local
scp 'hamilton8.dur.ac.uk:/nobackup/pjjg18/mkp-study/summary/mk_kp2_*.rds' dev/pilots/2026-05-12-prior-validation/summary/
scp 'hamilton8.dur.ac.uk:/nobackup/pjjg18/mkp-study/summary/mk_k9_*.rds'  dev/pilots/2026-05-12-prior-validation/summary/
scp 'hamilton8.dur.ac.uk:/nobackup/pjjg18/mkp-study/summary/mkp_geo_*.rds' dev/pilots/2026-05-12-prior-validation/summary/
```

Then adapt `u_post_eg.R` to `u_post_geo.R` (one-line change to the regex/labels) and extend `cid_four_prior.R` to `cid_six_prior.R` with mk_kp2 + mk_k9 arms.
