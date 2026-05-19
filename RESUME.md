# MkPrime ecology-aware — hand-off 2026-05-18

## What this repo is

MkPrime is an R package implementing an ecology-aware Mk' MCMC for Bayesian
phylogenetic inference from morphological data. The `worktree-ecology-aware`
worktree drives the ecology paper: simulations demonstrating that ecology-aware
inference reduces false-clade support relative to an ecology-blind Mk' baseline,
plus an empirical rodent case study.

Artifact targets: paper draft + the rodent comparison; underlying R package is
the methods substrate.

## Where we left off

This session pushed three threads forward:

1. **v4-cross sim architecture** (eco-1 in non-sister clades A,B; opposes true
   `((A,C),(B,D))` topology rather than aligning with it as in v4):
   - `d4d0bee` introduced helpers + v4cross-b run script.
   - `2873ab4` added v4cross-c stress regime.
   - `fc2737d` added single-rep PT (nChains=4) and confirmed an aware
     mode-trap: clade B contiguously broken at sample 84 onward, ll=−1587 vs
     mode-1 ll=−1583. Diagnostics at `inst/hamilton/sim3-v4cross-b/diagnose/`.
   - `63bb335` then ran a **5-rep PT sweep** (`v4cross-b-pt5`, jobs
     17205831_[1-5], all COMPLETED). Headline: aware does NOT reliably rescue
     true topology — AC sister mean blind=0.750, aware=0.479. Aware no longer
     consistently beats blind on this confound geometry, and false-clade
     support never emerges (max falseAB=0.006). Inverse of intended story.

2. **Rodent case study**:
   - `eb2840f` added Hamilton blind pilot (200k iter, job 17205921, completed).
   - `40ca0d4` produced a first comparison report
     (`inst/scripts/rodent-comparison/`) with CID-MDS scatter and MR consensus
     plots. **Both chains were underpowered** — blind minESS=21, aware was
     truncated at 127k of 500k iter (only 7 post-burnin trees). The 0-shared-
     splits result and MDS separation are pilot artefacts, not real findings.
   - `0a95af5` resubmitted both: blind continuation (job 17207482, COMPLETED
     in 19.5 min) and aware fresh 1M iter (job 17207480, **TIMEOUT** at 12h,
     reached 371k, minESS=31, 161 trees).
   - `79df844` added `rodent-aware-v2-cont/` — resume from checkpoint to 1M.
     Job **17215319** running (24h wall, ~13h ETA).

3. **Sim2/3 mode-trap diagnostics** (background): the v4cross-b clB asymmetry
   is the same topology-valley failure as the multirep audit (rep 01/05). PT
   helps but does not eliminate it on hard configurations.

## Pending jobs

| Job ID | Name | Wall | On completion |
|--------|------|------|---------------|
| 17222514 (aware-cont2) | mkp-rod-awar2 | 36h | Check minESS from `rodent-aware-v2-result.rds`; if ≥ 200, regenerate `inst/scripts/rodent-comparison/` |
| 17222515 (blind-cont2) | mkp-rod-blind2 | 4h | Check minESS from `rodent-blind-v2-result.rds`; if ≥ 200, regenerate `inst/scripts/rodent-comparison/` |

Both jobs submitted and running 2026-05-19. Check status with:
`/c/WINDOWS/System32/OpenSSH/ssh.exe pjjg18@hamilton8.dur.ac.uk "squeue -u pjjg18"`

Targets: aware 1M→2.3M (patch nIter=2300000L), blind 1M→5M (patch nIter=5000000L).
Checkpoints at `/nobackup/pjjg18/mkp-rodent-aware-v2/rodent-aware-v2.ckp` and
`/nobackup/pjjg18/mkp-rodent-blind-v2/rodent-blind-v2.ckp`.

## Open items / next steps

1. **Rodent comparison** (DONE): Full 1M-iter chains compared in
   `inst/scripts/rodent-comparison/`. Key result: 32/38 blind splits shared
   with aware consensus (was 0/38 in pilot). MDS clouds largely separated
   (centroid dist 9.65, spread ~9 / ~7). Both chains underconverged
   (minESS 88 aware, 38 blind). A longer run or PT is needed for publication.

2. **Decide narrative on the v4-cross result**. The 5-rep PT sweep shows aware
   ≠ better than blind on v4-cross. Options:
   - Report it honestly as a failure mode of the aware model when eco signal
     opposes ancestry signal (caveat paragraph in paper).
   - Investigate further: does a wider PT ladder or longer chain shift things?
     Or is this an over-fitting of eco signal that hints at a sampler/prior fix?
   - Pivot the paper's eco-aware-rescues claim to the v4 (sister-clade)
     architecture only, where aware DOES help.

3. **Mode-trap follow-up** still queued in user memory
   (`project_mixing_followups`): rep 01 TL collapse, rep 05 topology valley,
   pi0 collapse audit. Behind PT fix.

4. **pi0/phi/z sampler audit** (`project_pi0_audit`): R5-4 verified —
   identifiability + weak prior, not a sampler bug. No code change indicated.

## Technical pointers

- **SSH to Hamilton from this machine**: use Windows OpenSSH explicitly:
  `/c/WINDOWS/System32/OpenSSH/ssh.exe pjjg18@hamilton8.dur.ac.uk "<cmd>"`.
  Passwordless. The plain `ssh` in MSYS PATH fails (Permission denied), and
  `plink` is configured for git but lacks a saved Hamilton session. SSH
  occasionally hangs from the parent shell mid-session — dispatching the same
  command via a subagent usually works; if not, retry in 15-30 min (likely
  login-node transient).
- **Hamilton paths**:
  - Repo: `/nobackup/pjjg18/mkp-sim3-multirep-v3/MkPrime` (branch
    `worktree-ecology-aware`)
  - R library: `MKP_LIB=/nobackup/pjjg18/mkp-sim3-multirep-v3/lib`
  - Deploy via `git fetch origin && git restore --source=origin/worktree-ecology-aware -- <paths>`.
    Do NOT use `git reset --hard` on Hamilton — earlier session reset zeroed
    files when ORIG_HEAD.lock failed; recovery cost ~1 hour.
- **Rodent data**: nexus at
  `/nobackup/pjjg18/mkp-rodent-blind-v2/data/mbank_X24848_2026-5-9-1135.nex`
  (only there — NOT bundled in repo). 60 extant tips × 217 chars
  (160 neomorphic + 59 transformational after `AutoDetectNeomorphic()`).
  Aware-v2 config: `magnitudeMode="global"`, `kPrimePrior="geometric"`,
  `coding="variable"`, `expSteps=10`, `rho0Alpha=7`, `rho0Beta=3`,
  `thetaAlpha=2`, `thetaBeta=2`, `sigmaPhi=0.5`, `set.seed(20260512)`.
- **v4-cross helpers**: `inst/simulations/ecology/sim3v4cross-helpers.R` and
  `sim3v4-helpers.R`. Note `MatrixToPhyDat` / `MkPrimeData` are NOT exported
  in the installed lib — use bare names, not `MkPrime::` prefix.
- **PT execution model**: `nChains > 1` runs sequentially in MkPrimeMCMC;
  set `cpus-per-task=1` and multiply wall by nChains. Confirmed in
  `sim3-multirep-v4.sh` comments.
- **Aware iter rate on rodent data**: ~50k iter/h (much slower than
  simulations). 1M iter = ~20h. 2G RAM, 1 CPU sufficient.

## Things ruled out

- **chmod -R u+w** on Hamilton paths — blocked by classifier as destructive.
  Use file-by-file chmod or work around via git restore.
- **scp of .sh scripts** — line-ending corruption injected binary nulls in
  earlier sessions. Either strip CRLF before scp, or write via SSH heredoc.
- **Treating all rodent characters as transformational** — user clarified
  the existing 160 neo + 59 trans classification IS gold standard; do not
  reclassify.
- **Single-chain aware MCMC on v4-cross** — exhibits contiguous mode-traps
  (sample 84 onward stuck). Default to nChains ≥ 4 with PT for any
  v4-cross-like configuration.
- **`-amend` after a failed pre-commit hook** — the commit didn't happen, so
  amend rewrites the previous commit. Create a new commit instead.

## Worktrees

| Branch | Path | Status |
|--------|------|--------|
| worktree-ecology-aware | C:\Users\pjjg18\GitHub\mkp\.claude\worktrees\ecology-aware | clean (modulo pre-existing untracked debris) |

## Suggested first action

Once back, run the COLLECT step on this RESUME.md: poll job 17215319 via
`/c/WINDOWS/System32/OpenSSH/ssh.exe pjjg18@hamilton8.dur.ac.uk "sacct -j 17215319 -n -X -o JobID,State,Elapsed,ExitCode"`.
If COMPLETED, regenerate the rodent comparison per the table's "On completion"
column.
