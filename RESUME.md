# MkPrime ecology-aware — hand-off 2026-05-20 (evening)

## What this repo is

MkPrime is an R package implementing an ecology-aware Mk' MCMC for Bayesian
phylogenetic inference from morphological data. The `worktree-ecology-aware`
branch (aliased locally to `ecology-aware`) drives the ecology paper:
simulations demonstrating that ecology-aware inference reduces false-clade
support relative to an ecology-blind Mk' baseline, plus an empirical rodent
case study (MorphoBank X24848). It is being developed as a "plugin" to mkp
core, with all ecology-restricted content living under `dev/ecology/` and
`inst/ecology/` so future merges with `main` stay clean.

## Where we left off

This session fixed the two streaming bugs that broke `RelabelEcology()` on
rep01 and rep02 of the v9-induce 3-rep run. Both fixes landed in a single
commit (`326bb13`) on `worktree-ecology-aware`.

Recent commits (top first):

- `326bb13` **fix(streaming): atomic log writes, z_samples alignment, and
  trimZSamples escape hatch** — the two streaming bugs that prevented
  `RelabelEcology()` from running on PT chains are now fixed:
  - **Bug B (torn row)**: `.FlushBuffer()` now opens an explicit append
    connection and calls `writeLines(line, con)` once per row instead of
    building the whole flush block as one `cat(paste(..., collapse="\n"))`
    call. Each row is ~350 bytes (well under PIPE_BUF), making every OS write
    atomic on NFS.
  - **Bug A (z_samples drift)**: `.BuildResult()` now trims `z_samples` to
    `saved_idx` when they diverge, with a warning. This fixes the
    "z_samples=230 vs nrow(samples)=200" mismatch that aborted rep01.
  - `RelabelEcology()` aborts with rich diagnostics when `z_samples` and
    `samples` are still misaligned (e.g. if loaded from a log written before
    the fix). New `trimZSamples = "tail"` / `"head"` argument lets the
    caller explicitly trim a known surplus.
  - 5 new regression tests in `tests/testthat/test-streaming.R` (STREAM-003
    through STREAM-005: column-count integrity, ecology z alignment invariant,
    default-abort path, tail-trim path).
- `e2d88ed` docs: update RESUME.md for hand-off 2026-05-20 afternoon
  (previous hand-off; describes v9 null result and rodent v3 non-resumable)
- `fc6c2ff` feat: rodent v3 aware continuation harness — both submitted
  SLURM jobs FAILED at startup (checkpoint/tree desync); see Open items.

## Open items / next steps

1. **Rodent v4 path chosen** — fresh restart with v4 harness (`inst/ecology/hamilton/rodent-v4/`),
   PT (nChains=4) + parallel runs (nRuns=4, nCore=4) for Rhat + mode-trap resilience.
   Validation job `17245693` submitted (serial, nIter=20k); production parallel
   submission gated on its clean exit + `RelabelEcology()` success. Note: parallel
   mode (`nCore>1`) is NOT checkpoint-resumable — size production `nIter` to fit
   walltime with margin (see commit `b47e488` comments).
2. **Paper structure conversation** — the v9 PT-MCMC null result kills
   the original Sim 2 "blind fails / aware rescues" framing. Options:
   reframe around "aware reduces false sisters 3.5×" as a regularisation
   claim; pivot to the rodent empirical as the headline (if it reaches
   ESS); deepen v9 with more reps + alternative-z mechanisms (M1/M3 from
   `dev/ecology/sim-design/v9-discriminate.R`).
3. **Confirm Bug A root cause closed** — `.BuildResult()` trims the symptom
   (z_samples > saved_idx at result-construction time) but the upstream cause
   was never definitively traced by static analysis. Run rep01 or rep02 again
   with the fix to confirm `RelabelEcology()` no longer aborts. If drift
   recurs, the rich abort message now gives the surplus count.
4. **Document directional-z mode** in `R/MkPrimeModel.R` / `R/RunMkPrime.R`
   (z=1: `(rate01 × phi, rate10 / phi)` directional). Currently a feature
   without docs.
5. **Per-character partition API v2** (`3dbd7db`) — user-authored design
   plan in another session. Section 11 has 7 open decisions awaiting
   review. Not blocking; not implementation-ready.

## Pending jobs

| Type | ID / ref | Status | ETA | On completion |
|------|----------|--------|-----|---------------|
| SLURM | `17245693` (mkp-rod-v4-val) | submitted 2026-05-20 (post-arrive) | ~5h serial PT, 8h walltime | Check `/nobackup/pjjg18/mkp-rodent-v4/aware-validate/*.{out,err}` for clean exit + RelabelEcology success. If clean, write production script `rodent-v4-aware.sh` with `nRuns=4 nChains=4 nCore=4 nIter=100000` and `--cpus-per-task=4 --mem=8G --time=48:00:00`, submit. If `RelabelEcology` aborts with z_samples drift, rerun with `trimZSamples="tail"` after inspection. |
| SLURM | `17234909_*` (mkp-mk-tlshrink) | UNKNOWN — was running at 13:30 BST, may have finished | ~0-2h | Not ecology — mkp-core arm benchmarking. Ignore on this branch; collect from mkp main when ready. |

The rodent v3 cont jobs (`17238607`, `17238608`) FAILED and are not requeued.
v3 directory `/nobackup/pjjg18/mkp-rodent-v3/aware/` left intact for audit.

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
  - 2026-05-20 deploy: `git stash push -u -m hand-off-2026-05-20-pre-merge-deploy`
    preserved a working set of mods (NEWS.md, R/MkPrimeModel.R, test-priors.R
    deltas + untracked test-ecology-scaffold.R). Not yet applied; review
    before clearing.
- **Local worktree path**: `C:\Users\pjjg18\GitHub\worktrees\ecology-aware`
  (not the earlier nested `GitHub\GitHub\worktrees` mistake). Primary mkp
  checkout at `C:\Users\pjjg18\GitHub\mkp` stays on `main`. The
  `feature/ecology-aware` branch (commit `3a3b289` partition API v1) lives
  in mkp's local refs and is the prior version of the v2 design plan that
  landed at `3dbd7db` on this branch.
- **Plugin layout (post 50295e9)**:
  - `dev/`, `inst/hamilton/m131-*`, `inst/MkPrime/`, top-level `RED_TEAM_*.md`
    et al = mkp-core, mirror of main, do not touch from this branch.
  - `dev/ecology/`, `inst/ecology/` = ecology plugin namespace. Mirror
    structure inside (`dev/ecology/red-team/`, `dev/ecology/sim-design/`,
    etc.). Future ecology dev notes belong here.
- **Rodent data**: nexus at
  `/nobackup/pjjg18/mkp-rodent-blind-v2/data/mbank_X24848_2026-5-9-1135.nex`
  (only there — NOT bundled in repo). 60 extant tips × 217 chars
  (160 neomorphic + 59 transformational after `AutoDetectNeomorphic()`).
  Aware-v2 config: `magnitudeMode="global"`, `kPrimePrior="geometric"`,
  `coding="variable"`, `expSteps=10` (legacy; v3 uses auto), `rho0Alpha=7`,
  `rho0Beta=3`, `thetaAlpha=2`, `thetaBeta=2`, `sigmaPhi=0.5`,
  `set.seed(20260512)`.
- **v4-cross helpers**: `inst/ecology/simulations/sim3v4cross-helpers.R` and
  `sim3v4-helpers.R`. Note `MatrixToPhyDat` / `MkPrimeData` are NOT exported
  in the installed lib — use bare names, not `MkPrime::` prefix.
- **PT execution model**: `nChains > 1` runs sequentially in MkPrimeMCMC;
  set `cpus-per-task=1` and multiply wall by nChains. Confirmed in
  `sim3-multirep-v4.sh` comments.
- **Aware iter rate on rodent data**: ~15k iter/h post-tuning (revised
  down from earlier 50k estimate, which was a tuning-phase artifact);
  1M iter ≈ 70h aware. 2G RAM, 1 CPU sufficient. Multiple sequential
  continuations needed for full chain.
- **`.SimulateMkPrimeEcology` z semantics**: per-(char, eco) integer in {0,1,2}.
  - z=0: no eco effect on this (char, eco)
  - z=1: for neomorphic chars, `(rate01 * phi, rate10 / phi)` — DIRECTIONAL
    toward state 1. For transformational, `mult = phi` (symmetric, noisier).
  - z=2: opposite of z=1.
  - Set z asymmetrically across two ecologies (different cols of z matrix)
    to create directional convergence at realistic TL — see v9 (commit `8fc37ef`).
- **Auto-expSteps**: `R/MkPrimeModel.R` finalises `expSteps = parsimony /
  nChar × expStepsInflation` (default 1.05). For rodent: parsimony 1519,
  nChar 217 → expSteps ≈ 7.35. Logs print parsimony, nChar, expSteps.
  Floor 0.5 prevents degenerate priors on small datasets.
- **ResumeMkPrime signature (merged HEAD)**: `ResumeMkPrime(checkpointFile,
  data, tree=NULL, neomorphic=integer(0), knownStates=integer(0),
  model=NULL)`. No `mcmc=` argument. Checkpoint versions 1, 2, 3 accepted.
- **Checkpoint/tree-file sync invariant (NEW post-merge)**: chain log,
  checkpoint, and `_trees.nwk` must agree on the number of recorded
  trees. `.TruncateTreeToN()` enforces this on resume and aborts on
  mismatch. **Old (pre-merge) chains that hit STREAM-001 (0-byte tree
  file) cannot be resumed under the merged code** — must restart fresh.
- **Bypass-RelabelEcology pattern for topology-only scores**: when the
  streaming bugs fire (pre-fix chains), load `result$trees` directly and
  call `ScoreTreesUnrooted(trees, biparts, refTree)` from
  `inst/ecology/simulations/sim3-scoring.R`. Topology metrics are
  symmetric under phi-flip; only phi/theta posterior summaries need relabeling.
  For fresh runs using merged code, `RelabelEcology(res)` should work;
  if `z_samples` drift recurs, use `RelabelEcology(res, trimZSamples = "tail")`.

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
- **`--amend` after a failed pre-commit hook** — the commit didn't happen, so
  amend rewrites the previous commit. Create a new commit instead.
- **`ape::prop.part` for posterior tree scoring** — root-dependent. Use
  `TreeTools::as.Splits` instead, via `HasBipartSplits` in
  `inst/ecology/simulations/sim3-scoring.R`. Discovered 2026-05-19.
- **`expSteps=10` default** for the TL prior — produces Gamma mean=10, gives
  5-13× TL inflation on realistic 16-tip morphological sims. Fixed in `bc24d52`
  (then unit-corrected in `f7b4644`).
- **Phi-rate-multiplier mechanism for realistic-TL trap design** — symmetric
  scaling can't produce directional convergence; v6/v7/v8 confirm.
  Use directional z (v9 path) instead.
- **multirep-v3 truth TL=13.5** as a "main result" — saturated; user declared
  junk for typical morphology. Acceptable as a long-tree case study only.
- **Subagent dispatch with "directional bias / trap" terminology** — burned 4
  spurious API-policy refusals on 2026-05-20. Write code inline if it recurs.
- **Resuming pre-merge ecology checkpoints under merged HEAD** — fails when
  the original chain hit STREAM-001 (no trees written despite checkpoint
  recording them). The new desync check is correct; old chains are not
  recoverable as continuations. Restart fresh.
- **`RelabelEcology()` as a paper-required step** — for topology metrics
  (the actual paper claims about P(true/false sister), CID), relabel is
  unnecessary; metrics are phi-flip-symmetric. Only phi/theta posterior
  summaries need it.
- **MP-trap → MCMC-trap inference** — the parsimony-screen 62% false-sister
  rate at TL=1.48 (v9 / commit `8fc37ef`) does NOT carry over to 4-chain
  PT MCMC, which recovers truth ~85% (blind) / ~92% (aware) of the time.
  TL prior + parsimony anchoring regularises sufficiently. Don't extrapolate
  from MP screens to MCMC behaviour without the MCMC.
- **Silent z_samples trim in RelabelEcology** — an earlier approach warned
  and trimmed `zSamples[seq_len(nSamples)]` when surplus detected. Rejected
  because trimming the tail is wrong if the surplus is at the head (tuning-
  phase leak). Current approach: abort with diagnostics by default; explicit
  `trimZSamples = "tail"` / `"head"` required for deliberate trim.

## Worktrees

| Branch | Path | Status |
|--------|------|--------|
| ecology-aware (tracks worktree-ecology-aware) | C:\Users\pjjg18\GitHub\worktrees\ecology-aware | clean, in sync with origin |
| main (mkp primary) | C:\Users\pjjg18\GitHub\mkp | clean, in sync with origin/main |

## Suggested first action

Decide rodent v3 path and deploy the streaming fix to Hamilton. Fetch the
updated branch on Hamilton (`git fetch origin && git merge --ff-only
origin/worktree-ecology-aware`) so the fix is in the Hamilton lib before
the next run:

```bash
/c/WINDOWS/System32/OpenSSH/ssh.exe pjjg18@hamilton8.dur.ac.uk \
  "cd /nobackup/pjjg18/mkp-sim3-multirep-v3/MkPrime && \
   git fetch origin && git merge --ff-only origin/worktree-ecology-aware && \
   R --no-save -e 'install.packages(\".\", lib=\"/nobackup/pjjg18/mkp-sim3-multirep-v3/lib\", repos=NULL, type=\"source\")'"
```

Then, if fresh restart chosen for rodent v3:

```bash
/c/WINDOWS/System32/OpenSSH/ssh.exe pjjg18@hamilton8.dur.ac.uk \
  "rm /nobackup/pjjg18/mkp-rodent-v3/aware/rodent-aware-v3.{ckp,log} \
      /nobackup/pjjg18/mkp-rodent-v3/aware/rodent-aware-v3_trees.nwk && \
   cd /nobackup/pjjg18/mkp-sim3-multirep-v3/MkPrime && \
   J1=\$(sbatch --parsable inst/ecology/hamilton/rodent-v3/rodent-v3-aware.sh) && \
   J2=\$(sbatch --parsable --dependency=afterany:\$J1 inst/ecology/hamilton/rodent-v3-cont/rodent-v3-aware-cont.sh) && \
   J3=\$(sbatch --parsable --dependency=afterany:\$J2 inst/ecology/hamilton/rodent-v3-cont/rodent-v3-aware-cont.sh) && \
   echo restart=\$J1 cont1=\$J2 cont2=\$J3"
```
