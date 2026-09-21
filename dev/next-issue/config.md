# `/next-issue` config — MkPrime

Read by `~/.claude/skills/next-issue/SKILL.md`. Only what varies per repo lives
here; the doctrine is in the skill.

| Key | Value |
|-----|-------|
| `base_branch` | `main` on `origin` = `agent-issues/MkPrime`. **Never `upstream` (`Mk-prime/r`)** — its push URL is `no-push-use-gha` and a direct commit there turns every later fork sync into a conflicted merge. |
| `grouping` | `labels`, prefix `area:` (1–12). Fall back to `hot_files` for an unlabelled issue. |
| `exclude_label` | *None yet.* `needs-decision` is not created on this repo — exclude by judgment per the skill, and say so in the report. |
| `identity` | `agent` — `ms609-agent`, token env var `CLAUDE_GH_TOKEN`. Table in `~/.claude/CLAUDE.md`. |

## `hot_files` — never split across parallel chips

```
src/mcmc.cpp              # hot loop: propose → evaluate → accept/reject
src/mcmc_likelihood.cpp   # Felsenstein pruning, per-partition lik, Gibbs helpers
src/mcmc_state.h          # McmcData / McmcState structs + declarations
src/node_cl_cache.h       # partial-CL cache, dirty-flag invalidation
src/tree_moves.cpp        # SPR / NNI / TBR / pSPR
R/RunMkPrime.R            # MCMC engine, move dispatch, adaptation
```

`src/RcppExports.cpp` and `src/init.cpp` are **append-only shared files** — a
touch there is a coordination cost, not a dependency, and still branches from
`base_branch`.

## `pr_command`

```bash
GH_TOKEN=$CLAUDE_GH_TOKEN gh pr create --base main --head <branch> --reviewer ms609 --body-file <file>
```

Reads need no token prefix; `gh repo set-default` already points at the fork.
Never `export GH_TOKEN` — it silently re-identifies anything the human runs in
the same shell.

## `build` / `test`

From the parent `GitHub/` directory, with a per-chip agent id:

```bash
ID="a$(printf %s "$BRANCH" | sha1sum | cut -c1-5)"
bash build-agent.sh mkp "$ID"
bash test-agent.sh  mkp "$ID" [filter]
```

Builds a renamed package (`MkPrime.$ID`) into `.builds/MkPrime-$ID/`, avoiding
DLL lock conflicts. Targeted filter always — never a full suite.

The id **must** be derived from the branch, not the session: a resumed chip gets
a new session id and would silently build a second tree. There is no concurrency
cap — CPU oversubscription degrades gracefully, whereas a shared build id
corrupts, and unique ids are what prevent that.

## `branch_rule`

**Never `git checkout`, `switch`, `stash` or `reset --hard` in
`C:/Users/pjjg18/GitHub/mkp`.** Other agents share that checkout and it must stay
on `main`. All branch work happens in a named worktree alongside it (`mkp-gibbs/`,
`mkp-parallel/`); see `../AGENTS.md` → *Worktree discipline*.

Coordination files (`completed-tasks.md`, `coordination.md`,
`remote-jobs.md`, `dev/plans/`) are always read and written in `../mkp/`, never
in the worktree.

## `extra`

- **`RunMkPrime()` runs must be triple-guarded** — `maxTime` + `setTimeLimit()`
  + bash `timeout`. Template in `.AGENTS/memory/subprocess-runmkprime.md`.
  Convergence-criterion runs (`minEss`, `minTreeEss`, `maxRhat`) **must** also
  set `maxTime`.
- **Validation is GHA-first**; local builds are for iteration only.
- `src/Makevars.win` must never be committed — profiling flags live there.
- Seed reproducibility is not a feature: only the shape of a *converged* run is
  expected to reproduce. Never cite "changes RNG streams" as a cost.
- MCMC state and log-column names are **snake_case** (`tree_length`,
  `rate_loss`, `k_prime`) — domain terminology, deliberate divergence from the
  `r-conventions` default.
