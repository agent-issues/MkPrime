# Dispatcher Agent Brief

You are agent **{{AGENT_ID}}**, assigned to task **{{TASK_ID}}**: `{{TASK_ROW}}`

## Budget & Model

- **Budget**: {{BUDGET_MINUTES}} minutes (stay within this slice; if work won't fit, do a sub-step and check in)
- **Model assigned**: {{MODEL}}
- **Effort level**: {{EFFORT}}
- **Resume action** (if parked): {{RESUME_HINT}}

## Workflow

1. **Startup intake** (before claiming work):
   - Triage any new user reports (`a.*` and `u.*` files in project root)
   - Check `remote-jobs.md` for pending async results
   - See AGENTS.md for full protocols

2. **Read conventions**: See `AGENTS.md` (and `../AGENTS.md` for cross-package
   protocols) for:
   - Build/test/branch rules (GHA-first validation, `build-agent.sh` /
     `test-agent.sh`, `.builds/MkPrime-{{AGENT_ID}}/` isolation)
   - Shared-file coordination (append-only for `RcppExports.cpp`, `init.cpp`)
   - Feature branch lifecycle and mandatory pre-commit checks
   - Multi-agent workflow (worktree reserved tasks, user-report claim protocol)
   - Triple-guard for `RunMkPrime()` subprocesses
     (`.AGENTS/memory/subprocess-runmkprime.md`)

3. **Worktree rule**: MkPrime feature work happens in named worktrees
   alongside `mkp/` (e.g. `mkp-gibbs/`, `mkp-parallel/`). **Never** switch
   the main `C:/Users/pjjg18/GitHub/mkp` checkout to a different branch — it
   must stay on `main`. See `../AGENTS.md` → **Worktree discipline**.

4. **Build isolation**: From the parent `GitHub/` directory:
   ```bash
   bash build-agent.sh mkp {{AGENT_ID}}
   bash test-agent.sh  mkp {{AGENT_ID}} [filter]
   ```
   This builds a renamed package (`MkPrime.{{AGENT_ID}}`) into
   `.builds/MkPrime-{{AGENT_ID}}/`, avoiding DLL lock conflicts.

5. **Validation via GHA** (never run full test suites or R CMD check locally):
   - Push your branch: `git push -u origin feature/<name>`
   - Dispatch checks: `bash gha-dispatch.sh agent-check.yml feature/<name>`
   - Poll results: `bash gha-poll.sh <run_id>` (from another agent slice; don't block)

6. **Exit protocol**: run `bash dispatch.sh` from the repo root (it delegates
   to the global skill script).

   **When blocking on external wait** (GHA, Hamilton, human review):
   `bash dispatch.sh checkin {{AGENT_ID}} --kind=<gha|hamilton|human|other> --ref=<id> --eta=<iso-8601> --resume="<next action>"`

   Exit cleanly. The dispatcher will park this task and resume when the ETA passes.

   **When complete**:
   - Update `to-do.md` (delete task row; create new sections if needed)
   - Append summary row to `completed-tasks.md` under today's date
   - Run `bash dispatch.sh checkin {{AGENT_ID}} --done`
   - The dispatcher will mark the agent slot as free.

## Budget discipline

If the work won't fit in {{BUDGET_MINUTES}} minutes:
1. Do a **meaningful sub-step** (fix one bug, implement one small feature, resolve one blocker)
2. Check in with a resume action:
   `bash dispatch.sh checkin {{AGENT_ID}} --kind=other --eta=<next> --resume="<next step>"`
3. Exit cleanly rather than blowing the budget

## Tools

- `.AGENTS/memory/` — technical references (architecture, testing, performance, conventions)
- `bash dispatch.sh` — dispatcher CLI (locks, checkin, reap, kill)
- `gha-dispatch.sh` / `gha-poll.sh` — GitHub Actions integration (parent `GitHub/` directory)
- Claude Code skills — `Skill(skill: "hamilton-hpc")` for Hamilton SLURM,
  `Skill(skill: "r-package-profiling")` for VTune profiling
