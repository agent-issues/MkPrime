# dev/red-team

Working directory for the `/red-team` skill (project-level adversarial review
rotation). The skill definition holds the doctrine; this README records what
lives here and where a finding's life actually happens.

**The live work list is
[GitHub Issues on `agent-issues/MkPrime`](https://github.com/agent-issues/MkPrime/issues)**
(label `red-team` + `sev:*` + `area:N`). The markdown tables in this directory are
**frozen history**, kept so a finder can tell a fresh bug from a re-found one.
Nothing here carries a status column: `Fixes #N` closes a finding on merge, so
status cannot drift from merge state, and concurrent branches have no shared file
to conflict over.

## Which repo, and why a fork

`Mk-prime/r` is private, so the skill's default would be to file directly on it.
We deliberately use the `agent-issues/MkPrime` mirror instead, for the reason
`../TreeSearch` does: `ms609-agent` has write on the fork and **no access at all**
to `Mk-prime/r`, so every PR is reviewable by the maintainer (GitHub will not let
an account approve its own PR) and nothing reaches the release repo until the
human syncs the fork. See the mirror table in `../../AGENTS.md`.

Token for every write here is `CLAUDE_GH_TOKEN`. It is fine-grained and
repository-*selected*, so a newly created repo in the org is invisible to it
until it is both added to the selection **and** approved by an org owner — see
`~/.claude/CLAUDE.md`. The prefix
goes on `gh api graphql` for **Discussions** too, not just issues and PRs: a
discussion posted under the wrong account must be deleted and rewritten, there is
no re-attribution.

## Layout

```
dev/red-team/
  README.md                 # this file
  focus-areas.md            # rotation table (11 areas) — start_tier per area, edited rarely
  discussion-categories.md  # round records: category list, posting recipe, identity rules
  log.md                    # round records pre-Discussions + the model-version legend + last_focus
  findings-archive.md       # FROZEN 2026-09-18 — terminal-state findings; anti-duplication memory only
  migration-map.tsv         # old <PREFIX>-nnn id -> issue number / archive entry
  proofs/                   # analytic derivations backing specific findings
  patches/                  # standalone patches (toggles/experiments) referenced by findings
  numerical/                # numerical-auditor notes + stress drivers
  heavy-tests/              # standalone R harnesses (SBC/coverage/detailed-balance) + fixtures
  repros/                   # minimal reproducers cited from issue bodies
  reviews/                  # per-branch external-review reports
  notes/                    # investigation write-ups too long for an issue body
```

Legacy per-round investigation files (`MARGINAL-K-*.md`, `harness-001-investigation.md`,
`campaign-2026-05-26.md`) stay where they are; issue bodies cite them by path.

## Finding ids

Historical ids (`EG-001`, `LIKE-001`, `TREEMOVE-001`, `STREAM-003`, …) are
**frozen, not retired** — they persist in source comments, `log.md`, plan files
and memory. `migration-map.tsv` resolves each to its issue number or archive
entry. **Never mint a new prefixed id**: the GitHub issue number is the id from
now on. When searching for duplicates, do **not** filter on a legacy prefix —
new issues carry none.

## Escalation

`needs-escalation` is a GitHub label, not a file: routing checks
`gh issue list --label needs-escalation,area:N --state open` first. Add
`escalation-backlog.md` only if a soft signal appears that names no single issue.

## Tiers vs. model versions

`opus` / `fable` are **rungs on a capability ladder, not models** — the Agent tool
accepts only the rung alias, so which version a rung resolves to changes as new
models ship. Every round therefore stamps the version that actually ran, and a
rung version bump is its own revisit trigger: a seam that ran dry at `opus-4.8` is
re-visited at `opus-5` *before* escalating to `fable` or being left dormant.
Backward-looking verdicts (`ran dry`, `dormant`) are version-scoped; forward-looking
routing text ("escalate to opus") stays unversioned so it resolves to the current
version at dispatch. The alias-to-version legend is at the top of `log.md` —
reconcile it at the start of every round.

Run outputs under `heavy-tests/` and `numerical/` (`*-results/`, `*.log`, experiment
`*.rds`) are regenerable and git-ignored; committed regression fixtures are the
exception.
