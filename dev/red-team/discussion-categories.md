# Discussion categories — one per focus area

Round records live in GitHub Discussions on `agent-issues/MkPrime` (private repo
=> private discussions). Convention matches `agent-issues/TreeSearch` and
`StratoBayes/StratoBayes-R-package`: the name is `NN · <area name>`, the separator
is a middle dot (U+00B7), and **GitHub derives the slug automatically** — dropping
`&`, `:`, `+` and parentheses. **The leading `NN` is the only part the tooling
depends on.** The slugs in the table below were **read back from the live board on
2026-09-18** and are confirmed, not predicted; re-read them with the query at the
bottom of this file if a category is ever renamed.

**Format must be "Open-ended discussion"** for every area category. *Announcement*
format restricts posting to maintainers; `ms609-agent` is a collaborator, so its
post is rejected with a bare `Resource not accessible by personal access token`,
indistinguishable from a token-scope problem.

There is **no API for any of this** — no `createDiscussionCategory`, no delete, no
rename. Categories are hand-managed in Settings -> Discussions, so every change
costs manual UI work and cannot be scripted or reviewed in a diff.

## ⚠️ GitHub caps a repository at 25 discussion categories

**Live as of 2026-09-18: 16 of 25.** 12 area categories plus four retained defaults
(`Announcements`, `General`, `Ideas`, `Q&A`); `Polls` and `Show and tell` were
deleted at setup. That leaves **9 free slots**.

The cap bites late and hurts: GitHub makes you re-home every post in a category
before it can be deleted, so a spare default is free to remove while empty and
expensive later. StratoBayes hit 25/25 and can no longer add a focus area at all.
If the headroom here drops below about three, drop `Ideas` and `Q&A` before they
accumulate posts.

Adding a row to `focus-areas.md` later means creating its `area:N` label **and** its
category — the label half is scriptable, the category half is not, and an issue
filed against a missing label simply fails.

## The categories

| # | Category name | Slug (confirmed) |
|---|---|---|
| 1 | `01 · k′ prior families & normalisation` | `01-k-prior-families-normalisation` |
| 2 | `02 · k′ marginalisation & the p-samplers` | `02-k-marginalisation-the-p-samplers` |
| 3 | `03 · Logseries prior & LogPrior semantics` | `03-logseries-prior-logprior-semantics` |
| 4 | `04 · Hamilton harness: combine, recovery, ETAs` | `04-hamilton-harness-combine-recovery-etas` |
| 5 | `05 · Move scheduler & adaptive weights` | `05-move-scheduler-adaptive-weights` |
| 6 | `06 · Streaming, checkpoint & interrupt recovery` | `06-streaming-checkpoint-interrupt-recovery` |
| 7 | `07 · Tree-move proposals & Rcpp surface` | `07-tree-move-proposals-rcpp-surface` |
| 8 | `08 · Likelihood & partial CLs` | `08-likelihood-partial-cls` |
| 9 | `09 · Convergence diagnostics & ESS` | `09-convergence-diagnostics-ess` |
| 10 | `10 · Test-suite health` | `10-test-suite-health` |
| 11 | `11 · RB-oracle equivalence` | `11-rb-oracle-equivalence` |
| 12 | `12 · Red-team process meta-review` | `12-red-team-process-meta-review` |

**Categories 1 and 2 are pending a rename on the live board.** They still read
`01 · Empirical-geometric prior math & normalisation` and `02 · EG p-sampler: case 30
logit-MH`, from before those areas were rescoped on 2026-09-18. The table above is the
target. Nothing is broken meanwhile — only the leading `NN` is read by tooling, and a
rename does not move existing posts (discussion #28 sits in category 2).

Profiling does **not** get a category: the `/profile` skill keeps round records in
`dev/profiling/log.md` by design (one serial rotation, and `log.md` doubles as the
context `baselines.md` regressions are read against).

## Reading the category ids

Reads need no token prefix, but PowerShell 5.1 mangles an inline GraphQL query
(it strips the inner double quotes, and `owner: agent-issues` then fails to parse).
Write the query to a file and pass it with `-F query=@<file>`:

```bash
cat > /tmp/cats.graphql <<'Q'
{
  repository(owner: "agent-issues", name: "MkPrime") {
    id
    discussionCategories(first: 30) { nodes { id name slug } }
  }
}
Q
gh api graphql -F query=@/tmp/cats.graphql
```

## Posting a round record

The `createDiscussion` mutation needs the repository id and the category id from
the query above, and **must** carry the agent token — a discussion posted under the
maintainer's account cannot be re-attributed, only deleted and rewritten:

```bash
GH_TOKEN=$CLAUDE_GH_TOKEN gh api graphql -F query=@/tmp/post.graphql
```

Title format stamps the rung **and the version that ran**, per the skill:
`Round N — area #M (<name>) — <rung>-<version> — <date>`.
