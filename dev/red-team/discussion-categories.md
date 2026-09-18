# Discussion categories — one per focus area

Round records live in GitHub Discussions on `agent-issues/MkPrime` (private repo
=> private discussions). Convention matches `agent-issues/TreeSearch` and
`StratoBayes/StratoBayes-R-package`: the name is `NN · <area name>`, the separator
is a middle dot (U+00B7), and **GitHub derives the slug automatically** — dropping
`&`, `:`, `+` and parentheses. **The leading `NN` is the only part the tooling
depends on**; treat the slugs below as predicted, and read the real ones back with
the query at the bottom of this file before the first post.

**Format must be "Open-ended discussion"** for every area category. *Announcement*
format restricts posting to maintainers; `ms609-agent` is a collaborator, so its
post is rejected with a bare `Resource not accessible by personal access token`,
indistinguishable from a token-scope problem.

There is **no API for any of this** — no `createDiscussionCategory`, no delete, no
rename. Categories are hand-managed in Settings -> Discussions, so every change
costs manual UI work and cannot be scripted or reviewed in a diff.

## ⚠️ GitHub caps a repository at 25 discussion categories

12 area categories + 6 GitHub defaults = 18/25. **Delete the spare defaults while
they are still empty**: GitHub makes you re-home every post in a category before it
can be deleted, so an empty `General`, `Ideas`, `Polls`, `Q&A` or `Show and tell` is
free to remove today and expensive to remove later. `Announcements` is the one
default worth keeping (maintainer notices). Doing that leaves 12 free slots.

Adding a row to `focus-areas.md` later means creating its `area:N` label **and** its
category — the label half is scriptable, the category half is not, and an issue
filed against a missing label simply fails.

## The categories

| # | Category name | Predicted slug |
|---|---|---|
| 1 | `01 · Empirical-geometric prior math & normalisation` | `01-empirical-geometric-prior-math-normalisation` |
| 2 | `02 · EG p-sampler: case 30 logit-MH` | `02-eg-p-sampler-case-30-logit-mh` |
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
GH_TOKEN=$CLAUDE_GH_TOKEN_MS609 gh api graphql -F query=@/tmp/post.graphql
```

Title format stamps the rung **and the version that ran**, per the skill:
`Round N — area #M (<name>) — <rung>-<version> — <date>`.
