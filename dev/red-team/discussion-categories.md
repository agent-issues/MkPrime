# Discussion categories — one per focus area

Round records live in GitHub Discussions on `agent-issues/MkPrime` (private repo
=> private discussions). Convention matches `agent-issues/TreeSearch` and
`StratoBayes/StratoBayes-R-package`: the name is `NN · <area name>`, the separator
is a middle dot (U+00B7), and **GitHub derives the slug automatically** — dropping
`&`, `:`, `+` and parentheses. **The leading `NN` is the only part the tooling
depends on.** The slugs in the table below were **read back from the live board on
2026-09-25** and are confirmed, not predicted; re-read them with the query at the
bottom of this file if a category is ever renamed.

**Format must be "Open-ended discussion"** for every area category. *Announcement*
format restricts posting to maintainers; `ms609-agent` is a collaborator, so its
post is rejected with a bare `Resource not accessible by personal access token`,
indistinguishable from a token-scope problem.

There is **no API for any of this** — no `createDiscussionCategory`, no delete, no
rename. Categories are hand-managed in Settings -> Discussions, so every change
costs manual UI work and cannot be scripted or reviewed in a diff.

## ⚠️ GitHub caps a repository at 25 discussion categories

**Live as of 2026-09-25: 17 of 25.** 13 area categories plus four retained defaults
(`Announcements`, `General`, `Ideas`, `Q&A`); `Polls` and `Show and tell` were
deleted at setup. That leaves **8 free slots**.

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
| 3 | `03 · Model finalisation & stepping stone` | `03-model-finalisation-stepping-stone` |
| 4 | `04 · Research harnesses & result provenance` | `04-research-harnesses-result-provenance` |
| 5 | `05 · Move scheduler, tuning & tempering` | `05-move-scheduler-tuning-tempering` |
| 6 | `06 · Streaming, checkpoint, resume & GUI` | `06-streaming-checkpoint-resume-gui` |
| 7 | `07 · Tree-move proposals & Rcpp surface` | `07-tree-move-proposals-rcpp-surface` |
| 8 | `08 · Likelihood kernel & partial CLs` | `08-likelihood-kernel-partial-cls` |
| 9 | `09 · Convergence diagnostics & ESS` | `09-convergence-diagnostics-ess` |
| 10 | `10 · Test infrastructure` | `10-test-infrastructure` |
| 11 | `11 · RB-oracle equivalence` | `11-rb-oracle-equivalence` |
| 12 | `12 · Red-team process meta-review` | `12-red-team-process-meta-review` |
| 13 | `13 · Partition API & hierarchical hyperprior` | `13-partition-api-hierarchical-hyperprior` |

Profiling does **not** get a category: the `/profile` skill keeps round records in
`dev/profiling/log.md` by design (one serial rotation, and `log.md` doubles as the
context `baselines.md` regressions are read against).

## Reading the category ids

Desktop sessions only. A cloud session cannot reach GraphQL and does not need the
ids: see [From a cloud session](#from-a-cloud-session).

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

**In a cloud session, skip to [From a cloud session](#from-a-cloud-session)** —
the command below cannot run there, and a round record filed as an issue instead
lands under `ms609`.

The `createDiscussion` mutation needs the repository id and the category id from
the query above, and **must** carry the agent token — a discussion posted under the
maintainer's account cannot be re-attributed, only deleted and rewritten:

```bash
GH_TOKEN=$CLAUDE_GH_TOKEN gh api graphql -F query=@/tmp/post.graphql
```

Title format stamps the rung **and the version that ran**, per the skill:
`RT <date> - area <N> - <rung> (<Version>) - yield <n> (<h>h/<m>m/<l>l)`.
`yield` counts confirmed candidates, and the `h/m/l` split is by confirmed
candidate too, not by filed issue; say in the body how many issues they became.

### From a cloud session

The cloud proxy blocks GraphQL and replaces any token, `CLAUDE_GH_TOKEN` included,
with the maintainer's own, so nothing above runs and anything created directly
over REST is authored by `ms609`. Two workflows post as `ms609-agent` instead,
triggered by a REST `workflow_dispatch`:

```bash
gh api repos/agent-issues/MkPrime/actions/workflows/post-discussion.yml/dispatches --method POST --input payload.json
```

| Workflow | `inputs` |
|----------|----------|
| `post-discussion.yml` | `category` (a slug from the table above), `title`, `body`. To comment on an existing record: `discussion` (its number) and `body` |
| `post-as-agent.yml` | `action` (`issue`, `comment` or `pr`), `body`; `title` and comma-separated `labels` for an issue; `number` for a comment on an issue or PR; `title`, `head`, `base`, `draft` for a PR |

The workflow resolves the slug to the category id itself. The payload is
`{"ref":"main","inputs":{...}}` with every value a string; build it with `jq`
so the Markdown body is escaped for you:

```bash
jq -n --arg c 09-convergence-diagnostics-ess --arg t "RT <date> - area 9 - ..." --rawfile b body.md \
  '{ref: "main", inputs: {category: $c, title: $t, body: $b}}' > payload.json
```

**Payload limit.** A `workflow_dispatch` payload must stay under about 65 KB; a
50 KB body worked. Put the overflow of a longer record in a comment on it.

**Confirm the result.** The dispatch returns 204 and no run id, so read the
newest run of the workflow until it is `completed` with conclusion `success`:

```bash
gh api "repos/agent-issues/MkPrime/actions/workflows/post-discussion.yml/runs?per_page=1" \
  --jq '.workflow_runs[0] | {status, conclusion, html_url}'
```

Then read the record back and check that its author is `ms609-agent`. The REST
list is oldest first, so the new discussion is the last line:

```bash
gh api --paginate "repos/agent-issues/MkPrime/discussions?per_page=100" \
  --jq '.[] | "\(.number) \(.user.login) \(.title)"' | tail -1
```

For an issue or a comment posted through `post-as-agent.yml`:

```bash
gh api "repos/agent-issues/MkPrime/issues?creator=ms609-agent&state=all&per_page=1" --jq '.[0] | {number, title}'
gh api "repos/agent-issues/MkPrime/issues/<n>/comments?per_page=100" --jq '.[-1] | {user: .user.login, html_url}'
```
