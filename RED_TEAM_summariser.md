# Red-team audit — summariser & CID-to-truth (2026-05-19)

Subject: paradox that `mk_k40` beats truth-matching `mk` on mean CID, when
k_true ~ Geom(0.4)+2 (≈ floor of 2, decaying tail; k=40 is far above the
generative support).

Files audited
- `data-raw/hamilton/summarize_streamed.R` (committed canonical;
  `tmp_summarize_streamed.R` not present in worktree — references in
  `RESUME.md` describe an earlier transient file)
- `data-raw/hamilton/run_one.R` (NB: missing `mk_k40` branch on disk; patched
  in via `tmp_add_k40.py` on Hamilton — see WARNING below)
- `dev/pilots/2026-05-12-prior-validation/analysis/cid_eight_prior.R`
- `dev/pilots/2026-05-12-prior-validation/analysis/cid_nine_prior.R`
- `R/RunMkPrime.R` tree-write paths (L378, L970)
- `R/proposals.R` / `src/proposals.cpp` (rooting convention)
- `dev/red-team/findings.md` (prior rounds)

---

## Bottom line

**No correctness bug is biasing the CID comparison in favour of mk_k40.**

The CID computation itself is correct. The two real defects I found
(BURNIN-1, MK-K40-THIN) bias in the *opposite* direction — they
should disadvantage `mk_k40` relative to `mk`, so they cannot explain the
paradox.

The most likely explanation is a **property of mean CID under wide
posteriors** (POSTWIDTH-1 below): when posterior support is broad and
CID is bounded above, the *mean* CID is lowered by occasional near-truth
hits. mk_k40 (k far above the data's natural support) has a flatter
posterior surface and so a wider posterior than mk (k=kObs floor, very
constrained). The mean-CID-to-truth ranking is then **not** a clean proxy
for posterior concentration on the truth.

If the paradox is taken at face value, switch the headline statistic to
*median* CID, or to P(min_posterior_CID < ε), and re-run. The CID
*values* are fine; the *aggregation* is what's misleading.

---

## Findings table

| ID            | Severity | Status      | What                                                                                       |
|---------------|----------|-------------|--------------------------------------------------------------------------------------------|
| ROOT-1        | NOTE     | RULED OUT   | True tree NOT unrooted in streamed summariser; CID is rooting-invariant — no effect        |
| BURNIN-1      | WARNING  | CONFIRMED   | Tail-5000 thinning gives long chains a more burnin-pure window than short chains           |
| MK-K40-THIN   | WARNING  | CONFIRMED   | `mk` uses `thin=10`, `mk_k40` uses `thin=100` → 10× more trees per iter for mk             |
| MK-K40-MISSING| CRITICAL | PROCESS BUG | `mk_k40` branch absent from committed `run_one.R`; patched on Hamilton only via tmp script |
| POSTWIDTH-1   | WARNING  | LIKELY CAUSE| mean-CID under bounded metric favours wide posteriors                                      |
| AGG-1         | NOTE     | OK          | Mean-of-means across tasks is correct; no Reduce/intersect bug                             |
| TIPLABEL-1    | NOTE     | OK          | Tip labels are identical across arms for a given task (same NJ start tree)                 |
| CIDCALL-1     | NOTE     | OK          | `ClusteringInfoDistance(..., normalize=TRUE)` called identically across arms               |
| PARSER-1      | NOTE     | OK          | `count.fields(sep="\n")` + `readLines` correctly counts and reads tree lines               |
| KP-COL-1      | NOTE     | OK          | kPrime_ column alignment fix from commit c112167 is in place                               |

Detailed entries below.

---

### ROOT-1 (NOTE, RULED OUT) — true tree not unrooted, but CID is rooting-invariant

`summarize_streamed.R:142–143` reads the truth tree and only reorders it:
```r
true_tree <- read.tree(true_tree_file)
true_tree <- ape::reorder.phylo(true_tree, "cladewise")
```
By contrast, the legacy `run_one.R` "combine" mode at L83 *does* call
`UnrootTree(true_tree)`. The streamed summariser omits this.

MCMC trees are written by `R/RunMkPrime.R:970` with
`Nnode = length(tipLabels) - 1L` — this is the rooted count. The phylo
written is therefore rooted (a single root in the C++ proposal impl, see
`src/proposals.cpp:26` and `dev/red-team/findings.md` TREEMOVE-004).

So we are comparing rooted MCMC trees against the (rooted) true tree.
That's fine, **but** even if rooting were inconsistent between
arms, `TreeDist::ClusteringInfoDistance` operates on splits and is
provably invariant to rooting:

```r
> ClusteringInfoDistance(t1, t2, normalize=TRUE)       # rooted-rooted: 0.6754
> ClusteringInfoDistance(unroot(t1), unroot(t2), …)    # unrooted-unrooted: 0.6754
> ClusteringInfoDistance(t1, unroot(t2), …)            # mixed: 0.6754
```

User's previously-documented `ape::prop.part` root-dependency bug
(`memory: project_scoring_bug.md`) does NOT propagate to
ClusteringInfoDistance, which uses TreeTools' canonical Splits, not
`prop.part`. **Not a source of bias.**

Cosmetic: add `true_tree <- TreeTools::UnrootTree(true_tree)` to
`summarize_streamed.R:142` for consistency with `run_one.R` combine mode.

---

### BURNIN-1 (WARNING, CONFIRMED) — tail-5000 window favours long chains

`summarize_streamed.R:118–124`:
```r
tail_n   <- min(n_trees, max(n_thin, 5L * n_thin))     # = min(n_trees, 5000)
start_at <- n_trees - tail_n + 1L
keep_idx <- if (tail_n <= n_thin) seq.int(start_at, n_trees)
            else round(seq.int(start_at, n_trees, length.out = n_thin))
```

Behaviour by chain length (default `n_thin = 1000`):

| Trees written | tail_n | Effect                                          |
|---------------|--------|-------------------------------------------------|
| 800           | 800    | Keep ALL 800 trees (early ones included)        |
| 2,500         | 2,500  | Keep ALL 2,500, uniformly subsample to 1,000    |
| 5,000         | 5,000  | Keep last 5,000, subsample to 1,000             |
| 50,000        | 5,000  | Keep ONLY last 5,000, subsample to 1,000        |

A chain that ran longer gets a window starting at iter (n−5000); a chain
that ran shorter gets a window starting at iter 1. Although `_trees.nwk`
only receives post-warmup samples (`R/RunMkPrime.R:932`), the chain
continues to mix and converge after warmup — so the early post-warmup
trees are systematically further from truth than late ones.

**Reproduction** (slow-convergence multi-NNI chain, 15 tips):
```
long  chain (50000 iters → n_kept=1000): mean CID = 0.0673  (window: last 5000)
short chain  (2500 iters → n_kept=1000): mean CID = 0.1873  (window: ALL 2500)
```
A 2.8× CID inflation for the chain that ran fewer iters — purely from
the windowing rule.

**Direction**: this disadvantages arms that wrote fewer trees. Since
`mk_k40` is slower per iter (larger CL matrices) AND samples are
thinned 10× more (see MK-K40-THIN), `mk_k40` writes fewer trees and so
should be **penalised** by this windowing. The paradox is that
`mk_k40` wins anyway → this bug doesn't explain it, but it does mean
the *true* gap between mk_k40 and mk is even bigger than the summariser
reports.

**Fix**: discard a fraction-based burnin (`tail_n <- floor(n_trees / 2)`)
or, better, threshold-based on a convergence diagnostic. Make it
explicit so it's identical across arms.

---

### MK-K40-THIN (WARNING, CONFIRMED) — `thin` differs by 10× across arms

`data-raw/hamilton/run_one.R:161` defines:
```r
make_mcmc <- function(prefix, thin_iters = 10L)
```
and the `mk` arm at L218 takes the default `thin_iters = 10L`. Every
other arm (`mk_kp2`, `mk_k9`, `mk_k15`, `mk_k24`, `mkp_geo`, plus
`mk_k40` patched in by `tmp_add_k40.py`) explicitly passes
`thin_iters = 100L`.

This means, in 8h walltime:
- `mk` writes ~10× as many `_trees.nwk` lines per iter as the others
- `mk_k40`'s 2251–5601 trees (per task brief) correspond to
  ~225k–560k iterations
- `mk`'s tree count, if any, will be ~50k–500k trees

Combined with BURNIN-1, the CID comparison is being made over
non-comparable subsets of the chain. Whether the *number of MCMC iters*
or the *number of stored trees* is the right denominator is a design
question; right now the summariser uses neither — it uses "last 5000
stored trees" which conflates the two.

Note: in `R/MkPrimeMCMC.R`, scalar thinning and tree thinning may differ
via `treeThin`; for the production runs `treeThin = NULL` (default) so
tree-write frequency = scalar-write frequency = `thin`. Confirmed at
`tests/testthat/test-tree-thin.R:11`.

**Fix**: use the same `thin` for the mk arm as for the others, OR change
the summariser to thin by iteration index, not by stored-tree index.

---

### MK-K40-MISSING (CRITICAL, PROCESS BUG) — mk_k40 not in committed run_one.R

Commit `b71f52b` titled "feat: add mk_k40 arm — slurm script, summariser,
and arm list updates" actually does NOT update `data-raw/hamilton/run_one.R`.
Files changed:
```
data-raw/hamilton/mk_k40_array.slurm    +23
data-raw/hamilton/summarize_array.slurm  +-
tmp_add_k40.py                           +106
tmp_summarize_streamed.R                 +-
```
`run_one.R`'s `match.arg(args[5], c(…))` list ends with `"mk_k24",
"mkp", "mkp_eg", "mkp_geo", "combine"` — no `mk_k40`. Running the
committed file with `arm=mk_k40` would `stop()` immediately.

The production Hamilton install has been patched in-place by
`tmp_add_k40.py` (writes to `/nobackup/pjjg18/mkp-study/run_one.R`).
This patch is **not under version control** and not auditable
post-hoc. Anyone re-running this study from a fresh git checkout will
get a different program. Sibling concern in `dev/red-team/findings.md`
TREEMOVE-001: lost fixes due to "claimed fixed but never written".

This is a process bug, not a correctness bug per se, but it makes the
`mk_k40` results essentially impossible to verify from the repo. It's
also indistinguishable, from the commit log, whether `thin_iters = 100L`
was actually used for `mk_k40` (the python patcher template shows it
uses 100L, which matches `mk_k24`).

**Fix**: commit the run_one.R changes (or rerun the patch and commit).
Remove `tmp_add_k40.py` after.

---

### POSTWIDTH-1 (WARNING, LIKELY CAUSE OF PARADOX) — mean CID rewards wide posteriors

When the metric (CID) is bounded above (random pairs give ≈0.7–0.8) and
the *mean* is taken over a multi-modal/wide posterior, posteriors that
spend any meaningful mass near the truth have their mean dragged down by
the lower tail, regardless of where the posterior MODE sits.

Reproduction:
```
wide   posterior (uniform 0..12 NNI moves from truth): mean CID = 0.211
narrow posterior (always 6 NNI moves from truth):       mean CID = 0.239
```
The wide posterior wins on the mean despite having a less informative
shape relative to truth.

With k=40 vs k=kObs floor, the likelihood is *much* flatter for
`mk_k40` — k=40 lets every column accept a wide range of state-rate
configurations, so partial likelihoods are nearly indifferent to many
tree topologies. With k tight to kObs, the model is over-constrained for
some columns (especially those where the generative k_true ≥ kObs+2) and
*forced* into a more peaked but possibly biased posterior. Under a
mean-CID metric, the looser model can win even when its posterior MODE
is no closer to truth.

**Diagnostic to confirm**: re-compute
- `median(cid)` per task, then aggregate
- `min(cid)` per task, then aggregate
- `mean(cid >  threshold)` per task

If mk_k40's advantage shrinks or reverses under `median(cid)` or
`min(cid)`, POSTWIDTH-1 is the explanation. The cid_eight_prior.R and
cid_nine_prior.R scripts already compute `median_cid` (line 33–34 of
cid_eight_prior.R) — running those alongside the means would settle it.

---

### Other items verified clean

- **AGG-1 (OK)** — `cid_nine_prior.R:33–34` computes `mean_cid` per
  task; `cid_nine_prior.R:50` aggregates across tasks with simple
  `mean`. No nested-list flattening bug, no `Reduce`/intersect bug
  in pair-finding.
- **TIPLABEL-1 (OK)** — Both arms use `start_tree <- NJTree(pd)`
  (`run_one.R:158`), so `tip.label` is identical across arms for the
  same task. MCMC tree writes use `tipLabels` (`R/RunMkPrime.R:965`),
  preserving this. True tree from `tree_NN/tree.nwk` shares the
  simulated taxa.
- **CIDCALL-1 (OK)** — Single call site
  `ClusteringInfoDistance(reorder_multiPhylo(trees), true_tree,
  normalize = TRUE)` at `summarize_streamed.R:153`, identical for all
  arms (arm is just an output-file prefix at this point).
- **PARSER-1 (OK)** — `cat(write.tree(tr), "\n", file=…)` produces
  `tree;␣\n` with a trailing space-newline; `read.tree(text=…)` parses
  this correctly. `count.fields(sep="\n", quote="")` returns one entry
  per line.
- **KP-COL-1 (OK)** — kPrime_ column grep at L89 is correct; commit
  `c112167` ("Hamilton mkp study: mk_kp2/mk_k9 arms, fixed kp alignment
  in summariser") fixed a previous alignment bug.

---

## Recommendation

Before pivoting on the mk_k40 paradox, run two checks:

1. **Replace mean with median in `cid_nine_prior.R`**. The columns are
   already in the RDS files. If mk_k40's lead disappears under median,
   POSTWIDTH-1 is the story. Two lines to edit.

2. **Equalise thin across arms.** Either rerun `mk` with `thin = 100L`,
   or have the summariser thin to a fixed *iteration spacing* not a
   fixed *stored-tree count*. Otherwise BURNIN-1 + MK-K40-THIN make
   the comparison meaningless.

If the paradox survives both, then mk_k40 genuinely fits the trees
better than mk under the data's actual k distribution — which is a
substantive scientific finding (the geometric prior with floor=kObs is
too sharp; flat-but-wide priors win), not a bug.

3. **Commit the patched `run_one.R`** so the result is reproducible.
