# Lane L4 — Mk' Relabelling Correction

> **Lane.** L4 — first-principles derivation of `mk_prime_relabel_log(kPrime, kObs)`.
>
> **Files in scope.** `src/corrections.cpp`, `src/mcmc_likelihood.cpp`,
> `src/mcmc.cpp`, `src/node_cl_cache.h`, `R/MkPrimeModel.R`, `R/likelihood.R`,
> `tests/testthat/test-corrections.R`.
>
> **Status.** A repository-level proof of the same identity already lives at
> `./relabelling-correction-proof.md` (the file that motivated the original
> bugfix). This document is the red-team's independent rederivation; we cite
> the prior proof where useful but do not rely on it.

## Theorem (informal)

For a transformational character with `kObs` distinct observed states modelled
as Mk' with `kPrime = k' ≥ kObs` latent states under Jukes–Cantor symmetry, the
correction added to Felsenstein's log-likelihood to obtain the marginal log
data probability is the log of the falling factorial

$$
\log P(k',\,k_\text{obs}) \;=\; \log\!\frac{k'!}{(k' - k_\text{obs})!}
\;=\; \mathrm{lgamma}(k'+1) - \mathrm{lgamma}(k'-k_\text{obs}+1).
$$

The implementation at `src/corrections.cpp:44` is correct, the prior in
`R/MkPrimeModel.R::LogPrior` does **not** redundantly include the correction,
and every call site of `mk_prime_relabel_log` adds it to the **log-likelihood**
(never the log-prior).

## Assumptions

1. **Mk' substitution kernel is Jukes–Cantor with `k'` states.** All
   off-diagonal rates are equal, every state has stationary frequency
   `π_s = 1/k'`. Confirmed by inspection of `pruning_jc*` (the dispatch in
   `src/mcmc_likelihood.cpp:2587-2595` calls `pruning_jc[*flat]` whenever
   `useHet` is false, and the root frequencies are constructed as
   `NumericVector(kp, 1.0/kp)`).
2. **Observed labels are an ordered injection into the latent state set.**
   The Felsenstein pruner consumes integer tip codes
   `0, 1, …, kObs−1` directly as model-state indices (see `pinfo.tipStates`
   construction in `prepare_mcmc_data`, `src/mcmc_likelihood.cpp:2914`). No
   relabelling is done at runtime; the canonical assignment
   "observed label `i` ↔ latent state `i`" is implicit.
3. **The observer is exchangeable over latent labels.** No prior knowledge
   distinguishes one latent state from another. Under Mk' this is automatic
   from JC symmetry; under F81-Het it is *not* (see Verdict).
4. **`kObs ≥ 1`** (every character has at least one observed state) and
   **`k' ≥ kObs`**. Enforced at `src/corrections.cpp:36-41` and by the
   prior boundary in `R/MkPrimeModel.R:408`.
5. **The MCMC samples `k'` (not the latent-state allocation).** That is, the
   sampler integrates analytically over which `kObs` of the `k'` latent
   states are "the observed ones". This is what motivates the correction at
   all.

## Statement (formal)

Let `T` be a tree topology with positive branch lengths, let
$\mathcal{X} \in \{0,\dots,k_\text{obs}-1\}^n$ be the observed tip-state vector
for one character (`n` tips, exactly `kObs` distinct values appear), and let
$L_\text{fels}(\mathcal{X} \mid T, k')$ denote Felsenstein's pruning likelihood
computed under JC(`k'`) with the canonical embedding
"`observed label i ↔ latent state i`, latent states `kObs..k'−1` unused at tips".

Then the marginal probability of $\mathcal{X}$ under Mk', integrated over the
unobserved injection of observed labels into the latent state space, is

$$
P(\mathcal{X} \mid T, k')
\;=\; \frac{k'!}{(k'-k_\text{obs})!}\;\cdot\;L_\text{fels}(\mathcal{X} \mid T, k')
\quad\Longleftrightarrow\quad
\log P = \log L_\text{fels} + \underbrace{\log\frac{k'!}{(k'-k_\text{obs})!}}_{\text{relabel correction}}.
$$

## Proof

### Step 1 — Felsenstein's pruning sums over one labelling.

By construction, $L_\text{fels}(\mathcal{X} \mid T, k')$ marginalises over all
internal-node assignments and all root states, weighted by `π_s = 1/k'`. It
treats the observed tip codes `0..kObs−1` as **fixed model-state indices**
matching latent states `0..kObs−1` (Felsenstein 2004, §16.3). Latent states
`kObs..k'−1` receive zero weight at every tip in this canonical assignment.

### Step 2 — Counting equivalent labellings.

Suppose instead we draw a labelling: an injective map
$\sigma : \{0,\dots,k_\text{obs}-1\} \hookrightarrow \{0,\dots,k'-1\}$ from
observed labels to latent states. We then form the tip vector
$\sigma(\mathcal{X})$ (apply $\sigma$ element-wise) and evaluate the
Felsenstein likelihood under JC(`k'`):

$$
L_\text{fels}(\sigma(\mathcal{X}) \mid T, k').
$$

Under JC(`k'`), the kernel is invariant under any permutation
$\pi \in S_{k'}$ of latent labels (every state has identical off-diagonal
rate and identical stationary frequency `1/k'`). Therefore, for every
$\sigma$,

$$
L_\text{fels}(\sigma(\mathcal{X}) \mid T, k') \;=\; L_\text{fels}(\mathcal{X} \mid T, k')
\quad\text{(JC symmetry).}
$$

### Step 3 — The number of distinct injections.

The set of injections $\sigma : \{0,\dots,k_\text{obs}-1\} \to \{0,\dots,k'-1\}$
has cardinality

$$
\#\{\sigma\} \;=\; k' (k'-1)(k'-2)\cdots(k'-k_\text{obs}+1)
\;=\; \frac{k'!}{(k'-k_\text{obs})!}.
$$

This is the **falling factorial** $(k')_{k_\text{obs}}$, also written
$P(k', k_\text{obs})$ (Riordan 1958, §1.2; Lange 2010, eq. 1.4).

### Step 4 — Marginal data probability.

Since the latent labels are exchangeable (Assumption 3) and the data could
equally plausibly have been generated under any of the
$k'!/(k'-k_\text{obs})!$ injections,

$$
P(\mathcal{X} \mid T, k')
\;=\; \sum_{\sigma}\; L_\text{fels}(\sigma(\mathcal{X}) \mid T, k')
\;=\; \frac{k'!}{(k'-k_\text{obs})!}\;\cdot\; L_\text{fels}(\mathcal{X} \mid T, k').
$$

Taking logs gives the claimed correction. $\blacksquare$

### Step 5 — Numerical sanity check (independent of the repository's existing proof).

For a 3-tip star tree, branch length `t = 0.3`, observed pattern `(0, 0, 1)`
(so `kObs = 2`), we computed $L_\text{fels}$ for every injective map and
summed:

| k' | n_maps | sum L  | per-map L  | ratio sum/canonical |
|----|--------|--------|------------|---------------------|
| 2  | 2      | 0.17470145 | 0.08735072 | 2  |
| 3  | 6      | 0.16461328 | 0.02743555 | 6  |
| 4  | 12     | 0.15880135 | 0.01323345 | 12 |
| 5  | 20     | 0.15532411 | 0.00776621 | 20 |
| 6  | 30     | 0.15303815 | 0.00510127 | 30 |

The ratio column is exactly $k'(k'-1) = (k')_2$ to machine precision, and
every map produces the same Felsenstein value, confirming both steps 2 and 3
empirically.

## Implementation cross-check

### `mk_prime_relabel_log` itself

`src/corrections.cpp:35-45` defines

```cpp
double mk_prime_relabel_log(int kPrime, int kObs) {
  if (kPrime < kObs) {
    Rcpp::stop("kPrime (%d) must be >= kObs (%d)", kPrime, kObs);
  }
  if (kObs < 1) {
    Rcpp::stop("kObs must be >= 1");
  }
  return std::lgamma(kPrime + 1.0) - std::lgamma(kPrime - kObs + 1.0);
}
```

This is exactly $\log(k'!) - \log((k'-k_\text{obs})!) = \log P(k',\,k_\text{obs})$.
The batch form at `src/corrections.cpp:58-72` repeats the same expression in a
loop. Both match the theorem.

Tests at `tests/testthat/test-corrections.R:4-87` cover (i) the
$k'=k_\text{obs}$ identity $\log(k!)$, (ii) monotonic increase in $k'$ at fixed
$k_\text{obs}$, (iii) the $k' < k_\text{obs}$ error, (iv) batch == scalar
agreement, (v) the closed form via `lgamma`, (vi) brute-force enumeration
against $k'(k'-1)$ for `kObs = 2` and $k'(k'-1)(k'-2)$ for `kObs = 3`, and
(vii) the MH-delta property $c(3,2)-c(2,2) = \log 3$. All consistent with the
derivation.

### Call sites (likelihood, not prior)

Every call site adds the correction to a log-likelihood accumulator:

- `src/mcmc_likelihood.cpp:2607-2610` — homogeneous-`k'` partition path,
  added inside the per-partition likelihood loop after the `pruning_jc*` call.
- `src/mcmc_likelihood.cpp:2726-2729` — heterogeneous-`k'` sub-group path,
  added once per character with that character's own `kPrimePart[ci]`.
- `src/mcmc.cpp:1105-1118` — added to `candLL[ci]` (likelihood of candidate
  topologies in a topological move).
- `src/mcmc.cpp:1422`, `1815`, `2077` — same pattern in three other topological
  proposals (`grep` confirmed identical structure `relabelCorr +=
  mk_prime_relabel_log(...)`).
- `src/mcmc.cpp:3830-3831` — Gibbs `k'` sweep: `ll += mk_prime_relabel_log(k,
  tp.kObs)` *before* `w = beta * ll + logPrior_k` (line 3833). The likelihood
  carries the correction; the prior is added separately.
- `src/mcmc.cpp:5480-5481` — partitioned topological-move analogue.
- `src/node_cl_cache.h:663-666` — node-cache replay; correction reattached
  during root-likelihood reduction.
- `R/likelihood.R:240-244` — the R-side reference likelihood: `ll <- ll +
  sum(relabelCorr)`, where `relabelCorr` is built by
  `mk_prime_relabel_log_batch`.

In every site the correction lives in the **likelihood** branch.

### Prior does NOT include the correction

`R/MkPrimeModel.R::LogPrior` (lines 393–572) is the canonical log-prior.
`grep -n relabel R/MkPrimeModel.R` yields only `@param`, the constructor
field, and the print method; the body of `LogPrior` never references
`mk_prime_relabel_log*`. The per-character prior is either Geometric, Beta-
Geometric, Logseries, or Empirical-Geometric on the **integer state count**
`k'` itself (lines 467–520). None of these terms contain a falling-factorial
piece in `k'`, so the correction is not double-counted.

### `kObs > k'` is rejected upstream

- `LogPrior` returns `-Inf` (line 408) if any
  `state$kPrime[transIdx] < mkd$kObs[transIdx]`.
- The Gibbs k' sweep uses `k = tp.kObs + ko` with `ko ≥ 0`
  (`src/mcmc.cpp:3694`), so `k ≥ kObs` is structural.
- If `mk_prime_relabel_log` is ever called with `kPrime < kObs` directly,
  `Rcpp::stop` aborts (`src/corrections.cpp:37`). There is no path in normal
  operation where the correction would silently return `-Inf` or NaN.

## Edge cases

### `k' = 1`, `kObs = 1`

$\log P(1, 1) = \log(1!/0!) = \log 1 = 0$. The correction is zero — there is
exactly one labelling. Consistent with the theorem and the implementation
(`lgamma(2) - lgamma(1) = 0 - 0 = 0`).

### `k' = kObs`

$\log P(k, k) = \log(k!/0!) = \log(k!)$. Constant in `k'` only when `kObs` is
also held fixed; otherwise (e.g. in a Gibbs sweep that varies `k'` for a
character with fixed `kObs`) it is **not** constant — moving `k' = kObs` to
`k' = kObs + 1` changes the correction by
$\log((k_\text{obs}+1)!/1!) - \log(k_\text{obs}!) = \log(k_\text{obs}+1)$.
Test `test-corrections.R:79-87` checks the specific case
$c(3,2) - c(2,2) = \log 3$ ✓.

Within MH ratios for moves that hold `(k', kObs)` fixed at the focal
character (e.g. topology moves, branch-length moves), the correction is
shared between numerator and denominator and cancels — so contributing
`log(k!)` to the log-posterior is harmless. It is retained so absolute
log-posterior values (e.g. for stepping-stone marginal likelihood
estimation) are well-defined.

### `kObs = 0`

The Rcpp guard at `src/corrections.cpp:39` rejects `kObs < 1`. This is the
right behaviour: a character with zero observed states is not part of the
data (it's a column of NAs and is filtered before assembly). If the user
constructs an `MkPrimeData` with a column where every entry is `NA`, the
data validator should drop it; this is outside the proof's scope but worth
flagging as an upstream invariant.

### `k' = 0`

Out of support — `kPrime ≥ kObs ≥ 1`, so `k' = 0` would already fail the
`kPrime < kObs` check (`src/corrections.cpp:36`). No defence in depth needed.

### Very large `k'`

`lgamma` is asymptotically $k' \log k' - k' + O(\log k')$ (Stirling).
$\log P(k', k_\text{obs}) \approx k_\text{obs} \log k'$ for large $k'$.
This grows without bound, but the prior on `k'` (Geometric / Logseries /
Empirical-Geometric) decays at least geometrically in `k'`, so the joint
posterior remains proper. No numerical concern: `lgamma(int)` is exact-ish
for `k' < 1e6` and Mk' state counts never exceed ~50 in practice.

## Connection to the k'-prior — explicit double-count check

The k'-prior `P(k')` (Geometric / Logseries / etc.) is defined on the
**integer state count itself**, not on labelled allocations. To verify there
is no double-count:

- **Geometric.** `lp <- lp + length(transIdx) * log(state$p) + sum(u) *
  log1p(-state$p)` (`R/MkPrimeModel.R:471`). This is
  $\log p + u \log(1-p)$ per character with $u = k' - k_\text{obs}$. Contains
  no factorial of `k'`.
- **Beta-Geometric.** `lp <- lp + sum(lbeta(alpha + 1, beta_ + u) -
  lbeta(alpha, beta_))` (`R/MkPrimeModel.R:506`). The `lbeta` terms are in
  `alpha`, `beta`, `u`; they expand as `lgamma(alpha+1) + lgamma(beta + u) -
  lgamma(alpha + beta + 1 + u) - [lgamma(alpha) + lgamma(beta) -
  lgamma(alpha+beta)]`. No `lgamma(k')` or `lgamma(k' - kObs)` term.
- **Logseries.** `lp <- lp + sum(kp * log(c_ls) - log(kp)) - length(transIdx)
  * log(-log1p(-c_ls))` (`R/MkPrimeModel.R:518`). The `- log(kp)` is
  $-\log k'$, **not** $-\log(k'!)$ or $\log((k'-k_\text{obs})!)$.
- **Empirical-Geometric.** `.LogPriorEmpiricalGeometric` (lines 350–379)
  computes a convolution over `i = 2..m` of `logEmp[i-1] + log p +
  (m-i)*log(1-p)`; the empirical pmf `logEmp` is built by
  `.LogPemp` (lines 318–334) and contains no factorial term.

None of the four priors contains $\log(k'!) - \log((k'-k_\text{obs})!)$ or
any algebraically equivalent expression. The relabelling correction lives
**only** in the likelihood. No double-count.

## Verdict

**Watertight with caveats.**

The derivation closes cleanly under Assumptions 1–5. The implementation at
`src/corrections.cpp:44` matches the theorem exactly; every call site adds it
to the log-likelihood, and the prior in `R/MkPrimeModel.R::LogPrior` does not
duplicate the term. Numerical brute-force enumeration on a star tree confirms
the falling factorial form to machine precision.

**Caveats and open issues:**

1. **JC symmetry is essential.** Step 2 of the proof relies on Assumption 1
   (JC). Under `qHeterogeneity = TRUE` the per-character kernel is an F81
   mixture (`pruning_f81_het_acrv_*` at `src/mcmc_likelihood.cpp:2656-2668`),
   and the relabelling correction is still applied (lines 2728). The two
   F81-Het branches in the code (`src/mcmc_likelihood.cpp:2607-2610` and
   `:2726-2729`) call `mk_prime_relabel_log` unconditionally — but F81 is
   **not** invariant under arbitrary latent-label permutations; only
   permutations that map the equilibrium-frequency vector
   $(\pi_0, \dots, \pi_{k'-1})$ to itself preserve the likelihood. If the
   Het bins are exchangeable in distribution (each bin is a symmetric
   marginal of `Dirichlet(α, …, α)`, see lines 78–95 of
   `R/MkPrimeModel.R`), then averaging over bins restores the exchangeability
   *in expectation*. This is a subtle point and the present proof does not
   discharge it; it is a candidate for a follow-up lane covering F81-Het.
2. **`kObs` itself is data-dependent.** The correction conditions on
   `kObs`. The full Mk' likelihood is the conditional
   $P(\mathcal{X} \mid T, k', k_\text{obs}) \times P(k_\text{obs} \mid k',
   \dots)$, but `kObs` is a deterministic function of `X`, so this collapses
   to $P(\mathcal{X} \mid T, k')$ as derived. No issue, but worth noting.
3. **Ascertainment correction.** The constant-site (and singleton-site)
   ascertainment correction (`coding != 0`) is applied **after**
   `pruning_jc*` returns and **before** the relabelling term (e.g.
   `src/mcmc_likelihood.cpp:2596-2605` then `:2607-2610`). The proof assumes
   no interaction; verifying that ascertainment and relabelling factorise
   correctly (i.e. that the relabelling correction does not need to be
   re-derived under the ascertainment-truncated sample space) is a separate
   exercise. Heuristically, since the truncation factor $1/(1 - p_\text{const})$
   depends only on the *pattern of state-occupancies* and not on the
   labelling, it is invariant under the same `k'!/(k'-k_\text{obs})!`
   relabellings; this means the correction should commute with ascertainment.
   A rigorous statement is left to a follow-up.

**No patch attached.** The implementation already encodes the correct
formula. (The historical bug it replaced is documented in
`./relabelling-correction-proof.md` and was already fixed in the worktree
base.)

**Worktree note.** No edits to `R/`, `src/`, or `tests/testthat/`. Only this
proof file was created; no commit, no merge, no push, per the
math-prover/red-team protocol.

## References

- Felsenstein, J. (2004) *Inferring Phylogenies*. Sinauer. §§16.3, 18.1.
- Lange, K. (2010) *Applied Probability* (2nd ed.). Springer. Eq. 1.4
  (falling factorial).
- Lewis, P.O. (2001) A likelihood approach to estimating phylogeny from
  discrete morphological character data. *Systematic Biology* 50: 913–925.
- Riordan, J. (1958) *An Introduction to Combinatorial Analysis*. Wiley.
- Tuffley, C. & Steel, M. (1997) Links between maximum likelihood and
  maximum parsimony under a simple model of site substitution.
  *Bulletin of Mathematical Biology* 59: 581–607.
- Yang, Z. (2014) *Molecular Evolution: A Statistical Approach*. OUP.
- Repository: `./relabelling-correction-proof.md` (the proof that motivated
  the original bugfix; this document is an independent rederivation,
  consistent with it on every check we ran).
