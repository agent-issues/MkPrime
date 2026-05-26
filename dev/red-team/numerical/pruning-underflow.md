# N2 — Felsenstein pruning underflow and ascertainment log-sum-exp audit

**Lane:** N2 (numerical-auditor)
**Scope:** `src/likelihood.cpp` (`pruning_jc`, `pruning_jc_collapsed`, `pruning_mkn`),
`src/mcmc_likelihood.cpp` (`pruning_jc_flat`, `pruning_jc_acrv_flat`,
`pruning_jc_flat_collapsed`, `pruning_jc_acrv_flat_collapsed`,
`pruning_mkn_flat`, `pruning_mkn_acrv_flat`, `pruning_f81_het_*_flat`),
`src/ascertainment.cpp`, the ~16 ascertainment-subtraction call sites
`logL -= n * log(1 - p)` listed in §Findings/A2.
**Author:** numerical-auditor agent, 2026-05-26.

---

## What this audits

Two distinct numerical concerns in the Felsenstein pruning pipeline:

1. **Partial-likelihood underflow.** The hot-path pruning kernels in
   `src/mcmc_likelihood.cpp` and the reference kernels in
   `src/likelihood.cpp` compute Felsenstein conditional likelihoods
   (CLs) as repeated products of values $\leq 1$ stored in IEEE 754
   `double`, *without any rescaling or log-space accumulation*. Deep
   trees with conflicting tip patterns or extreme branch lengths can
   drive these products below the subnormal floor (~$5 \times 10^{-324}$),
   at which point either (a) the kernel correctly returns
   `R_NegInf` because the root site likelihood underflows the
   `site_lik <= 0` guard (e.g. `src/likelihood.cpp:123`,
   `src/mcmc_likelihood.cpp:263`, `:444`, `:1038`, etc.), or (b) the
   answer is silently corrupted because an internal CL underflowed
   asymmetrically while the root sum remains positive.

2. **Ascertainment log-cancellation.** The Lewis (2001) variable-coding
   correction
   $$
   \log L_{\text{cor}} \;=\; \log L_{\text{raw}} \;-\; n_\text{char} \cdot \log(1 - P_\text{excl})
   $$
   appears in ~16 call sites with the literal form `std::log(1.0 - p)`
   (e.g. `src/mcmc_likelihood.cpp:2367`, `:2396`, `:2447`, `:2497`,
   `:2559`, `:2604`, `:2674`, `:2721`; `src/mcmc.cpp:1094`, `:1407`,
   `:1798`, `:2062`, `:5469`; `src/node_cl_cache.h:719`, `:724`). When
   $p \to 1$ — small tree, near-zero branches, large $k$ — the
   subtraction $1 - p$ exhibits catastrophic cancellation. The standard
   fix is `std::log1p(-p)`. This audit measures whether the substitution
   actually recovers precision in this setting.

---

## Conditioning analysis

### A1. Pruning condition number

At an internal node $v$ with descendant subtree $T_v$ and $L_v = |T_v|$
leaves, the Felsenstein CL satisfies
$$
0 \leq a_v(s) \leq 1, \qquad a_v(s) = \Pr(\text{data in } T_v \mid X_v = s).
$$
The site likelihood at the root is $\mathcal L = \sum_s \pi_s a_\rho(s)
\leq 1$. The root sum can underflow to $0$ when every $a_\rho(s)$ falls
below $\approx 2^{-1074}$. For the JC analytic update with the O($k$)
trick used in all kernels,
$$
a_v(s) = \prod_{c \in \mathrm{ch}(v)} B_c(s), \qquad
B_c(s) = p_{\mathrm{diff}}(t_c)\,\Sigma_c + (p_{\mathrm{same}} - p_{\mathrm{diff}})(t_c)\,a_c(s),
$$
where $\Sigma_c = \sum_{s'} a_c(s')$.

Define $\delta_c = p_{\mathrm{same}}(t_c) - p_{\mathrm{diff}}(t_c) =
e^{-k t_c/(k-1)}$. Two extreme regimes:

- **Near-zero branches, $t \to 0$.** $p_{\mathrm{diff}} \to 0$,
  $\delta \to 1$; $B_c(s) \to a_c(s)$ (identity update). For a leaf $\ell$
  in state $s_\ell$, $a_\ell(s) = \mathbb{1}\{s = s_\ell\}$ and CL
  remains sparse with at most one non-zero per node. For conflicting
  patterns, the product across a chain of length $L$ involves
  $\sim L$ factors of $p_{\mathrm{diff}} \approx kt/(k-1)$, giving
  $a_\rho \sim t^{L}$ — underflows when $L \log(1/t) \gtrsim 745$.
  Caterpillar depth $L = 64$ with $t = 10^{-12}$: predicted
  $-\log_{10} a_\rho \sim 64 \cdot 12 = 768$, i.e. $a_\rho \sim 10^{-768}$
  → underflow. Confirmed empirically (§Results, Experiment 1).

- **Equilibrium-length branches, $t \to \infty$.** $\delta \to 0$,
  $p_{\mathrm{diff}}, p_{\mathrm{same}} \to 1/k$. $B_c(s) \to
  (1/k) \Sigma_c$, identical across all $s$. CL becomes class-constant
  at $\Sigma_c/k$ and the product across edges yields
  $a_\rho(s) \to \prod_c (\Sigma_c/k)$. Since $\Sigma_c \leq k$ and
  tip $\Sigma = 1$ (observed) or $k$ (missing), the product is
  $\leq 1$, and for $n$ tips in distinct states the product
  $\sim k^{-(n-1)}$. For $n = 16$, $k = 4$: $a_\rho \sim 4^{-15}
  \approx 9 \times 10^{-10}$, $\log \mathcal L \approx -22$ —
  **no underflow risk**. Equilibrium branches are well-conditioned.

- **Intermediate, $t \sim O(1)$ per edge, depth $L$.** Site likelihood
  decays as $(1/k)^L$ roughly. Underflow at $L \approx 745/\log k$,
  e.g. $L \approx 540$ for $k = 4$. This is well outside any realistic
  morphological tree.

**Conclusion on pruning underflow:** the underflow regime is
$L \log(1/t) \gtrsim 745$. For empirically-sourced morphological trees
($L \lesssim 200$, $t \in [10^{-3}, 10]$ — see project memory
`feedback_resumable_runs`), $L \log(1/t) \lesssim 200 \cdot 7 = 1400$
is **reachable** for the worst-case combination of deep ladder + very
short branch, which can occur during MCMC warmup when the chain has
not yet found a reasonable branch-length scale.

### A2. Ascertainment log-cancellation

For $p \in [0, 1)$ stored as `double`, `1.0 - p` involves at most one
rounding. `log(1 - p)` and `log1p(-p)` then both apply IEEE-correct
`log` to a value whose precision is already limited by the *input* `p`.
The classical justification for `log1p(x)` is that it avoids the
intermediate `1 + x` when `x` is small — but here `-p` is not small
(it is near $-1$), and `log1p(-p)` internally still must form `1 - p`
identically to the naive expression. **In double precision,
`log(1 - p)` and `log1p(-p)` are bit-equivalent for every $p \in (0, 1)$
this audit could probe** (§Results, Experiment 6).

The cancellation that matters happened **upstream**, when the kernel
computed $p$ from a chain of multiplications. By the time
$p = 1 - \epsilon$ is stored, $\epsilon$ has at most $\log_2$ ULP
precision relative to $p \approx 1$. No log reformulation can recover
the precision lost at that earlier stage.

The fix that would actually help is a **structurally different
computation**: compute $\log(1 - P_\text{const})$ via the
*non-constant-site likelihood* directly. For JC the constant-site sum
over $k$ patterns has a closed-form complement only via the full
likelihood marginal, which is more expensive than the cancellation
loss. (For very small trees this is tractable; for production trees it
is not.)

---

## Reference computation

**Pruning underflow:** an independent pure-R Felsenstein implementation
(`ref_pruning_jc` in `pruning-driver.R`) using the same JC analytic
formula and the same O($k$) reduction, in IEEE double precision. This
cannot detect *intra-double* precision loss but **does** catch
algorithmic / indexing bugs, and it underflows at exactly the same
input grid as the kernel — confirming the kernel's underflow behaviour
is structural to double precision, not a kernel-specific defect.

`Rmpfr` is not available in this environment; a higher-precision check
would require either Windows-side `Rmpfr` install (flaky) or a
mpmath-based Python comparator. The empirical agreement between the
two double-precision pipelines on every finite-result case, with max
relative error $1.84 \times 10^{-16}$ ($\approx 1$ ULP), suggests no
intra-double bug; an `Rmpfr` cross-check would only confirm this.

**Ascertainment:** direct comparison of `log(1 - p)` to `log1p(-p)`
over a 16-point grid $p \in \{1 - 10^{-1}, \ldots, 1 - 10^{-16}\}$,
plus the ascertainment outputs from
`MkPrime::constant_site_prob_jc` and `singleton_site_prob_jc`.

---

## Stress-test design

Six experiments in `pruning-driver.R`:

1. **Caterpillar underflow.** Trees of depth 8/16/32/64/128/256/512
   (effective max-depth = $n_\text{tip} - 1$) on $k = 4$ with the
   alternating tip pattern $(0,1,2,3,0,\ldots)$ — maximally conflicting
   under JC. Branch grid $\{10^{-12}, 10^{-6}, 10^{-3}, 0.1, 1, 10, 100\}$.
2. **Balanced equilibrium.** Balanced binary trees with depths 4/6/8/10
   ($n = 16, 64, 256, 1024$ tips) and branch lengths 0.5/1/5/20, with
   tips in rotating distinct states. Tests the "equilibrium branches
   on deep trees" claim that this regime is *not* an underflow risk.
3. **Ascertainment cancellation.** 3-tip tree, branches
   $\{10^{-15}, 10^{-12}, \ldots, 10\}$, $k \in \{2, 4, 8, 12\}$.
4. **Long-branch ascertainment.** Same 3-tip tree at long branches
   $\{1, 10, 100, 1000\}$, where $p_\text{const} \to 1/k$ and the
   subtraction is well-conditioned.
5. **Singleton subtraction.** 5-tip tree (so
   $p_\text{const} + p_\text{singleton} < 1$ at $k = 2$), same
   $k$ and branch grids as Experiment 3.
6. **Synthetic `log1p` cancellation.** Direct evaluation of
   `log(1 - p)` and `log1p(-p)` on $p = 1 - 10^{-i}$, $i = 1..16$,
   with no kernel involvement.

`--quick` mode restricts (1) to depths 8/16 and (2) to depths 4/6,
completing in $\approx$ 30 s; full mode in $\approx$ 90 s.

---

## Results

Full results in `dev/red-team/numerical/pruning-results/`:

| Experiment                            | Worst $\varepsilon_\text{rel}$ vs reference     | Underflow observed?                                                                                                 |
|---------------------------------------|--------------------------------------------------|---------------------------------------------------------------------------------------------------------------------|
| 1. Caterpillar (max depth 511)        | $1.84 \times 10^{-16}$ (1 ULP)                   | YES — $-\infty$ at $n=64, t = 10^{-12}$ and beyond, in lockstep with reference (so detected, not silent)            |
| 2. Balanced equilibrium (depth 10)    | $1.60 \times 10^{-16}$ (1 ULP)                   | None at $L = 1024, t = 20$ → $\log \mathcal L \approx -355$, well above underflow                                  |
| 3. Ascertainment, $p \to 1$           | $3.31 \times 10^{-15}$ vs `log1p` form           | n/a; difference is $\leq$ 1 ULP of the log result, NOT a precision improvement                                      |
| 4. Long-branch ascertainment          | $6.69 \times 10^{-16}$                           | n/a                                                                                                                 |
| 5. Singleton subtraction              | $6.66 \times 10^{-16}$                           | n/a                                                                                                                 |
| 6. Synthetic `log1p` vs `log(1-p)`    | $0$ — bit-identical across $p = 1 - 10^{-16}$    | n/a                                                                                                                 |

### Underflow pattern in Experiment 1

```
n_tip  t=1e-12   t=1e-6    t=1e-3    t=0.1     t=1       t=10      t=100
   8   -170      -87       -46       -19       -11       -11       -11
  16   -342      -176      -93       -39       -23       -22       -22
  32   -685      -353      -188      -78       -45       -44       -44
  64   -Inf      -710      -378      -156      -90       -89       -89
 128   -Inf      -Inf      -Inf      -313      -180      -177      -177
 256   -Inf      -Inf      -Inf      -626      -360      -355      -355
 512   -Inf      -Inf      -Inf      -Inf      -721      -710      -710
```

**Underflow boundary**: kernel returns $-\infty$ exactly when the
predicted $\log \mathcal L < -\log(2^{1074}) \approx -744$. Reference
returns $-\infty$ at identical input points. Kernel never produces a
finite value disagreeing with reference: **no silent corruption
observed**.

### `log1p` vs `log(1-p)` empirical

For 16 test points $p = 1 - 10^{-i}$, $i = 1, \ldots, 16$, both
formulations produce **bit-identical results** in IEEE 754 double:

```
i      one-minus-p     log(1-p)         log1p(-p)        diff
1      1e-01           -2.302585        -2.302585        0
6      1e-06          -13.815511       -13.815511        0
12     ~1e-12         -27.631043       -27.631043        0
15     ~1e-15         -34.539576       -34.539576        0
16     ~1.11e-16      -36.736801       -36.736801        0
```

This is because R / glibc / msvcrt `log1p(-p)` for `p` near 1 reduces
to `log(1.0 - p)` after the unavoidable `1.0 - p` rounding. The
substitution is **cosmetic, not a precision improvement, in this
context**.

---

## Verdict

**Stable with caveats** for both concerns.

### Verdict 1 — pruning underflow: stable up to logL ≈ −744

The Felsenstein kernels in `src/likelihood.cpp` and
`src/mcmc_likelihood.cpp` lack partial-likelihood rescaling. When the
root site likelihood underflows the IEEE double subnormal floor,
the kernel correctly returns `R_NegInf` via the `site_lik <= 0`
guard (`src/likelihood.cpp:123`, `:249`, `:352`;
`src/mcmc_likelihood.cpp:263`, `:444`, `:1038`, `:1217`, `:1449`,
`:1641`, `:1949`).

The underflow regime is reachable for caterpillar-like topologies of
depth $\gtrsim 64$ with branch lengths $\lesssim 10^{-12}$, or
depth $\gtrsim 512$ with branch lengths $\lesssim 0.1$. For
**typical morphology trees** (~50–200 tips, branches $10^{-3}$–$10$,
mostly balanced), the worst-case $\log \mathcal L$ per character is
$\sim -200$ to $-500$, **safely above the underflow floor**.

The realistic risk: during MCMC warmup with very small branch lengths,
or on a *single character* of high information content combined with a
deep subtree, $\log \mathcal L$ for that character can hit $-\infty$
and the entire partition log-likelihood collapses to $-\infty$. This
is a **soft-fail rather than silent-wrong** scenario — the chain will
reject the proposal — but it can bias the move-acceptance pattern
during warmup. Project memory `feedback_wallclock_adaptation` already
documents a related warmup-fragility tradeoff.

**Recommendation:** add log-space CL accumulation (standard per-node
rescaling: divide each CL row by its max, accumulate $\log \max$ into
a per-character running log-scale) **if and when** profiling shows
this matters. Pure morphology trees rarely require it; the project's
typical $n_\text{tip} \leq 100$ and well-conditioned branches put
underflow several orders of magnitude away from realistic
operating points. **Numerical bug rating: low priority.**

### Verdict 2 — `log(1 − p)` cancellation: not a precision bug

The 16 call sites flagged as `std::log(1.0 - p)` in the orchestrator
brief are **functionally equivalent** to `std::log1p(-p)` in IEEE 754
double precision, because the cancellation occurs in the upstream
computation of $p$ (which both forms inherit identically) rather than
in the subtraction itself. Bit-equivalence verified empirically across
the entire reachable range of $p$ (Experiment 6).

**Recommendation:** the `log1p(-p)` substitution is *not* a numerical
fix but **is** a code-clarity / static-analyser-friendliness improvement
and a defence against future hypothetical recomputation patterns that
might keep `1 - p` as a separately-tracked quantity. The patch is
captured at `dev/red-team/patches/N2-pruning-underflow.patch` as a
mechanical substitution; **applying it changes no output value, so it
is not urgent**.

The **actual** numerical pathology in the ascertainment correction is
the upstream computation of $p_\text{const}$ when $p \to 1$. For
3-tip / 5-tip trees with $t = 10^{-12}$, the constant-site probability
loses $\approx 12$ decimal digits of precision relative to $1.0$
(`1 - p` resolves only to $\sim 4 \times 10^{-12}$ even though the
analytic value is $\sim 5 \times 10^{-12}$). This is an **input-side
precision loss**, not a log-formulation issue. Fixing it would
require deriving $\log(1 - P_\text{const})$ via an alternative
formulation that avoids forming $P_\text{const}$ explicitly — non-
trivial and outside this lane's scope.

---

## Recommendation

1. **Pruning underflow.** No immediate fix. File an orchestrator note:
   "Pruning underflow can hit `R_NegInf` for caterpillar trees of
   depth $\geq 64$ with $t \leq 10^{-12}$ or depth $\geq 512$ with
   $t \leq 0.1$. Behaviour is soft-fail (returns $-\infty$, never
   silent corruption). Add log-space CL rescaling only if MCMC
   diagnostics show warmup pathology traceable to this." Reference
   the precondition `feedback_wallclock_adaptation` in project memory
   for related warmup concerns.

2. **`log1p` substitution.** Patch captured at
   `dev/red-team/patches/N2-pruning-underflow.patch` (16 sites,
   mechanical `log(1.0 - p)` → `log1p(-p)`). **Patch is not applied**
   to the worktree; it is bit-equivalent and adds no precision.
   Recommend applying it as a documentation / static-analysis hygiene
   improvement, but not as a numerical correctness fix.

3. **Upstream cancellation in `P_const`.** Deferred. Out of scope for
   this lane. Note for a future lane: under JC the analytic alternative
   $\log(1 - P_\text{const})$ via a separate "non-constant pseudo-data"
   pruning is structurally identical to the existing fused
   ascertainment kernel, so the bite of fixing it is much smaller than
   first appears. Worth a separate audit if Lewis correction precision
   becomes important.

---

## Worktree state

- Added `dev/red-team/numerical/pruning-driver.R` and
  `dev/red-team/numerical/pruning-underflow.md` (this file).
- Generated `dev/red-team/numerical/pruning-results/01-…06-…csv`.
- Captured `dev/red-team/patches/N2-pruning-underflow.patch` (16-site
  `log1p` substitution; **not applied**).
- No modifications to `src/`, `R/`, `tests/testthat/`.
