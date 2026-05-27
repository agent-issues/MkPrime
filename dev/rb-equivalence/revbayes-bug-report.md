## [BUG] `srMinESS` triggers SIGFPE at iter 0 of MCMC simulation phase (non-MPI `mcmcmc`)

### Summary

Adding `srMinESS(N, "<log>.log", FALSE)` to the stopping-rules list of a non-MPI `mcmcmc(... nruns = 2, nchains = 4, tuneHeat = TRUE)` analysis causes the binary to die with `SIGFPE` (exit code 136 = 128 + signal 8) immediately after the burn-in completes — the `Iter 0` row is printed and the process aborts before generation 1. Removing only the `srMinESS` rule (keeping `srMaxTime`) makes the same script run to completion on identical data. We have reproduced this with two independent datasets on RevBayes 1.3.2 / Rocky Linux 8.10.

This looks closely related to closed issue [#1023](https://github.com/revbayes/revbayes/issues/1023) (segfault from `MinEssStoppingRule::getStatistic` at iter 0 on a cluster), which was labelled MPI-specific. The crash we observe is in the **non-MPI** binary, so we believe the underlying defect is broader than the MPI build: `MinEssStoppingRule` does not safely handle being queried before any post-burn-in samples have been written to the log file.

### Environment

- **RevBayes:** `1.4.0-dev` (single-thread `rb`, not `rb-mpi`)
- **Source:** branch `power-posterior-checkpoint`, commit `8f71e5f0532819bb7d28a7458c34e355407c5204` ("No more free infinite asterisks", 2026-05-27)
- **Build:** clean rebuild on 2026-05-27 16:09 UTC via `./build.sh clean && ./build.sh -j 8` with `module load gcc/11.2 boost`. Binary at `/nobackup/pjjg18/revbayes/projects/cmake/build/rb` (52 MB).
- **Runtime modules:** `module load gcc/13.2 boost`
- **OS:** Rocky Linux 8.10 (Green Obsidian), HPC cluster (Durham Hamilton8)
- **Job scheduler:** SLURM (single-cpu, 8 G RAM, `shared` partition)

### Reproduction

#### Data — `project950.neo.nex` (12 taxa × 2 neomorphic chars)

```nexus
#NEXUS
BEGIN DATA;
  DIMENSIONS NTAX=12 NCHAR=2;
  FORMAT DATATYPE=STANDARD MISSING=? GAP=- INTERLEAVE=NO symbols="0123456789";
  MATRIX
    Asterias_asteriformis             01
    Anemonia_sulcata_acanthoneme      11
    Aiptasia_sp.                      11
    Cereus_pedunculatus               11
    Diadumene_leucolena               11
    Haliplanella_lineata_pB2d         11
    Urticina_felina                   11
    Metridium_senile_pB1              11
    Palythoa_sp                       00
    Rhodactis_sp.                     00
    Stomphia_sp.                      11
    Actinostola_sp.                   01
  ;
END;
```

#### Data — `project950.trans.nex` (12 taxa × 6 transformational chars)

```nexus
#NEXUS
BEGIN DATA;
  DIMENSIONS NTAX=12 NCHAR=6;
  FORMAT DATATYPE=STANDARD MISSING=? GAP=- INTERLEAVE=NO symbols="0123456789";
  MATRIX
    Asterias_asteriformis             100010
    Anemonia_sulcata_acanthoneme      ?11110
    Aiptasia_sp.                      000001
    Cereus_pedunculatus               001101
    Diadumene_leucolena               101101
    Haliplanella_lineata_pB2d         000001
    Urticina_felina                   001101
    Metridium_senile_pB1              011101
    Palythoa_sp                       1102-?
    Rhodactis_sp.                     1102-0
    Stomphia_sp.                      010010
    Actinostola_sp.                   000001
  ;
END;
```

#### Model — `by_nt_9v.Rev`

Asymmetric neomorphic Q (free `fnFreeK` with separate gain/loss rates) plus a single `fnJC(9)` transformational partition under variable coding, both with discretised log-normal among-site rate variation.

```rev
# Separate Neomorphic partition
partitioned[1] <- neo
partitioned[2] <- trans
nChar <- v(neo.nchar(), trans.nchar())

# Dataset properties
taxa <- neo.names()
nTaxa <- neo.size()
nEdge <- 2 * nTaxa - 3

moves = VectorMoves()

# Uniform prior on tree topologies
topology ~ dnUniformTopology(taxa)
moves.append( mvNNI(topology, weight = nEdge / 2.0) )
moves.append( mvSPR(topology, weight = nEdge / 8.0) )

# Compound dirichlet prior on edge length variability
gamma_shape <- 2
exp_steps <- 1
tree_length ~ dnGamma(shape = gamma_shape, rate = gamma_shape / exp_steps)
moves.append( mvScale(tree_length, weight = 1.0) )

rel_br_lengths ~ dnDirichlet( rep(1.0, nEdge) )
moves.append( mvBetaSimplex(rel_br_lengths, weight = nEdge / 3.0) )
moves.append( mvDirichletSimplex(rel_br_lengths, weight = nEdge / 20.0) )
br_lengths := rel_br_lengths * tree_length

phylogeny := treeAssembly(topology, br_lengths)

# Log-normal distributed rate variation
rate_log_sd ~ dnGamma( 1, 1 )
moves.append( mvScale(rate_log_sd, weight = 1.5) )
rate_categories := fnDiscretizeDistribution( dnLognormal( 0, rate_log_sd ), 6)

# Neomorphic partition may have different rate to transformational
rate_neo ~ dnLognormal(mean = 0, sd = 2)
moves.append( mvScale(rate_neo, weight = 0.5) )
partition_rate := [ rate_neo / (1 + rate_neo), (1 / (1 + rate_neo)) ] / nChar * sum(nChar)

# Forward (gain) and reverse (loss) rates differ for neomorphic chars
rate_loss ~ dnLognormal(mean = 0, sd = 2)
moves.append( mvScale(rate_loss, lambda = 1, weight = 0.5) )

rate01 := 2 / (1 + rate_loss)
rate10 := 2 * rate_loss / (1 + rate_loss)
rates := [ [0.0, rate01],
           [rate10, 0.0] ]
neoQ := fnFreeK(rates)
stationaryDist := Simplex( rate10, rate01 )

m_morph[1] ~ dnPhyloCTMC(
  tree = phylogeny,
  branchRates = partition_rate[1],
  siteRates = rate_categories,
  Q = neoQ,
  rootFrequencies = stationaryDist,
  type = "Standard",
  coding = "variable"
)
m_morph[1].clamp(neo)

# Transformational partition: all chars together under fnJC(9)
transQ := fnJC(9)
m_morph[2] ~ dnPhyloCTMC(
  tree = phylogeny,
  branchRates = partition_rate[2],
  siteRates = rate_categories,
  Q = transQ,
  type = "Standard",
  coding = "variable"
)
m_morph[2].clamp(trans)

mymodel = model(phylogeny)

monitors = VectorMonitors()
monitors.append( mnScreen(printgen = 1000, prior = FALSE, posterior = FALSE, rate_neo, rate_loss) )
```

#### Driver — `long_by_nt_9v.Rev` (crashing version)

```rev
print("RB cell starting at", time("year"), time("day"), time("seconds") / 60 / 60)
seed(0)

neo   <- readDiscreteCharacterData("project950.neo.nex")
trans <- readDiscreteCharacterData("project950.trans.nex")

source("by_nt_9v.Rev")

monitors.append( mnModel(filename = "by_nt_9v.log",   printgen = 18) )
monitors.append( mnModel(filename = "by_nt_9v.p.log", printgen = 6, stochasticOnly = TRUE, exclude = ["rel_br_lengths"]) )
monitors.append( mnFile(filename  = "by_nt_9v.trees", printgen = 12, phylogeny) )

mymc3 = mcmcmc(mymodel, monitors, moves,
               nruns      = 2,
               nchains    = 4,
               tuneHeat   = TRUE,
               swapMethod = "both",
               swapMode   = "multiple",
               combine    = "none")

mymc3.burnin(generations = 200, tuningInterval = 25)

stopping_rules[1] = srMaxTime(30, "minutes")
stopping_rules[2] = srMinESS(128, "by_nt_9v.log", FALSE)   # <-- remove this line to fix

mymc3.run(generations       = 1000000,
          rules             = stopping_rules,
          checkpointFile    = "by_nt_9v.ckp",
          checkpointInterval = 1000)

q()
```

#### Invocation

```bash
module load gcc/13.2 boost
/path/to/revbayes/projects/cmake/build/rb long_by_nt_9v.Rev
```

### Observed behaviour

The script reads the two matrices, processes the model, runs the 200-iter burn-in to completion, prints the post-burn-in `Iter 0` row, and then the process dies with `SIGFPE`:

```
   Running MCMC simulation
   This simulation runs 2 independent replicates.
   The MCMCMC simulator runs 1 cold chain and 3 heated chains.
   The simulator uses 8 different moves in a random move schedule with 24.675 moves per iteration

   Stopping rules:
       Maximum allowed time: 00:30:00
       Target value for effective sample size (ESS): > 128


Iter        |     Likelihood   |      rate_loss   |       rate_neo   |    elapsed   |        ETA   |
----------------------------------------------------------------------------------------------------
0           |       -85.0238   |       0.120873   |        0.14511   |   00:00:00   |   --:--:--   |
/var/spool/slurmd/job17301940/slurm_script: line 18: 2260074 Floating point exception ... rb long_by_nt_9v.Rev
```

No core dump is produced; SLURM reports clean termination with the SIGFPE signal.

#### SLURM jobs

| Job ID    | Binary               | Dataset       | Tips × chars        | Stopping rules           | State     | Elapsed   | ExitCode |
|-----------|----------------------|---------------|---------------------|--------------------------|-----------|-----------|----------|
| **17301940** | **1.4.0-dev @ 8f71e5f (rebuilt today)** | pid 950 (Hexacorallia) | 12 × 8 (6 trans + 2 neo) | `srMaxTime` + `srMinESS` | **COMPLETED with SIGFPE** | 00:13 | **136**  |
| 17296889  | 1.3.2 (stale binary from older source) | pid 950 (Hexacorallia) | 12 × 8 (6 trans + 2 neo) | `srMaxTime` + `srMinESS` | FAILED | 22:21 | 8:0 |
| 17297082  | 1.3.2 (stale)        | pid 635 (Podalyria)    | 19 × 17 (12 trans + 5 neo) | `srMaxTime` + `srMinESS` | FAILED | 00:25 | 8:0 |
| 17297095  | 1.3.2 (stale)        | pid 635 (Podalyria) — same cell with **only `srMinESS` removed** | 19 × 17 | `srMaxTime` only | COMPLETED | 30:26 | 0:0 |

Job 17301940 is the **definitive reproducer against current HEAD** (`8f71e5f`): the serial `rb` was clean-rebuilt today from that commit, and the run dies with SIGFPE 13 seconds in, immediately after printing the post-burn-in `Iter 0` row. The slurm wrapper captured exit 136 (= 128 + signal 8).

Job 17297095 is the controlled comparator (against the older 1.3.2 binary): identical model body, identical priors, identical data, identical `mcmcmc` configuration. The single edit of removing the `srMinESS` line from the stopping-rules list lets the run reach the 30-minute time cap normally. We did not separately re-run the comparator against the fresh `1.4.0-dev` binary, but the same cell sans `srMinESS` was run against the fresh binary as part of an unrelated MkPrime-equivalence test on 2026-05-27 16:04 and ran to completion without incident.

### Expected behaviour

`srMinESS` should either

1. silently skip its ESS check until enough post-burn-in samples have been written for ESS to be defined, **or**
2. raise a Rev-level error message ("not enough samples in `<log>.log` to compute ESS") and let the user decide.

`SIGFPE` is not an acceptable outcome under any input combination.

### Hypothesis

The `Iter 0` row is printed by `mnScreen` from inside the MCMC loop *before* any periodic stopping-rule check fires; the very first stopping-rule pass then calls `MinEssStoppingRule::getStatistic`, which reads `by_nt_9v.log`. At that moment the log file contains only the header line written by `mnModel` (the first data row is not flushed until `printgen = 18` is reached). `coda`-style ESS computation on a zero-row or one-row trace divides by zero (or by `var = 0`), raising `SIGFPE`.

The MPI-flavoured crash reported in [#1023](https://github.com/revbayes/revbayes/issues/1023) has an almost identical signature — backtrace through `MinEssStoppingRule::getStatistic -> TraceContinuousReader` at `Iter 0` — but on the MPI build the read of an empty/half-written log file produces a `SIGSEGV` rather than a `SIGFPE`. The two are almost certainly the same defect surfacing differently depending on what the empty trace lands on inside the ESS code path.

### Workaround

Drop `srMinESS` from the stopping-rules list. Use `srMaxTime` (or a generation cap) only, and compute ESS post-hoc from the log files with `coda::effectiveSize` / `tracerer` / Tracer. This is what we have done in our pipeline.

### Searched prior reports

- [#1023](https://github.com/revbayes/revbayes/issues/1023) — closed, labelled `MPI`. Same call-site (`MinEssStoppingRule::getStatistic` at iter 0). Closed without a code fix landing on master as far as we can see. This report documents that the same defect occurs in the **non-MPI** binary, so the MPI label is too narrow.
- [#846](https://github.com/revbayes/revbayes/issues/846) (`[Feature Request] srMinESS updates`), [#848](https://github.com/revbayes/revbayes/issues/848) (missing docs) — related but not duplicates.
- [#642](https://github.com/revbayes/revbayes/issues/642), [#543](https://github.com/revbayes/revbayes/issues/543), [#422](https://github.com/revbayes/revbayes/issues/422) — other stopping-rule bugs, not the same failure mode.
