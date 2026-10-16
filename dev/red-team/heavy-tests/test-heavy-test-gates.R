#!/usr/bin/env Rscript
# Checks that the heavy-test gates can fail, and that their verdicts can be
# regenerated from committed config (#291, #292). Everything here runs on
# fabricated fixtures; nothing needs Hamilton.
#
# Run from the repository root:
#   Rscript dev/red-team/heavy-tests/test-heavy-test-gates.R
#
# With MkPrime installed (or loadable by pkgload), the final section also runs a
# few quick sbc.R simulations; set MKP_SKIP_SBC_SMOKE=1 to skip it.

HT <- "dev/red-team/heavy-tests"
MK <- file.path(HT, "marginal-k")
stopifnot(dir.exists(MK))

ok <- function(label, cond) {
  cat(sprintf("%-66s %s\n", label, if (isTRUE(cond)) "PASS" else "*** FAIL"))
  if (!isTRUE(cond)) stop(label)
}

# Evaluates one top-level definition from a driver without running the driver,
# whose first act is to load and compile the package.
Definition <- function(path, name, env = new.env(parent = globalenv())) {
  exprs <- parse(path, keep.source = FALSE)
  hit <- Filter(function(e) {
    is.call(e) && (identical(e[[1L]], quote(`<-`)) ||
                     identical(e[[1L]], quote(`=`))) &&
      identical(as.character(e[[2L]]), name)
  }, as.list(exprs))
  if (length(hit) != 1L) stop(name, " is not defined exactly once in ", path)
  eval(hit[[1L]], env)
  get(name, env)
}

# Every call to `fn` anywhere in a file.
CallsTo <- function(path, fn) {
  Walk <- function(e) {
    if (!is.call(e)) return(list())
    args <- Filter(function(a) !identical(a, quote(expr = )), as.list(e)[-1L])
    inner <- unlist(lapply(args, Walk), recursive = FALSE)
    if (identical(e[[1L]], as.name(fn))) c(list(e), inner) else inner
  }
  unlist(lapply(as.list(parse(path, keep.source = FALSE)), Walk),
         recursive = FALSE)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

## ---- HARN2-008: priorVariant pinned in every marginal-k driver -------------
for (f in list.files(MK, "\\.R$", full.names = TRUE)) {
  calls <- CallsTo(f, "MkPrimeModel")
  if (!length(calls)) next
  ok(sprintf("%s pins priorVariant", basename(f)),
     all(vapply(calls, function(cl) "priorVariant" %in% names(cl), logical(1L))))
}
ovlDriver <- file.path(MK, "T-OVL-sampled-vs-marginal.R")
WithPriorVariant <- function(value, code) {
  old <- Sys.getenv("MARGINAL_K_PRIOR_VARIANT", unset = NA)
  on.exit(if (is.na(old)) Sys.unsetenv("MARGINAL_K_PRIOR_VARIANT")
          else Sys.setenv(MARGINAL_K_PRIOR_VARIANT = old))
  if (is.na(value)) Sys.unsetenv("MARGINAL_K_PRIOR_VARIANT")
  else Sys.setenv(MARGINAL_K_PRIOR_VARIANT = value)
  code
}
ok("T-OVL defaults to the shipped Model A",
   identical(WithPriorVariant(NA, Definition(ovlDriver, "PRIOR_VARIANT")), "unconditional"))
ok("T-OVL prior variant is overridable",
   identical(WithPriorVariant("conditional", Definition(ovlDriver, "PRIOR_VARIANT")),
             "conditional"))
ok("the Phase-2 overnight script reproduces its recorded Model B runs",
   any(grepl("^export MARGINAL_K_PRIOR_VARIANT=conditional$",
             readLines(file.path(MK, "run-phase2-overnight.sh")))))

## ---- HARN2-012: marginal SBC output dir is per batch ----------------------
sbcDriver <- file.path(MK, "T-SBC-marginal-geometric.R")
override <- file.path(tempdir(), "batch-2")
Sys.setenv(MARGINAL_K_SBC_OUTDIR = override)
ok("T-SBC-marginal honours MARGINAL_K_SBC_OUTDIR",
   identical(Definition(sbcDriver, "OUT_DIR"), override))
Sys.unsetenv("MARGINAL_K_SBC_OUTDIR")
ok("T-SBC-marginal default output dir is unchanged",
   identical(Definition(sbcDriver, "OUT_DIR"),
             "dev/red-team/heavy-tests/marginal-k/sbc-results"))

BatchMap <- function(path) {
  grep("^\\s*[0-9]+\\) SEEDBASE=", readLines(path), value = TRUE)
}
arrMap <- BatchMap(file.path(HT, "submit-marginal-k-sbc.sh"))
aggMap <- BatchMap(file.path(HT, "submit-marginal-k-agg.sh"))
outSub <- sub(".*OUTSUB=([^ ;]+).*", "\\1", arrMap)
ok("marginal SBC array maps each batch to its own dir",
   length(arrMap) >= 2L && !anyDuplicated(outSub))
ok("marginal SBC aggregate reads the batch the array wrote",
   identical(arrMap, aggMap))
ok("marginal SBC array exports the dir to the driver",
   any(grepl("^export MARGINAL_K_SBC_OUTDIR=",
             readLines(file.path(HT, "submit-marginal-k-sbc.sh")))))

## ---- HARN2-010: SBC arms seeded by name, not position ----------------------
sbcR <- file.path(HT, "sbc.R")
ArmSeedOffset <- Definition(sbcR, "ArmSeedOffset")
arms <- vapply(Definition(sbcR, "ALL_ARMS"), `[[`, "", "name")
offsets <- vapply(arms, ArmSeedOffset, integer(1L))
ok("each SBC arm has its own seed block",
   !anyDuplicated(offsets) && min(diff(sort(offsets))) >= 1000L)
# Pinned, so an edit to the hash that would silently re-seed an arm fails here.
ok("arm seed offsets are unchanged",
   identical(unname(offsets), c(10935000L, 9647000L, 29083000L, 9765000L,
                                11067000L)))
ok("seeds stay within integer range",
   all(20260528 + offsets + 1000 < .Machine$integer.max))
seedLine <- grep("--seed", readLines(file.path(HT, "submit-sbc.sh")), value = TRUE)
ok("submit-sbc.sh passes one seedBase to every task",
   length(seedLine) && !any(grepl("SLURM_ARRAY_TASK_ID", seedLine)))

## ---- HARN2-013: submit-marginal-k-ovl.sh reads the verdict the driver writes
tagLines <- grep("^TAG=", readLines(file.path(HT, "submit-marginal-k-ovl.sh")),
                 value = TRUE)
ShellTag <- function(extra) {
  script <- paste(c(tagLines, 'printf %s "$TAG"'), collapse = "\n")
  env <- if (is.na(extra)) character(0L) else sprintf("MARGINAL_K_OVL_EXTRA=%s", extra)
  system2("env", c(env, "bash", "-c", shQuote(script)), stdout = TRUE)
}
ok("OVL submit computes the gated tag", identical(ShellTag(NA), "gated"))
ok("OVL submit computes a move-schedule tag like the driver",
   identical(ShellTag("weightedSpr,blockGibbsBranch"),
             paste(c("weightedSpr", "blockGibbsBranch"), collapse = "+")))
ok("OVL submit cats that tag's verdict file",
   any(grepl("T-OVL-${TAG}-verdict.txt",
             readLines(file.path(HT, "submit-marginal-k-ovl.sh")), fixed = TRUE)))

## ---- HARN2-014: T-OVL needs a floor of conclusive parameters ---------------
OvlVerdict <- Definition(ovlDriver, "OvlVerdict", {
  e <- new.env(); e$MIN_CONCLUSIVE <- 0.5; e
})
ok("one PASS among 47 INCONCLUSIVE is not a PASS",
   identical(OvlVerdict(c("PASS", rep("INCONCLUSIVE", 47L))), "INCONCLUSIVE"))
ok("half conclusive and passing is a PASS",
   identical(OvlVerdict(rep(c("PASS", "INCONCLUSIVE"), 24L)), "PASS"))
ok("any FAIL is a FAIL",
   identical(OvlVerdict(c(rep("PASS", 47L), "FAIL")), "FAIL"))
ok("nothing conclusive is INCONCLUSIVE",
   identical(OvlVerdict(rep("INCONCLUSIVE", 3L)), "INCONCLUSIVE"))

## ---- HARN2-016: pool-sampled-batches.R fails an NA Anderson–Darling --------
poolScript <- normalizePath(file.path(MK, "pool-sampled-batches.R"))
FakeBatch <- function(dir, n, rls, L = 100L) {
  dir.create(dir, recursive = TRUE)
  ranks <- round(seq(1, L - 1L, length.out = n))
  sims <- lapply(seq_len(n), function(i) {
    r <- ranks[(i * 7L) %% n + 1L]
    list(skipped = FALSE, L = L,
         p_true = if (r >= 30 && r <= 70 && i %% 3L == 0L) 0.05 else 0.5,
         ranks = list(tree_length = ranks[i], rate_log_sd = rls(i, ranks),
                      p = r))
  })
  saveRDS(sims, file.path(dir, "sims.rds"))
}
RunPool <- function(rls) {
  wd <- tempfile("pool"); dir.create(file.path(wd, MK), recursive = TRUE)
  FakeBatch(file.path(wd, "b1"), 60L, rls)
  FakeBatch(file.path(wd, "b2"), 60L, rls)
  old <- setwd(wd); on.exit(setwd(old))
  out <- suppressWarnings(system2("Rscript", c(shQuote(poolScript), "b1", "b2"),
                                  stdout = TRUE, stderr = TRUE))
  list(status = attr(out, "status") %||% 0L, text = paste(out, collapse = "\n"))
}
finite <- RunPool(function(i, ranks) ranks[(i * 11L) %% length(ranks) + 1L])
ok("pooling fixture passes when every AD is finite",
   finite$status == 0L && grepl("VERDICT: PASS", finite$text))
naRls <- RunPool(function(i, ranks) NA_real_)
ok("an NA AD fails clause (a)",
   naRls$status != 0L && grepl("\\(a\\) global AD\\s+: FAIL", naRls$text))

## ---- HARN2-009: sbc.R simulates from the drawn k and ranks the k' prior ----
canLoad <- requireNamespace("MkPrime", quietly = TRUE) ||
  requireNamespace("pkgload", quietly = TRUE)
if (!nzchar(Sys.getenv("MKP_SKIP_SBC_SMOKE")) && canLoad) {
  smokeOut <- tempfile("sbc-smoke")
  for (arm in c("Mkp_geometric", "MkNT_geometric")) {
    status <- system2("Rscript", c(shQuote(sbcR), "--quick", "--arm", arm,
                                   "--out", shQuote(smokeOut), "--run-id", "smoke"),
                      stdout = FALSE, stderr = FALSE)
    ok(sprintf("sbc.R --quick runs %s", arm), identical(status, 0L))
  }
  Good <- function(arm) {
    sims <- readRDS(file.path(smokeOut, arm, "summary.rds"))$sims
    Filter(function(x) !isTRUE(x$skipped), sims)
  }
  mkp <- Good("Mkp_geometric")
  ok("Mk' characters are simulated at their drawn k, not all at 2",
     length(mkp) >= 3L && any(unlist(lapply(mkp, `[[`, "kObs")) > 2L))
  ok("Mk' k' draws respect the truncation cap",
     all(unlist(lapply(mkp, `[[`, "kTrue")) <= Definition(sbcR, "K_MAX_PRIOR")))
  ok("Mk' ranks p and k'",
     all(vapply(mkp, function(x) is.finite(x$ranks$p) &&
                  is.finite(x$ranks$kPrime_sum), logical(1L))))
  mknt <- Good("MkNT_geometric")
  ok("MkNT characters are simulated at their drawn k",
     length(mknt) >= 3L && any(unlist(lapply(mknt, `[[`, "kObs")) > 2L))
  ok("the verdict is stamped with the run id",
     any(readLines(file.path(smokeOut, "verdict-MkNT_geometric.txt")) ==
           "run_id:   smoke"))
} else {
  cat("(sbc.R smoke skipped: MkPrime not loadable or MKP_SKIP_SBC_SMOKE set)\n")
}

## ---- HARN3-003/004/005: sbc-mixed simulates the package's MkN and gates hard
mixedR <- file.path(HT, "sbc-mixed.R")
SimNeo <- Definition(mixedR, ".simNeoChar")
# P(tip1 = 0, tip2 = 1) on a two-tip tree, root at stationarity.
JointPattern <- function(r, t, lambda) {
  pi <- c(r / (1 + r), 1 / (1 + r))
  Pr <- function(a, b) {
    move <- pi[b + 1L] * (1 - exp(-lambda * t))
    if (a == b) 1 - move else move
  }
  sum(vapply(0:1, function(rt) pi[rt + 1L] * Pr(rt, 0L) * Pr(rt, 1L), 0))
}
# Package rates (src/mkn_rates.h): rate01 = (1+r)/(2r), rate10 = (1+r)/2.
PackageP01 <- function(r, t) {
  Q <- matrix(c(-(1 + r) / (2 * r), (1 + r) / (2 * r),
                (1 + r) / 2, -(1 + r) / 2), 2L, byrow = TRUE)
  e <- eigen(Q)
  (e$vectors %*% diag(exp(e$values * t)) %*% solve(e$vectors))[1L, 2L]
}
rTrue <- 4
tBranch <- 0.3
lambdaPkg <- (1 + rTrue)^2 / (2 * rTrue)
ok("lambda = (1+r)^2/(2r) reproduces the package P01",
   abs(PackageP01(rTrue, tBranch) -
         (1 / (1 + rTrue)) * (1 - exp(-lambdaPkg * tBranch))) < 1e-9)
set.seed(1)
twoTip <- ape::read.tree(text = sprintf("(t1:%g,t2:%g);", tBranch, tBranch))
draws <- replicate(40000L, SimNeo(twoTip, rTrue, 1))
pHat <- mean(draws[1L, ] == 0L & draws[2L, ] == 1L)
ok("sbc-mixed MkN forward sim matches the package rate",
   abs(pHat - JointPattern(rTrue, tBranch, lambdaPkg)) < 0.008)
ok("... and would have rejected the pre-eecee51 lambda = 2",
   abs(pHat - JointPattern(rTrue, tBranch, 2)) > 0.015)

MixedVerdict <- Definition(mixedR, "MixedVerdict")
all3 <- c("tree_length", "rate_neo", "rate_loss")
FakeSims <- function(n, params = all3, nSkip = 0L, L = 100L) {
  lapply(seq_len(n), function(i) {
    if (i <= nSkip) return(list(skipped = TRUE, reason = "x"))
    list(skipped = FALSE, L = L,
         ranks = stats::setNames(rep(list((i * 37L) %% L), length(params)),
                                 params))
  })
}
Gate <- function(sims) {
  MixedVerdict(sims, all3, 100L, 0.001, 0.9, 50L,
               function(x) suppressWarnings(stats::ks.test(x, "punif"))$p.value)$verdict
}
ok("sbc-mixed passes a complete, uniform fixture", identical(Gate(FakeSims(200L)), "PASS"))
ok("sbc-mixed fails when an expected parameter is missing",
   identical(Gate(FakeSims(200L, c("tree_length", "rate_neo"))), "FAIL"))
ok("sbc-mixed fails when ranks are INSUFFICIENT",
   identical(Gate(FakeSims(40L)), "FAIL"))
ok("sbc-mixed fails when > 10% of sims are skipped",
   identical(Gate(FakeSims(200L, nSkip = 30L)), "FAIL"))
ok("sbc-mixed tolerates a few skips",
   identical(Gate(FakeSims(200L, nSkip = 10L)), "PASS"))
ok("sbc-mixed fails on a non-uniform rank histogram",
   identical(Gate(lapply(FakeSims(200L), function(s) {
     s$ranks$rate_loss <- 2L
     s
   })), "FAIL"))
mixedLines <- readLines(mixedR)
ok("sbc-mixed ranks no inert rate_log_sd",
   !any(grepl("ranks$rate_log_sd", mixedLines, fixed = TRUE)))
ok("sbc-mixed doc and Bonferroni constant agree",
   any(grepl("FWER 0.001 / K_arm", mixedLines, fixed = TRUE)) &&
     identical(Definition(mixedR, "FWER"), 0.001))

## ---- HARN3-008/010/013: other gates ----------------------------------------
spr <- readLines(file.path(HT, "gibbs-spr-db.R"))
ok("gibbs-spr-db FAIL also covers the composite br_* arms (q2Fail)",
   any(grepl("q1Fail || q2Fail", spr, fixed = TRUE)))
ok("gibbs-spr-db exits non-zero on WARN and INCONCLUSIVE",
   any(grepl("WARN = 2L, INCONCLUSIVE = 3L", spr, fixed = TRUE)))
for (f in c("gibbs-spr-db.R", "marginal-k/T-SBC-marginal-geometric.R",
            "marginal-k/T-SBC-sampled-geometric.R",
            "marginal-k/T-OVL-sampled-vs-marginal.R")) {
  src <- readLines(file.path(HT, f))
  ok(sprintf("%s loads the package through load-mkprime.R", basename(f)),
     !any(grepl("pkgload::load_all", src)) &&
       any(grepl("LoadMkPrime(", src, fixed = TRUE)))
}
for (f in c("T-SBC-marginal-geometric.R", "T-SBC-sampled-geometric.R")) {
  ok(sprintf("%s fails when > 10%% of sims are missing", f),
     any(grepl("ceiling(0.9 * N_SIM)", readLines(file.path(MK, f)),
               fixed = TRUE)))
}
for (f in c("submit-marginal-k-sbc.sh", "submit-marginal-k-sbc-sampled.sh")) {
  tl <- sub(".*--time=", "",
            grep("--time=", readLines(file.path(HT, f)), value = TRUE))
  mins <- sum(as.integer(strsplit(tl, ":")[[1L]]) * c(60L, 1L, 0L))
  ok(sprintf("%s walltime covers 5 sims at ~3 min plus start-up", f),
     mins >= 30L)
}
funnel <- readLines(file.path(HT, "funnel-stress-hyperprior-sigma.R"))
ok("funnel-stress pins nCat, expSteps and priorVariant",
   any(grepl("MODEL_PIN <- list(nCat", funnel, fixed = TRUE)) &&
     any(grepl("priorVariant = MODEL_PIN", funnel, fixed = TRUE)))
ok("funnel-stress prints no NEWS.md headline", !any(grepl("NEWS", funnel)))
ok("subtree-swap-db is marked retired",
   any(grepl("^# RETIRED", readLines(file.path(HT, "subtree-swap-db.R"))[1:6])))

cat("\nALL CHECKS PASSED\n")
