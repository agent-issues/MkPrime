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
    is.call(e) && as.character(e[[1L]]) %in% c("<-", "=") &&
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
ok("an arm's seeds do not depend on the other arms",
   identical(ArmSeedOffset("Mkp_logseries"),
             vapply(rev(arms), ArmSeedOffset, integer(1L))[["Mkp_logseries"]]))
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

cat("\nALL CHECKS PASSED\n")
