#!/usr/bin/env Rscript
## Unit checks for the helpers run_one.R uses to guard its inputs and outputs.
##
## Run from the repository root:
##   Rscript data-raw/hamilton/test_run_one_helpers.R
##
## The helpers are lifted out of run_one.R by name rather than re-stated here,
## so this cannot drift from what the harness actually runs.

suppressPackageStartupMessages(library(MkPrime))

src   <- parse("data-raw/hamilton/run_one.R")
wanted <- c(".kObsFromMatrix", ".UPostMeans", ".BURNIN_FRAC", ".arm_own_files",
            ".validate_ckp", ".IsCorruptCheckpoint", ".run_arm",
            ".CheckCharOrder", ".VariableCharIdx", ".SaveResult", ".BuildSha",
            ".BUILD_SHA", ".CHAR_ORDER_STAMP", ".CheckNotFinished",
            ".PinnedTreeLengthPrior", ".PinnedModel")
env   <- new.env(parent = globalenv())
## The helpers close over the script's ckp_dir.
env$ckp_dir <- tempfile("ckp"); dir.create(env$ckp_dir)
found <- character(0L)
for (e in src) {
  if (is.call(e) && identical(as.character(e[[1L]]), "<-") &&
      is.name(e[[2L]]) && as.character(e[[2L]]) %in% wanted) {
    eval(e, envir = env)
    found <- c(found, as.character(e[[2L]]))
  }
}
stopifnot(setequal(found, wanted))

ok <- function(label, cond) {
  cat(sprintf("%-58s %s\n", label, if (isTRUE(cond)) "PASS" else "*** FAIL"))
  if (!isTRUE(cond)) stop(label)
}

## ---- .kObsFromMatrix --------------------------------------------------------
mat <- cbind(
  plain      = c("0", "1", "0", "1"),
  with_gaps  = c("0", "1", "?", "-"),
  polymorph  = c("0", "1", "{01}", "?"),
  paren_poly = c("0", "1", "(01)", "-"),
  invariant  = c("1", "1", "1", "?")
)
k <- env$.kObsFromMatrix(mat)
ok("plain binary counts 2", k[["plain"]] == 2L)
ok("? and - are not states", k[["with_gaps"]] == 2L)
ok("{01} is not a third state", k[["polymorph"]] == 2L)
ok("(01) is not a third state", k[["paren_poly"]] == 2L)
ok("invariant character counts 1", k[["invariant"]] == 1L)

## The old expression, for contrast: it is the inflation #104 describes.
old <- apply(mat, 2L, function(col) length(unique(col[!col %in% c("?", "-")])))
ok("old expression did inflate the polymorphic columns",
   old[["polymorph"]] == 3L && old[["paren_poly"]] == 3L)

## ---- .arm_own_files ---------------------------------------------------------
## The shared per-(tree, rep) directory as every arm leaves it.
siblings <- c(
  "mk_checkpoint.rds", "mk_run_1.log", "mk_run_2.log",
  "mk_trees_1.nwk", "mk_trees_2.nwk", "mk_trees.nwk", "mk_run.log",
  "mk_run_1.log.gz", "mk_char_order.csv",
  "mk_k9_checkpoint.rds", "mk_k9_run_1.log.gz", "mk_k9_run_1.log", "mk_k9_trees_1.nwk",
  "mk_k15_run_1.log", "mk_k24_run_1.log", "mk_k40_trees_2.nwk",
  "mk_kp1_checkpoint.rds", "mk_kp2_run_1.log", "mk_ktrue_trees_1.nwk",
  "mk_tlshrink_run_2.log",
  "mkp_checkpoint.rds", "mkp_run_1.log",
  "mkp_eg_run_1.log", "mkp_geo_checkpoint.rds", "mkp_highk_trees_2.nwk",
  "mkp_logs_run_1.log",
  ".slurm_job_id"
)
invisible(file.create(file.path(env$ckp_dir, siblings)))

own_mk  <- basename(env$.arm_own_files("mk"))
own_mkp <- basename(env$.arm_own_files("mkp"))

ok("mk claims its own eight files, gzipped log included", setequal(own_mk, c(
  "mk_checkpoint.rds", "mk_run_1.log", "mk_run_2.log", "mk_run_1.log.gz",
  "mk_trees_1.nwk", "mk_trees_2.nwk", "mk_trees.nwk", "mk_run.log")))
ok("mk claims nothing of mk_k9 / mk_kp1 / mk_ktrue / mk_tlshrink",
   !any(grepl("^mk_(k[0-9]|kp[0-9]|ktrue|tlshrink)", own_mk)))
ok("per-run checkpoints are the arm's own",
   identical(basename(env$.arm_own_files("mkp_geo")), "mkp_geo_checkpoint.rds") &&
     {
       file.create(file.path(env$ckp_dir, "mkp_geo_checkpoint_2.rds"))
       "mkp_geo_checkpoint_2.rds" %in% basename(env$.arm_own_files("mkp_geo"))
     })
ok("mkp claims its own two files",
   setequal(own_mkp, c("mkp_checkpoint.rds", "mkp_run_1.log")))
ok("mkp claims nothing of mkp_eg / mkp_geo / mkp_highk / mkp_logs",
   !any(grepl("^mkp_(eg|geo|highk|logs)", own_mkp)))
ok("the job sentinel is never purged",
   !(".slurm_job_id" %in% c(own_mk, own_mkp)))

## The old prefix glob, for contrast: this is the blast radius in #98.
old_mk <- grep("^mk_", siblings, value = TRUE)
ok("old prefix glob did reach eleven foreign files",
   sum(!old_mk %in% own_mk & old_mk != "mk_char_order.csv") == 11L)

## ---- .UPostMeans ------------------------------------------------------------
## Two runs whose k' means differ sharply, so reading one run or skipping the
## burn-in both give visibly different answers from reading both correctly.
nChar <- 3L
mkd   <- list(nChar = nChar, kObs = rep(2L, nChar))

write_log <- function(path, n, kp_value) {
  hdr <- paste(c("Sample", "tree_length", paste0("kPrime_", seq_len(nChar))),
               collapse = "\t")
  rows <- vapply(seq_len(n), function(i) {
    ## First half of each run sits at 99 so a missing burn-in cut is obvious.
    kp <- if (i <= n / 2) 99 else kp_value
    paste(c(i, 1.0, rep(kp, nChar)), collapse = "\t")
  }, character(1L))
  writeLines(c(hdr, rows), path)
}

dir <- tempfile("upost"); dir.create(dir)
f1 <- file.path(dir, "a_run_1.log")
f2 <- file.path(dir, "a_run_2.log")
write_log(f1, 100L, 10)
write_log(f2, 100L, 20)

res <- list(logFile = c(f1, f2), samples = NULL)
up  <- env$.UPostMeans(res, mkd, burninFrac = 0.5)

ok("both runs read", up$nRuns == 2L)
ok("post-burn-in row count is 50 + 50", up$n == 100L)
ok("burn-in excludes the 99 block", all(abs(up$means - (15 - 2)) < 1e-9))
ok("burninFrac is reported back", identical(up$burninFrac, 0.5))

## Run 1 alone, no burn-in -- what the old code computed -- differs.
old_means <- colMeans(ReadMkLog(f1)[, paste0("kPrime_", seq_len(nChar)),
                                    drop = FALSE]) - mkd$kObs
ok("old run-1-only, no-burn-in estimate differs",
   all(abs(old_means - up$means) > 1))

## Degrades to res$samples when no log file is on disk.
res_nolog <- list(logFile = character(0L),
                  samples = ReadMkLog(f2))
ok("falls back to res$samples", env$.UPostMeans(res_nolog, mkd)$nRuns == 1L)

## No k' columns at all (the fixed-k arms): NA of the right length.
res_nokp <- list(logFile = character(0L),
                 samples = matrix(1, nrow = 4L, ncol = 1L,
                                  dimnames = list(NULL, "tree_length")))
ok("no kPrime_ columns gives NA per character",
   identical(env$.UPostMeans(res_nokp, mkd)$means, rep(NA_real_, nChar)))

ok("default burn-in fraction is 25%", identical(env$.BURNIN_FRAC, 0.25))

## ---- .IsCorruptCheckpoint / .run_arm ----------------------------------------
env$ckp_dir <- tempfile("purge"); dir.create(env$ckp_dir)
validCkp <- list(runs = list(1), mcmc = list(nRuns = 1L), iter = 10L)
ArmFiles <- function() basename(env$.arm_own_files("mkp"))
Populate <- function() {
  saveRDS(validCkp, file.path(env$ckp_dir, "mkp_checkpoint.rds"))
  writeLines("Sample", file.path(env$ckp_dir, "mkp_run_1.log"))
}
FailOnce <- function(msg) {
  calls <- new.env()
  calls$n <- 0L
  function() {
    calls$n <- calls$n + 1L
    if (calls$n == 1L) stop(msg) else "resumed"
  }
}

## Real messages from R/RunMkPrime.R that name the checkpoint but are not
## corruption: the old regex purged on every one of them (#289).
for (msg in c("Unsupported checkpoint version 3.",
              "Checkpoint holds state for 1 of 2 runs.",
              "Resumed move set differs from the one this checkpoint was written with.",
              "Checkpointed move counts do not match the rebuilt move schedule.")) {
  Populate()
  res <- tryCatch(env$.run_arm(FailOnce(msg), "mkp"), error = identity)
  ok(paste0("not purged: ", substr(msg, 1L, 40L)),
     inherits(res, "error") &&
       setequal(ArmFiles(), c("mkp_checkpoint.rds", "mkp_run_1.log")))
}
ok("old regex would have purged those",
   grepl("reading from connection|checkpoint",
         "Unsupported checkpoint version 3.", ignore.case = TRUE))

Populate()
ok("a connection read error purges and retries",
   identical(env$.run_arm(FailOnce("error reading from connection"), "mkp"),
             "resumed") && !length(ArmFiles()))

Populate()
writeLines("not an rds", file.path(env$ckp_dir, "mkp_checkpoint.rds"))
ok("a checkpoint failing .validate_ckp purges whatever the message",
   identical(env$.run_arm(FailOnce("anything at all"), "mkp"), "resumed") &&
     !length(ArmFiles()))

## ---- .CheckCharOrder --------------------------------------------------------
env$ckp_dir <- tempfile("order"); dir.create(env$ckp_dir)
numeric <- file.path("rep", sprintf("chr%d.nex", 1:11))
lexical <- sort(numeric)
orderCsv <- file.path(env$ckp_dir, "mkp_eg_char_order.csv")
ckpFile  <- file.path(env$ckp_dir, "mkp_eg_checkpoint.rds")

env$.CheckCharOrder("mkp_eg", numeric)
rec <- read.csv(orderCsv)
ok("first start records the order", identical(rec$file, basename(numeric)) &&
     identical(rec$char_idx, 1:11))

saveRDS(validCkp, ckpFile)
ok("resume under the recorded order proceeds",
   !inherits(tryCatch(env$.CheckCharOrder("mkp_eg", numeric), error = identity),
             "error"))

refusal <- tryCatch(env$.CheckCharOrder("mkp_eg", lexical), error = identity)
ok("resume under another order is refused", inherits(refusal, "error"))
ok("refusal leaves the saved state in place", file.exists(ckpFile) &&
     identical(read.csv(orderCsv)$file, basename(numeric)))
ok("refusal names this arm's files, never a prefix glob",
   grepl("mkp_eg_checkpoint.rds", conditionMessage(refusal)) &&
     !grepl("mkp_eg_\\*", conditionMessage(refusal)))
ok("refusal names both orders",
   grepl("chr1.nex chr2.nex chr3.nex", conditionMessage(refusal)) &&
     grepl("chr1.nex chr10.nex chr11.nex", conditionMessage(refusal)))
ok("refusal would not trigger the purge",
   !env$.IsCorruptCheckpoint(refusal, "mkp_eg"))

invisible(file.remove(orderCsv))
ok("saved state with no recorded order is refused",
   inherits(tryCatch(env$.CheckCharOrder("mkp_eg", numeric), error = identity),
            "error") && !file.exists(orderCsv))

invisible(file.remove(ckpFile))
saveRDS(validCkp, file.path(env$ckp_dir, "mkp_eg_checkpoint_2.rds"))
ok("a per-run checkpoint alone counts as saved state",
   inherits(tryCatch(env$.CheckCharOrder("mkp_eg", numeric), error = identity),
            "error"))
invisible(file.remove(file.path(env$ckp_dir, "mkp_eg_checkpoint_2.rds")))
env$.CheckCharOrder("mkp_eg", lexical)
ok("with no saved state, the new order is recorded",
   identical(read.csv(orderCsv)$file, basename(lexical)))

## ---- Birth stamp on the order record (#388) --------------------------------
ok("a fresh start stamps the record",
   all(c("lifecycle", "build_sha") %in% names(rec)))

## A record rewritten by a pre-5a890df resume: right order, no stamp.
unstamped <- rec[, c("position", "file", "char_idx")]
write.csv(unstamped, orderCsv, row.names = FALSE)
saveRDS(validCkp, ckpFile)
Sys.unsetenv("MKP_ACCEPT_UNSTAMPED_ORDER")
refusal <- tryCatch(env$.CheckCharOrder("mkp_eg", numeric), error = identity)
ok("resume against an unstamped record is refused even in the same order",
   inherits(refusal, "error") && grepl("birth stamp", conditionMessage(refusal)))
ok("the refusal names check_kprime_order.R",
   grepl("check_kprime_order.R", conditionMessage(refusal)))
ok("refusal leaves the record and the checkpoint alone",
   file.exists(ckpFile) && !("lifecycle" %in% names(read.csv(orderCsv))))
Sys.setenv(MKP_ACCEPT_UNSTAMPED_ORDER = "1")
ok("MKP_ACCEPT_UNSTAMPED_ORDER=1 lets it through in the same order",
   !inherits(tryCatch(env$.CheckCharOrder("mkp_eg", numeric), error = identity),
             "error"))
ok("but never under another order",
   inherits(tryCatch(env$.CheckCharOrder("mkp_eg", lexical), error = identity),
            "error"))
Sys.unsetenv("MKP_ACCEPT_UNSTAMPED_ORDER")
invisible(file.remove(ckpFile))
env$.CheckCharOrder("mkp_eg", numeric)
ok("a fresh start over an unstamped record re-stamps it",
   "lifecycle" %in% names(read.csv(orderCsv)))

## ---- .CheckNotFinished (#388) ------------------------------------------------
fin <- tempfile("fin"); dir.create(fin)
ckpD <- file.path(fin, "t01_r01"); dir.create(ckpD)
Refused <- function(arm = "mkp")
  inherits(tryCatch(env$.CheckNotFinished(arm, "t01_r01", fin, ckpD),
                    error = identity), "error")
Sys.unsetenv("MKP_FORCE_RESTART")
ok("no final result: a fresh start is allowed", !Refused())
saveRDS(list(stop_reason = "max_time"), file.path(fin, "mkp_t01_r01.rds"))
ok("final result and no checkpoint: a fresh start is refused", Refused())
msg <- conditionMessage(tryCatch(env$.CheckNotFinished("mkp", "t01_r01", fin, ckpD),
                                 error = identity))
ok("the refusal names MKP_FORCE_RESTART and the files to delete",
   grepl("MKP_FORCE_RESTART=1", msg, fixed = TRUE) &&
     grepl("delete", msg, fixed = TRUE))
ok("another arm's finished result is not this arm's", !Refused("mkp_eg"))
invisible(file.create(file.path(ckpD, "mkp_checkpoint_2.rds")))
ok("final result with a checkpoint: it resumes", !Refused())
invisible(file.remove(file.path(ckpD, "mkp_checkpoint_2.rds")))
Sys.setenv(MKP_FORCE_RESTART = "1")
ok("MKP_FORCE_RESTART=1 permits the restart", !Refused())
Sys.unsetenv("MKP_FORCE_RESTART")

## ---- .VariableCharIdx -------------------------------------------------------
## 11 characters, character 5 invariant: the variable subset is not 1:10.
taxa <- paste0("t", 1:6)
mat  <- vapply(1:11, function(i) {
  if (i == 5L) rep("1", 6L) else as.character(c(0, 1, 0, 1, i %% 3, 0))
}, character(6L))
rownames(mat) <- taxa
tmpNex <- tempfile(fileext = ".nex")
ape::write.nexus.data(setNames(lapply(seq_along(taxa), function(i) mat[i, ]),
                               taxa),
                      file = tmpNex, format = "standard")
pd <- TreeTools::ReadAsPhyDat(tmpNex)
ok("invariant character dropped from char_idx",
   identical(env$.VariableCharIdx(pd, 1:11), c(1:4, 6:11)))
ok("result aligns with MkPrimeData()'s kObs",
   length(env$.VariableCharIdx(pd, 1:11)) ==
     suppressWarnings(MkPrimeData(pd))$nChar)
ok("wrong-length index is refused",
   inherits(tryCatch(env$.VariableCharIdx(pd, 1:10), error = identity), "error"))

## ---- .SaveResult ------------------------------------------------------------
env$.PIN <- list()
env$out_dir <- tempfile("out"); dir.create(env$out_dir)
env$arm <- "mkp_eg"
env$tag <- "t01_r01"
env$var_char_idx <- c(1:4, 6:11)
env$.SaveResult(list(stop_reason = "max_time"))
saved <- readRDS(file.path(env$out_dir, "mkp_eg_t01_r01.rds"))
ok("result records char_idx", identical(saved$char_idx, c(1:4, 6:11)))
ok("result records the build sha", identical(saved$build_sha, env$.BUILD_SHA))

## ---- .BuildSha -------------------------------------------------------------
FakeLib <- function(extra) {
  lib <- tempfile("lib")
  dir.create(file.path(lib, "MkPrime"), recursive = TRUE)
  writeLines(c("Package: MkPrime", "Version: 0.0.0.9000", extra),
             file.path(lib, "MkPrime", "DESCRIPTION"))
  lib
}
ok("a stamped build reports its sha",
   identical(env$.BuildSha(FakeLib("RemoteSha: 5a890df")), "5a890df"))
ok("an unstamped build reports NA",
   identical(env$.BuildSha(FakeLib(character(0L))), NA_character_))
ok("result keeps its own fields", identical(saved$stop_reason, "max_time"))

## ---- Pinned tree-length prior (#389) ----------------------------------------
tree <- TreeTools::NJTree(pd, edgeLengths = TRUE)
mkd  <- suppressWarnings(MkPrimeData(pd))
oldDefault <- max(1, MkPrime:::.FitchScore(tree, mkd))      # before 3c091f9
headDefault <- MkPrime:::.FinalizeModel(MkPrimeModel(coding = "variable"),
                                        NULL, mkd)$expSteps
ok("fixture discriminates the old default from HEAD's",
   !isTRUE(all.equal(oldDefault, headDefault)))
pin <- env$.PinnedTreeLengthPrior(pd, tree)
ok("the pin is the pre-3c091f9 expSteps", identical(pin$expSteps, oldDefault))
ok("the pin is not HEAD's default", !isTRUE(all.equal(pin$expSteps, headDefault)))
ok("the pin's rate keeps the prior mean at expSteps",
   isTRUE(all.equal(pin$treeLengthRate, 2 / oldDefault)))
env$.PIN <- pin
fin <- MkPrime:::.FinalizeModel(env$.PinnedModel(), NULL, mkd)
ok("an arm's resolved model carries the pin, not the new default",
   identical(fin$expSteps, oldDefault) &&
     isTRUE(all.equal(fin$treeLengthRate, 2 / oldDefault)))
fin <- MkPrime:::.FinalizeModel(
  env$.PinnedModel(treeLengthShape = 20, treeLengthRate = 20 / 0.7), NULL, mkd)
ok("an arm with its own tree-length prior keeps it",
   fin$treeLengthShape == 20 && isTRUE(all.equal(fin$treeLengthRate, 20 / 0.7)))
ok("no arm calls MkPrimeModel() unpinned",
   !any(grepl("MkPrimeModel(", readLines("data-raw/hamilton/run_one.R"),
              fixed = TRUE)))
env$arm <- "mkp"; env$var_char_idx <- 1:3
env$.SaveResult(list(stop_reason = "max_time"))
ok("result records the resolved tree-length prior",
   identical(readRDS(file.path(env$out_dir, "mkp_t01_r01.rds"))$tree_length_prior,
             pin))

cat("\nALL CHECKS PASSED\n")
