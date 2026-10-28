# Synthetic-fixture checks for arm labelling in cid_three_prior.R and
# redteam_tl_hypothesis.R.
# Run: Rscript dev/pilots/2026-05-12-prior-validation/analysis/test_arm_labels.R
library(data.table)

Check <- function(label, ok) {
  cat(sprintf("%s %s\n", if (isTRUE(ok)) "PASS" else "FAIL", label))
  if (!isTRUE(ok)) quit(status = 1L)
}

# Evaluate one top-level function definition from a script without running the script.
LoadFunction <- function(script, name, env) {
  for (e in parse(script)) {
    if (is.call(e) && identical(as.character(e[[1]]), "<-") && identical(as.character(e[[2]]), name)) {
      eval(e, env)
      return(invisible())
    }
  }
  stop("no ", name, " in ", script)
}

dir <- "dev/pilots/2026-05-12-prior-validation/analysis"

# --- cid_three_prior: the mk glob must not pool mk_k15, mk_kp1, ... ---
tmp <- tempfile(); dir.create(tmp)
Write <- function(arm, tag, cid) {
  saveRDS(list(tag = tag, tree_idx = 1L, rep_idx = 1L, cid = rep(cid, 4)),
          file.path(tmp, sprintf("%s_%s.rds", arm, tag)))
}
for (arm in c("mk", "mk_k15", "mk_kp1", "mk_tlshrink")) {
  Write(arm, "t01_r01", if (arm == "mk") 0.1 else 0.9)
}
src <- readLines(file.path(dir, "cid_three_prior.R"))
perTask <- new.env()
perTask$SUMMARY_DIR <- tmp
LoadFunction(file.path(dir, "cid_three_prior.R"), "per_task_cid", perTask)
callSite <- grep('per_task_cid\\(.*"mk"\\)', src, value = TRUE)
Check("cid_three_prior has exactly one mk call", length(callSite) == 1L)
mk <- eval(parse(text = sub("^[^<]*<-\\s*", "", callSite)), perTask)
Check("mk arm is one file, not nine pooled arms", nrow(mk) == 1L)
Check("mk arm is the mk arm's value, not mk_k15's", isTRUE(all.equal(mk$cid_mean, 0.1)))

# --- redteam_tl_hypothesis: mk's k comes from the data ---
e <- new.env()
LoadFunction(file.path(dir, "redteam_tl_hypothesis.R"), "AssignFixedK", e)
dt <- data.table(arm = c("mk", "mk", "mk_k9", "mkp_eg"), mean_kobs = c(3.5, 2.25, NA, NA))
armK <- c(mk = NA, mk_k9 = 9L, mkp_eg = NA)
out <- e$AssignFixedK(dt, armK)
Check("mk k follows each task's kObs", identical(out$k_fixed[1:2], c(3.5, 2.25)))
Check("fixed-k arms keep their k", out$k_fixed[3] == 9)
Check("mk with no kObs errors rather than anchoring at 2",
      inherits(try(e$AssignFixedK(data.table(arm = "mk", mean_kobs = NA_real_), armK), silent = TRUE),
               "try-error"))
