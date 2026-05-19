#!/usr/bin/env python3
"""Add mk_k40 arm to run_one.R, summarize_streamed.R, and summarize_array.slurm."""

import re

# ------------------------------------------------------------------
# 1. run_one.R  — match.arg
# ------------------------------------------------------------------
path = "/nobackup/pjjg18/mkp-study/run_one.R"
with open(path) as f:
    txt = f.read()

old = '"mk_k15", "mk_k24",'
new = '"mk_k15", "mk_k24", "mk_k40",'
if old not in txt:
    print(f"WARN: pattern not found in {path}")
else:
    txt = txt.replace(old, new, 1)
    print(f"run_one.R  match.arg: updated")

# Insert mk_k40 arm block just before the final cat("  Done.\n")
k40_block = r"""} else if (arm == "mk_k40") {
  # Mk with knownStates = 40 across all variable characters. Extended
  # endpoint; tests whether the k-ramp continues past k=24.
  kobs_raw <- apply(combined_mat, 2L, function(col) {
    length(unique(col[!col %in% c("?", "-")]))
  })
  var_orig <- which(kobs_raw > 1L)
  k40_for_mk <- setNames(rep(40L, length(var_orig)),
                          as.character(var_orig))

  res <- .run_arm(function() {
    RunMkPrime(
      pd, start_tree,
      knownStates = k40_for_mk,
      model = MkPrimeModel(coding = "variable"),
      mcmc  = make_mcmc("mk_k40", thin_iters = 100L)
    )
  }, "mk_k40")

  cat(sprintf("  Mk(40) done: %d trees, stop=%s\n",
              length(res$trees), res$stop_reason))

  partial <- list(trees = res$trees, stop_reason = res$stop_reason,
                  acceptance = res$acceptance)
  saveRDS(partial, file.path(out_dir, sprintf("mk_k40_%s.rds", tag)))

"""

marker = '\ncat("  Done.\\n")'
if marker not in txt:
    print("WARN: Done marker not found — mk_k40 block NOT inserted")
elif 'arm == "mk_k40"' in txt:
    print("run_one.R  mk_k40 block: already present")
else:
    txt = txt.replace(marker, "\n" + k40_block + marker, 1)
    print("run_one.R  mk_k40 block: inserted")

with open(path, "w") as f:
    f.write(txt)

# ------------------------------------------------------------------
# 2. summarize_streamed.R  — match.arg
# ------------------------------------------------------------------
path = "/nobackup/pjjg18/mkp-study/summarize_streamed.R"
with open(path) as f:
    txt = f.read()

old = '"mk_k15", "mk_k24")'
new = '"mk_k15", "mk_k24", "mk_k40")'
if old not in txt:
    print(f"WARN: pattern not found in {path}")
else:
    txt = txt.replace(old, new, 1)
    print("summarize_streamed.R  match.arg: updated")

with open(path, "w") as f:
    f.write(txt)

# ------------------------------------------------------------------
# 3. summarize_array.slurm  — case statement
# ------------------------------------------------------------------
path = "/nobackup/pjjg18/mkp-study/summarize_array.slurm"
with open(path) as f:
    txt = f.read()

old = "mk|mkp|mkp_eg|mk_kp1|mk_kp2|mk_k9|mk_k15|mk_k24)"
new = "mk|mkp|mkp_eg|mk_kp1|mk_kp2|mk_k9|mk_k15|mk_k24|mk_k40)"
if old not in txt:
    print(f"WARN: case pattern not found in {path} — trying alternate")
    # try without mk_k24 already there
    old2 = "mk_k15|mk_k24)"
    new2 = "mk_k15|mk_k24|mk_k40)"
    if old2 in txt:
        txt = txt.replace(old2, new2, 1)
        print("summarize_array.slurm  case: updated (alternate match)")
    else:
        print("WARN: could not update summarize_array.slurm — check manually")
else:
    txt = txt.replace(old, new, 1)
    print("summarize_array.slurm  case: updated")

with open(path, "w") as f:
    f.write(txt)

print("Done.")
