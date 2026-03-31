# M-131: Upload to Hamilton and submit SLURM array job
#
# Run interactively in RStudio console:
#   source("inst/hamilton/m131-upload-and-submit.R")
#
# Requires: ssh package, passphrase for id_ed25519

library(ssh)

# --- Connect (will prompt for passphrase) ---
session <- ssh::ssh_connect("pjjg18@hamilton8.dur.ac.uk",
                             keyfile = "C:/Users/pjjg18/id_ed25519")

BASEDIR <- "/nobackup/pjjg18/m131"

# --- Create directories ---
ssh::ssh_exec_wait(session,
  sprintf("mkdir -p %s/{lib,logs,results,data-raw}", BASEDIR))

# --- Upload MkPrime tarball ---
cat("Uploading MkPrime tarball...\n")
ssh::scp_upload(session, "MkPrime_0.0.0.9000.tar.gz", BASEDIR)

# --- Upload validation scripts ---
cat("Uploading scripts...\n")
for (f in list.files("inst/hamilton", full.names = TRUE)) {
  ssh::scp_upload(session, f, BASEDIR)
}

# --- Upload nexus data files ---
cat("Uploading datasets...\n")
datasets <- c("Longrich2010", "Vinther2008", "DeAssis2011", "Wortley2006",
               "Sun2018", "Wilson2003", "Zhu2013", "Dikow2009")
for (ds in datasets) {
  f <- file.path("../TreeSearch-a/data-raw", paste0(ds, ".nex"))
  ssh::scp_upload(session, f, file.path(BASEDIR, "data-raw"))
}

# --- Fix line endings ---
cat("Fixing line endings...\n")
ssh::ssh_exec_wait(session,
  sprintf("sed -i 's/\r$//' %s/*.R %s/*.slurm %s/*.sh", BASEDIR, BASEDIR, BASEDIR))

# --- Install dependencies + MkPrime ---
cat("Installing R packages on Hamilton (this takes a few minutes)...\n")
ssh::ssh_exec_wait(session,
  sprintf("cd %s && bash m131-setup.sh", BASEDIR))

# --- Submit SLURM array job ---
cat("Submitting SLURM job...\n")
out <- ssh::ssh_exec_internal(session,
  sprintf("cd %s && sbatch m131-warmup-validation.slurm", BASEDIR))
cat(rawToChar(out$stdout))
if (out$status != 0) cat("stderr:", rawToChar(out$stderr), "\n")

# The output will contain the job ID (e.g. "Submitted batch job 12345678")
# Record this in remote-jobs.md

cat("\nDone! Check job status with:\n")
cat('  ssh::ssh_exec_wait(session, "squeue -u $USER")\n')
