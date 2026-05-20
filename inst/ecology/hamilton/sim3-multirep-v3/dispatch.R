# dispatch.R — orchestrate Sim 3 v3 multirep dispatch from local laptop
#
# Connects to Hamilton, syncs setup script + run_rep + sbatch script,
# triggers setup, submits the array job.
#
# Run after SSH env vars are set (sshLogin, sshKey, sshPass) or the
# user has authenticated interactively.

suppressPackageStartupMessages(library(ssh))

# Reuse a top-level session if present
if (!exists("session") || !inherits(session, "ssh_session")) {
  session <- ssh::ssh_connect(
    host = Sys.getenv("sshLogin", "pjjg18@hamilton8.dur.ac.uk"),
    keyfile = Sys.getenv("sshKey", "C:/Users/pjjg18/id_ed25519"),
    passwd = Sys.getenv("sshPass", "")
  )
}

PROJECT_REMOTE <- "/nobackup/pjjg18/mkp-sim3-multirep-v3"
HPC_DIR <- "/nobackup/pjjg18/mkp-sim3-multirep-v3/MkPrime/inst/ecology/hamilton/sim3-multirep-v3"

# 1. Run one-time project setup (clones + builds package).
cat("Running setup_project.sh on Hamilton...\n")
out <- ssh::ssh_exec_internal(session,
  sprintf("mkdir -p %s/logs && cd %s && bash %s/setup_project.sh 2>&1 | tail -40",
          PROJECT_REMOTE, PROJECT_REMOTE, HPC_DIR),
  error = FALSE)
cat(rawToChar(out$stdout))
if (out$status != 0) {
  # First-time: HPC_DIR doesn't exist yet because we haven't cloned.
  cat("setup_project.sh not yet on Hamilton (or first-time clone). Bootstrapping...\n")
  # Bootstrap: upload setup_project.sh first
  setup_local <- "inst/ecology/hamilton/sim3-multirep-v3/setup_project.sh"
  ssh::ssh_exec_internal(session,
    sprintf("mkdir -p %s/bootstrap", PROJECT_REMOTE), error = FALSE)
  ssh::scp_upload(session, setup_local,
                  to = file.path(PROJECT_REMOTE, "bootstrap/"))
  ssh::ssh_exec_internal(session,
    sprintf("sed -i 's/\\r$//' %s/bootstrap/setup_project.sh",
            PROJECT_REMOTE), error = FALSE)
  out2 <- ssh::ssh_exec_internal(session,
    sprintf("bash %s/bootstrap/setup_project.sh 2>&1 | tail -40", PROJECT_REMOTE),
    error = FALSE)
  cat(rawToChar(out2$stdout))
  if (out2$status != 0) {
    cat("\nstderr:\n"); cat(rawToChar(out2$stderr))
    stop("Setup failed.")
  }
}

# 2. Sync sbatch script (fix DOS line endings)
cat("\nUploading sbatch script and run_rep.R...\n")
ssh::scp_upload(session, "inst/ecology/hamilton/sim3-multirep-v3/sim3-multirep-v3.sh",
                to = paste0(HPC_DIR, "/"))
ssh::scp_upload(session, "inst/ecology/hamilton/sim3-multirep-v3/run_rep.R",
                to = paste0(HPC_DIR, "/"))
ssh::ssh_exec_internal(session,
  sprintf("sed -i 's/\\r$//' %s/sim3-multirep-v3.sh %s/run_rep.R",
          HPC_DIR, HPC_DIR),
  error = FALSE)

# 3. Submit the array job
cat("\nSubmitting array job (8 reps)...\n")
out <- ssh::ssh_exec_internal(session,
  sprintf("cd %s && sbatch sim3-multirep-v3.sh", HPC_DIR), error = FALSE)
cat(rawToChar(out$stdout))
if (out$status != 0) {
  cat("stderr:\n"); cat(rawToChar(out$stderr))
}
