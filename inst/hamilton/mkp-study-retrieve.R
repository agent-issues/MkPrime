# Retrieve a representative SAMPLE of Hamilton job 17140607 / 17141180
# results.  Pulls the LAST 50k samples of each mkp_eg_run_{1,2}.log from
# a diversified subset of stable subdirs (not currently being written),
# plus all 260 mkp_eg_checkpoint.rds files (tiny, metadata only).
#
# The full per-task .log files are ~400 MB each; we tail to ~50 MB per
# file before transfer, which is enough for logP decomposition, ESS,
# and topology-hash autocorrelation analyses.
#
# Pulls trees only for ONE subdir as a smoke test (~180 MB).
#
# Usage:
#   readRenviron("~/Documents/.Renviron")     # if not autoloaded
#   source("inst/hamilton/mkp-study-retrieve.R")

library(ssh)

if (!exists("session") || !inherits(session, "ssh_session")) {
  if (!nzchar(Sys.getenv("sshPass"))) {
    readRenviron("~/Documents/.Renviron")
  }
  session <- ssh::ssh_connect(
    host    = Sys.getenv("sshLogin", "pjjg18@hamilton8.dur.ac.uk"),
    keyfile = Sys.getenv("sshKey", "C:/Users/pjjg18/id_ed25519"),
    passwd  = Sys.getenv("sshPass", "")
  )
}

JOB_ID <- "17140607-17141180"
REMOTE <- "/nobackup/pjjg18/mkp-study/results"
LOCAL  <- file.path("inst", "hamilton", paste0("mkp-study-", JOB_ID))
dir.create(LOCAL, showWarnings = FALSE, recursive = TRUE)

# Diversified sample: 1 rep from each of 5 evenly-spaced trees, chosen
# from the "stable" (mtime > 1 h old) pool seen in the inventory.
SAMPLE_DIRS <- c("t01_r01", "t05_r05", "t08_r05", "t04_r09",
                 "t03_r09", "t06_r03")
TAIL_LINES  <- 50000L

cat(sprintf("Pulling tail(%d) of mkp_eg_run_{1,2}.log + checkpoint for %d subdirs:\n",
            TAIL_LINES, length(SAMPLE_DIRS)))
for (sd in SAMPLE_DIRS) cat("  ", sd, "\n")

# 1. Pull tail of each log to a tmp file on the remote, then scp_download.
remote_tmp <- "/tmp/mkp_eg_sample.tar.gz"
build_cmd <- c(
  "set -e",
  "rm -rf /tmp/mkp_eg_sample && mkdir -p /tmp/mkp_eg_sample"
)
for (sd in SAMPLE_DIRS) {
  build_cmd <- c(build_cmd,
    sprintf("mkdir -p /tmp/mkp_eg_sample/%s", sd),
    sprintf("cp %s/%s/mkp_eg_checkpoint.rds /tmp/mkp_eg_sample/%s/ 2>/dev/null || true",
            REMOTE, sd, sd),
    sprintf("head -1 %s/%s/mkp_eg_run_1.log > /tmp/mkp_eg_sample/%s/mkp_eg_run_1.log",
            REMOTE, sd, sd),
    sprintf("tail -n %d %s/%s/mkp_eg_run_1.log | grep -v '^#' >> /tmp/mkp_eg_sample/%s/mkp_eg_run_1.log",
            TAIL_LINES, REMOTE, sd, sd),
    sprintf("head -1 %s/%s/mkp_eg_run_2.log > /tmp/mkp_eg_sample/%s/mkp_eg_run_2.log 2>/dev/null || true",
            REMOTE, sd, sd),
    sprintf("tail -n %d %s/%s/mkp_eg_run_2.log 2>/dev/null | grep -v '^#' >> /tmp/mkp_eg_sample/%s/mkp_eg_run_2.log",
            TAIL_LINES, REMOTE, sd, sd)
  )
}
build_cmd <- c(build_cmd,
  sprintf("cd /tmp && tar czf %s mkp_eg_sample/", remote_tmp))
cat("\nPacking on remote...\n")
out <- ssh::ssh_exec_internal(session,
  paste(build_cmd, collapse = " && "), error = FALSE)
if (out$status != 0) {
  cat("Remote packing stderr:\n", rawToChar(out$stderr), "\n")
}

cat("Downloading...\n")
ssh::scp_download(session, remote_tmp, LOCAL)
ssh::ssh_exec_internal(session,
  sprintf("rm -f %s && rm -rf /tmp/mkp_eg_sample", remote_tmp))

cat("Extracting locally...\n")
old_wd <- getwd(); on.exit(setwd(old_wd), add = TRUE)
setwd(LOCAL)
untar("mkp_eg_sample.tar.gz")
setwd(old_wd)

cat("\nLocal files retrieved:\n")
print(list.files(LOCAL, recursive = TRUE))
