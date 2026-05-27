# Inventory of /nobackup/pjjg18/mkp-study/results/ — figure out which
# array tasks are stable (done) vs in-flight (still being appended).
# Sources Documents/.Renviron explicitly because Rscript HOME != RStudio R_USER.

readRenviron("~/Documents/.Renviron")
session <- ssh::ssh_connect(
  host    = Sys.getenv("sshLogin"),
  keyfile = Sys.getenv("sshKey"),
  passwd  = Sys.getenv("sshPass")
)

cat("=== mtime of mkp_eg_run_2.log per subdir (sorted) ===\n")
cmd <- paste(
  "for d in /nobackup/pjjg18/mkp-study/results/t*_r*/; do",
  "  if [ -f \"$d/mkp_eg_run_2.log\" ]; then",
  "    printf \"%s\\t%s\\t%s\\n\" \"$(stat -c %Y \"$d/mkp_eg_run_2.log\")\"",
  "      \"$(stat -c %y \"$d/mkp_eg_run_2.log\" | cut -d. -f1)\"",
  "      \"$(basename \"$d\")\"",
  "  fi",
  "done | sort -n"
)
out <- ssh::ssh_exec_internal(session, cmd, error = FALSE)
all_lines <- strsplit(rawToChar(out$stdout), "\n")[[1]]
all_lines <- all_lines[nzchar(all_lines)]
cat(sprintf("Total subdirs with run_2.log: %d\n", length(all_lines)))
cat("\nOldest 10 (most likely stable):\n")
for (ln in head(all_lines, 10)) cat("  ", ln, "\n")
cat("\nNewest 10 (potentially being written):\n")
for (ln in tail(all_lines, 10)) cat("  ", ln, "\n")

cat("\n=== sacct array 17141180 state distribution ===\n")
out <- ssh::ssh_exec_internal(session,
  "sacct -j 17141180 --format=JobID,State -X | tail -n +3 | awk '{print $2}' | sort | uniq -c",
  error = FALSE)
cat(rawToChar(out$stdout))

cat("\n=== file sizes for one representative subdir (t01_r01) ===\n")
out <- ssh::ssh_exec_internal(session,
  "du -sh /nobackup/pjjg18/mkp-study/results/t01_r01/*", error = FALSE)
cat(rawToChar(out$stdout))

cat("\n=== count fully-thinned samples in one log (lines - header) ===\n")
out <- ssh::ssh_exec_internal(session,
  paste("for f in /nobackup/pjjg18/mkp-study/results/t01_r01/mkp_eg_run_1.log",
        "          /nobackup/pjjg18/mkp-study/results/t01_r01/mkp_eg_run_2.log; do",
        "  printf \"%s\\t%d lines\\n\" \"$(basename \"$f\")\" \"$(wc -l < \"$f\")\"",
        "done"),
  error = FALSE)
cat(rawToChar(out$stdout))

cat("\n=== ascii head of one log (column names + first 2 rows) ===\n")
out <- ssh::ssh_exec_internal(session,
  "head -3 /nobackup/pjjg18/mkp-study/results/t01_r01/mkp_eg_run_1.log",
  error = FALSE)
cat(rawToChar(out$stdout))

ssh::ssh_disconnect(session)
