.libPaths(c("/nobackup/pjjg18/mkp-study/lib", .libPaths()))
if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table",
                   lib = "/nobackup/pjjg18/mkp-study/lib",
                   repos = "https://cloud.r-project.org")
}
cat("data.table version: ", as.character(packageVersion("data.table")), "\n")
