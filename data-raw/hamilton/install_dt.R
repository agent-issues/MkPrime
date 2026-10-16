.libPaths(c("/nobackup/pjjg18/mkp-study/lib", .libPaths()))
# fread(fill = Inf) in summarize_streamed.R needs data.table >= 1.16.0.
if (!requireNamespace("data.table", quietly = TRUE) ||
    packageVersion("data.table") < "1.16.0") {
  install.packages("data.table",
                   lib = "/nobackup/pjjg18/mkp-study/lib",
                   repos = "https://cloud.r-project.org")
}
cat("data.table version: ", as.character(packageVersion("data.table")), "\n")
