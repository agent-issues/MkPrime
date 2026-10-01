# Synthetic-fixture checks for the m9-pilot, m131 and install-stamp scripts.
# Run from the repo root:  Rscript dev/m9-pilot/test_harness_scripts.R
# The m131 end-to-end check is skipped when MkPrime is not installed.

Check <- function(label, ok) {
  cat(sprintf("%s %s\n", if (isTRUE(ok)) "PASS" else "FAIL", label))
  if (!isTRUE(ok)) quit(status = 1L)
}
Skip <- function(label) cat("SKIP", label, "\n")
root <- normalizePath(".")

# --- Every shipped script parses (m131-analyse-warmup.R had an invalid "\." escape) ---
for (f in c(Sys.glob("inst/hamilton/m131-*.R"), Sys.glob("dev/m9-pilot/*.R"))) {
  Check(paste("parses:", f), !inherits(try(parse(f), silent = TRUE), "try-error"))
}

# --- m131 result files are found by the analyser's pattern ---
src <- readLines("inst/hamilton/m131-analyse-warmup.R")
patLine <- grep('pattern = ', src, value = TRUE)[1]
pat <- eval(parse(text = sub('.*pattern = ("[^"]*").*', "\\1", patLine)))
Check("analyser pattern matches m131_X_seedN.rds only",
      identical(grepl(pat, c("m131_Sun2018_seed1.rds", "m131xrds", "m131_analysis.pdf")),
                c(TRUE, FALSE, FALSE)))

# --- process_pilot: congeneric tips are kept distinct until selection ---
LoadFunction <- function(script, name, env) {
  for (e in parse(script)) {
    if (is.call(e) && identical(as.character(e[[1]]), "<-") && identical(as.character(e[[2]]), name)) {
      eval(e, env)
      return(invisible())
    }
  }
  stop("no ", name, " in ", script)
}
if (requireNamespace("TreeTools", quietly = TRUE)) {
  suppressPackageStartupMessages({library(ape); library(TreeTools)})
  env <- new.env()
  for (fn in c("GenusOf", "ToGenusTree", "de_zz")) LoadFunction("dev/m9-pilot/process_pilot.R", fn, env)
  env$KeepTip <- KeepTip
  tr <- read.tree(text = "(((zzCanis_lupus,Canis_dirus),Felis_catus),(Ursus_arctos,Vulpes_vulpes));")
  tr <- env$de_zz(tr)
  Check("de_zz keeps the species suffix", "Canis_lupus" %in% tr$tip.label)
  out <- env$ToGenusTree(tr, c("Canis", "Felis", "Ursus"))
  Check("one tip per genus, relabelled", setequal(out$tip.label, c("Canis", "Felis", "Ursus")))
  Check("genus outside the WCT dropped", !"Vulpes" %in% out$tip.label)
  Check("collapse-first would have produced duplicates",
        anyDuplicated(env$GenusOf(tr$tip.label)) > 0L)
} else {
  Skip("ToGenusTree (TreeTools not installed)")
}

# --- build_mkp_tarball.sh stamps the commit into the tarball ---
tmp <- tempfile(); dir.create(tmp); setwd(tmp)
system2("git", c("init", "-q", "."))
dir.create("R"); writeLines("f <- function() 1", "R/f.R")
writeLines(c("Package: toypkg", "Version: 0.0.1", "Title: Toy", "Description: Toy.",
             "License: GPL-3", "Author: A", "Maintainer: A <a@b.c>"), "DESCRIPTION")
writeLines("export(f)", "NAMESPACE")
system2("git", c("add", "-A"))
system2("git", c("-c", "user.name=t", "-c", "user.email=t@t", "commit", "-q", "-m", "x"))
sha <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
outdir <- file.path(tmp, "out"); dir.create(outdir)
res <- system2("bash", c(file.path(root, "data-raw/hamilton/build_mkp_tarball.sh"), outdir),
               stdout = TRUE, stderr = TRUE)
tb <- Sys.glob(file.path(outdir, "toypkg_*.tar.gz"))
Check("build script produces a tarball", length(tb) == 1L)
d <- tempfile(); dir.create(d)
untar(tb, files = "toypkg/DESCRIPTION", exdir = d)
desc <- read.dcf(file.path(d, "toypkg/DESCRIPTION"))
Check("tarball DESCRIPTION carries RemoteSha == HEAD", identical(unname(desc[1, "RemoteSha"]), sha))
Check("working tree DESCRIPTION untouched", !any(grepl("RemoteSha", readLines("DESCRIPTION"))))
setwd(root)

# --- m131-warmup-validation.R runs end to end on a toy matrix ---
if (requireNamespace("MkPrime", quietly = TRUE) && requireNamespace("TreeTools", quietly = TRUE)) {
  work <- tempfile(); dir.create(file.path(work, "data-raw"), recursive = TRUE)
  set.seed(1)
  m <- matrix(sample(0:1, 6 * 12, TRUE), 6, 12, dimnames = list(paste0("t", 1:6), NULL))
  ape::write.nexus.data(setNames(lapply(seq_len(6), function(i) m[i, ]), rownames(m)),
                        file = file.path(work, "data-raw", "Toy.nex"), format = "standard")
  setwd(work)
  status <- system2("timeout", c("600", "Rscript", file.path(root, "inst/hamilton/m131-warmup-validation.R")),
                    stdout = FALSE, stderr = FALSE,
                    env = c("DATASET=Toy", "SEED=1", "NITER=3000", "OUTDIR=res"))
  setwd(root)
  Check("m131-warmup-validation.R exits cleanly", identical(status, 0L))
  rds <- readRDS(file.path(work, "res", "m131_Toy_seed1.rds"))
  Check("m131 result records taxa and characters", rds$nTip == 6L && rds$nChar == 12L)
  Check("m131 result holds a warmup trace", length(rds$warmup_trace) > 0L)
} else {
  Skip("m131-warmup-validation.R end to end (MkPrime not installed)")
}
