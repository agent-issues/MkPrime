#!/usr/bin/env Rscript
# Build a RevBayes fixed-state probe for one (pid, model) cell: K posterior
# states from a finished RB cell, re-evaluated under each model toggle.
#
# Usage:
#   Rscript make_probe.R <pid> <model> --rb-dir=<cell dir with <model>_run_1.{log,trees}>
#     [--out-dir=probe] [--remote-dir=/nobackup/$USER/rt11/probe]
#     [--neotrans=/nobackup/$USER/neotrans] [--local-neotrans=../../../neotrans]
#     [--k=8] [--seed=1]
#
# Writes probe_<pid>_<model>.Rev, trees_<pid>_<model>.nex,
# sample_<pid>_<model>.rds and, where the matrix has polymorphisms,
# project<pid>.{trans,neo}.depoly.nex (polymorphisms -> '?').
# Copy <out-dir> to <remote-dir>, run run_probe.sh there, copy out_*.txt back,
# then run check_probe.R.

suppressPackageStartupMessages(library("ape"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("Usage: make_probe.R <pid> <model> --rb-dir=DIR [opts]")
pid <- args[[1]]
model <- args[[2]]
stopifnot(model %in% c("by_nt_9v", "by_nt_kv"))
user <- Sys.getenv("USER", "pjjg18")
opt <- list(rb_dir = NA_character_, out_dir = "probe",
            remote_dir = sprintf("/nobackup/%s/rt11/probe", user),
            neotrans = sprintf("/nobackup/%s/neotrans", user),
            local_neotrans = "../../../neotrans", k = "8", seed = "1")
for (a in args[-(1:2)]) {
  kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
  if (length(kv) != 2L) stop("Bad option: ", a)
  opt[[gsub("-", "_", kv[[1]])]] <- kv[[2]]
}
if (is.na(opt$rb_dir)) stop("--rb-dir is required")
nState <- as.integer(opt$k)
set.seed(as.integer(opt$seed))
dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)
cell <- sprintf("%s_%s", pid, model)
localProj <- file.path(opt$local_neotrans, "inst", "projects")
remoteProj <- file.path(opt$neotrans, "inst", "projects")

# --- Sample K post-burn-in states where log and trees coincide -------------
trTab <- read.delim(file.path(opt$rb_dir, sprintf("%s_run_1.trees", model)),
                    stringsAsFactors = FALSE)
lgTab <- read.delim(file.path(opt$rb_dir, sprintf("%s_run_1.log", model)),
                    stringsAsFactors = FALSE, check.names = FALSE)
common <- intersect(trTab$Iteration, lgTab$Iteration)
common <- common[common >= stats::quantile(common, 0.25)]
pick <- sort(sample(common, nState))
trees <- structure(lapply(pick, function(it) {
  read.tree(text = gsub("\\[[^]]*\\]", "", trTab$phylogeny[trTab$Iteration == it]))
}), class = "multiPhylo")
par <- lgTab[match(pick, lgTab$Iteration),
             c("Iteration", "Likelihood", "rate_loss", "rate_neo", "rate_log_sd")]
write.nexus(trees, file = file.path(opt$out_dir, sprintf("trees_%s.nex", cell)))
saveRDS(list(trees = trees, par = par),
        file.path(opt$out_dir, sprintf("sample_%s.rds", cell)))

# --- Polymorphism -> '?' copies, so data handling can be isolated -----------
depolyPath <- function(split) {
  src <- readLines(file.path(localProj, sprintf("project%s.%s.nex", pid, split)))
  dst <- gsub("\\([0-9,]+\\)|\\{[0-9,]+\\}", "?", src)
  if (identical(src, dst)) return(file.path(remoteProj, sprintf("project%s.%s.nex", pid, split)))
  name <- sprintf("project%s.%s.depoly.nex", pid, split)
  writeLines(dst, file.path(opt$out_dir, name))
  file.path(opt$remote_dir, name)
}
neoD <- depolyPath("neo")
transD <- depolyPath("trans")

# --- Rev script --------------------------------------------------------------
transBlock <- function(tag, dataVar, coding) {
  if (model == "by_nt_9v") {
    sprintf(paste0('t_%1$s ~ dnPhyloCTMC(tree=tr, branchRates=pr[2], siteRates=rc, ',
                   'Q=fnJC(9), type="Standard", coding="%3$s")\n',
                   't_%1$s.clamp(%2$s)\nlt_%1$s <- t_%1$s.lnProbability()'),
            tag, dataVar, coding)
  } else {
    sprintf(paste0('lt_%1$s <- 0\nfor (k in 2:10) {\n',
                   '  tbs_%1$s[k-1] <- %2$s\n  tbs_%1$s[k-1].setNumStatesPartition(k)\n',
                   '  tbs_%1$s[k-1].removeExcludedCharacters()\n',
                   '  if (tbs_%1$s[k-1].nchar() > 0) {\n',
                   '    tk_%1$s[k-1] ~ dnPhyloCTMC(tree=tr, branchRates=pr[2], siteRates=rc, ',
                   'Q=fnJC(k), type="Standard", coding="%3$s")\n',
                   '    tk_%1$s[k-1].clamp(tbs_%1$s[k-1])\n',
                   '    lt_%1$s <- lt_%1$s + tk_%1$s[k-1].lnProbability()\n  }\n}'),
            tag, dataVar, coding)
  }
}
# Neo variants: rescaled (T) or not (F); coding variable (V) or all (A);
# original data or polymorphisms -> ? (_D).
neoVariants <- list(c("TV", "qT", "variable", "neo"), c("TA", "qT", "all", "neo"),
                    c("FV", "qF", "variable", "neo"), c("FA", "qF", "all", "neo"),
                    c("TV_D", "qT", "variable", "neoD"), c("TA_D", "qT", "all", "neoD"),
                    c("FV_D", "qF", "variable", "neoD"), c("FA_D", "qF", "all", "neoD"))
transVariants <- list(c("V", "trans", "variable"), c("A", "trans", "all"),
                      c("V_D", "transD", "variable"), c("A_D", "transD", "all"))
cols <- c(paste0("n", vapply(neoVariants, `[`, "", 1)),
          paste0("t", vapply(transVariants, `[`, "", 1)))

rev <- c(
  sprintf('neo <- readDiscreteCharacterData("%s/project%s.neo.nex")', remoteProj, pid),
  sprintf('trans <- readDiscreteCharacterData("%s/project%s.trans.nex")', remoteProj, pid),
  sprintf('neoD <- readDiscreteCharacterData("%s")', neoD),
  sprintf('transD <- readDiscreteCharacterData("%s")', transD),
  "nChar <- v(neo.nchar(), trans.nchar())",
  sprintf('trees <- readTrees("%s/trees_%s.nex")', opt$remote_dir, cell),
  sprintf('print("HEADER i iter %s")', paste(cols, collapse = " "))
)
for (i in seq_len(nState)) {
  p <- par[i, ]
  rev <- c(rev,
    sprintf("tr <- trees[%d]", i),
    sprintf("rate_log_sd <- %.10g", p$rate_log_sd),
    sprintf("rate_loss <- %.10g", p$rate_loss),
    sprintf("rate_neo <- %.10g", p$rate_neo),
    "rc_raw := fnDiscretizeDistribution(dnLognormal(-rate_log_sd*rate_log_sd/2, rate_log_sd), 6)",
    "rc := rc_raw / mean(rc_raw)",
    "pr := [rate_neo/(1+rate_neo), 1/(1+rate_neo)] / nChar * sum(nChar)",
    "rates := [[0.0, 2/(1+rate_loss)], [2*rate_loss/(1+rate_loss), 0.0]]",
    "qT := fnFreeK(rates)", "qF := fnFreeK(rates, rescaled=FALSE)",
    "sd := Simplex(2*rate_loss/(1+rate_loss), 2/(1+rate_loss))",
    vapply(neoVariants, function(v) sprintf(paste0(
      'n_%1$s ~ dnPhyloCTMC(tree=tr, branchRates=pr[1], siteRates=rc, Q=%2$s, ',
      'rootFrequencies=sd, type="Standard", coding="%3$s")\n',
      'n_%1$s.clamp(%4$s)\nln_%1$s <- n_%1$s.lnProbability()'),
      v[[1]], v[[2]], v[[3]], v[[4]]), ""),
    vapply(transVariants, function(v) transBlock(v[[1]], v[[2]], v[[3]]), ""),
    sprintf('print("ROW", %d, %d, %s, sep=" ")', i, p$Iteration,
            paste(c(paste0("ln_", vapply(neoVariants, `[`, "", 1)),
                    paste0("lt_", vapply(transVariants, `[`, "", 1))), collapse = ", "))
  )
}
writeLines(c(rev, "q()"), file.path(opt$out_dir, sprintf("probe_%s.Rev", cell)))
cat(sprintf("[make_probe] Wrote %s/probe_%s.Rev (%d states)\n", opt$out_dir, cell, nState))
