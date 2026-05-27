#!/bin/bash
set -e
module load r/4.5.1 >/dev/null 2>&1
export R_LIBS=/nobackup/pjjg18/mkp-study/red-team/lib:/nobackup/pjjg18/mkp-study/lib
Rscript -e 'library(MkPrime); cat("MkPrime loaded OK\n"); cat("RunMkPrime kPrimePrior:\n"); fr <- formals(MkPrime::RunMkPrime); if (!is.null(fr$kPrimePrior)) print(fr$kPrimePrior) else cat("(default not eval-printed; checking partition-api)\n"); src <- tryCatch(readLines(system.file("R", "RunMkPrime.R", package = "MkPrime"), warn = FALSE), error=function(e) NULL); cat("priors found in shipped sources:", length(grep("empirical_geometric|beta_geometric|logseries", src)), "matches\n")'
