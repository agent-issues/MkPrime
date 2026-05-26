#!/usr/bin/env Rscript
# Lane N1 — fast_exp + closed-form P(t) stability driver
#
# Usage:
#   Rscript dev/red-team/numerical/fast-exp-driver.R           # full grid
#   Rscript dev/red-team/numerical/fast-exp-driver.R --quick   # smoke test (<2 min)
#
# Strategy:
#   There is no eigendecomposition in this codebase — all transition
#   probabilities are closed-form (JC / MkN / F81; see src/rate_matrix.cpp,
#   src/likelihood.cpp, src/gibbs_partial_cl.h, src/mcmc_likelihood.cpp).
#   So this audit reduces to:
#     (a) verify that src/fast_exp.h::fast_neg_exp matches std::exp within
#         the claimed <1e-15 relative error;
#     (b) verify the closed-form formulas as a whole (including any
#         catastrophic cancellation in 1 - E) across the stress grid.
#
#   Rmpfr is NOT installed on this machine; mpmath via Python is the
#   reference oracle. This R driver shells out to the Python driver and
#   relays its CSV output. Re-implementing fast_neg_exp bit-exactly in R
#   would be slower and less faithful than in Python (which exposes
#   IEEE-754 reinterpretation via struct.pack).
#
# Lane: numerical-auditor / N1

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

args <- commandArgs(trailingOnly = TRUE)
quick <- "--quick" %in% args

# Resolve script directory robustly: --file=... when Rscripted, fallback to cwd
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg)) {
  script_path <- sub("^--file=", "", script_arg[1])
  here <- normalizePath(dirname(script_path), winslash = "/", mustWork = FALSE)
} else {
  # Sourced from a session
  ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  here <- if (!is.null(ofile)) normalizePath(dirname(ofile), winslash = "/", mustWork = FALSE) else getwd()
}
if (!nzchar(here) || here == ".") {
  here <- normalizePath(".", winslash = "/")
}

# Locate python
python_paths <- c(
  Sys.getenv("PYTHON", unset = NA),
  "C:/Python312/python.exe",
  "python3", "python"
)
python_paths <- python_paths[!is.na(python_paths) & nzchar(python_paths)]
python <- NULL
for (p in python_paths) {
  status <- tryCatch(
    suppressWarnings(system2(p, c("-c", shQuote("import mpmath")),
                             stdout = FALSE, stderr = FALSE)),
    error = function(e) 127L,
    warning = function(w) 127L
  )
  if (is.numeric(status) && length(status) == 1 && status == 0) {
    python <- p
    break
  }
}
if (is.null(python)) {
  stop("Need a python with mpmath available. Tried: ",
       paste(python_paths, collapse = ", "),
       "\nInstall: <python> -m pip install --user mpmath")
}
cat("Using python at:", python, "\n")

py_driver <- file.path(here, "fast-exp-driver.py")
if (!file.exists(py_driver)) {
  stop("Python driver missing: ", py_driver)
}

cmd_args <- c(py_driver)
if (quick) cmd_args <- c(cmd_args, "--quick")

status <- system2(python, cmd_args, stdout = "", stderr = "")
if (status != 0) {
  stop("Python driver failed with status ", status)
}

cat("\nDriver finished. Results CSV under dev/red-team/numerical/fast-exp-results/\n")
