# Code to be run with
#   R -d "valgrind --tool=memcheck --leak-check=full --error-exitcode=1" --vanilla < memcheck/thisfile.R
# Tests the installed package: load_package = "source" would recompile src/
# in place via pkgload and leave the installed build untested.
testthat::test_local(load_package = "installed")
