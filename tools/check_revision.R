# Usage: Rscript buRst/tools/check_revision.R buRst [--integration]
args <- commandArgs(trailingOnly = TRUE)
package_dir <- if (length(args)) args[[1L]] else "."
if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Install devtools before running this development check.")
}
cat("Focused likelihood, initialization, and numerical regression tests\n")
devtools::test(package_dir, filter = "likelihood|initialization|numeric")
if ("--integration" %in% args) {
  Sys.setenv(NOT_CRAN = "true")
  cat("Full simulation / HDF5 / fitting integration test\n")
  devtools::test(package_dir, filter = "full-circle")
}
print(sessionInfo())
