args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
root <- dirname(dirname(normalizePath(script)))
if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("Running the registry tests requires an existing testthat installation.")
}
testthat::test_dir(file.path(root, "tests", "testthat"), stop_on_failure = TRUE)
