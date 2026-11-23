runner_module <- file.path("..", "..", "analysis", "R", "analysis_runner.R")
if (file.exists(runner_module)) source(runner_module)

runner_fixture <- function() {
  testthat::skip_if_not_installed("SingleCellExperiment")
  testthat::skip_if_not_installed("digest")
  root <- tempfile("analysis-runner-")
  dir.create(root)
  counts <- Matrix::sparseMatrix(i = c(1, 2, 2), j = c(1, 2, 3), x = c(1, 4, 8),
                                 dims = c(2, 3),
                                 dimnames = list(c("A1BG", "B2M"), c("a", "b", "c")))
  object <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
  input <- file.path(root, "existing_sce.rds")
  saveRDS(object, input)
  code <- file.path(root, "method.R")
  writeLines("# Declared method source for the experiment", code)
  list(root = root, directory = file.path(root, "decisions"), object = object,
       input = input, code = code)
}

runner_qc_step <- function(object, parameters) {
  object$nUMI <- Matrix::colSums(SummarizedExperiment::assay(object, "counts"))
  object$qc_keep <- object$nUMI >= parameters$min_umi
  list(object = object, diagnostics = list(retained = sum(object$qc_keep), cells = ncol(object)))
}
