testthat::test_that("initial filtering keeps exactly 100 UMIs for single and multiple samples", {
  testthat::skip_if_not_installed("SingleCellExperiment")
  source(file.path("..", "..", "site_data", "R", "init_filter.R"), local = TRUE)
  counts <- Matrix::Matrix(matrix(c(98, 1, 99, 1, 100, 1), nrow = 2,
    dimnames = list(c("gene_a", "gene_b"), c("below", "boundary", "above"))), sparse = TRUE)
  metadata <- list(cells = colnames(counts), nUMI = c(99, 100, 101),
    nGene = rep(2, 3), log10GenesPerUMI = log10(2) / log10(c(99, 100, 101)),
    sample = rep("sample", 3), mtUMI = rep(0, 3), mitoRatio = rep(0, 3),
    rbcUMI = rep(0, 3), rbcRatio = rep(0, 3), counts = counts)
  input <- list(sample = metadata)
  for (multiple in c(FALSE, TRUE)) {
    samples <- if (multiple) c(input, list(second = metadata)) else input
    result <- initial_filt(samples, multiple)
    testthat::expect_identical(names(result), names(samples))
    for (sample in result) {
      object <- sample$se
      testthat::expect_identical(colnames(object), c("boundary", "above"))
      testthat::expect_equal(object$nUMI, c(100, 101))
      testthat::expect_identical(SummarizedExperiment::assay(object, "counts"), counts[, 2:3])
    }
  }
  testthat::expect_identical(input$sample$counts, counts)
  testthat::expect_equal(input$sample$nUMI, c(99, 100, 101))
})

testthat::test_that("QC review reports the inclusive initial UMI boundary without filtering", {
  source(file.path("..", "..", "analysis", "R", "legacy_qc_adapter.R"), local = TRUE)
  f <- runner_fixture()
  on.exit(unlink(f$root, recursive = TRUE))
  object <- f$object
  object$nUMI <- c(99, 100, 101)
  object$mitoRatio <- rep(0, 3)
  object$sample <- rep("sample", 3)
  proposal <- list(object = object, run_id = "boundary-review",
    diagnostics = list(plots = data.frame(function_name = character(), status = character(),
      warnings = character(), error = character())), artifacts = data.frame(path = character()))
  path <- file.path(f$root, "review.html")
  write_legacy_qc_review(object, proposal, path)
  html <- paste(readLines(path), collapse = "\n")
  testthat::expect_match(html, "nUMI &gt;= 100", fixed = TRUE)
  testthat::expect_match(html, "would remove <strong>1 cells</strong>", fixed = TRUE)
  testthat::expect_identical(colnames(proposal$object), colnames(object))
})
