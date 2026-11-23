testthat::test_that("existing metadata calculations run through a proposal while preserving the SCE", {
  source(file.path("..", "..", "analysis", "R", "legacy_qc_adapter.R"), local = TRUE)
  f <- runner_fixture()
  on.exit(unlink(f$root, recursive = TRUE))
  counts <- Matrix::Matrix(matrix(c(10,0,1,1,0,0, 20,2,2,2,1,1, 2,0,0,0,1,1), nrow = 6,
    dimnames = list(c("A1BG","B2M","MT-ND1","MT-ND2","HBB","HBD"), c("a","b","c"))), sparse = TRUE)
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
  sce$geo_accession <- rep("GSM101", 3)
  S4Vectors::metadata(sce)$geo <- list(source = "original")
  saveRDS(sce, f$input)
  steps <- make_legacy_qc_steps(file.path("..", "..", "site_data"))
  proposal <- propose_analysis_step(sce, f$directory, "existing-qc-metadata", steps$metadata,
    list(organism = "Homo sapiens", sample_column = "geo_accession"), input_path = f$input,
    source_files = steps$source_files)
  testthat::expect_equal(proposal$object$nUMI, c(12,28,4))
  testthat::expect_equal(proposal$object$nGene, c(3,6,3))
  testthat::expect_equal(proposal$object$mitoRatio, c(1/6,1/7,0))
  testthat::expect_equal(proposal$object$rbcRatio, c(0,1/14,1/2))
  testthat::expect_identical(SummarizedExperiment::assay(proposal$object, "counts"), counts)
  testthat::expect_identical(S4Vectors::metadata(proposal$object)$geo, list(source = "original"))
  testthat::expect_null(sce$nUMI)
  testthat::expect_identical(read_analysis_history(f$directory)$status, "proposed")
  audit <- proposal$diagnostics$gene_audit
  testthat::expect_true(all(audit$checks$max_absolute_difference == 0))
  testthat::expect_true(all(c("MT-ND1","MT-ND2") %in% audit$genes$gene_name[audit$genes$set == "mitochondrial" & audit$genes$present]))
  testthat::expect_true(all(c("HBB","HBD") %in% audit$genes$gene_name[audit$genes$set == "hemoglobin" & audit$genes$present]))
  testthat::expect_false(audit$genes$present[audit$genes$gene_name == "HBBP1"])
})

testthat::test_that("an existing sample identity is never silently overwritten", {
  source(file.path("..", "..", "analysis", "R", "legacy_qc_adapter.R"),local=TRUE)
  f <- runner_fixture()
  on.exit(unlink(f$root,recursive=TRUE))
  f$object$geo_accession <- rep("GSM101",3)
  f$object$sample <- rep("original_sample",3)
  saveRDS(f$object,f$input)
  steps <- make_legacy_qc_steps(file.path("..","..","site_data"))
  testthat::expect_error(propose_analysis_step(f$object,f$directory,"metadata",steps$metadata,
    list(organism="Homo sapiens",sample_column="geo_accession"),input_path=f$input),"conflicts with existing sample metadata")
  testthat::expect_identical(f$object$sample,rep("original_sample",3))
  testthat::expect_null(f$object$nUMI)
})

testthat::test_that("existing filter_cells uses explicit thresholds after metadata acceptance", {
  source(file.path("..", "..", "analysis", "R", "legacy_qc_adapter.R"), local = TRUE)
  f <- runner_fixture()
  on.exit(unlink(f$root, recursive = TRUE))
  counts <- Matrix::Matrix(matrix(c(10,0,1,1,0,0, 20,2,2,2,1,1, 2,0,0,0,1,1), nrow = 6,
    dimnames = list(c("A1BG","B2M","MT-ND1","MT-ND2","HBB","HBD"), c("a","b","c"))), sparse = TRUE)
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
  sce$geo_accession <- rep("GSM101", 3)
  saveRDS(sce, f$input)
  steps <- make_legacy_qc_steps(file.path("..", "..", "site_data"))
  proposal <- propose_analysis_step(sce, f$directory, "metadata", steps$metadata,
    list(organism = "Homo sapiens", sample_column = "geo_accession"), input_path = f$input)
  thresholds <- list(minUMIsPerCell = 10, maxUMIsPerCell = 28, minGenesPerCell = 0,
    maxGenesPerCell = 6, noveltyPerCell = 0, minMitoPerCell = 0, maxMitoPerCell = 1,
    minRbcPerCell = 0, maxRbcPerCell = 1)
  testthat::expect_error(propose_analysis_step(proposal$object, f$directory, "filter", steps$filter,
    list(thresholds = thresholds), parent_run = proposal$run_id), "accepted")
  accepted <- decide_analysis_step(proposal, "accept", "Fixture metadata matches known counts")
  filtered <- propose_analysis_step(accepted, f$directory, "filter", steps$filter,
    list(thresholds = thresholds), parent_run = proposal$run_id)
  testthat::expect_identical(colnames(filtered$object), c("a","b"))
  testthat::expect_identical(filtered$object$geo_accession, rep("GSM101", 2))
  testthat::expect_identical(SummarizedExperiment::assay(filtered$object, "counts"), counts[,1:2])
  testthat::expect_equal(filtered$diagnostics$excluded_cells, 1)
  testthat::expect_identical(ncol(accepted), 3L)
})
