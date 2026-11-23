testthat::test_that("registry joins GEO metadata by accession and preserves its evidence", {
  f <- registry_fixture()
  on.exit(unlink(f$root, recursive = TRUE))
  registry <- build_gse176031_registry(f$data_root, f$metadata)

  testthat::expect_identical(registry$samples$accession, c("GSM102", "GSM101"))
  testthat::expect_identical(registry$samples$patient, c("P1", "P1"))
  testthat::expect_identical(registry$samples$source_name,
                           c("Organoid", "Radical Prostatectomy Paired Normal"))
  testthat::expect_identical(registry$samples$protocol, rep("Seq-Well S3", 2))
  testthat::expect_true(all(is.na(registry$samples$capture_id)))
  testthat::expect_true(all(is.na(registry$samples$biological_sample_id)))
  testthat::expect_identical(registry$samples$dge_file,
                           c("GSM102_organoid_dge.txt.gz", "GSM101_tissue_dge.txt.gz"))
  testthat::expect_identical(registry$samples$n_cells, c(2L, 2L))
  testthat::expect_identical(registry$barcodes[["GSM101_tissue_dge.txt.gz"]],
                           c("shared", "cell-b"))
  testthat::expect_true(any(registry$metadata$field == "patient"))
  testthat::expect_false(registry$capabilities$background_droplets_confirmed)
  testthat::expect_false(registry$capabilities$independent_quality_labels_confirmed)
})

testthat::test_that("ambiguous and missing DGE evidence is reported without selecting a file", {
  f <- registry_fixture()
  on.exit(unlink(f$root, recursive = TRUE))
  unlink(file.path(f$data_root, "GSM102_organoid_dge.txt.gz"))
  f$write_dge("GSM101_other_dge.txt.gz")
  f$write_dge("GSM999_unmatched_dge.txt.gz")
  registry <- build_gse176031_registry(f$data_root, f$metadata)

  testthat::expect_true(all(is.na(registry$samples$dge_file)))
  testthat::expect_true(all(c("missing_dge", "duplicate_dge_accession", "unmatched_dge") %in%
                            registry$issues$code))
  testthat::expect_equal(nrow(registry$inputs), 3L)
})

testthat::test_that("conflicting reference declarations remain explicit instead of being resolved by inference", {
  f <- registry_fixture()
  on.exit(unlink(f$root, recursive = TRUE))
  f$write_metadata(c(
    '!Sample_data_processing\t"Aligned to mm8 whole genome"\t"Aligned to mm8 whole genome"',
    '!Sample_data_processing\t"Genome_build: hg19"\t"Genome_build: hg19"'
  ))
  registry <- build_gse176031_registry(f$data_root, f$metadata)

  testthat::expect_identical(registry$samples$genome_build, c("hg19", "hg19"))
  testthat::expect_identical(registry$samples$donor_id, c("P1", "P1"))
  testthat::expect_equal(sum(registry$issues$code == "reference_conflict"), 2L)
  testthat::expect_true(all(registry$metadata$source_field[
    registry$metadata$field == "data_processing"] == "!Sample_data_processing"))
  testthat::expect_true("capture_mapping_unresolved" %in% registry$issues$code)
  testthat::expect_true("repeated_barcodes_across_inputs" %in% registry$issues$code)
})

testthat::test_that("later tickets resolve a relocated registered input and reject stale headers", {
  f <- registry_fixture()
  on.exit(unlink(f$root, recursive = TRUE))
  registry <- build_gse176031_registry(f$data_root, f$metadata)
  moved <- file.path(f$root, "moved-dge")
  testthat::expect_true(file.rename(f$data_root, moved))
  input <- resolve_gse176031_input(registry, "GSM101", data_root = moved)
  testthat::expect_identical(input$dge_path,
                           normalizePath(file.path(moved, "GSM101_tissue_dge.txt.gz")))
  testthat::expect_identical(input$barcodes, c("shared", "cell-b"))
  testthat::expect_identical(input$sample$patient, "P1")
  con <- gzfile(input$dge_path, "wt")
  writeLines("GENE\tchanged\tcell-b", con)
  close(con)
  testthat::expect_error(resolve_gse176031_input(registry, "GSM101", data_root = moved),
                         "changed")
  testthat::expect_error(resolve_gse176031_input(registry, "GSM999", data_root = moved),
                         "exactly one")
})

testthat::test_that("duplicate GEO accessions are diagnosed and cannot select an input", {
  f <- registry_fixture()
  on.exit(unlink(f$root, recursive = TRUE))
  rows <- readLines(gzfile(f$metadata), warn = FALSE)
  rows <- sub('"GSM102"\t"GSM101"', '"GSM101"\t"GSM101"', rows, fixed = TRUE)
  con <- gzfile(f$metadata, "wt")
  writeLines(rows, con)
  close(con)
  registry <- build_gse176031_registry(f$data_root, f$metadata)
  testthat::expect_true("duplicate_metadata_accession" %in% registry$issues$code)
  testthat::expect_true(all(is.na(registry$samples$dge_file)))
  testthat::expect_error(resolve_gse176031_input(registry, "GSM101"), "exactly one")
})

testthat::test_that("existing R artifacts are preferred candidates without deserialization", {
  f <- registry_fixture()
  on.exit(unlink(f$root, recursive = TRUE))
  artifact_root <- file.path(f$root, "artifacts")
  dir.create(artifact_root)
  saveRDS(list(note = "Not a validated expression object"),
          file.path(artifact_root, "GSM101_existing.rds"))
  writeLines("not a readable HDF5 file", file.path(artifact_root, "GSM101_filtered.h5ad"))
  registry <- build_gse176031_registry(f$data_root, f$metadata,
                                      artifact_roots = artifact_root)
  testthat::expect_identical(registry$artifacts$format, "rds")
  testthat::expect_identical(registry$artifacts$inspection_status[1], "R_serialization_header_only")
  testthat::expect_false(any(registry$artifacts$counts_verified))
  testthat::expect_false(any(registry$artifacts$compatibility_verified))
  testthat::expect_equal(nrow(resolve_gse176031_input(registry, "GSM101")$artifacts), 1L)
})

testthat::test_that("confirmed hg19 resolves the reference conflict without replacing source text", {
  f <- registry_fixture()
  on.exit(unlink(f$root, recursive = TRUE))
  f$write_metadata(c(
    '!Sample_data_processing\t"Aligned to mm8 whole genome"\t"Aligned to mm8 whole genome"',
    '!Sample_data_processing\t"Genome_build: hg19"\t"Genome_build: hg19"'
  ))
  registry <- build_gse176031_registry(f$data_root, f$metadata, confirmed_genome_build = "hg19")
  testthat::expect_identical(registry$samples$confirmed_genome_build, c("hg19", "hg19"))
  testthat::expect_true(all(grepl("mm8", registry$samples$data_processing)))
  testthat::expect_true("reference_conflict_resolved" %in% registry$issues$code)
  testthat::expect_false("reference_conflict" %in% registry$issues$code)
  testthat::expect_identical(registry$provenance$reference_confirmation, "explicit user confirmation: hg19")
})

testthat::test_that("GEO sample metadata attaches to existing SCE and Seurat objects by cell identity", {
  testthat::skip_if_not_installed("SingleCellExperiment")
  testthat::skip_if_not_installed("SeuratObject")
  f <- registry_fixture()
  on.exit(unlink(f$root, recursive = TRUE))
  registry <- build_gse176031_registry(f$data_root, f$metadata, confirmed_genome_build = "hg19")
  counts <- Matrix::sparseMatrix(i = c(1, 2, 3), j = c(1, 2, 1), x = c(3, 5, 2),
                                 dims = c(3, 2),
                                 dimnames = list(c("A1BG", "B2M", "MT-CO1"), c("a", "b")))
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
  sce$sample <- c("local-1", "local-2")
  S4Vectors::metadata(sce)$existing_note <- "retain"
  mapping <- c(b = "GSM102", a = "GSM101")
  attached <- attach_geo_sample_metadata(sce, registry, mapping)
  testthat::expect_identical(SummarizedExperiment::assay(attached, "counts"), counts)
  testthat::expect_identical(as.character(attached$geo_source_name),
                           c("Radical Prostatectomy Paired Normal", "Organoid"))
  testthat::expect_identical(attached$sample, c("local-1", "local-2"))
  testthat::expect_identical(S4Vectors::metadata(attached)$existing_note, "retain")
  testthat::expect_equal(nrow(S4Vectors::metadata(attached)$geo$sample_metadata), nrow(registry$metadata))
  testthat::expect_identical(as.character(attached$geo_confirmed_genome_build), c("hg19", "hg19"))
  testthat::expect_identical(attach_geo_sample_metadata(attached, registry, mapping), attached)

  seu <- SeuratObject::CreateSeuratObject(counts = counts)
  seu@misc$existing_note <- "retain"
  seu <- attach_geo_sample_metadata(seu, registry, mapping)
  testthat::expect_identical(as.character(seu$geo_accession), c("GSM101", "GSM102"))
  testthat::expect_identical(seu@misc$existing_note, "retain")
  testthat::expect_equal(nrow(seu@misc$geo$sample_metadata), nrow(registry$metadata))
  testthat::expect_equal(as.matrix(SeuratObject::LayerData(seu, layer = "counts")), as.matrix(counts))
  testthat::expect_error(attach_geo_sample_metadata(seu, registry, "GSM999"), "unregistered")
})

testthat::test_that("R command produces reusable local metadata outputs without overwriting existing files", {
  f <- registry_fixture()
  on.exit(unlink(f$root, recursive = TRUE))
  script <- normalizePath(file.path("..", "..", "scripts", "register_gse176031.R"), mustWork = FALSE)
  output <- file.path(f$root, "registry")
  args <- c("--vanilla", shQuote(script), "--data-root", shQuote(f$data_root),
            "--series-matrix", shQuote(f$metadata), "--output-dir", shQuote(output),
            "--confirmed-genome-build", "hg19")
  log <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), args, stdout = TRUE, stderr = TRUE))
  testthat::expect_null(attr(log, "status"), info = paste(log, collapse = "\n"))
  testthat::expect_true(all(file.exists(file.path(output,
    c("registry.rds", "samples.csv", "source_metadata.csv", "inputs.csv", "issues.csv", "artifacts.csv")))))
  if (!file.exists(file.path(output, "registry.rds"))) return(invisible(NULL))
  registry <- readRDS(file.path(output, "registry.rds"))
  testthat::expect_identical(resolve_gse176031_input(registry, "GSM101")$barcodes,
                           c("shared", "cell-b"))
  before <- unname(tools::md5sum(file.path(output, "registry.rds")))
  repeated <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), args, stdout = TRUE, stderr = TRUE))
  testthat::expect_identical(attr(repeated, "status"), 1L)
  testthat::expect_match(paste(repeated, collapse = "\n"), "not empty")
  testthat::expect_identical(unname(tools::md5sum(file.path(output, "registry.rds"))), before)
})

testthat::test_that("source literals and reserved accession fields cannot corrupt cell mapping", {
  testthat::skip_if_not_installed("SingleCellExperiment")
  f <- registry_fixture()
  on.exit(unlink(f$root, recursive = TRUE))
  f$write_metadata(c(
    '!Sample_characteristics_ch1\t"accession: external-a"\t"accession: external-b"',
    '!Sample_characteristics_ch1\t"note: NA"\t"note: NA"'
  ))
  registry <- build_gse176031_registry(f$data_root, f$metadata)
  counts <- Matrix::sparseMatrix(i = 1, j = 1, x = 2, dims = c(1, 1),
                                 dimnames = list("A1BG", "cell"))
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
  sce <- attach_geo_sample_metadata(sce, registry, "GSM101")
  testthat::expect_identical(as.character(sce$geo_accession), "GSM101")
  testthat::expect_true("NA" %in% S4Vectors::metadata(sce)$geo$sample_metadata$value)
  testthat::expect_true("external-b" %in% S4Vectors::metadata(sce)$geo$sample_metadata$value)
})

testthat::test_that("an empty final barcode is reported rather than silently dropped", {
  f <- registry_fixture()
  on.exit(unlink(f$root, recursive = TRUE))
  f$write_dge("GSM101_tissue_dge.txt.gz", cells = c("cell-a", ""))
  registry <- build_gse176031_registry(f$data_root, f$metadata)
  testthat::expect_true("invalid_dge_header" %in% registry$issues$code)
  testthat::expect_true(is.na(registry$samples$dge_file[registry$samples$accession == "GSM101"]))
  testthat::expect_error(resolve_gse176031_input(registry, "GSM101"), "unambiguous")
})
