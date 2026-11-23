testthat::test_that("a proposal records evidence without accepting it or altering the input object", {
  f <- runner_fixture()
  on.exit(unlink(f$root, recursive = TRUE))
  before <- SummarizedExperiment::assay(f$object, "counts")
  proposed <- propose_analysis_step(f$object, f$directory, "qc", runner_qc_step,
                                    parameters = list(min_umi = 3), seed = 42,
                                    input_path = f$input, source_files = f$code)
  testthat::expect_s4_class(proposed$object, "SingleCellExperiment")
  testthat::expect_identical(SummarizedExperiment::assay(proposed$object, "counts"), before)
  testthat::expect_null(f$object$qc_keep)
  testthat::expect_identical(as.logical(proposed$object$qc_keep), c(FALSE, TRUE, TRUE))
  testthat::expect_identical(proposed$diagnostics$retained, 2L)
  history <- read_analysis_history(f$directory)
  testthat::expect_identical(history$status, "proposed")
  testthat::expect_identical(history$parameters[[1]], list(min_umi = 3))
  testthat::expect_identical(history$seed, 42L)
  testthat::expect_true(nzchar(history$code[[1]]$function_body))
  testthat::expect_true(nzchar(history$code[[1]]$files[[1]]$checksum))
  testthat::expect_true(nzchar(history$input[[1]]$checksum))
})

testthat::test_that("explicit decisions gate downstream work and survive reloading the log", {
  f <- runner_fixture()
  on.exit(unlink(f$root, recursive = TRUE))
  pending <- propose_analysis_step(f$object, f$directory, "qc", runner_qc_step,
                                   list(min_umi = 3), input_path = f$input)
  testthat::expect_error(propose_analysis_step(pending$object, f$directory, "next", runner_qc_step,
    list(min_umi = 8), parent_run = pending$run_id), "accepted")
  testthat::expect_error(decide_analysis_step(pending, "accept", ""), "reason")
  accepted <- decide_analysis_step(pending, "accept", "Retains two cells; reviewed for pilot testing.")
  testthat::expect_s4_class(accepted, "SingleCellExperiment")
  child <- propose_analysis_step(accepted, f$directory, "stricter-qc", runner_qc_step,
                                 list(min_umi = 8), parent_run = pending$run_id)
  testthat::expect_null(decide_analysis_step(child, "reject", "Only one retained cell; keep the earlier setting."))
  history <- read_analysis_history(f$directory)
  testthat::expect_identical(history$status[match(c(pending$run_id, child$run_id), history$run_id)],
                           c("accepted", "rejected"))
  testthat::expect_true(all(nzchar(history$reason)))
  testthat::expect_error(propose_analysis_step(child$object, f$directory, "next", runner_qc_step,
    list(min_umi = 3), parent_run = child$run_id), "accepted")
  testthat::expect_error(decide_analysis_step(pending, "reject", "Change decision"), "already")
})

testthat::test_that("revising an upstream choice preserves earlier decisions and supersedes descendants", {
  f <- runner_fixture()
  on.exit(unlink(f$root, recursive = TRUE))
  first <- propose_analysis_step(f$object, f$directory, "qc", runner_qc_step,
                                 list(min_umi = 3), input_path = f$input)
  accepted <- decide_analysis_step(first, "accept", "Pilot choice A")
  child <- propose_analysis_step(accepted, f$directory, "downstream", runner_qc_step,
                                 list(min_umi = 5), parent_run = first$run_id)
  old_head <- decide_analysis_step(child, "accept", "Pilot downstream result")
  revision <- propose_analysis_step(f$object, f$directory, "qc", runner_qc_step,
                                    list(min_umi = 1), input_path = f$input)
  revised <- decide_analysis_step(revision, "accept", "Keep all cells pending further evidence")
  history <- read_analysis_history(f$directory)
  rows <- match(c(first$run_id, child$run_id, revision$run_id), history$run_id)
  testthat::expect_identical(history$status[rows], c("superseded", "superseded", "accepted"))
  testthat::expect_identical(history$decision[rows], rep("accept", 3))
  testthat::expect_identical(history$reason[rows],
    c("Pilot choice A", "Pilot downstream result", "Keep all cells pending further evidence"))
  testthat::expect_identical(history$run_id[history$current], revision$run_id)
  testthat::expect_error(propose_analysis_step(old_head, f$directory, "stale", runner_qc_step,
    list(min_umi = 3), parent_run = child$run_id), "current analysis lineage")
  testthat::expect_identical(as.logical(revised$qc_keep), c(TRUE, TRUE, TRUE))
})

testthat::test_that("failed and altered proposals cannot be accepted and random state is restored", {
  f <- runner_fixture()
  on.exit(unlink(f$root, recursive = TRUE))
  set.seed(31)
  random_state <- .Random.seed
  failing <- function(object, parameters) { stats::runif(1); stop("method failed") }
  testthat::expect_error(propose_analysis_step(f$object, f$directory, "failure", failing,
    input_path = f$input, seed = 1), "method failed")
  testthat::expect_identical(.Random.seed, random_state)
  testthat::expect_identical(read_analysis_history(f$directory)$status, "failed")
  testthat::expect_match(read_analysis_history(f$directory)$reason, "method failed")
  proposed <- propose_analysis_step(f$object, f$directory, "qc", runner_qc_step,
                                    list(min_umi = 3), input_path = f$input)
  proposed$object$qc_keep <- rep(TRUE, 3)
  testthat::expect_error(decide_analysis_step(proposed, "accept", "Try modified result"), "changed")
  testthat::expect_null(decide_analysis_step(proposed, "reject", "Result was modified outside its recorded step"))
})

testthat::test_that("conversion preserves GEO evidence and existing Seurat objects use the same runner", {
  testthat::skip_if_not_installed("SeuratObject")
  f <- runner_fixture()
  on.exit(unlink(f$root, recursive = TRUE))
  f$object$geo_accession <- rep("GSM101", 3)
  S4Vectors::metadata(f$object)$geo <- list(note = "original GEO source")
  saveRDS(f$object, f$input)
  convert <- function(object, parameters) {
    SeuratObject::CreateSeuratObject(counts = SummarizedExperiment::assay(object, "counts"))
  }
  proposed <- propose_analysis_step(f$object, f$directory, "representation", convert, input_path = f$input)
  testthat::expect_s4_class(proposed$object, "Seurat")
  testthat::expect_identical(as.character(proposed$object$geo_accession), rep("GSM101", 3))
  testthat::expect_identical(proposed$object@misc$geo, list(note = "original GEO source"))
  accepted <- decide_analysis_step(proposed, "accept", "Confirmed cell identity and metadata transfer")
  testthat::expect_s4_class(accepted, "Seurat")
  testthat::expect_identical(as.matrix(SeuratObject::LayerData(accepted, layer = "counts")),
                           as.matrix(SummarizedExperiment::assay(f$object, "counts")))
  next_step <- propose_analysis_step(accepted, f$directory, "seurat-check", function(object, parameters) object,
                                    parent_run = proposed$run_id)
  testthat::expect_identical(next_step$object@misc$geo, accepted@misc$geo)
})

testthat::test_that("diagnostic files are recorded as evidence and changed evidence blocks acceptance", {
  f <- runner_fixture()
  on.exit(unlink(f$root, recursive = TRUE))
  evidence <- file.path(f$root, "qc-summary.csv")
  utils::write.csv(data.frame(retained = 2), evidence, row.names = FALSE)
  with_evidence <- function(object, parameters) {
    result <- runner_qc_step(object, parameters)
    result$artifacts <- parameters$evidence
    result
  }
  proposed <- propose_analysis_step(f$object, f$directory, "qc", with_evidence,
                                    list(min_umi = 3, evidence = evidence), input_path = f$input)
  history <- read_analysis_history(f$directory)
  testthat::expect_identical(history$artifacts[[1]]$path, normalizePath(evidence))
  testthat::expect_match(history$artifacts[[1]]$checksum, "^[a-f0-9]{32}$")
  utils::write.csv(data.frame(retained = 1), evidence, row.names = FALSE)
  testthat::expect_error(decide_analysis_step(proposed, "accept", "Plots reviewed"), "evidence.*changed")
  testthat::expect_null(decide_analysis_step(proposed, "reject", "Evidence was edited; rerun for review"))
})

testthat::test_that("only the recorded input and current accepted lineage can advance the analysis", {
  f <- runner_fixture()
  on.exit(unlink(f$root, recursive = TRUE))
  option_a <- propose_analysis_step(f$object, f$directory, "qc", runner_qc_step,
                                    list(min_umi = 3), input_path = f$input)
  option_b <- propose_analysis_step(f$object, f$directory, "qc", runner_qc_step,
                                    list(min_umi = 5), input_path = f$input)
  accepted <- decide_analysis_step(option_a, "accept", "Compare B before continuing")
  testthat::expect_error(decide_analysis_step(option_b, "accept", "Use stale comparison"), "advanced")
  testthat::expect_null(decide_analysis_step(option_b, "reject", "Regenerate comparison from the accepted state"))
  altered <- accepted
  altered$nUMI <- rep(0, 3)
  testthat::expect_error(propose_analysis_step(altered, f$directory, "next", runner_qc_step,
    list(min_umi = 5), parent_run = option_a$run_id), "does not match")
  next_step <- propose_analysis_step(accepted, f$directory, "next", runner_qc_step,
                                     list(min_umi = 5), parent_run = option_a$run_id)
  saveRDS(altered, f$input)
  testthat::expect_error(decide_analysis_step(next_step, "accept", "Input was replaced"), "checkpoint changed")
  testthat::expect_error(propose_analysis_step(accepted, f$directory, "next", runner_qc_step,
    list(min_umi = 5), parent_run = option_a$run_id), "checkpoint changed")
})

testthat::test_that("the QC example produces inspectable evidence without excluding cells", {
  source(file.path("..", "..", "analysis", "examples", "qc_metrics.R"), local = TRUE)
  f <- runner_fixture()
  on.exit(unlink(f$root, recursive = TRUE))
  evidence <- file.path(f$directory, "qc.pdf")
  proposal <- propose_analysis_step(f$object, f$directory, "inspect-metrics", qc_metric_step,
                                    list(plot_path = evidence), input_path = f$input)
  testthat::expect_equal(proposal$diagnostics,
    list(cells = 3L, features = 2L, total_umis = 13, median_umis = 4, median_genes = 1, excluded_cells = 0L))
  testthat::expect_identical(SummarizedExperiment::assay(proposal$object, "counts"),
                            SummarizedExperiment::assay(f$object, "counts"))
  testthat::expect_identical(colnames(proposal$object), c("a", "b", "c"))
  testthat::expect_identical(readChar(evidence, 4L), "%PDF")
  testthat::expect_s4_class(decide_analysis_step(proposal, "accept", "Fixture evidence reviewed"), "SingleCellExperiment")
})

testthat::test_that("subsetting preserves source metadata but changed cell identities or GEO fields fail", {
  f <- runner_fixture()
  on.exit(unlink(f$root, recursive = TRUE))
  f$object$geo_accession <- c("GSM101", "GSM102", "GSM103")
  saveRDS(f$object, f$input)
  subset_step <- function(object, parameters) object[, c(3, 1)]
  proposal <- propose_analysis_step(f$object, f$directory, "subset", subset_step, input_path = f$input)
  testthat::expect_identical(colnames(proposal$object), c("c", "a"))
  testthat::expect_identical(proposal$object$geo_accession, c("GSM103", "GSM101"))
  changed_geo <- function(object, parameters) { object$geo_accession <- rep("GSM999", 3); object }
  testthat::expect_error(propose_analysis_step(f$object, f$directory, "bad-geo", changed_geo), "GEO cell metadata")
  changed_ids <- function(object, parameters) { colnames(object) <- c("x", "y", "z"); object }
  testthat::expect_error(propose_analysis_step(f$object, f$directory, "bad-ids", changed_ids), "cell identifiers")
})

testthat::test_that("a baseline must reference a readable checkpoint file", {
  f <- runner_fixture()
  on.exit(unlink(f$root, recursive = TRUE))
  testthat::expect_error(propose_analysis_step(f$object, f$directory, "qc", runner_qc_step,
    list(min_umi = 3), input_path = f$root), "checkpoint file")
  testthat::expect_false(dir.exists(f$directory))
})
