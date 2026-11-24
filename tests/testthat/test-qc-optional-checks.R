testthat::test_that("cell-cycle scoring adds diagnostic fields without filtering or regression", {
  testthat::skip_if_not_installed("Seurat")
  source(file.path("..", "..", "analysis", "R", "qc_optional_checks.R"), local = TRUE)
  f <- runner_fixture()
  on.exit(unlink(f$root, recursive = TRUE))
  counts <- Matrix::Matrix(matrix(rep(1:11, length.out=160), nrow=20,
    dimnames=list(c("PCNA","MCM5","TOP2A","MKI67","Gene_X","Gene-X",paste0("Gene",1:14)),paste0("cell",1:8))), sparse=TRUE)
  sce <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=counts))
  sce$geo_accession <- rep("GSM101",8)
  saveRDS(sce,f$input)
  reference <- list(name="test marker reference",species="Homo sapiens",
                    s_genes=c("PCNA","MCM5"),g2m_genes=c("TOP2A","MKI67"))
  proposal <- propose_analysis_step(sce,f$directory,"cell-cycle",cell_cycle_qc_step,
    list(organism="Homo sapiens",reference=reference,annotation=cell_cycle_fixture_annotation(reference),
      nbin=2L,ctrl=2L,method_seed=1L),input_path=f$input)
  testthat::expect_identical(SummarizedExperiment::assay(proposal$object,"counts"),counts)
  testthat::expect_identical(colnames(proposal$object),colnames(sce))
  testthat::expect_true(all(is.finite(proposal$object$cell_cycle_S_score)))
  testthat::expect_true(all(is.finite(proposal$object$cell_cycle_G2M_score)))
  testthat::expect_true(all(proposal$object$cell_cycle_phase %in% c("G1","S","G2M","Undecided")))
  testthat::expect_false(proposal$diagnostics$regression_applied)
  testthat::expect_identical(proposal$diagnostics$excluded_cells,0L)
  testthat::expect_identical(proposal$diagnostics$coverage$matched,c(2L,2L))
  testthat::expect_identical(proposal$object$geo_accession,sce$geo_accession)
})

testthat::test_that("unresolved capture information never becomes a doublet prediction", {
  source(file.path("..", "..", "analysis", "R", "qc_optional_checks.R"), local=TRUE)
  f <- runner_fixture()
  on.exit(unlink(f$root,recursive=TRUE))
  f$object$geo_capture_id <- rep(NA_character_,3)
  saveRDS(f$object,f$input)
  proposal <- propose_analysis_step(f$object,f$directory,"doublets",doublet_qc_step,
    list(capture_column="geo_capture_id",technology="Seq-Well S3"),input_path=f$input)
  testthat::expect_identical(proposal$diagnostics$status,"not_assessed")
  testthat::expect_match(paste(proposal$diagnostics$reasons,collapse=" "),"Capture identities are unresolved")
  testthat::expect_null(proposal$object$scDblFinder.score)
  testthat::expect_null(proposal$object$scDblFinder.class)
  testthat::expect_identical(SummarizedExperiment::assay(proposal$object,"counts"),
                            SummarizedExperiment::assay(f$object,"counts"))
})

testthat::test_that("doublet detection flags known-capture input without removing cells or counts", {
  testthat::skip_if_not_installed("scDblFinder")
  source(file.path("..", "..", "analysis", "R", "qc_optional_checks.R"),local=TRUE)
  f <- runner_fixture()
  on.exit(unlink(f$root,recursive=TRUE))
  set.seed(14)
  sce <- scDblFinder::mockDoubletSCE(ncells=c(60,60),ngenes=300,dbl.rate=0.1)
  sce$capture <- rep("simulated",ncol(sce))
  saveRDS(sce,f$input)
  original <- SummarizedExperiment::assay(sce,"counts")
  rates <- data.frame(capture="simulated",expected_rate=0.1,uncertainty=0.02,source="Simulation input parameter")
  proposal <- propose_analysis_step(sce,f$directory,"doublets",doublet_qc_step,
    list(capture_column="capture",capture_evidence="Known simulated capture",technology="simulation",
         rate_table=rates,method_seed=1L),input_path=f$input)
  testthat::expect_identical(proposal$diagnostics$status,"scored")
  testthat::expect_identical(SummarizedExperiment::assay(proposal$object,"counts"),original)
  testthat::expect_identical(colnames(proposal$object),colnames(sce))
  testthat::expect_true(all(is.finite(proposal$object$scDblFinder.score)))
  testthat::expect_true(all(proposal$object$scDblFinder.class %in% c("singlet","doublet")))
  testthat::expect_identical(proposal$diagnostics$excluded_cells,0L)
})

testthat::test_that("missing cell-cycle markers do not become zero scores", {
  source(file.path("..", "..", "analysis", "R", "qc_optional_checks.R"), local = TRUE)
  f <- runner_fixture()
  on.exit(unlink(f$root,recursive=TRUE))
  reference <- list(name="missing reference",species="Homo sapiens",s_genes="AbsentS",g2m_genes="B2M")
  testthat::expect_error(propose_analysis_step(f$object,f$directory,"cell-cycle",cell_cycle_qc_step,
    list(organism="Homo sapiens",reference=reference,annotation=cell_cycle_fixture_annotation(reference),
      method_seed=1L),input_path=f$input),"no matching features")
  testthat::expect_null(f$object$cell_cycle_S_score)
})

testthat::test_that("cell-cycle scoring on Seurat preserves its active assay and RNA counts", {
  testthat::skip_if_not_installed("Seurat")
  source(file.path("..", "..", "analysis", "R", "qc_optional_checks.R"),local=TRUE)
  f <- runner_fixture()
  on.exit(unlink(f$root,recursive=TRUE))
  counts <- Matrix::Matrix(matrix(rep(1:11,length.out=160),nrow=20,
    dimnames=list(c("PCNA","MCM5","TOP2A","MKI67",paste0("Gene",1:16)),paste0("cell",1:8))),sparse=TRUE)
  object <- SeuratObject::CreateSeuratObject(counts=counts)
  object[["ADT"]] <- SeuratObject::CreateAssay5Object(counts=counts[1:2,])
  SeuratObject::DefaultAssay(object) <- "ADT"
  object@misc$geo <- list(source="original")
  saveRDS(object,f$input)
  reference <- list(name="test marker reference",species="Homo sapiens",s_genes=c("PCNA","MCM5"),g2m_genes=c("TOP2A","MKI67"))
  proposal <- propose_analysis_step(object,f$directory,"cell-cycle",cell_cycle_qc_step,
    list(organism="Homo sapiens",reference=reference,annotation=cell_cycle_fixture_annotation(reference),
      nbin=2L,ctrl=2L,method_seed=1L),input_path=f$input)
  testthat::expect_identical(SeuratObject::DefaultAssay(proposal$object),"ADT")
  testthat::expect_identical(SeuratObject::LayerData(proposal$object,assay="RNA",layer="counts"),counts)
  testthat::expect_identical(proposal$object@misc$geo,list(source="original"))
  testthat::expect_null(object[[]]$cell_cycle_S_score)
  testthat::expect_true(all(is.finite(proposal$object$cell_cycle_S_score)))
})

testthat::test_that("ambient RNA without suitable background inputs remains not assessed", {
  source(file.path("..", "..", "analysis", "R", "qc_optional_checks.R"), local = TRUE)
  f <- runner_fixture()
  on.exit(unlink(f$root,recursive=TRUE))
  proposal <- propose_analysis_step(f$object,f$directory,"ambient-readiness",ambient_rna_readiness_step,
    list(technology="Seq-Well S3",droplet_compatible=FALSE),input_path=f$input)
  testthat::expect_identical(proposal$diagnostics$status,"not_assessed")
  testthat::expect_match(paste(proposal$diagnostics$reasons,collapse=" "),"background")
  testthat::expect_null(proposal$object$ambient_contamination)
  testthat::expect_identical(SummarizedExperiment::assay(proposal$object,"counts"),
                            SummarizedExperiment::assay(f$object,"counts"))
})

testthat::test_that("pooled DecontX analysis needs an explicit exploratory choice", {
  source(file.path("..", "..", "analysis", "R", "qc_optional_checks.R"),local=TRUE)
  f <- runner_fixture()
  on.exit(unlink(f$root,recursive=TRUE))
  testthat::expect_error(propose_analysis_step(f$object,f$directory,"ambient-decontx",decontx_qc_step,
    list(method_seed=1L),input_path=f$input),"batch_column or explicitly enable exploratory_pooled")
})

testthat::test_that("DecontX estimates retain original counts and separate preliminary clusters", {
  testthat::skip_if_not_installed("decontX")
  source(file.path("..", "..", "analysis", "R", "qc_optional_checks.R"),local=TRUE)
  f <- runner_fixture()
  on.exit(unlink(f$root,recursive=TRUE))
  first <- c(rep(100,50),rep(1,50))
  second <- rev(first)
  counts <- Matrix::Matrix(cbind(matrix(first,100,20),matrix(second,100,20)),sparse=TRUE)
  dimnames(counts) <- list(paste0("gene",1:100),paste0("cell",1:40))
  sce <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=counts))
  sce$known_cluster <- rep(c("A","B"),each=20)
  sce$capture <- rep("simulation",40)
  saveRDS(sce,f$input)
  proposal <- propose_analysis_step(sce,f$directory,"decontx",decontx_qc_step,
    list(batch_column="capture",clusters_column="known_cluster",method_seed=1L,
         model_args=list(maxIter=50L)),input_path=f$input)
  testthat::expect_identical(proposal$diagnostics$status,"estimated")
  testthat::expect_identical(SummarizedExperiment::assay(proposal$object,"counts"),counts)
  testthat::expect_identical(SummarizedExperiment::assayNames(proposal$object),"counts")
  testthat::expect_true(all(is.finite(proposal$object$ambient_decontX_fraction)))
  testthat::expect_true(all(proposal$object$ambient_decontX_fraction>=0 & proposal$object$ambient_decontX_fraction<=1))
  testthat::expect_identical(proposal$object$known_cluster,sce$known_cluster)
  testthat::expect_false(proposal$diagnostics$counts_corrected)
  testthat::expect_true(length(proposal$diagnostics$fit_diagnostics)>0L)
  testthat::expect_identical(proposal$diagnostics$fit_diagnostics[[1]]$iteration_limit,50L)
})
