testthat::test_that("AnnotationHub mappings distinguish missing, ambiguous and duplicated records", {
  source(file.path("..","..","analysis","R","qc_optional_checks.R"),local=TRUE)
  annotation <- list(table=data.frame(
    gene_name=c("Sok","Sok","Samb","Samb","Gok","Gcollision","Alias","OtherAlias"),
    gene_id=c("ENSG1","ENSG1","ENSG2","ENSG3","ENSG4","ENSG5","ENSG6","ENSG6")),
    species="Homo sapiens",feature_id_type="gene_id",provenance=list(source="Synthetic annotation"))
  requested <- list(S=c("Sok","Samb","AbsentAnnotation","Alias"),G2M=c("Gok","Gcollision","AbsentFeature"))
  features <- c("ENSG1.1","ENSG2.1","ENSG4.2","ENSG5.1","ENSG5.2","ENSG6.1","AbsentAnnotation")
  result <- .qc_cell_cycle_mapping(features,requested,annotation,"Homo sapiens",TRUE)
  testthat::expect_identical(result$matched,list(S="ENSG1.1",G2M="ENSG4.2"))
  testthat::expect_identical(result$mapping$status,c("matched","ambiguous_annotation","missing_annotation",
    "ambiguous_annotation","matched","ambiguous_features","missing_annotation"))
  # An annotated gene absent from the matrix is different from an absent annotation.
  annotation$table <- rbind(annotation$table,data.frame(gene_name="AbsentFeature",gene_id="ENSG7"))
  mapped <- .qc_cell_cycle_mapping(features,requested,annotation,"Homo sapiens",TRUE)
  testthat::expect_identical(tail(mapped$mapping$status,1),"missing_feature")
  exact <- .qc_cell_cycle_mapping(features,requested,annotation,"Homo sapiens",FALSE)
  testthat::expect_length(exact$matched$S,0)
  annotation$table <- annotation$table[nrow(annotation$table):1,]
  reversed <- .qc_cell_cycle_mapping(features,requested,annotation,"Homo sapiens",TRUE)
  testthat::expect_identical(reversed$matched,mapped$matched)
  testthat::expect_error(.qc_cell_cycle_mapping(features,requested,annotation,"Mus musculus"),"species")
  testthat::expect_error(.qc_cell_cycle_mapping(features,requested,NULL,"Homo sapiens"),"Supply annotation")
  both <- .qc_cell_cycle_mapping(features,list(S="Sok",G2M="Sok"),annotation,"Homo sapiens",TRUE)
  testthat::expect_true(all(both$mapping$status=="ambiguous_program"))
  testthat::expect_true(all(is.na(both$mapping$matched_feature)))
})

testthat::test_that("symbol and Ensembl scoring agree through annotations without rewriting counts", {
  testthat::skip_if_not_installed("Seurat")
  source(file.path("..","..","analysis","R","qc_optional_checks.R"),local=TRUE)
  reference <- list(name="test programs",species="Homo sapiens",s_genes=c("PCNA","MCM5"),g2m_genes=c("TOP2A","MKI67"))
  symbols <- c(reference$s_genes,reference$g2m_genes,paste0("Background",1:16))
  counts <- Matrix::Matrix(matrix(rep(1:11,length.out=160),nrow=20,
    dimnames=list(symbols,paste0("cell",1:8))),sparse=TRUE)
  sce <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=counts))
  sce$geo_accession <- rep("GSM101",8)
  a <- cell_cycle_fixture_annotation(reference)
  p <- list(organism="Homo sapiens",reference=reference,annotation=a,nbin=2L,ctrl=2L,method_seed=1L)
  symbol_result <- cell_cycle_qc_step(sce,p)
  ids <- c(paste0(a$table$gene_id,".5"),paste0("ENSG",100:115,".2"))
  ensembl <- sce
  rownames(ensembl) <- ids
  original <- SummarizedExperiment::assay(ensembl,"counts")
  p$annotation$feature_id_type <- "gene_id"
  p$strip_ensembl_version <- TRUE
  id_result <- cell_cycle_qc_step(ensembl,p)
  testthat::expect_identical(id_result$object$cell_cycle_S_score,symbol_result$object$cell_cycle_S_score)
  testthat::expect_identical(id_result$object$cell_cycle_G2M_score,symbol_result$object$cell_cycle_G2M_score)
  testthat::expect_identical(id_result$object$cell_cycle_phase,symbol_result$object$cell_cycle_phase)
  testthat::expect_identical(SummarizedExperiment::assay(id_result$object,"counts"),original)
  testthat::expect_identical(rownames(id_result$object),ids)
  testthat::expect_identical(id_result$object$geo_accession,sce$geo_accession)
  testthat::expect_null(ensembl$cell_cycle_phase)
  testthat::expect_identical(id_result$diagnostics$matched_genes$S,ids[1:2])
  testthat::expect_identical(id_result$diagnostics$marker_mapping$matched_feature,ids[1:4])
  testthat::expect_identical(id_result$diagnostics$annotation$table_hash,digest::digest(a$table,algo="sha256"))

  # Exercise the combined QC stage and its exported mapping, not just the scorer.
  source(file.path("..","..","analysis","R","legacy_qc_adapter.R"),local=TRUE)
  source(file.path("..","..","analysis","examples","preprocessing_qc_step.R"),local=TRUE)
  directory <- tempfile("annotation-qc-export-")
  on.exit(unlink(directory,recursive=TRUE))
  p$annotation <- a
  result <- preprocessing_qc_step(sce,list(app_root=file.path("..","..","site_data"),
    legacy=list(organism="Homo sapiens",sample_column="geo_accession"),
    cell_cycle=p,diagnostic_dir=directory))
  exported <- read.csv(file.path(directory,"cell-cycle-marker-coverage.csv"))
  testthat::expect_identical(exported$matched_feature,symbols[1:4])
  testthat::expect_true(all(exported$in_matrix))
  testthat::expect_true(file.exists(file.path(directory,"cell-cycle-annotation-provenance.txt")))
  testthat::expect_identical(SummarizedExperiment::assay(result$object,"counts"),counts)
  testthat::expect_identical(result$diagnostics$checks$doublets$status,"disabled")
})

testthat::test_that("the app annotation cache is reused with honest provenance", {
  source(file.path("..","..","analysis","R","qc_optional_checks.R"),local=TRUE)
  a <- load_cell_cycle_annotation("Homo sapiens","gene_name",file.path("..","..","site_data"))
  testthat::expect_true(all(c("gene_id","gene_name") %in% names(a$table)))
  testthat::expect_identical(a$provenance$md5,unname(tools::md5sum(a$provenance$path)))
  testthat::expect_true(is.na(a$provenance$genome_build))
  testthat::expect_true(is.na(a$provenance$resource_id))
  testthat::expect_error(load_cell_cycle_annotation("Homo sapiens","auto"),"feature_id_type")
})
