testthat::test_that("sample labels change atomically without changing data or stable identities", {
  f <- input_fixture(); on.exit(unlink(f$root,recursive=TRUE))
  sce <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=f$counts))
  sce$sample <- "original metadata"
  original <- list(A=new_count_import(f$counts,"A",object=sce),B=new_count_import(f$counts,"B"))
  renamed <- rename_import_samples(original,c(" Patient 1 ","Control"))
  testthat::expect_identical(names(renamed),c("Patient 1","Control"))
  testthat::expect_identical(renamed[[1]]$sample_id,"Patient 1")
  testthat::expect_identical(renamed[[1]]$source$sample_id_at_import,"A")
  for (i in 1:2) {
    testthat::expect_identical(renamed[[i]]$counts,original[[i]]$counts)
    testthat::expect_identical(renamed[[i]]$cell_identity,original[[i]]$cell_identity)
    testthat::expect_identical(renamed[[i]]$source_object,original[[i]]$source_object)
  }
  testthat::expect_identical(renamed[[1]]$sample_name_history[[1]]$from,"A")
  testthat::expect_identical(rename_import_samples(renamed,names(renamed)),renamed)
  for (bad in list(c("same","same"),c("same"," same "),c("","B"),c(NA,"B"),"A",c("../A","B")))
    testthat::expect_error(rename_import_samples(original,bad))
  testthat::expect_identical(names(original),c("A","B"))
})

testthat::test_that("app sample naming refreshes QC labels and invalidates prepared results only on Apply", {
  f <- input_fixture(); on.exit(unlink(f$root,recursive=TRUE))
  server <- function(input,output,session) {
    importer <- input_import_server(input,output,session,function(...) list(A=new_count_import(f$counts,"A"),B=new_count_import(f$counts,"B")))
    qc <- annotation_qc_server(input,output,session,importer$ready,normalizePath(file.path("..","..","site_data")))
  }
  shiny::testServer(server,{
    session$setInputs(startQC=0,applySampleNames=0,addMetadata=0,initFiltering=0,
      typeOrganism="Homo sapiens",annotationKey_1="gene_name",annotationKey_2="gene_name")
    session$setInputs(startQC=1); session$flushReact()
    session$setInputs(addMetadata=1); session$flushReact()
    session$setInputs(initFiltering=1)
    testthat::expect_false(is.null(qc$baseline()))
    identity <- importer$ready()$A$cell_identity
    session$setInputs(importSampleName_1="Tumour",importSampleName_2="Normal")
    testthat::expect_identical(names(importer$ready()),c("A","B"))
    testthat::expect_false(is.null(qc$baseline()))
    session$setInputs(applySampleNames=1)
    testthat::expect_identical(names(importer$ready()),c("Tumour","Normal"))
    testthat::expect_identical(importer$ready()$Tumour$cell_identity,identity)
    testthat::expect_null(qc$baseline())
    session$setInputs(addMetadata=2); session$flushReact(); session$setInputs(initFiltering=2)
    testthat::expect_true(all(qc$baseline()$Tumour$qc_sample=="Tumour"))
    testthat::expect_identical(S4Vectors::metadata(qc$baseline()$Tumour)$scrna_explorer_qc$import_sample_names[[1]]$from,"A")
    session$setInputs(importSampleName_2="Tumour",applySampleNames=2)
    testthat::expect_match(importer$status(),"unique")
    testthat::expect_identical(names(importer$ready()),c("Tumour","Normal"))
    testthat::expect_false(is.null(qc$baseline()))
  })
})
