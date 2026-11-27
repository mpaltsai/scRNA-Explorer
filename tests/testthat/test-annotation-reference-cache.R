testthat::test_that("bundled provenance remains honest and its loader never checks online", {
  f <- reference_cache_fixture(); on.exit(unlink(f$root,recursive=TRUE))
  files_before <- list.files(f$root,recursive=TRUE)
  for (sp in c("Homo sapiens","Mus musculus")) {
    ref <- load_gene_annotation_cache(sp,f$app_root,include_ribosomal=FALSE,cache_dir=f$root)
    testthat::expect_match(ref$provenance$source,"origin confirmed by the user",fixed=TRUE)
    testthat::expect_true(is.na(ref$provenance$resource_id))
    testthat::expect_true(is.na(ref$provenance$genome_build))
    testthat::expect_identical(ref$provenance$reference_id,"bundled")
    testthat::expect_true(nrow(ref$table)>50000)
  }
  testthat::expect_identical(list.files(f$root,recursive=TRUE),files_before)
})

testthat::test_that("monthly checks are build-specific, cached and do not retrieve resources", {
  f <- reference_cache_fixture(); on.exit(unlink(f$root,recursive=TRUE))
  calls <- 0L
  catalog <- function(cache_dir) {calls<<-calls+1L;f$fetch(cache_dir)}
  check <- function(day,force=FALSE,fetch=catalog) f$api$check_annotation_updates("Homo sapiens","hg19",
    f$app_root,f$root,force=force,today=as.Date(day),catalog_fetch=fetch)
  first <- check("2026-01-31")
  testthat::expect_identical(first$candidates$resource_id,c("AH102","AH101"))
  testthat::expect_identical(first$next_check,as.Date("2026-02-28"))
  testthat::expect_false(first$due)
  testthat::expect_identical(calls,1L)
  testthat::expect_false(dir.exists(file.path(f$root,"references")))
  invisible(check("2026-02-27"));testthat::expect_identical(calls,1L)
  second <- check("2026-02-28");testthat::expect_identical(calls,2L)
  testthat::expect_identical(second$next_check,as.Date("2026-03-28"))
  failed <- check("2026-03-28",fetch=function(...) stop("offline"))
  testthat::expect_match(failed$error,"offline")
  testthat::expect_identical(failed$last_attempt,as.Date("2026-03-28"))
  testthat::expect_identical(failed$last_success,as.Date("2026-02-28"))
  testthat::expect_identical(failed$candidates,second$candidates)
  testthat::expect_true(failed$due)
  testthat::expect_error(f$api$check_annotation_updates("Homo sapiens","",f$app_root,f$root),"species-matched genome build")
  testthat::expect_error(f$api$check_annotation_updates("Mus musculus","hg19",f$app_root,f$root),"species-matched genome build")
  other <- f$api$annotation_update_state("Homo sapiens","GRCh38",f$root)
  testthat::expect_true(other$due)
  testthat::expect_null(other$candidates)
})

testthat::test_that("downloads are reviewed, verified and immutable while bundled data stay intact", {
  f <- reference_cache_fixture(); on.exit(unlink(f$root,recursive=TRUE))
  candidate <- f$rows[1,,drop=FALSE]
  paths <- file.path(f$app_root,c("annotations_hsap.Rda","annotations_mmus.Rda"))
  hashes <- tools::md5sum(paths)
  calls <- 0L
  fetch <- function(...) {calls<<-calls+1L;f$payload(...)}
  add <- function(reason="Test compatible reference",resource_fetch=fetch) f$api$download_annotation_reference(
    candidate,"2026-04-23",reason,f$app_root,f$root,resource_fetch)
  testthat::expect_error(add("  "),"Provide a reason")
  testthat::expect_identical(calls,0L)
  wrong <- function(...) {x<-f$payload(...);x$genome_build<-"GRCh38";x}
  testthat::expect_error(add(resource_fetch=wrong),"species/build does not match")
  testthat::expect_false(dir.exists(file.path(f$root,"references")))
  entry <- add()
  testthat::expect_match(entry$reference_id,"^AH101_2026-04-23_")
  testthat::expect_identical(entry$provenance$reason,"Test compatible reference")
  ref <- load_gene_annotation_cache("Homo sapiens",f$app_root,FALSE,entry$reference_id,f$root)
  testthat::expect_identical(ref$table,f$payload()$table)
  testthat::expect_identical(ref$provenance$annotation_release,"105")
  testthat::expect_identical(ref$provenance$genome_build,"GRCh37")
  testthat::expect_identical(ref$provenance$resource_id,"AH101")
  testthat::expect_identical(ref$provenance$snapshot_date,"2026-04-23")
  source(file.path("..","..","analysis","R","qc_optional_checks.R"),local=TRUE)
  cc_annotation <- load_cell_cycle_annotation("Homo sapiens","gene_name",f$app_root,
    reference_id=entry$reference_id,cache_dir=f$root)
  testthat::expect_identical(cc_annotation$table,ref$table)
  testthat::expect_identical(cc_annotation$provenance$reference_id,entry$reference_id)
  testthat::expect_true(file.exists(entry$path))
  before <- tools::md5sum(entry$path)
  repeated <- add()
  testthat::expect_identical(repeated$reference_id,entry$reference_id)
  testthat::expect_identical(tools::md5sum(entry$path),before)
  # Two different references coexist and earlier analyses retain their original
  # table, provenance and counts without another expression checkpoint.
  counts <- Matrix::Matrix(matrix(c(9,10,80,0,10,20,70,0,0,1,90,10),nrow=4,
    dimnames=list(c("MT-A","HBA","GeneA","Shared"),c("below","boundary","above"))),sparse=TRUE)
  r <- new_count_import(counts,"A")
  p <- prepare_sample_qc(r,ref,"gene_name")
  sce <- make_baseline_sce(r,p)
  original <- S4Vectors::metadata(sce)$scrna_explorer_qc$config$reference
  next_fetch <- function(...) {x<-f$payload(...);x$table$description[3]<-"new evidence";x}
  next_entry <- add(resource_fetch=next_fetch)
  testthat::expect_false(identical(next_entry$reference_id,entry$reference_id))
  testthat::expect_identical(S4Vectors::metadata(sce)$scrna_explorer_qc$config$reference,original)
  testthat::expect_identical(SummarizedExperiment::assay(sce,"counts"),r$counts[,2:3])
  testthat::expect_identical(tools::md5sum(paths),hashes)
  # A changed file is rejected rather than used with yesterday's provenance.
  saveRDS(data.frame(changed=TRUE),entry$path)
  testthat::expect_error(load_gene_annotation_cache("Homo sapiens",f$app_root,FALSE,entry$reference_id,f$root),"checksum mismatch")
  testthat::expect_error(add(),"failed verification")
})

testthat::test_that("Shiny checks and saves references without accepting or invalidating prepared data", {
  f <- reference_cache_fixture(); on.exit(unlink(f$root,recursive=TRUE))
  calls <- 0L
  server <- function(input,output,session) {
    ref <- annotation_reference_server(input,output,session,f$app_root,f$root,
      catalog_fetch=function(...) {calls<<-calls+1L;f$fetch(...)},resource_fetch=f$payload)
  }
  shiny::testServer(server,{
    session$setInputs(annotationCheckUpdates=0,annotationDownload=0,typeOrganism="Homo sapiens",annotationDataBuild="GRCh37",annotationReference="bundled")
    testthat::expect_match(output$annotationUpdateStatus,"Monthly annotation check is due")
    testthat::expect_null(output$annotationDownloadControls)
    testthat::expect_identical(calls,0L)
    session$setInputs(annotationCheckUpdates=1);session$flushReact()
    testthat::expect_identical(calls,1L)
    testthat::expect_identical(ref$checked()$candidates$resource_id,c("AH102","AH101"))
    testthat::expect_identical(ref$configuration()$reference_id,"bundled")
    session$setInputs(annotationCandidate="AH101",annotationDownload=1);session$flushReact()
    testthat::expect_match(ref$status(),"Confirm the reference download")
    testthat::expect_false(dir.exists(file.path(f$root,"references")))
    session$setInputs(annotationDownloadConsent=TRUE,annotationUpdateReason="Review hg19 mappings",annotationDownload=2);session$flushReact()
    testthat::expect_match(ref$status(),"Saved GRCh37.*AH101")
    testthat::expect_identical(ref$configuration()$reference_id,"bundled")
    session$setInputs(annotationDataBuild="GRCh38")
    testthat::expect_null(ref$checked()$candidates)
    testthat::expect_identical(calls,1L)
  })
})

testthat::test_that("actual QC stays pinned through maintenance and changes only after reference selection", {
  f <- reference_cache_fixture(); on.exit(unlink(f$root,recursive=TRUE))
  counts <- Matrix::Matrix(matrix(c(9,10,80,0,10,20,70,0,0,1,90,10),nrow=4,
    dimnames=list(c("MT-A","HBA","GeneA","Shared"),c("below","boundary","above"))),sparse=TRUE)
  server <- function(input,output,session) {
    qc <- annotation_qc_server(input,output,session,function() list(A=new_count_import(counts,"A")),f$app_root,
      reference_services=list(cache_dir=f$root,catalog_fetch=f$fetch,resource_fetch=f$payload))
  }
  shiny::testServer(server,{
    session$setInputs(addMetadata=0,initFiltering=0,annotationCheckUpdates=0,annotationDownload=0,
      typeOrganism="Homo sapiens",annotationKey_1="gene_name",annotationDataBuild="GRCh37",annotationReference="bundled")
    session$setInputs(addMetadata=1);session$flushReact()
    session$setInputs(initFiltering=1)
    original <- qc$baseline()$A
    prepared_hash <- digest::digest(qc$prepared())
    baseline_hash <- digest::digest(qc$baseline())
    testthat::expect_identical(S4Vectors::metadata(original)$scrna_explorer_qc$config$reference$reference_id,"bundled")
    session$setInputs(annotationCheckUpdates=1);session$flushReact()
    testthat::expect_identical(digest::digest(qc$prepared()),prepared_hash)
    testthat::expect_identical(digest::digest(qc$baseline()),baseline_hash)
    session$setInputs(annotationCandidate="AH101",annotationDownloadConsent=TRUE,annotationUpdateReason="Compare compatible reference",annotationDownload=1)
    session$flushReact()
    testthat::expect_identical(digest::digest(qc$prepared()),prepared_hash)
    testthat::expect_identical(digest::digest(qc$baseline()),baseline_hash)
    entries <- f$api$list_annotation_references("Homo sapiens",f$app_root,f$root)
    id <- setdiff(names(entries),"bundled")
    testthat::expect_length(id,1)
    session$setInputs(annotationReference=id)
    testthat::expect_null(qc$prepared())
    testthat::expect_null(qc$baseline())
    session$setInputs(addMetadata=2);session$flushReact()
    session$setInputs(initFiltering=2)
    testthat::expect_identical(qc$prepared()$A$config$reference$reference_id,id)
    testthat::expect_identical(SummarizedExperiment::assay(qc$baseline()$A,"counts"),counts[,2:3])
    session$setInputs(annotationDataBuild="GRCh38")
    session$setInputs(addMetadata=3);session$flushReact()
    testthat::expect_null(qc$prepared())
    testthat::expect_null(qc$baseline())
    testthat::expect_match(output$paramsAddMetadata,"differs from the declared counts build")
    testthat::expect_identical(S4Vectors::metadata(original)$scrna_explorer_qc$config$reference$reference_id,"bundled")
    testthat::expect_identical(SummarizedExperiment::assay(original,"counts"),counts[,2:3])
  })
})
