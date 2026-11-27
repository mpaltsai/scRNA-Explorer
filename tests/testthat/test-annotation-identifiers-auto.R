identifier_test_record <- function(ids) {
  new_count_import(Matrix::Matrix(matrix(10,nrow=length(ids),ncol=2,
    dimnames=list(ids,c("cellA","cellB"))),sparse=TRUE))
}

testthat::test_that("automatic suggestions distinguish gene IDs, symbols and review cases", {
  ref <- annotation_fixture()
  ensembl <- identifier_test_record(c("ENSG00000141510.18","ENSG00000139618.12"))
  before <- ensembl
  s <- suggest_feature_identifiers(ensembl,ref)
  testthat::expect_identical(s$suggested_type,"gene_id")
  testthat::expect_identical(s$versioned_features,2L)
  testthat::expect_true(s$strip_version_suggested)
  testthat::expect_identical(ensembl,before)
  testthat::expect_identical(suggest_feature_identifiers(identifier_test_record(c("GeneA","Unknown")),ref)$suggested_type,"gene_name")
  testthat::expect_identical(suggest_feature_identifiers(identifier_test_record("GeneA"))$suggested_type,"")
  for (ids in list(c("ENSG1","GeneA"),c("ENSG1","ENSMUSG2"),c("ENST1","ENST2"),c("ENSP1","ENSP2"),c("10","20"),"Unknown"))
    testthat::expect_identical(suggest_feature_identifiers(identifier_test_record(ids),ref)$suggested_type,"")
  testthat::expect_identical(resolve_feature_identifier_choice(s)$mode,"automatic")
  manual <- resolve_feature_identifier_choice(s,paste0("gene_name:",s$token))
  testthat::expect_identical(manual,list(key="gene_name",mode="manual"))
  # A manual browser choice for a different feature list cannot override auto.
  stale <- resolve_feature_identifier_choice(s,"gene_name:old-upload-token")
  testthat::expect_identical(stale,list(key="gene_id",mode="automatic"))
})

testthat::test_that("automatic annotation handles symbols and versioned IDs together without renaming counts", {
  app_root <- normalizePath(file.path("..","..","site_data"))
  ref <- load_gene_annotation_cache("Homo sapiens",app_root,include_ribosomal=FALSE)
  counts <- Matrix::Matrix(matrix(c(9,10,80,10,20,70,0,1,100),nrow=3,
    dimnames=list(c("MT-ND1","HBB","A1BG"),c("below","boundary","above"))),sparse=TRUE)
  id_counts <- counts
  ids <- ref$table$gene_id[match(rownames(counts),ref$table$gene_name)]
  rownames(id_counts) <- paste0(ids,".77")
  records <- list(A=new_count_import(counts,"A"),B=new_count_import(id_counts,"B"))
  server <- function(input,output,session) {
    controller <- annotation_qc_server(input,output,session,function() records,app_root)
  }
  shiny::testServer(server,{
    session$setInputs(addMetadata=0,initFiltering=0,typeOrganism="Homo sapiens",annotationKey_1="auto",annotationKey_2="auto")
    testthat::expect_identical(unname(controller$identifier_options()$keys),c("gene_name","gene_id"))
    testthat::expect_true(controller$identifier_options()$strip_ensembl_version)
    testthat::expect_match(output$annotationVersionControl$html,"checked")
    session$setInputs(addMetadata=1);session$flushReact()
    before <- digest::digest(controller$prepared())
    # Browser acknowledgement of the suggested checkbox value is not a new
    # scientific decision and must not clear an already prepared result.
    do.call(session$setInputs,stats::setNames(list(TRUE),controller$version_control_id()))
    testthat::expect_identical(digest::digest(controller$prepared()),before)
    session$setInputs(unrelatedPlotOption="changed")
    testthat::expect_identical(digest::digest(controller$prepared()),before)
    session$setInputs(initFiltering=1)
    b <- controller$baseline()
    testthat::expect_identical(rownames(b$B),rownames(id_counts))
    testthat::expect_identical(SummarizedExperiment::assay(b$B,"counts"),id_counts[,2:3])
    fields <- setdiff(.qc_metric_names,"sample")
    testthat::expect_equal(baseline_qc_metrics(b$A)[,fields],baseline_qc_metrics(b$B)[,fields],tolerance=0)
    testthat::expect_identical(S4Vectors::metadata(b$B)$scrna_explorer_qc$config$identifier_decision$mode,"automatic")
  })
})

testthat::test_that("suffix opt-out and explicit identifier choices remain under user control", {
  app_root <- normalizePath(file.path("..","..","site_data"))
  record <- identifier_test_record(c("ENSG00000141510.18","ENSG00000139618.12"))
  server <- function(input,output,session) {
    imported <- shiny::reactiveVal(list(A=record))
    controller <- annotation_qc_server(input,output,session,imported,app_root)
  }
  shiny::testServer(server,{
    session$setInputs(typeOrganism="Homo sapiens",addMetadata=0,annotationKey_1="auto")
    do.call(session$setInputs,stats::setNames(list(FALSE),controller$version_control_id()))
    testthat::expect_false(controller$identifier_options()$strip_ensembl_version)
    session$setInputs(addMetadata=1);session$flushReact()
    testthat::expect_true(all(controller$prepared()$A$feature_annotation$annotation_status=="unmapped"))
    s <- controller$identifier_suggestions()$A
    session$setInputs(annotationKey_1=paste0("gene_name:",s$token))
    testthat::expect_identical(controller$identifier_options()$keys[[1]],"gene_name")
    testthat::expect_identical(controller$identifier_options()$identifier_modes[[1]],"manual")
    testthat::expect_null(controller$prepared())
    # A new feature list gets fresh defaults and ignores an old qualified choice.
    imported(list(A=identifier_test_record("ENSMUSG00000000001.2")))
    testthat::expect_identical(controller$identifier_options()$keys[[1]],"gene_id")
    testthat::expect_true(controller$identifier_options()$strip_ensembl_version)
  })
})

testthat::test_that("ribosomal feature counts include annotated zero-count rows and retain unknowns", {
  ref <- annotation_fixture();ref$table$entrezid <- as.list(as.character(1:5))
  ref$ribosomal <- list(status="available",species="Homo sapiens",groups=list(
    ribosomal_protein="3",cytosolic_ribosomal_protein="3",mitochondrial_ribosomal_protein=character()))
  counts <- Matrix::Matrix(matrix(c(0,100,10,0,200,20),nrow=3,
    dimnames=list(c("GeneA","MT-A","Unmapped"),c("cellA","cellB"))),sparse=TRUE)
  p <- prepare_sample_qc(new_count_import(counts),ref,"gene_name")
  r <- p$ribosomal_summary[p$ribosomal_summary$group=="ribosomal_protein",]
  testthat::expect_identical(r$annotated_features,1L)
  testthat::expect_identical(r$unresolved_features,1L)
  testthat::expect_true(p$feature_annotation$ribosomal_protein[1])
  testthat::expect_true(is.na(p$feature_annotation$ribosomal_protein[3]))
  testthat::expect_equal(p$metrics$nUMI,c(110,220))
})

testthat::test_that("mixed identifiers require an explicit user choice and preserve unmapped rows", {
  app_root <- normalizePath(file.path("..","..","site_data"))
  record <- identifier_test_record(c("ENSG1","GeneA"))
  server <- function(input,output,session) {
    controller <- annotation_qc_server(input,output,session,function() list(A=record),app_root)
  }
  shiny::testServer(server,{
    session$setInputs(typeOrganism="Homo sapiens",annotationKey_1="auto",addMetadata=0)
    session$setInputs(addMetadata=1);session$flushReact()
    testthat::expect_null(controller$prepared())
    testthat::expect_match(output$annotationIdentifierSummary,"Choose explicitly")
    token <- controller$identifier_suggestions()$A$token
    session$setInputs(annotationKey_1=paste0("gene_id:",token),addMetadata=2);session$flushReact()
    testthat::expect_equal(nrow(controller$prepared()$A$feature_annotation),2)
    testthat::expect_identical(controller$prepared()$A$feature_annotation$source_feature_id,c("ENSG1","GeneA"))
    testthat::expect_identical(controller$prepared()$A$config$identifier_decision$mode,"manual")
  })
})
