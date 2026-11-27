biomart_supplement_fixture <- function() {
  f <- reference_cache_fixture()
  candidate <- f$rows[1,,drop=FALSE]
  candidate$resource_id <- "AH10684";candidate$title <- "Homo_sapiens.GRCh37.75.gtf"
  candidate$rdataclass <- "GRanges";candidate$sourcetype <- "GTF"
  table <- annotation_fixture()$table
  table$description <- c(NA,NA,NA,"Original description",NA)
  table$entrezid <- c(NA,NA,NA,"99",NA)
  parent <- f$api$download_annotation_reference(candidate,"2026-04-23","Fixture GRCh37 parent",f$app_root,f$root,
    resource_fetch=function(...) list(table=table,species="Homo sapiens",genome_build="GRCh37",annotation_release="75",source_metadata=list()))
  payload <- list(metadata=list(source="Ensembl GRCh37 BioMart",species="Homo sapiens",genome_build="GRCh37",
    annotation_release="75",host="https://grch37.ensembl.org:443/biomart/martservice",
    dataset="hsapiens_gene_ensembl",dataset_version="GRCh37.p13",mart_version="Ensembl Genes 75",
    release_evidence_url="https://jun2026.archive.ensembl.org/info/website/tutorials/grch37.html",
    retrieved_at="2026-11-27 12:00:00 UTC",biomart_version="test",attributes=c("ensembl_gene_id","description","entrezgene_id")),
    rows=data.frame(gene_id=c("ENSG1","ENSG1","ENSG2","ENSG3","ENSG3","ENSG4","ENSG9"),
      description=c("Mitochondrial gene","Mitochondrial gene","hemoglobin alpha","First","Second","Changed description","Unused gene"),
      entrezid=c("10","11","20","30","30","100","90")))
  f$parent <- parent;f$payload <- payload;f$table <- table
  f
}

testthat::test_that("exact gene-ID supplements fill gaps and retain conflicts and multiple Entrez links", {
  f <- biomart_supplement_fixture();on.exit(unlink(f$root,recursive=TRUE))
  ref <- f$api$load_versioned_annotation(f$parent$reference_id,"Homo sapiens",f$app_root,f$root)
  before <- ref
  result <- f$api$merge_grch37_biomart(ref,f$payload)
  testthat::expect_identical(ref,before)
  testthat::expect_identical(result$table[,1:4],f$table[,1:4])
  testthat::expect_equal(nrow(result$table),5)
  testthat::expect_identical(result$table$description,c("Mitochondrial gene","hemoglobin alpha",NA,"Original description",NA))
  testthat::expect_identical(result$table$entrezid[[1]],c("10","11"))
  testthat::expect_identical(result$table$entrezid[[4]],"99")
  testthat::expect_length(result$table$entrezid[[5]],0)
  testthat::expect_identical(result$review$description_status[3],"conflicting_unresolved")
  testthat::expect_identical(result$review$description_status[4],"existing_conflict_preserved")
  testthat::expect_identical(result$review$entrezid_status[1],"filled_multiple")
  testthat::expect_identical(result$summary$descriptions_filled,2L)
  bad <- f$payload;bad$metadata$genome_build <- "GRCh38"
  testthat::expect_error(f$api$merge_grch37_biomart(ref,bad),"metadata must verify")
  bad <- f$payload;bad$metadata$annotation_release <- "116"
  testthat::expect_error(f$api$merge_grch37_biomart(ref,bad),"metadata must verify")
  updated <- f$payload;updated$metadata$mart_version <- "Ensembl Genes 116"
  testthat::expect_identical(f$api$merge_grch37_biomart(ref,updated)$table,result$table)
  bad <- f$payload;bad$rows$entrezid[1] <- "not-an-id"
  testthat::expect_error(f$api$merge_grch37_biomart(ref,bad),"invalid Entrez")
})

testthat::test_that("supplements are immutable separate references with verified audit downloads", {
  f <- biomart_supplement_fixture();on.exit(unlink(f$root,recursive=TRUE))
  parent_hash <- tools::md5sum(f$parent$path)
  add <- function(fetch=function() f$payload) f$api$supplement_grch37_reference(
    f$parent$reference_id,"Restore missing GRCh37 fields",f$app_root,f$root,fetch)
  original <- names(f$api$list_annotation_references("Homo sapiens",f$app_root,f$root))
  testthat::expect_error(add(function() stop("HTTP 429")),"HTTP 429")
  testthat::expect_identical(names(f$api$list_annotation_references("Homo sapiens",f$app_root,f$root)),original)
  entry <- add()
  testthat::expect_false(identical(entry$reference_id,f$parent$reference_id))
  testthat::expect_match(entry$label,"GRCh37.*GTF \\+ BioMart")
  testthat::expect_identical(tools::md5sum(f$parent$path),parent_hash)
  ref <- f$api$load_versioned_annotation(entry$reference_id,"Homo sapiens",f$app_root,f$root)
  testthat::expect_identical(ref$provenance$parent_reference$reference_id,f$parent$reference_id)
  testthat::expect_identical(ref$provenance$supplement$annotation_release,"75")
  testthat::expect_identical(ref$provenance$supplement$summary$entrez_multiple,1L)
  review <- f$api$load_annotation_supplement_review(entry$reference_id,f$app_root,f$root)
  testthat::expect_identical(review$description_status[3],"conflicting_unresolved")
  repeated <- add(function() {p<-f$payload;p$metadata$retrieved_at<-"2026-12-01 00:00:00 UTC";p})
  testthat::expect_identical(repeated$reference_id,entry$reference_id)
  testthat::expect_error(f$api$supplement_grch37_reference(entry$reference_id,"Nested",f$app_root,f$root),"unsupplemented")
  writeLines("changed",file.path(dirname(entry$path),"supplement-review.csv"))
  testthat::expect_error(f$api$load_annotation_supplement_review(entry$reference_id,f$app_root,f$root),"checksum mismatch")
  testthat::expect_error(add(),"checksum mismatch")
})

testthat::test_that("Shiny requires supplement consent and retains pinned QC until explicit selection", {
  f <- biomart_supplement_fixture();on.exit(unlink(f$root,recursive=TRUE))
  counts <- Matrix::Matrix(matrix(c(50,10,40,0,0,80,20,100,0,0),nrow=5,
    dimnames=list(f$table$gene_name,c("cellA","cellB"))),sparse=TRUE)
  # Shared source symbols need the normal explicit import decision first.
  record <- resolve_import_duplicates(new_count_import(counts),"keep_separate","Preserve every source feature")
  calls <- 0L
  fail_fetch <- FALSE
  server <- function(input,output,session) {
    qc <- annotation_qc_server(input,output,session,function() list(sample=record),f$app_root,
      reference_services=list(cache_dir=f$root,supplement_fetch=function() {
        calls<<-calls+1L
        if (fail_fetch) stop("HTTP 429")
        f$payload
      }))
  }
  shiny::testServer(server,{
    session$setInputs(addMetadata=0,initFiltering=0,annotationSupplement=0,typeOrganism="Homo sapiens",
      annotationDataBuild="GRCh37",annotationReference=f$parent$reference_id,annotationKey_1="gene_name")
    session$setInputs(addMetadata=1);session$flushReact();session$setInputs(initFiltering=1)
    baseline <- digest::digest(qc$baseline())
    session$setInputs(annotationSupplementParent=f$parent$reference_id,annotationSupplement=1)
    session$flushReact();testthat::expect_identical(calls,0L)
    testthat::expect_match(output$annotationUpdateStatus,"confirm BioMart retrieval")
    session$setInputs(annotationSupplementConsent=TRUE,annotationSupplementReason="Restore hemoglobin and Entrez evidence",annotationSupplement=2)
    session$flushReact();testthat::expect_identical(calls,1L)
    testthat::expect_identical(digest::digest(qc$baseline()),baseline)
    entries <- f$api$list_annotation_references("Homo sapiens",f$app_root,f$root)
    id <- names(entries)[vapply(entries,function(e) !is.null(e$provenance$supplement),logical(1))]
    testthat::expect_length(id,1)
    session$setInputs(annotationReference=id)
    testthat::expect_null(qc$baseline())
    session$setInputs(addMetadata=2);session$flushReact();session$setInputs(initFiltering=2)
    testthat::expect_true(all(is.finite(qc$prepared()$sample$metrics$rbcRatio)))
    testthat::expect_match(output$annotationProvenance,"Ensembl GRCh37 BioMart")
    testthat::expect_identical(SummarizedExperiment::assay(qc$baseline()$sample,"counts"),record$counts)
    testthat::expect_match(output$annotationSupplementReviewControls$html,"Download BioMart supplement review")
    current_baseline <- digest::digest(qc$baseline())
    fail_fetch <<- TRUE
    session$setInputs(annotationSupplement=3);session$flushReact()
    testthat::expect_match(output$annotationUpdateStatus,"HTTP 429")
    testthat::expect_identical(digest::digest(qc$baseline()),current_baseline)
    testthat::expect_null(output$annotationUpdateBusy)
  })
})

testthat::test_that("BioMart XML retrieval verifies the dataset and refuses partial query responses", {
  f <- reference_cache_fixture();on.exit(unlink(f$root,recursive=TRUE))
  response <- function(parameters) {
    if (!is.null(parameters$query)) return("ENSG1\themoglobin alpha\t10\n[success]\n")
    switch(parameters$type,
      registry='<MartRegistry><MartURLLocation name="ENSEMBL_MART_ENSEMBL" host="grch37.ensembl.org" displayName="Ensembl Genes 116"/></MartRegistry>',
      datasets=" \nTableSet\thsapiens_gene_ensembl\tHuman genes (GRCh37.p13)\t1\tGRCh37.p13",
      attributes="ensembl_gene_id\tID\tLong ID description\ndescription\tDescription\nentrezgene_id\tEntrez")
  }
  p <- f$api$grch37_biomart_fetch(response)
  testthat::expect_identical(p$metadata$dataset_version,"GRCh37.p13")
  testthat::expect_identical(p$metadata$annotation_release,"75")
  testthat::expect_identical(p$metadata$mart_version,"Ensembl Genes 116")
  testthat::expect_identical(p$rows$description,"hemoglobin alpha")
  testthat::expect_identical(p$rows$entrezid,"10")
  testthat::expect_error(f$api$grch37_biomart_fetch(function(q) {
    if (!is.null(q$query)) return("ENSG1\tpartial\t10\n")
    response(q)
  }),"incomplete or unsuccessful")
  testthat::expect_error(f$api$grch37_biomart_fetch(function(q) {
    if (identical(q$type,"datasets")) return("TableSet\thsapiens_gene_ensembl\tHuman genes\t1\tGRCh38.p14")
    response(q)
  }),"metadata must verify")
  testthat::expect_error(f$api$grch37_biomart_fetch(function(q) {
    if (identical(q$type,"registry")) return(sub("grch37.ensembl.org","www.ensembl.org",response(q),fixed=TRUE))
    response(q)
  }),"registry host")
})
