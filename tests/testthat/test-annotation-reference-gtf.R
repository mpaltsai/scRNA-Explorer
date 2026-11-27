testthat::test_that("GRCh37 GTF and GRCh38 EnsDb references are discovered without mixing builds", {
  f <- reference_cache_fixture(); on.exit(unlink(f$root,recursive=TRUE))
  gtf <- f$rows[1,,drop=FALSE]
  gtf$resource_id <- "AH10684"; gtf$title <- "Homo_sapiens.GRCh37.75.gtf"
  gtf$rdataclass <- "GRanges"; gtf$sourcetype <- "GTF"
  other <- gtf; other$resource_id <- "AH999"; other$title <- "Unrelated peaks"; other$sourcetype <- "BED"
  rows <- f$rows; rows$sourcetype <- "Ensembl"
  rows <- rbind(rows,gtf,other)
  catalog <- function(...) list(rows=rows,snapshot_date="2026-04-23",annotationhub_version="test")
  human37 <- f$api$check_annotation_updates("Homo sapiens","hg19",f$app_root,f$root,TRUE,catalog_fetch=catalog)
  human38 <- f$api$check_annotation_updates("Homo sapiens","hg38",f$app_root,f$root,TRUE,catalog_fetch=catalog)
  testthat::expect_setequal(human37$candidates$resource_id,c("AH101","AH102","AH10684"))
  testthat::expect_identical(human38$candidates$resource_id,"AH103")
  testthat::expect_true(all(human37$candidates$genome=="GRCh37"))
  testthat::expect_false(dir.exists(file.path(f$root,"references")))
})

testthat::test_that("GTF export keeps conflicting gene mappings and records unavailable fields", {
  f <- reference_cache_fixture(); on.exit(unlink(f$root,recursive=TRUE))
  candidate <- f$rows[1,,drop=FALSE]
  candidate$resource_id <- "AH10684"; candidate$title <- "Homo_sapiens.GRCh37.75.gtf"
  candidate$rdataclass <- "GRanges"; candidate$sourcetype <- "GTF"
  gr <- GenomicRanges::GRanges(c("1","1","1","MT"),IRanges::IRanges(1:4,width=1))
  gr$gene_id <- c("ENSG1","ENSG1","ENSG1","ENSG2")
  gr$gene_name <- c("A","A","Different","MT-A")
  gr$gene_biotype <- c(rep("protein_coding",3),"Mt_rRNA")
  GenomeInfoDb::genome(gr) <- "GRCh37"
  payload <- f$api$annotation_gtf_table(gr,candidate)
  testthat::expect_equal(nrow(payload$table),3)
  testthat::expect_identical(payload$annotation_release,"75")
  testthat::expect_true(all(is.na(payload$table$entrezid)))
  testthat::expect_true(all(is.na(payload$table$description)))
  testthat::expect_setequal(payload$table$gene_name[payload$table$gene_id=="ENSG1"],c("A","Different"))
  entry <- f$api$download_annotation_reference(candidate,"2026-04-23","Review GRCh37",f$app_root,f$root,
    resource_fetch=function(...) payload)
  testthat::expect_match(entry$label,"GRCh37.*Ensembl 75.*AH10684.*GTF")
  ref <- load_gene_annotation_cache("Homo sapiens",f$app_root,FALSE,entry$reference_id,f$root)
  testthat::expect_identical(ref$table,payload$table)
  testthat::expect_identical(ref$provenance$source,"AnnotationHub Ensembl GTF")
  testthat::expect_match(ref$provenance$source_metadata$build_evidence,"agree")
  counts <- Matrix::Matrix(matrix(c(100,20,80,30),nrow=2,
    dimnames=list(c("ENSG1","ENSG2"),c("cellA","cellB"))),sparse=TRUE)
  prepared <- prepare_sample_qc(new_count_import(counts),ref,"gene_id")
  testthat::expect_identical(prepared$summary$hemoglobin_status,"missing_reference_definition")
  testthat::expect_true(all(is.na(prepared$metrics$rbcRatio)))
  testthat::expect_equal(prepared$metrics$nUMI,c(120,110))
  server <- function(input,output,session) {
    qc <- annotation_qc_server(input,output,session,function() list(sample=new_count_import(counts)),f$app_root,
      reference_services=list(cache_dir=f$root))
  }
  shiny::testServer(server,{
    session$setInputs(addMetadata=0,typeOrganism="Homo sapiens",annotationKey_1="gene_id",
      annotationDataBuild="GRCh37",annotationReference=entry$reference_id)
    session$setInputs(addMetadata=1);session$flushReact()
    testthat::expect_false(is.null(qc$prepared()))
    testthat::expect_match(output$paramsAddMetadata,"Hemoglobin QC is unavailable")
    testthat::expect_match(output$paramsAddMetadata,"QC filtering remains blocked")
  })
  GenomeInfoDb::genome(gr) <- "GRCh38"
  testthat::expect_error(f$api$annotation_gtf_table(gr,candidate),"differs")
  GenomeInfoDb::genome(gr) <- NA_character_
  testthat::expect_match(f$api$annotation_gtf_table(gr,candidate)$source_metadata$build_evidence,"unavailable")
  gr$gene_biotype <- NULL
  testthat::expect_error(f$api$annotation_gtf_table(gr,candidate),"lacks")
})
