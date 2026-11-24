testthat::test_that("ribosomal protein flags use reference membership rather than broad names", {
  ref <- annotation_fixture()
  ref$table <- data.frame(gene_id=paste0("ID",1:8),
    gene_name=c("RPL3","FAU","MRPL12","RPS6KA1","RNA18S","MT-RNR1","RNA18SP","Unknown"),
    gene_biotype=c(rep("protein_coding",4),"rRNA","Mt_rRNA","rRNA_pseudogene","protein_coding"),
    seq_name=c(rep("1",5),"MT","1","1"),description=rep("test annotation",8))
  ref$table$entrezid <- list("1","2","3","4","5","6","7",NA_character_)
  ref$ribosomal <- list(status="available",species="Homo sapiens",groups=list(
    ribosomal_protein=c("1","2","3"),cytosolic_ribosomal_protein=c("1","2"),mitochondrial_ribosomal_protein="3"))
  counts <- Matrix::Matrix(matrix(1:16,nrow=8,dimnames=list(ref$table$gene_name,c("cell1","cell2"))),sparse=TRUE)
  r <- new_count_import(counts)
  p <- prepare_sample_qc(r,ref,"gene_name")
  a <- p$feature_annotation
  testthat::expect_identical(a$ribosomal_protein,c(TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,FALSE,NA))
  testthat::expect_true(a$cytosolic_ribosomal_protein[2]) # FAU lacks an RPL/RPS prefix.
  testthat::expect_true(a$mitochondrial_ribosomal_protein[3])
  testthat::expect_false(a$ribosomal_protein[4]) # RPS6 kinase is not a structural ribosomal protein.
  testthat::expect_true(a$nuclear_rRNA[5])
  testthat::expect_true(a$mitochondrial_rRNA[6])
  testthat::expect_true(a$rRNA_pseudogene[7])
  testthat::expect_false(a$ribosomal_protein[7])
  testthat::expect_identical(a$ribosomal_annotation_status[8],"unresolved")
  testthat::expect_equal(p$metrics$nUMI,as.numeric(Matrix::colSums(counts)))
  testthat::expect_identical(r$counts,counts)
  rownames(r$counts) <- ref$table$gene_id
  r$feature_identity$source_feature_id <- rownames(r$counts)
  r$feature_identity$internal_feature_id <- rownames(r$counts)
  q <- prepare_sample_qc(r,ref,"gene_id")
  testthat::expect_identical(q$feature_annotation[, .qc_ribosomal_fields],a[, .qc_ribosomal_fields])
  b <- make_baseline_sce(r,q)
  for (g in .qc_ribosomal_fields) testthat::expect_identical(SummarizedExperiment::rowData(b)[[paste0("qc_",g)]],a[[g]])
  testthat::expect_identical(rownames(b),rownames(r$counts))
})

testthat::test_that("mixed ribosomal candidates remain unknown and missing resources are explicit", {
  ref <- annotation_fixture()
  ref$table$entrezid <- list("1","2","3","4","5")
  ref$ribosomal <- list(status="available",species="Homo sapiens",groups=list(
    ribosomal_protein="4",cytosolic_ribosomal_protein="4",mitochondrial_ribosomal_protein=character()))
  counts <- Matrix::Matrix(matrix(c(100,100),nrow=1,dimnames=list("Shared",c("a","b"))),sparse=TRUE)
  p <- prepare_sample_qc(new_count_import(counts),ref,"gene_name")
  testthat::expect_true(is.na(p$feature_annotation$ribosomal_protein))
  testthat::expect_identical(p$feature_annotation$ribosomal_annotation_status,"unresolved")
  ref$ribosomal <- list(status="not_assessed",species="Homo sapiens",reason="Reference unavailable")
  p <- prepare_sample_qc(new_count_import(counts),ref,"gene_name")
  testthat::expect_true(is.na(p$feature_annotation$ribosomal_protein))
  testthat::expect_identical(p$config$ribosomal_reference$status,"not_assessed")
  ref$ribosomal$species <- "Mus musculus"
  testthat::expect_error(prepare_sample_qc(new_count_import(counts),ref,"gene_name"),"species")
})

testthat::test_that("existing human and mouse GO resources retain versions, sets and controls", {
  for (species in c("Homo sapiens","Mus musculus")) {
    package <- if (species=="Homo sapiens") "org.Hs.eg.db" else "org.Mm.eg.db"
    testthat::skip_if_not_installed(package)
    ref <- load_ribosomal_reference(species)
    testthat::expect_identical(ref$status,"available")
    testthat::expect_identical(ref$species,species)
    testthat::expect_true(length(ref$groups$ribosomal_protein)>50)
    testthat::expect_true(length(ref$groups$mitochondrial_ribosomal_protein)>20)
    testthat::expect_identical(ref$groups_hash,digest::digest(ref$groups,algo="sha256"))
    testthat::expect_true("GOSOURCEDATE" %in% ref$source_metadata$name)
    symbols <- if (species=="Homo sapiens") c("RPL3","FAU","MRPL12","RPS6KA1") else c("Rpl3","Fau","Mrpl12","Rps6ka1")
    db <- getExportedValue(package,package)
    ids <- suppressMessages(AnnotationDbi::select(db,keys=symbols,keytype="SYMBOL",columns="ENTREZID"))
    flag <- ids$ENTREZID %in% ref$groups$ribosomal_protein
    testthat::expect_true(all(flag[ids$SYMBOL %in% symbols[1:3]]))
    testthat::expect_false(any(flag[ids$SYMBOL==symbols[4]]))
  }
})
