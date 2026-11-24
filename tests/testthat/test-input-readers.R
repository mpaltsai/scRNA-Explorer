testthat::test_that("text readers preserve names, feature column and exact counts", {
  f <- input_fixture(); on.exit(unlink(f$root,recursive=TRUE))
  for (sep in c(",","\t")) {
    p <- file.path(f$root,"upload-without-extension")
    utils::write.table(data.frame(feature=rownames(f$counts),as.matrix(f$counts),check.names=FALSE),p,
      sep=sep,quote=TRUE,row.names=FALSE)
    r <- read_count_table(p)
    testthat::expect_equal(as.matrix(r$counts),as.matrix(f$counts))
    testthat::expect_identical(r$status,"ready")
    testthat::expect_false(any(r$cell_identity$generated))
  }
})

testthat::test_that("invalid UMI values and cell identifiers are not silently repaired", {
  f <- input_fixture(); on.exit(unlink(f$root,recursive=TRUE))
  for (bad in list(-1,NA_real_,Inf,0.5)) {
    x <- as.matrix(f$counts); x[1,1] <- bad
    testthat::expect_error(new_count_import(x),"UMI counts|Fractional")
  }
  x <- as.matrix(f$counts); x[1,1] <- "text"
  testthat::expect_error(new_count_import(x),"numeric")
  x <- f$counts; colnames(x) <- c("same","same")
  testthat::expect_error(new_count_import(x),"unique within")
  x <- f$counts; rownames(x)[1] <- ""
  testthat::expect_error(new_count_import(x),"primary identifier")
  x <- f$counts; x[,] <- 0
  testthat::expect_silent(new_count_import(x))
})

testthat::test_that("numeric-looking identifiers keep their original spelling", {
  f <- input_fixture(); on.exit(unlink(f$root,recursive=TRUE))
  p <- file.path(f$root,"numeric.tsv")
  writeLines(c("gene\t0001\t0002","001\t0\t100","002\t99\t1"),p)
  r <- read_count_table(p)
  testthat::expect_identical(rownames(r$counts),c("001","002"))
  testthat::expect_identical(colnames(r$counts),c("0001","0002"))
})

testthat::test_that("compressed 10X uses primary IDs and supplied barcodes", {
  f <- input_fixture(); on.exit(unlink(f$root,recursive=TRUE))
  for (compress in c(FALSE,TRUE)) {
    d <- file.path(f$root,paste0("sample-",compress)); write_10x_fixture(d,f$counts,compress)
    r <- read_count_directory(d)
    testthat::expect_equal(as.matrix(r$counts),as.matrix(f$counts))
    testthat::expect_identical(r$feature_data$gene_name,c("Shared","Shared","Other"))
    testthat::expect_identical(r$status,"ready")
    testthat::expect_false(any(r$cell_identity$generated))
  }
  # Shiny temp files lose their original suffix, including .gz.
  paths <- list.files(d,full.names=TRUE)
  renamed <- file.path(f$root,paste0("upload",seq_along(paths)))
  file.copy(paths,renamed)
  r <- read_count_10x(renamed[grepl("matrix",paths)],renamed[grepl("features",paths)],renamed[grepl("barcodes",paths)])
  testthat::expect_equal(as.matrix(r$counts),as.matrix(f$counts))
})

testthat::test_that("missing barcodes generate reproducible, reversible sample-qualified IDs", {
  f <- input_fixture(); on.exit(unlink(f$root,recursive=TRUE))
  d <- file.path(f$root,"sample"); write_10x_fixture(d,f$counts,barcodes=FALSE)
  a <- read_count_directory(d,"sample:A"); b <- read_count_directory(d,"sample:B")
  testthat::expect_true(all(a$cell_identity$generated))
  testthat::expect_identical(a$cell_identity,read_count_directory(d,"sample:A")$cell_identity)
  testthat::expect_length(intersect(a$cell_identity$internal_cell_id,b$cell_identity$internal_cell_id),0)
  testthat::expect_equal(unname(as.matrix(a$counts)),unname(as.matrix(f$counts)))
  for (id in a$cell_identity$internal_cell_id) {
    parts <- strsplit(id,"::",fixed=TRUE)[[1]]
    testthat::expect_identical(utils::URLdecode(parts[1]),"sample:A")
  }
})

testthat::test_that("10X rejects incompatible dimensions and mixed modalities", {
  f <- input_fixture(); on.exit(unlink(f$root,recursive=TRUE))
  d <- file.path(f$root,"sample"); write_10x_fixture(d,f$counts)
  writeLines("one",file.path(d,"barcodes.tsv"))
  testthat::expect_error(read_count_directory(d),"Barcode-file dimensions")
  writeLines(colnames(f$counts),file.path(d,"barcodes.tsv"))
  utils::write.table(data.frame(rownames(f$counts),letters[1:3],c("Gene Expression","Antibody Capture","Gene Expression")),
    file.path(d,"features.tsv"),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
  testthat::expect_error(read_count_directory(d),"Gene Expression features only")
})

testthat::test_that("SCE import preserves its object, assays and all metadata without densifying", {
  f <- input_fixture(); on.exit(unlink(f$root,recursive=TRUE))
  sce <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=f$counts,extra=f$counts))
  sce$geo_accession <- c("GSM1","GSM1")
  SummarizedExperiment::rowData(sce)$custom <- letters[1:3]
  S4Vectors::metadata(sce)$geo <- list(build="hg19",source="fixture")
  SingleCellExperiment::reducedDim(sce,"example") <- matrix(1:4,nrow=2)
  p <- file.path(f$root,"sce.rds"); saveRDS(sce,p)
  hash <- tools::md5sum(p); r <- read_count_rds(p)
  testthat::expect_identical(r$counts,f$counts)
  testthat::expect_identical(r$source_object,sce)
  testthat::expect_identical(r$cell_data,SummarizedExperiment::colData(sce))
  testthat::expect_identical(r$feature_data,SummarizedExperiment::rowData(sce))
  testthat::expect_identical(tools::md5sum(p),hash)
  saveRDS(list(sce),file.path(f$root,"list.rds"))
  testthat::expect_error(read_count_rds(file.path(f$root,"list.rds")),"single SingleCellExperiment")
})

testthat::test_that("duplicate primary IDs pause import and require a lossless explicit decision", {
  f <- input_fixture(); on.exit(unlink(f$root,recursive=TRUE))
  x <- f$counts; rownames(x) <- c("gene","gene","gene__import_row_1")
  r <- new_count_import(x)
  testthat::expect_identical(r$status,"needs_feature_decision")
  testthat::expect_equal(import_duplicate_rows(r)$source_row,c(1,2))
  testthat::expect_error(import_counts(list(sample=r)),"Resolve or cancel")
  testthat::expect_error(resolve_import_duplicates(r,"keep_separate"," "),"reason")
  testthat::expect_error(resolve_import_duplicates(r,"sum","same gene"),"Choose")
  resolved <- resolve_import_duplicates(r,"keep_separate","Retain the two source rows for review")
  testthat::expect_false(anyDuplicated(rownames(resolved$counts))>0)
  testthat::expect_identical(resolved$feature_identity$source_feature_id,rownames(x))
  testthat::expect_equal(unname(as.matrix(resolved$counts)),unname(as.matrix(x)))
  testthat::expect_identical(resolved$decisions[[1]]$reason,"Retain the two source rows for review")
  testthat::expect_identical(r$counts,x)
  cancelled <- resolve_import_duplicates(r,"cancel","I will correct the source file")
  testthat::expect_identical(cancelled$status,"cancelled")
  testthat::expect_error(import_counts(list(sample=cancelled)),"Resolve or cancel")
})

testthat::test_that("ZIP readers support modern 10X and tables without assuming nesting depth", {
  f <- input_fixture(); on.exit(unlink(f$root,recursive=TRUE))
  samples <- file.path(f$root,"source"); dir.create(samples)
  write_10x_fixture(file.path(samples,"nested","A"),f$counts,gzip=TRUE)
  dir.create(file.path(samples,"nested","B"))
  utils::write.table(data.frame(feature=rownames(f$counts),as.matrix(f$counts),check.names=FALSE),
    file.path(samples,"nested","B","counts.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
  archive <- make_input_zip(samples,file.path(f$root,"input.zip"))
  first <- file.path(f$root,"session1"); second <- file.path(f$root,"session2")
  dir.create(first); dir.create(second)
  a <- read_count_zip(archive,first); b <- read_count_zip(archive,second)
  testthat::expect_identical(names(a),c("A","B"))
  testthat::expect_equal(as.matrix(a$A$counts),as.matrix(f$counts))
  testthat::expect_equal(as.matrix(a$B$counts),as.matrix(f$counts))
  testthat::expect_false(identical(a$A$source$paths,b$A$source$paths))
  testthat::expect_true(all(startsWith(a$A$source$paths,first)))
  testthat::expect_true(all(startsWith(b$A$source$paths,second)))
})

testthat::test_that("ZIP input rejects traversal before extraction", {
  f <- input_fixture(); on.exit(unlink(f$root,recursive=TRUE))
  directory <- file.path(f$root,"inside"); dir.create(directory)
  writeLines("outside",file.path(f$root,"outside.tsv"))
  old <- setwd(directory); on.exit(setwd(old),add=TRUE)
  archive <- file.path(f$root,"unsafe.zip")
  utils::zip(archive,"../outside.tsv",flags="-q")
  testthat::expect_error(read_count_zip(archive,directory),"unsafe")
  testthat::expect_length(list.files(directory,all.files=TRUE,no..=TRUE),0)
})

testthat::test_that("legacy entry points retain named sample lists and provide import audits", {
  source(file.path("..","..","site_data","R","readData.R"),local=TRUE)
  source(file.path("..","..","site_data","R","my_read10X.R"),local=TRUE)
  source(file.path("..","..","site_data","R","read_zip_dir.R"),local=TRUE)
  f <- input_fixture(); on.exit(unlink(f$root,recursive=TRUE))
  d <- file.path(f$root,"A"); write_10x_fixture(d,f$counts)
  for (result in list(readData(FALSE,FALSE,NULL,d),my_read10X(file.path(d,"matrix.mtx"),file.path(d,"features.tsv"),file.path(d,"barcodes.tsv")))) {
    testthat::expect_identical(names(result),"sample")
    testthat::expect_equal(as.matrix(result$sample),as.matrix(f$counts))
    testthat::expect_true(!is.null(attr(result,"import_records")))
  }
  testthat::expect_identical(names(readZip(c(f$root,d))),"A")
})
