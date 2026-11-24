source(file.path("..","..","site_data","R","input_readers.R"))
source(file.path("..","..","site_data","R","input_import_server.R"))

input_fixture <- function() {
  root <- tempfile("input-fixture-"); dir.create(root)
  counts <- Matrix::Matrix(matrix(c(99,1,0,100,0,1),nrow=3,
    dimnames=list(c("ENSG1","ENSG2","ENSG3"),c("AAAC-1","TTTG-1"))),sparse=TRUE)
  list(root=root,counts=counts)
}

write_10x_fixture <- function(directory, counts, gzip=FALSE, barcodes=TRUE) {
  dir.create(directory,recursive=TRUE,showWarnings=FALSE)
  Matrix::writeMM(counts,file.path(directory,"matrix.mtx"))
  utils::write.table(data.frame(rownames(counts),c("Shared","Shared","Other")),
    file.path(directory,"features.tsv"),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
  if (barcodes) writeLines(colnames(counts),file.path(directory,"barcodes.tsv"))
  if (gzip) for (path in list.files(directory,full.names=TRUE)) {
    bytes <- readBin(path,"raw",n=file.info(path)$size)
    con <- gzfile(paste0(path,".gz"),"wb"); writeBin(bytes,con); close(con)
    unlink(path)
  }
}

make_input_zip <- function(directory, archive) {
  old <- setwd(directory); on.exit(setwd(old))
  utils::zip(archive,list.files(".",recursive=TRUE),flags="-q")
  archive
}
