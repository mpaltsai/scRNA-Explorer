# Retain the old signature, now respecting supplied barcodes and primary feature
# IDs. Symbols are annotations, not a reason to merge distinct feature rows.
my_read10X <- function(mtxPath,genesPath,barcodesPath=NULL,DGEList=FALSE) {
  records <- list(sample=read_count_10x(mtxPath,genesPath,barcodesPath))
  result <- import_counts(records)
  attr(result,"import_records") <- records
  result
}
