# Compatibility entry point for existing callers. New Shiny code retains full
# import records directly; callers of this wrapper can inspect import_records.
readData <- function(count.matrix,seurat,inputFile,inputDir) {
  record <- if (!count.matrix) read_count_directory(inputDir,"sample") else if (seurat) {
    read_count_rds(inputFile)
  } else read_count_table(inputFile)
  records <- list(sample=record)
  result <- import_counts(records)
  attr(result,"import_records") <- records
  result
}
