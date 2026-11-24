# Compatibility for already extracted directories. The Shiny path instead uses
# read_count_zip() with session-owned extraction and path validation.
readZip <- function(dirs) {
  candidates <- dirs[vapply(dirs,function(d) any(grepl("(matrix\\.mtx(\\.gz)?$|\\.(rds|csv|tsv)(\\.gz)?$)",list.files(d),ignore.case=TRUE)),logical(1))]
  if (!length(candidates) || anyDuplicated(basename(candidates))) stop("Provide unique sample directories containing supported counts.")
  records <- stats::setNames(lapply(candidates,read_count_directory),basename(candidates))
  result <- import_counts(records)
  attr(result,"import_records") <- records
  result
}
