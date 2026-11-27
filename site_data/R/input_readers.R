# Import counts without gene annotation, filtering, normalization or disk checkpoints.
# An import record carries the sparse counts plus source identifiers/metadata. The
# original SCE is retained when supplied; later stages can reuse it directly.
validate_umi_counts <- function(counts) {
  if (!is.matrix(counts) && !inherits(counts,"sparseMatrix")) stop("Supply a numeric matrix of UMI counts.")
  if (length(dim(counts))!=2L || any(dim(counts)==0L)) stop("The count matrix must contain features and cells.")
  values <- if (inherits(counts,"sparseMatrix")) methods::as(methods::as(counts,"CsparseMatrix"),"dgCMatrix")@x else counts
  if (!is.numeric(values) || anyNA(values) || any(!is.finite(values)) || any(values<0)) {
    stop("UMI counts must be numeric, finite and nonnegative, with no missing values.")
  }
  if (any(values!=floor(values))) stop("Fractional values are not raw UMI counts; no rounding is performed.")
  invisible(TRUE)
}

new_count_import <- function(counts, sample_id="sample", feature_data=NULL, cell_data=NULL,
                             object=NULL, source=list()) {
  validate_umi_counts(counts)
  if (length(sample_id)!=1L || is.na(sample_id) || !nzchar(sample_id)) stop("Supply a nonempty sample identifier.")
  features <- rownames(counts)
  if (is.null(features) || anyNA(features) || any(!nzchar(features))) stop("Every feature needs a nonempty primary identifier.")
  generated <- is.null(colnames(counts))
  cells <- if (generated) paste0("generated_cell_",seq_len(ncol(counts))) else colnames(counts)
  if (anyNA(cells) || any(!nzchar(cells)) || anyDuplicated(cells)) stop("Cell identifiers must be nonempty and unique within each sample.")
  if (is.null(feature_data)) feature_data <- data.frame(row.names=seq_len(nrow(counts)))
  if (is.null(cell_data)) cell_data <- data.frame(row.names=seq_len(ncol(counts)))
  if (nrow(feature_data)!=nrow(counts) || nrow(cell_data)!=ncol(counts)) stop("Imported metadata dimensions do not match counts.")
  # These audit fields live outside imported colData/rowData: conflicting source
  # fields are never overwritten. Internal IDs are reversible across samples.
  identity <- data.frame(source_cell_id=cells,
    internal_cell_id=paste(utils::URLencode(sample_id,reserved=TRUE),
      vapply(cells,utils::URLencode,character(1),reserved=TRUE),sep="::"),
    generated=rep(generated,length(cells)),stringsAsFactors=FALSE)
  counts <- methods::as(methods::as(counts,"CsparseMatrix"),"dgCMatrix")
  colnames(counts) <- cells
  duplicates <- duplicated(features) | duplicated(features,fromLast=TRUE)
  source$sample_id_at_import <- sample_id
  list(counts=counts,sample_id=sample_id,feature_data=feature_data,cell_data=cell_data,
    feature_identity=data.frame(source_row=seq_along(features),source_feature_id=features,
      internal_feature_id=features,duplicate_primary_id=duplicates),cell_identity=identity,
    source_object=object,source=source,decisions=list(),sample_name_history=list(),
    status=if (any(duplicates)) "needs_feature_decision" else "ready")
}

# Sample names are user-facing labels. Renaming must not reconstruct cell IDs,
# change source metadata or alter counts. Apply the complete set atomically so
# duplicate/empty labels cannot create ambiguous per-sample QC state.
rename_import_samples <- function(records,sample_names) {
  if (!length(records) || !is.character(sample_names) || length(sample_names)!=length(records) || anyNA(sample_names))
    stop("Provide a name for every imported sample.")
  sample_names <- trimws(sample_names)
  if (any(!nzchar(sample_names))) stop("Sample names cannot be empty.")
  if (anyDuplicated(sample_names)) stop("Sample names must be unique.")
  if (any(sample_names %in% c(".","..")) || any(grepl("[[:cntrl:]/\\\\]",sample_names)))
    stop("Sample names cannot contain path separators or control characters.")
  for (i in seq_along(records)) {
    old <- records[[i]]$sample_id
    if (identical(old,sample_names[[i]])) next
    records[[i]]$sample_name_history <- c(records[[i]]$sample_name_history,
      list(list(from=old,to=sample_names[[i]],time=format(Sys.time(),tz="UTC",usetz=TRUE))))
    records[[i]]$sample_id <- sample_names[[i]]
  }
  stats::setNames(records,sample_names)
}

# Plain text may be gzipped even when Shiny's temporary filename has no suffix.
# gzfile reads both compressed and uncompressed files. The first column contains
# feature IDs; column headers after it are cell IDs. Do not let read.table repair names.
.import_table <- function(path, sep="\t", header=FALSE) {
  con <- gzfile(path,open="rt"); on.exit(close(con))
  utils::read.table(con,sep=sep,header=header,check.names=FALSE,stringsAsFactors=FALSE,
    quote="\"",comment.char="",fill=FALSE,colClasses="character")
}

read_count_table <- function(path, sample_id="sample") {
  con <- gzfile(path,open="rt"); first <- readLines(con,n=1L,warn=FALSE); close(con)
  sep <- if (grepl("\t",first,fixed=TRUE)) "\t" else ","
  d <- .import_table(path,sep=sep,header=TRUE)
  if (ncol(d)<2L) stop("Count tables need a feature-ID column followed by cell columns.")
  # Preserve IDs such as 001 as text; only expression columns are numeric.
  columns <- lapply(d[,-1,drop=FALSE],function(x) suppressWarnings(as.numeric(x)))
  counts <- do.call(cbind,columns); rownames(counts) <- d[[1]]
  new_count_import(counts,sample_id,source=list(format="count_table",path=normalizePath(path),md5=unname(tools::md5sum(path))))
}

# Read a single SCE through public accessors. Retain its other assays, reductions,
# colData, rowData and metadata in source_object rather than rebuilding it here.
read_count_rds <- function(path, sample_id="sample") {
  object <- readRDS(path)
  if (!inherits(object,"SingleCellExperiment")) stop("RDS import currently requires a single SingleCellExperiment with a counts assay; arbitrary Seurat objects and sample lists are not supported.")
  if (!"counts" %in% SummarizedExperiment::assayNames(object)) stop("The SingleCellExperiment has no counts assay.")
  new_count_import(SummarizedExperiment::assay(object,"counts"),sample_id,
    feature_data=SummarizedExperiment::rowData(object),cell_data=SummarizedExperiment::colData(object),
    object=object,source=list(format="SingleCellExperiment",path=normalizePath(path),md5=unname(tools::md5sum(path))))
}

# Feature column 1 is the primary ID; column 2 is a symbol annotation. Identical
# symbols do not imply identical genes. All three file paths may be gzipped.
read_count_10x <- function(matrix_path, features_path, barcodes_path=NULL, sample_id="sample") {
  con <- gzfile(matrix_path,open="rt"); on.exit(close(con))
  counts <- Matrix::readMM(con)
  f <- .import_table(features_path)
  if (nrow(f)!=nrow(counts) || !ncol(f) %in% c(1L,2L,3L)) stop("Feature-file dimensions do not match the matrix.")
  # Multimodal inputs need an explicit assay-selection workflow, not silent mixing.
  if (ncol(f)==3L && any(f[[3]]!="Gene Expression")) stop("This reader accepts Gene Expression features only; separate other feature types before import.")
  rownames(counts) <- as.character(f[[1]])
  feature_data <- data.frame(gene_id=as.character(f[[1]]),gene_name=if(ncol(f)>=2L) as.character(f[[2]]) else NA_character_)
  if (!is.null(barcodes_path) && length(barcodes_path) && nzchar(barcodes_path)) {
    b <- .import_table(barcodes_path)
    if (ncol(b)!=1L || nrow(b)!=ncol(counts)) stop("Barcode-file dimensions do not match the matrix.")
    colnames(counts) <- as.character(b[[1]])
  }
  paths <- c(matrix_path,features_path,barcodes_path)
  new_count_import(counts,sample_id,feature_data=feature_data,
    source=list(format="10X",paths=normalizePath(paths),md5=unname(tools::md5sum(paths))))
}

read_count_directory <- function(directory, sample_id=basename(directory)) {
  files <- list.files(directory,full.names=TRUE)
  choose <- function(pattern) files[grepl(pattern,basename(files),ignore.case=TRUE)]
  m <- choose("^matrix\\.mtx(\\.gz)?$")
  if (length(m)) {
    f <- choose("^(features|genes)\\.tsv(\\.gz)?$"); b <- choose("^barcodes\\.tsv(\\.gz)?$")
    if (length(m)!=1L || length(f)!=1L || length(b)>1L) stop("Ambiguous or incomplete 10X directory: ",sample_id)
    return(read_count_10x(m,f,if(length(b)) b else NULL,sample_id))
  }
  candidates <- choose("\\.(rds|csv|tsv)(\\.gz)?$")
  if (length(candidates)!=1L) stop("Each sample directory needs exactly one supported count input: ",sample_id)
  if (grepl("\\.rds$",candidates,ignore.case=TRUE)) read_count_rds(candidates,sample_id) else read_count_table(candidates,sample_id)
}

# Extract into a fresh directory owned by this session. Reject traversal and
# repeated entries before extraction; never delete or reuse ./userFile.
read_count_zip <- function(path, session_directory) {
  entries <- utils::unzip(path,list=TRUE)$Name
  normalized <- gsub("\\\\","/",entries)
  if (!length(entries) || any(grepl("(^/|^[A-Za-z]:|(^|/)\\.\\.(/|$))",normalized)) || anyDuplicated(normalized)) stop("ZIP contains unsafe or repeated paths.")
  if (!dir.exists(session_directory)) stop("Missing session import directory.")
  destination <- tempfile("upload-",tmpdir=session_directory)
  dir.create(destination)
  utils::unzip(path,exdir=destination,unzip="internal")
  directories <- list.dirs(destination,recursive=TRUE)
  supported <- vapply(directories,function(d) any(grepl("(matrix\\.mtx(\\.gz)?$|\\.(rds|csv|tsv)(\\.gz)?$)",list.files(d),ignore.case=TRUE)),logical(1))
  directories <- directories[supported]
  if (!length(directories)) stop("ZIP contains no supported sample directories.")
  ids <- basename(directories)
  if (anyDuplicated(ids)) stop("Sample directory names must be unique within the ZIP.")
  stats::setNames(lapply(seq_along(directories),function(i) read_count_directory(directories[i],ids[i])),ids)
}

import_duplicate_rows <- function(record) {
  rows <- which(record$feature_identity$duplicate_primary_id)
  data.frame(sample=record$sample_id,record$feature_identity[rows,,drop=FALSE],
    total_umis=as.numeric(Matrix::rowSums(record$counts)[rows]),check.names=FALSE)
}

# The first UI supports a lossless choice or cancellation/corrected re-upload.
# Merging/excluding rows is deliberately not implemented as an implicit repair.
resolve_import_duplicates <- function(record, action, reason) {
  if (!identical(record$status,"needs_feature_decision")) stop("No unresolved feature IDs in this sample.")
  if (length(reason)!=1L || is.na(reason) || !nzchar(trimws(reason))) stop("Record a reason for the import decision.")
  if (!action %in% c("keep_separate","cancel")) stop("Choose keep_separate or cancel.")
  record$decisions <- c(record$decisions,list(list(action=action,reason=trimws(reason),
    affected_rows=record$feature_identity$source_row[record$feature_identity$duplicate_primary_id],time=format(Sys.time(),tz="UTC",usetz=TRUE))))
  if (action=="cancel") {record$status <- "cancelled"; return(record)}
  # make.unique avoids collisions with existing names too. Original identifiers
  # and row numbers remain in the identity table for annotation and audit.
  ids <- make.unique(record$feature_identity$source_feature_id,sep="__import_row_")
  rownames(record$counts) <- ids
  record$feature_identity$internal_feature_id <- ids
  record$status <- "ready"
  record
}

# Bridge to the current named-list-of-counts interface only after import review.
# Full records remain session state so the annotation step can reuse metadata.
import_counts <- function(records) {
  if (!length(records) || !all(vapply(records,function(x) identical(x$status,"ready"),logical(1)))) stop("Resolve or cancel pending import decisions before proceeding.")
  stats::setNames(lapply(records,`[[`,"counts"),names(records))
}

# Route the existing upload choices through the same readers used by tests.
# Shiny upload paths have no meaningful extension, so use input type/content.
read_app_input <- function(input, session_directory, app_root=".") {
  upload_path <- function(file,label) {
    if (is.null(file$datapath) || length(file$datapath)!=1L || !file.exists(file$datapath)) stop("Upload ",label," before reading data.")
    file$datapath
  }
  if (identical(as.character(input$testData),"1")) {
    return(list(sample=read_count_rds(file.path(app_root,"raw_se_S190.rds"))))
  }
  if (!identical(as.character(input$testData),"2")) stop("Select the test dataset or upload a dataset.")
  type <- as.character(input$typeData)
  if (identical(type,"4")) return(read_count_zip(upload_path(input$zipFile,"a ZIP file"),session_directory))
  record <- if (identical(type,"1")) {
    read_count_table(upload_path(input$csvFile,"a count table"))
  } else if (identical(type,"3")) {
    read_count_rds(upload_path(input$rdsFile,"a SingleCellExperiment RDS"))
  } else if (identical(type,"2")) {
    read_count_10x(upload_path(input$mtxInput,"a count matrix"),
      upload_path(input$genes10XInput,"a feature table"),input$barcodes10XInput$datapath)
  } else stop("Select a supported input format.")
  # Shiny datapath is a temporary filename; retain the uploaded filename too.
  upload_fields <- switch(type,"1"="csvFile","2"=c("mtxInput","genes10XInput","barcodes10XInput"),"3"="rdsFile")
  record$source$upload_names <- lapply(upload_fields,function(field) input[[field]]$name)
  names(record$source$upload_names) <- upload_fields
  list(sample=record)
}
