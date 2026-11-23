# Standalone input/metadata support. Building a registry never reads expression arrays.
.gse_read_series <- function(source) {
  local <- source
  if (grepl("^https://", source)) {
    local <- tempfile(fileext = ".txt.gz")
    on.exit(unlink(local), add = TRUE)
    utils::download.file(source, local, mode = "wb", quiet = TRUE, method = "libcurl",
                         headers = c("User-Agent" = "scRNA-Explorer/registry (GSE176031 metadata)",
                                     "Accept" = "application/octet-stream"))
  }
  if (!file.exists(local)) stop("Series matrix does not exist: ", source)
  con <- gzfile(local, "rt")
  on.exit(close(con), add = TRUE)
  rows <- list()
  line_number <- 0L
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line) || identical(line, "!series_matrix_table_begin")) break
    line_number <- line_number + 1L
    if (!startsWith(line, "!")) next
    fields <- scan(text = line, what = character(), sep = "\t", quote = '"',
                   comment.char = "", na.strings = character(), quiet = TRUE)
    rows[[length(rows) + 1L]] <- list(key = fields[1], values = fields[-1],
                                    row = line_number)
  }
  accessions <- Filter(function(x) x$key == "!Sample_geo_accession", rows)
  if (length(accessions) != 1L || !length(accessions[[1]]$values)) {
    stop("Series matrix must contain one nonempty !Sample_geo_accession row.")
  }
  n <- length(accessions[[1]]$values)
  sample_rows <- Filter(function(x) startsWith(x$key, "!Sample_"), rows)
  if (any(vapply(sample_rows, function(x) length(x$values) != n, logical(1)))) {
    stop("Series matrix sample rows have inconsistent column counts.")
  }
  list(rows = rows, sample_rows = sample_rows, accessions = accessions[[1]]$values)
}

resolve_gse176031_input <- function(registry, accession, data_root = registry$roots$data) {
  if (!identical(registry$schema_version, 1L)) stop("Unsupported registry schema version.")
  if (length(accession) != 1L || is.na(accession) || !nzchar(accession)) {
    stop("Supply exactly one nonempty GEO accession.")
  }
  sample <- registry$samples[registry$samples$accession == accession, , drop = FALSE]
  if (nrow(sample) != 1L) stop("Expected exactly one registered sample for ", accession)
  if (is.na(sample$dge_file)) stop("No unambiguous valid DGE input for ", accession)
  path <- file.path(data_root, sample$dge_file)
  if (!file.exists(path)) stop("Registered DGE input is unavailable: ", sample$dge_file)
  cells <- .gse_header(path)
  expected <- registry$inputs[registry$inputs$file == sample$dge_file, , drop = FALSE]
  if (!identical(cells, registry$barcodes[[sample$dge_file]]) ||
      nrow(expected) != 1L || file.info(path)$size != expected$size_bytes) {
    stop("DGE input changed since registry creation; rebuild the registry before reuse.")
  }
  artifacts <- registry$artifacts
  if (!is.null(artifacts)) {
    artifacts <- artifacts[is.na(artifacts$accession) | artifacts$accession == accession, , drop = FALSE]
  }
  list(sample = sample, dge_path = normalizePath(path), barcodes = cells,
       artifacts = artifacts)
}

attach_geo_sample_metadata <- function(object, registry, accessions) {
  if (!identical(registry$schema_version, 1L)) stop("Unsupported registry schema version.")
  sce <- inherits(object, "SingleCellExperiment")
  seurat <- inherits(object, "Seurat")
  if (!sce && !seurat) stop("Supply an existing SingleCellExperiment or Seurat object.")
  cells <- colnames(object)
  if (!length(cells) || anyNA(cells) || any(!nzchar(cells)) || anyDuplicated(cells)) {
    stop("The existing object must have nonempty, unique cell identifiers.")
  }
  if (!is.character(accessions) || anyNA(accessions) || any(!nzchar(accessions))) {
    stop("Supply GEO accessions as a nonempty character mapping.")
  }
  if (!is.null(names(accessions))) {
    if (anyDuplicated(names(accessions)) || !setequal(names(accessions), cells)) {
      stop("Named accession mappings must match the object's cells exactly.")
    }
    accessions <- unname(accessions[match(cells, names(accessions))])
  } else if (length(accessions) == 1L) {
    accessions <- rep(accessions, length(cells))
  } else if (length(accessions) != length(cells)) {
    stop("Supply one accession or one accession per cell, in object order.")
  }
  used <- unique(accessions)
  multiplicity <- vapply(used, function(id) sum(registry$samples$accession == id), integer(1))
  if (any(multiplicity == 0L)) stop("Accession mapping contains an unregistered sample.")
  if (any(multiplicity != 1L)) stop("Accession mapping contains ambiguous GEO sample columns.")
  samples <- registry$samples[match(used, registry$samples$accession), , drop = FALSE]
  original <- registry$metadata[registry$metadata$accession %in% used, , drop = FALSE]
  fields <- unique(original$field)
  derived_fields <- setdiff(names(samples), c("accession", "dge_file", "n_cells"))
  reserved <- unique(c("geo_accession", paste0("geo_", derived_fields)))
  field_columns <- tail(make.unique(c(reserved, paste0("geo_", make.names(fields)))), length(fields))
  wide <- data.frame(row.names = used)
  wide$geo_accession <- used
  for (i in seq_along(fields)) {
    wide[[field_columns[i]]] <- vapply(used, function(id) {
      values <- unique(original$value[original$accession == id & original$field == fields[i]])
      if (length(values)) paste(values, collapse = " | ") else NA_character_
    }, character(1), USE.NAMES = FALSE)
  }
  for (field in derived_fields) {
    column <- paste0("geo_", field)
    if (!column %in% names(wide)) wide[[column]] <- samples[[field]]
  }
  projected <- wide[match(accessions, rownames(wide)), , drop = FALSE]
  rownames(projected) <- cells
  record <- list(samples = samples, sample_metadata = original,
                 field_columns = data.frame(field = fields, column = field_columns),
                 series_metadata = Filter(function(row) startsWith(row$key, "!Series_"), registry$series_rows),
                 provenance = registry$provenance)
  current <- if (sce) as.data.frame(SummarizedExperiment::colData(object)) else object[[]]
  source_metadata <- if (sce) S4Vectors::metadata(object) else object@misc
  if (!is.null(source_metadata$geo) && !identical(source_metadata$geo, record)) {
    stop("Existing GEO provenance differs; it will not be overwritten.")
  }
  for (column in intersect(names(projected), names(current))) {
    if (!identical(as.character(current[[column]]), as.character(projected[[column]]))) {
      stop("Existing metadata would be overwritten: ", column)
    }
  }
  added <- projected[, setdiff(names(projected), names(current)), drop = FALSE]
  if (sce) {
    for (column in names(added)) SummarizedExperiment::colData(object)[[column]] <- added[[column]]
    S4Vectors::metadata(object)$geo <- record
  } else {
    if (ncol(added)) object <- SeuratObject::AddMetaData(object, added)
    object@misc$geo <- record
  }
  object
}

.gse_field <- function(rows, key, n) {
  selected <- Filter(function(x) x$key == key, rows)
  if (!length(selected)) return(rep(NA_character_, n))
  vapply(seq_len(n), function(i) {
    values <- unique(vapply(selected, function(x) x$values[i], character(1)))
    paste(values, collapse = " | ")
  }, character(1))
}

.gse_header <- function(path) {
  con <- gzfile(path, "rt")
  on.exit(close(con))
  header <- readLines(con, n = 1L, warn = FALSE)
  if (!length(header)) stop("Empty DGE file: ", basename(path))
  if (endsWith(header, "\t")) stop("DGE header has an empty final cell identifier: ", basename(path))
  fields <- strsplit(header, "\t", fixed = TRUE)[[1]]
  if (length(fields) < 2L || fields[1] != "GENE") {
    stop("DGE header must start with GENE and contain cell identifiers: ", basename(path))
  }
  cells <- fields[-1]
  if (any(!nzchar(trimws(cells))) || anyDuplicated(cells)) {
    stop("DGE header has empty or duplicate cell identifiers: ", basename(path))
  }
  cells
}

.gse_filename_accessions <- function(files) {
  names <- basename(files)
  as.character(ifelse(grepl("^GSM[0-9]+_", names), sub("_.*", "", names), NA_character_))
}

.gse_walk_counts <- function(input, consume, chunk_rows = 500L) {
  if (length(chunk_rows) != 1L || !is.numeric(chunk_rows) || is.na(chunk_rows) ||
      !is.finite(chunk_rows) || chunk_rows < 1 || chunk_rows != floor(chunk_rows)) {
    stop("chunk_rows must be a positive integer.")
  }
  con <- gzfile(input$dge_path, "rt")
  on.exit(close(con))
  header <- readLines(con, n = 1L, warn = FALSE)
  if (!identical(strsplit(header, "\t", fixed = TRUE)[[1]][-1], input$barcodes)) {
    stop("DGE header changed before count reading.")
  }
  seen <- character()
  repeat {
    lines <- readLines(con, n = chunk_rows, warn = FALSE)
    if (!length(lines)) break
    table <- utils::read.table(text = lines, sep = "\t", header = FALSE, quote = "",
                               comment.char = "", na.strings = character(),
                               colClasses = c("character", rep("numeric", length(input$barcodes))),
                               check.names = FALSE)
    if (ncol(table) != length(input$barcodes) + 1L) stop("DGE count columns do not match its header.")
    genes <- table[[1]]
    if (anyNA(genes) || any(!nzchar(trimws(genes))) || anyDuplicated(c(seen, genes))) {
      stop("DGE feature identifiers must be nonempty and unique for verified reuse.")
    }
    values <- as.matrix(table[-1])
    if (any(!is.finite(values)) || any(values < 0) || any(values != floor(values))) {
      stop("DGE values must be finite, nonnegative integer counts.")
    }
    dimnames(values) <- list(genes, input$barcodes)
    consume(genes, values)
    seen <- c(seen, genes)
  }
  if (!length(seen)) stop("DGE file contains no feature rows.")
  invisible(seen)
}

read_gse176031_counts <- function(registry, accession, data_root = registry$roots$data,
                                 chunk_rows = 500L) {
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Count import requires an existing Matrix installation.")
  input <- resolve_gse176031_input(registry, accession, data_root)
  blocks <- list()
  .gse_walk_counts(input, function(genes, values) {
    nonzero <- which(values != 0, arr.ind = TRUE)
    blocks[[length(blocks) + 1L]] <<- Matrix::sparseMatrix(
      i = nonzero[, 1], j = nonzero[, 2], x = values[nonzero], dims = dim(values),
      dimnames = list(genes, input$barcodes))
  }, chunk_rows)
  do.call(rbind, blocks)
}

inspect_gse176031_r_artifact <- function(registry, path, accessions, original_barcodes = NULL,
                                        data_root = registry$roots$data, chunk_rows = 500L) {
  path <- normalizePath(path, mustWork = TRUE)
  row <- which(registry$artifacts$path == path)
  if (length(row) != 1L || registry$artifacts$format[row] != "rds") {
    stop("Explicit inspection requires exactly one registered RDS artifact.")
  }
  # This operation is explicit: ordinary discovery never deserializes an R object.
  object <- readRDS(path)
  object <- attach_geo_sample_metadata(object, registry, accessions)
  cells <- colnames(object)
  accessions <- as.character(object$geo_accession)
  if (is.null(original_barcodes)) original_barcodes <- cells
  if (!is.character(original_barcodes) || anyNA(original_barcodes) ||
      length(original_barcodes) != length(cells) || any(!nzchar(original_barcodes))) {
    stop("Supply one original DGE barcode per object cell.")
  }
  if (!is.null(names(original_barcodes))) {
    if (anyDuplicated(names(original_barcodes)) || !setequal(names(original_barcodes), cells)) {
      stop("Original-barcode mappings must match the object's cells exactly.")
    }
    original_barcodes <- unname(original_barcodes[match(cells, names(original_barcodes))])
  }
  counts <- if (inherits(object, "SingleCellExperiment")) {
    SummarizedExperiment::assay(object, "counts")
  } else {
    layers <- SeuratObject::Layers(object[["RNA"]], search = "^counts($|\\.)")
    if (!identical(layers, "counts")) stop("Explicit inspection requires one RNA counts layer.")
    SeuratObject::LayerData(object, assay = "RNA", layer = "counts")
  }
  genes <- rownames(counts)
  if (!length(genes) || anyNA(genes) || any(!nzchar(genes)) || anyDuplicated(genes) ||
      !all(cells %in% colnames(counts))) {
    stop("Existing counts must have unique feature IDs and contain every mapped cell.")
  }
  for (id in unique(accessions)) {
    selected <- which(accessions == id)
    barcodes <- original_barcodes[selected]
    input <- resolve_gse176031_input(registry, id, data_root)
    if (anyDuplicated(barcodes) || !all(barcodes %in% input$barcodes)) {
      stop("Existing cell barcodes do not match the declared DGE source.")
    }
    found <- logical(length(genes))
    .gse_walk_counts(input, function(source_genes, values) {
      positions <- match(source_genes, genes)
      included <- which(!is.na(positions))
      if (!length(included)) return(invisible(NULL))
      observed <- as.matrix(counts[positions[included], cells[selected], drop = FALSE])
      expected <- values[included, match(barcodes, input$barcodes), drop = FALSE]
      if (anyNA(observed) || any(observed != expected)) {
        stop("Existing counts do not match the original DGE values for the declared cells/features.")
      }
      found[positions[included]] <<- TRUE
    }, chunk_rows)
    if (!all(found)) stop("Existing features do not match the original DGE feature identifiers.")
  }
  registry$artifacts$inspection_status[row] <- "existing_R_object_source_verified"
  registry$artifacts$counts_verified[row] <- TRUE
  registry$artifacts$compatibility_verified[row] <- TRUE
  registry$artifacts$detail[row] <- paste(
    class(object)[1], "counts verified against DGE for the declared feature/cell subset;",
    "GEO metadata attached without replacing counts. Filtering/biological quality is not established.")
  list(object = object, registry = registry)
}

.gse_artifacts <- function(roots, explicit_paths, add_issue) {
  paths <- explicit_paths
  for (root in roots) {
    if (!dir.exists(root)) {
      add_issue("artifact_root_missing", NA_character_, root)
      next
    }
    paths <- c(paths, list.files(root, pattern = "\\.(rds|rda|rdata)$",
                                 ignore.case = TRUE, recursive = TRUE, full.names = TRUE))
  }
  paths <- unique(paths)
  result <- data.frame(path = character(), accession = character(), format = character(),
                       stage_hint = character(), inspection_status = character(),
                       counts_verified = logical(), compatibility_verified = logical(),
                       detail = character(), stringsAsFactors = FALSE)
  for (path in paths) {
    name <- basename(path)
    format <- tolower(tools::file_ext(path))
    stage <- if (grepl("merged", name, ignore.case = TRUE)) "merged_name_hint" else
      if (grepl("filtered", name, ignore.case = TRUE)) "filtered_name_hint" else
        if (grepl("raw", name, ignore.case = TRUE)) "raw_name_hint" else "unspecified"
    row <- data.frame(path = path,
                      accession = .gse_filename_accessions(name),
                      format = format, stage_hint = stage, inspection_status = "uninspected",
                      counts_verified = FALSE, compatibility_verified = FALSE,
                      detail = "Count values, dataset association, and analysis compatibility are not verified.",
                      stringsAsFactors = FALSE)
    tryCatch({
      if (!file.exists(path)) stop("Artifact does not exist.")
      row$path <- normalizePath(path)
      if (format == "rds") {
        info <- infoRDS(path)
        row$inspection_status <- "R_serialization_header_only"
        row$detail <- paste(row$detail, "Minimum R reader:", info$min_reader_version)
      } else if (format %in% c("rda", "rdata")) {
        row$inspection_status <- "R_workspace_not_loaded"
      } else stop("Unsupported optional artifact format.")
    }, error = function(e) {
      row$inspection_status <<- "inspection_failed"
      row$detail <<- conditionMessage(e)
      add_issue("artifact_inspection_failed", row$accession, paste(name, conditionMessage(e)))
    })
    result <- rbind(result, row)
  }
  result[order(!result$format %in% c("rds", "rda", "rdata"), result$path), , drop = FALSE]
}

build_gse176031_registry <- function(data_root, series_matrix, artifact_roots = character(),
                                    r_artifacts = character(), confirmed_genome_build = NULL) {
  if (!is.null(confirmed_genome_build) &&
      (!is.character(confirmed_genome_build) || length(confirmed_genome_build) != 1L ||
       is.na(confirmed_genome_build) || !nzchar(confirmed_genome_build))) {
    stop("confirmed_genome_build must be one nonempty reference name.")
  }
  data_root <- normalizePath(data_root, mustWork = TRUE)
  if (!dir.exists(data_root)) stop("data_root must be a directory.")
  series <- .gse_read_series(series_matrix)
  ids <- series$accessions
  issues <- data.frame(severity = character(), code = character(), accession = character(),
                       detail = character(), stringsAsFactors = FALSE)
  add_issue <- function(code, accession, detail, severity = "warning") {
    issues[nrow(issues) + 1L, ] <<- list(severity, code, accession, detail)
  }
  duplicate_metadata <- unique(ids[duplicated(ids)])
  for (id in duplicate_metadata) {
    add_issue("duplicate_metadata_accession", id,
              "More than one GEO sample column uses this accession.", "error")
  }
  metadata <- do.call(rbind, lapply(series$sample_rows, function(row) {
    field <- rep(sub("^!Sample_", "", row$key), length(ids))
    value <- row$values
    if (row$key == "!Sample_characteristics_ch1") {
      split <- regexpr(":", value, fixed = TRUE)
      has_key <- split > 0L
      field[has_key] <- trimws(substr(value[has_key], 1L, split[has_key] - 1L))
      value[has_key] <- trimws(sub("^[^:]*:", "", value[has_key]))
    }
    data.frame(accession = ids, field = field, value = value, source_value = row$values,
               source_field = row$key, source_row = row$row, stringsAsFactors = FALSE)
  }))
  characteristic <- function(field) vapply(ids, function(id) {
    values <- unique(metadata$value[metadata$accession == id & metadata$field == field])
    if (length(values) == 1L && nzchar(values)) values else NA_character_
  }, character(1), USE.NAMES = FALSE)
  files <- sort(list.files(data_root, pattern = "_dge\\.txt(\\.gz)?$", ignore.case = TRUE))
  file_ids <- .gse_filename_accessions(files)
  barcodes <- setNames(lapply(seq_along(files), function(i) {
    tryCatch(.gse_header(file.path(data_root, files[i])), error = function(e) {
      add_issue("invalid_dge_header", file_ids[i], conditionMessage(e), "error")
      character()
    })
  }), files)
  info <- file.info(file.path(data_root, files))
  inputs <- data.frame(file = files, accession = file_ids, size_bytes = info$size,
                       n_cells = lengths(barcodes), stringsAsFactors = FALSE)
  for (id in unique(file_ids[!is.na(file_ids) & duplicated(file_ids)])) {
    add_issue("duplicate_dge_accession", id,
              paste(inputs$file[which(inputs$accession == id)], collapse = " | "), "error")
  }
  matched <- match(ids, inputs$accession)
  for (i in seq_along(ids)) {
    candidates <- which(inputs$accession == ids[i])
    if (!length(candidates)) {
      add_issue("missing_dge", ids[i], "No local DGE input matches this GEO accession.", "error")
    }
    if (length(candidates) != 1L || !inputs$n_cells[candidates[1]] ||
        ids[i] %in% duplicate_metadata) matched[i] <- NA_integer_
  }
  for (i in which(is.na(inputs$accession) | !inputs$accession %in% ids)) {
    add_issue("unmatched_dge", inputs$accession[i], inputs$file[i])
  }
  design <- .gse_field(series$rows, "!Series_overall_design", 1L)
  protocol <- if (!is.na(design) && grepl("Seq.Well", design, ignore.case = TRUE)) {
    "Seq-Well S3"
  } else NA_character_
  samples <- data.frame(
    accession = ids,
    title = .gse_field(series$sample_rows, "!Sample_title", length(ids)),
    patient = characteristic("patient"), tissue = characteristic("tissue"),
    source_name = .gse_field(series$sample_rows, "!Sample_source_name_ch1", length(ids)),
    organism = .gse_field(series$sample_rows, "!Sample_organism_ch1", length(ids)),
    platform = .gse_field(series$sample_rows, "!Sample_platform_id", length(ids)),
    protocol = rep(protocol, length(ids)), capture_id = NA_character_,
    biological_sample_id = NA_character_, sequencing_pool_id = NA_character_,
    dge_file = inputs$file[matched], n_cells = as.integer(inputs$n_cells[matched]),
    stringsAsFactors = FALSE
  )
  samples$donor_id <- samples$patient
  samples$data_processing <- .gse_field(series$sample_rows, "!Sample_data_processing", length(ids))
  samples$genome_build <- vapply(ids, function(id) {
    values <- metadata$value[metadata$accession == id & metadata$field == "data_processing"]
    values <- unique(trimws(sub("^Genome_build:", "", values[startsWith(values, "Genome_build:")])))
    if (length(values) == 1L) values else NA_character_
  }, character(1), USE.NAMES = FALSE)
  samples$confirmed_genome_build <- if (is.null(confirmed_genome_build)) NA_character_ else confirmed_genome_build
  for (i in seq_along(ids)) {
    text <- samples$data_processing[i]
    if (is.na(text)) next
    mouse <- grepl("\\b(mm[0-9]+|GRCm[0-9]+)\\b", text, ignore.case = TRUE)
    human <- grepl("\\b(hg[0-9]+|GRCh[0-9]+)\\b", text, ignore.case = TRUE)
    if (mouse && (human || isTRUE(samples$organism[i] == "Homo sapiens"))) {
      if (is.null(confirmed_genome_build)) {
        add_issue("reference_conflict", ids[i],
                  "GEO declares mouse-reference alignment alongside human organism/build evidence; unresolved.")
      } else {
        add_issue("reference_conflict_resolved", ids[i],
                  paste("Original GEO contradiction retained; explicit user confirmation:", confirmed_genome_build), "info")
      }
    }
  }
  add_issue("capture_mapping_unresolved", NA_character_,
            "GEO accession/patient fields do not establish capture, biological-sample, or sequencing-pool identity.", "info")
  repeated <- sum(duplicated(unlist(barcodes, use.names = FALSE)))
  if (repeated) {
    add_issue("repeated_barcodes_across_inputs", NA_character_,
              paste(repeated, "barcode occurrences repeat across inputs; preserve accession/source context."), "info")
  }
  artifacts <- .gse_artifacts(artifact_roots, r_artifacts, add_issue)
  list(schema_version = 1L, samples = samples, inputs = inputs, barcodes = barcodes,
       metadata = metadata, series_rows = series$rows, issues = issues, artifacts = artifacts,
       provenance = list(series_matrix = series_matrix, overall_design = design,
                         reference_confirmation = if (is.null(confirmed_genome_build)) "not provided" else
                           paste("explicit user confirmation:", confirmed_genome_build),
                         R_version = R.version.string, expression_values_read = FALSE),
       roots = list(data = data_root, artifacts = artifact_roots),
       capabilities = list(background_droplets_confirmed = FALSE,
                           independent_quality_labels_confirmed = FALSE))
}
