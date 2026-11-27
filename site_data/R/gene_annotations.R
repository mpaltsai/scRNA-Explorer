# Cached AnnotationHub-derived references and identifier-preserving QC preparation.
# No network lookup, ortholog conversion, gene removal or count correction occurs.
# Private sourcing also supports the standalone runner without Shiny autoloading.
.annotation_cache_api <- function(app_root) {
  e <- new.env(parent=baseenv())
  sys.source(file.path(app_root,"R","annotation_reference_cache.R"),envir=e)
  sys.source(file.path(app_root,"R","annotation_biomart_supplement.R"),envir=e)
  e
}

load_gene_annotation_cache <- function(organism, app_root=".", include_ribosomal=TRUE,
    reference_id="bundled",cache_dir=NULL) {
  api <- .annotation_cache_api(app_root)
  if (is.null(cache_dir)) cache_dir <- api$annotation_cache_dir(app_root)
  implementation <- normalizePath(file.path(app_root,"R","gene_annotations.R"),mustWork=TRUE)
  result <- api$load_versioned_annotation(reference_id,organism,app_root,cache_dir)
  result$provenance$loader_path <- implementation
  result$provenance$loader_md5 <- unname(tools::md5sum(implementation))
  result$provenance$cache_loader_md5 <- unname(tools::md5sum(file.path(app_root,"R","annotation_reference_cache.R")))
  if (include_ribosomal) result$ribosomal <- load_ribosomal_reference(organism)
  result
}

.qc_metric_names <- c("cells","nUMI","nGene","log10GenesPerUMI","sample","mtUMI","mitoRatio","rbcUMI","rbcRatio")
.qc_ribosomal_fields <- c("ribosomal_protein","cytosolic_ribosomal_protein","mitochondrial_ribosomal_protein","nuclear_rRNA","mitochondrial_rRNA","rRNA_pseudogene")
.qc_annotation_fields <- c("source_feature_id","gene_id","gene_name","gene_biotype","seq_name","description","annotation_status","annotation_candidates","mitochondrial","hemoglobin",.qc_ribosomal_fields,"ribosomal_annotation_status")
.ribosomal_reference_cache <- new.env(parent=emptyenv())

# Use existing, species-specific OrgDb snapshots, not symbol prefixes (which
# would also match RPS6 kinases). No package install or reference download occurs.
# The protein definition is explicitly GO-supported structural constituents;
# reference nonmembership does not prove that a gene has no ribosomal function.
load_ribosomal_reference <- function(organism) {
  packages <- c("Homo sapiens"="org.Hs.eg.db","Mus musculus"="org.Mm.eg.db")
  if (length(organism)!=1L || !organism %in% names(packages)) stop("Select a species-specific ribosomal reference.")
  package <- packages[[organism]]
  needed <- c("AnnotationDbi",package)
  if (!all(vapply(needed,requireNamespace,logical(1),quietly=TRUE))) return(list(
    status="not_assessed",species=organism,reason="Existing AnnotationDbi/OrgDb packages unavailable; no installation attempted."))
  tryCatch({
    db <- getExportedValue(package,package)
    info <- file.info(AnnotationDbi::dbfile(db))
    cache_key <- paste(organism,utils::packageVersion(package),utils::packageVersion("AnnotationDbi"),info$size,info$mtime,sep="|")
    if (exists(cache_key,envir=.ribosomal_reference_cache,inherits=FALSE)) return(get(cache_key,envir=.ribosomal_reference_cache))
    terms <- c(structural="GO:0003735",cytosolic="GO:0022626",mitochondrial="GO:0005761")
    rows <- suppressMessages(AnnotationDbi::select(db,keys=unname(terms),keytype="GOALL",columns="ENTREZID"))
    sets <- lapply(terms,function(term) sort(unique(rows$ENTREZID[!is.na(rows$ENTREZID) & rows$GOALL==term])))
    types <- suppressMessages(AnnotationDbi::select(db,keys=sets$structural,keytype="ENTREZID",columns="GENETYPE"))
    coding <- types$ENTREZID[!is.na(types$GENETYPE) & types$GENETYPE=="protein-coding"]
    proteins <- intersect(sets$structural,coding)
    groups <- list(ribosomal_protein=proteins,cytosolic_ribosomal_protein=intersect(proteins,sets$cytosolic),
      mitochondrial_ribosomal_protein=intersect(proteins,sets$mitochondrial))
    result <- list(status="available",species=organism,package=package,package_version=as.character(utils::packageVersion(package)),
      annotation_dbi_version=as.character(utils::packageVersion("AnnotationDbi")),
      source_metadata=AnnotationDbi::metadata(db),go_terms=terms,key_type="ENTREZID",
      go_key_type="GOALL",evidence_policy="All evidence codes in the installed OrgDb snapshot; not a completeness or quality claim",
      definition="Protein-coding Entrez genes annotated to structural constituent of ribosome; cytosolic/mitochondrial subgroups require the respective GOALL location",
      groups=groups,groups_hash=digest::digest(groups,algo="sha256"))
    assign(cache_key,result,envir=.ribosomal_reference_cache)
    result
  },error=function(e) list(status="not_assessed",species=organism,reason=conditionMessage(e)))
}

# Resolve only agreement across candidate annotation rows/Entrez IDs. Unknown or
# conflicting classifications remain NA. This annotation step never drops rows,
# changes count denominators, or excludes genes from variable-feature selection.
.ribosomal_feature_flags <- function(a, matches, reference) {
  n <- nrow(a)
  row_flags <- list()
  for (group in .qc_ribosomal_fields[1:3]) row_flags[[group]] <- vapply(seq_len(n),function(i) {
    biotype <- as.character(a$gene_biotype[i])
    if (is.na(biotype) || !nzchar(biotype)) return(NA)
    if (!biotype %in% c("protein_coding","LRG_gene")) return(FALSE)
    if (!identical(reference$status,"available") || !"entrezid" %in% names(a)) return(NA)
    ids <- unique(as.character(unlist(a$entrezid[i],use.names=FALSE)))
    ids <- ids[!is.na(ids) & nzchar(ids)]
    if (!length(ids)) return(NA)
    values <- ids %in% reference$groups[[group]]
    if (length(unique(values))==1L) values[1] else NA
  },logical(1))
  biotype <- tolower(as.character(a$gene_biotype))
  chromosome <- as.character(a$seq_name)
  row_flags$nuclear_rRNA <- biotype=="rrna" & chromosome!="MT"
  row_flags$mitochondrial_rRNA <- biotype=="mt_rrna" | (biotype=="rrna" & chromosome=="MT")
  row_flags$rRNA_pseudogene <- ifelse(is.na(biotype),NA,grepl("rrna.*pseudogene",biotype))
  consensus <- function(values) if (!length(values) || anyNA(values) || length(unique(values))!=1L) NA else values[1]
  mapped <- as.data.frame(lapply(row_flags,function(flag) unname(vapply(matches,function(rows) consensus(flag[rows]),logical(1)))))
  mapped$ribosomal_annotation_status <- apply(mapped,1,function(v) {
    if (anyNA(v)) "unresolved" else if (any(v)) "reference_member" else "not_annotated_member"
  })
  mapped
}

# All imported fields remain intact. Own current metrics use qc_* columns; the
# legacy names are supplied only when absent, for backward-compatible exports.
.check_qc_namespace <- function(record) {
  own <- if (is.null(record$source_object)) NULL else S4Vectors::metadata(record$source_object)[["scrna_explorer_qc"]]
  recognised <- is.list(own) && identical(own$producer,"scRNA-Explorer") && identical(own$schema_version,1L)
  reserved_col <- paste0("qc_",c(.qc_metric_names,"cell_id","generated_cell_id"))
  reserved_row <- paste0("qc_",.qc_annotation_fields)
  if (!recognised && (any(reserved_col %in% names(record$cell_data)) ||
      any(reserved_row %in% names(record$feature_data)) || !is.null(own))) {
    stop("Imported metadata uses the reserved QC namespace. Review/rename those fields before proceeding; none were overwritten.")
  }
}

# Reuse recorded source identifiers, including those retained by an earlier
# lossless duplicate-row decision. Suggestions and actual lookup must inspect
# the same IDs, not the repaired internal row names of an exported app SCE.
.annotation_source_features <- function(record) {
  original <- as.character(record$feature_identity$source_feature_id)
  prior <- if (is.null(record$source_object)) NULL else S4Vectors::metadata(record$source_object)[["scrna_explorer_qc"]]
  if (is.list(prior) && identical(prior$producer,"scRNA-Explorer") &&
      identical(prior$schema_version,1L) && "qc_source_feature_id" %in% names(record$feature_data))
    original <- as.character(record$feature_data$qc_source_feature_id)
  original
}

# Format recognition is a suggestion about lookup keys, not a species/build
# assignment or a biological mapping repair. Mixed/unsupported identifiers need
# an explicit choice. Symbols require at least one exact match in the selected
# local reference; unknown rows are retained and reported by normal annotation.
suggest_feature_identifiers <- function(record,reference=NULL) {
  ids <- .annotation_source_features(record)
  gene <- grepl("^ENS[A-Z]*G[0-9]+(\\.[0-9]+)?$",ids)
  versioned <- gene & grepl("\\.[0-9]+$",ids)
  other <- grepl("^ENS[A-Z]*[TPE][0-9]+(\\.[0-9]+)?$",ids) | grepl("^[0-9]+$",ids)
  human <- any(grepl("^ENSG[0-9]",ids));mouse <- any(grepl("^ENSMUSG[0-9]",ids))
  symbols <- if (is.null(reference)) rep(FALSE,length(ids)) else ids %in% reference$table$gene_name
  ref_ids <- if (is.null(reference)) rep(FALSE,length(ids)) else ids %in% reference$table$gene_id
  result <- list(suggested_type="",status="needs_review",reason="No supported identifier type established; choose the type explicitly.",
    features=length(ids),ensembl_features=sum(gene),versioned_features=sum(versioned),
    symbol_matches=sum(symbols),reference_id_matches=sum(ref_ids),
    strip_version_suggested=any(versioned),token=substr(digest::digest(ids,algo="sha256"),1,16))
  if (length(ids) && all(gene) && !(human && mouse)) {
    result$suggested_type <- "gene_id";result$status <- "suggested"
    result$reason <- "Ensembl gene-ID format detected. Confirm species and genome build separately."
  } else if (any(gene) || any(other)) {
    result$reason <- if (human && mouse) "Human and mouse Ensembl prefixes are mixed; review the input before choosing a lookup type." else
      "Mixed gene IDs, numeric IDs or transcript/protein/exon IDs need review; choose a lookup type explicitly or correct the input."
  } else if (any(symbols) && !any(ref_ids & !symbols)) {
    result$suggested_type <- "gene_name";result$status <- "suggested"
    result$reason <- "Gene symbols suggested from exact matches in the selected local reference; unmatched features remain for review."
  } else if (any(ref_ids) && !any(symbols)) {
    result$suggested_type <- "gene_id";result$status <- "suggested"
    result$reason <- "Gene IDs suggested from exact matches in the selected local reference."
  }
  result
}

# Browser choices carry a feature-list token. A manual choice left over from a
# different upload cannot silently control a new feature list. Plain key values
# remain accepted for existing programmatic callers/tests of this controller.
resolve_feature_identifier_choice <- function(suggestion,value=NULL) {
  if (length(value)==1L && !is.na(value)) {
    for (key in c("gene_id","gene_name"))
      if (value %in% c(key,paste(key,suggestion$token,sep=":"))) return(list(key=key,mode="manual"))
  }
  list(key=suggestion$suggested_type,mode="automatic")
}

# Return small annotation/metric tables, not another expression object. Lookups
# use source IDs even when the user retained duplicated rows with new internal IDs.
prepare_sample_qc <- function(record, reference, feature_id_type, strip_ensembl_version=FALSE) {
  if (!identical(record$status,"ready")) stop("Complete import review first.")
  if (length(feature_id_type)!=1L || !feature_id_type %in% c("gene_id","gene_name")) stop("Select gene symbols or Ensembl gene IDs for each sample.")
  if (!is.logical(strip_ensembl_version) || length(strip_ensembl_version)!=1L || is.na(strip_ensembl_version)) stop("Version stripping must be explicitly TRUE or FALSE.")
  a <- reference$table
  required <- c("gene_id","gene_name","gene_biotype","seq_name","description")
  if (!is.data.frame(a) || !all(required %in% names(a))) stop("The annotation reference lacks required gene fields.")
  if (!is.list(reference$provenance) || length(reference$provenance$source)!=1L ||
      is.na(reference$provenance$source) || !nzchar(reference$provenance$source)) stop("Record the annotation reference source.")
  if (length(reference$species)!=1L || !reference$species %in% c("Homo sapiens","Mus musculus")) stop("Reference species must be explicit.")
  if (!is.null(reference$ribosomal) && !identical(reference$ribosomal$species,reference$species)) stop("Ribosomal reference species must match gene annotations.")
  .check_qc_namespace(record)
  normalise <- function(x) if (feature_id_type=="gene_id" && strip_ensembl_version) {
    sub("^(ENS[A-Z]*G[0-9]+)\\.[0-9]+$","\\1",x)
  } else x
  keys <- normalise(as.character(a[[feature_id_type]]))
  original <- .annotation_source_features(record)
  lookup <- normalise(original)
  valid <- !is.na(keys) & nzchar(keys)
  groups <- split(which(valid),keys[valid])
  matches <- groups[match(lookup,names(groups))]
  one_value <- function(rows,field) {
    values <- unique(as.character(a[[field]][rows]))
    values <- values[!is.na(values) & nzchar(values)]
    if (length(values)==1L) values else NA_character_
  }
  mapped <- data.frame(source_feature_id=original,stringsAsFactors=FALSE)
  for (field in required) mapped[[field]] <- vapply(matches,one_value,character(1),field=field)
  mapped$annotation_status <- vapply(matches,function(rows) {
    if (!length(rows)) return("unmapped")
    ids <- unique(a$gene_id[rows]); ids <- ids[!is.na(ids) & nzchar(ids)]
    symbols <- unique(a$gene_name[rows]); symbols <- symbols[!is.na(symbols) & nzchar(symbols)]
    if (length(ids)>1L || length(symbols)>1L) "ambiguous" else if (!length(ids) || !length(symbols)) "partial" else "matched"
  },character(1))
  mapped$annotation_candidates <- vapply(matches,function(rows) paste(unique(a$gene_id[rows]),collapse=";"),character(1))
  # Preserve the user's existing gene definitions exactly, expressed in the
  # declared identifier type. This set membership does not resolve ambiguous IDs.
  mt_rows <- !is.na(a$seq_name) & a$seq_name=="MT"
  hb_rows <- !is.na(a$description) & grepl("^hemoglobin",a$description,ignore.case=TRUE)
  mt_set <- unique(keys[valid & mt_rows]); hb_set <- unique(keys[valid & hb_rows])
  mapped$mitochondrial <- lookup %in% mt_set
  mapped$hemoglobin <- lookup %in% hb_set
  mapped <- cbind(mapped,.ribosomal_feature_flags(a,matches,reference$ribosomal))
  counts <- record$counts
  umi <- as.numeric(Matrix::colSums(counts)); genes <- as.numeric(Matrix::colSums(counts>0))
  metric_status <- function(set,flag) if (!length(set)) "missing_reference_definition" else if (!any(flag)) "no_matching_features" else "available"
  mt_status <- metric_status(mt_set,mapped$mitochondrial)
  hb_status <- metric_status(hb_set,mapped$hemoglobin)
  # Missing annotation evidence is unavailable, not a biological zero. Known
  # genes present as zero-count rows do yield genuine observed zero fractions.
  totals <- function(flag,status) if (status=="available") as.numeric(Matrix::colSums(counts[flag,,drop=FALSE])) else rep(NA_real_,ncol(counts))
  mt <- totals(mapped$mitochondrial,mt_status); hb <- totals(mapped$hemoglobin,hb_status)
  fraction <- function(x) {y<-x/umi; y[!is.finite(y)]<-NA_real_; y}
  novelty <- log10(genes)/log10(umi); novelty[!is.finite(novelty)] <- NA_real_
  metrics <- data.frame(cells=colnames(counts),nUMI=umi,nGene=genes,log10GenesPerUMI=novelty,
    sample=rep(record$sample_id,ncol(counts)),mtUMI=mt,mitoRatio=fraction(mt),rbcUMI=hb,rbcRatio=fraction(hb),
    row.names=colnames(counts),stringsAsFactors=FALSE)
  conflicts <- .qc_metric_names[vapply(.qc_metric_names,function(field) {
    field %in% names(record$cell_data) && !isTRUE(all.equal(record$cell_data[[field]],metrics[[field]],check.attributes=FALSE))
  },logical(1))]
  list(metrics=metrics,feature_annotation=mapped,
    config=list(organism=reference$species,feature_id_type=feature_id_type,
      strip_ensembl_version=strip_ensembl_version,reference=reference$provenance,
      reference_table_hash=digest::digest(a,algo="sha256"),
      gene_definitions=list(mitochondrial="seq_name == MT",hemoglobin="description starts with hemoglobin (case insensitive)"),
      gene_sets=list(mitochondrial=mt_set,hemoglobin=hb_set),
      ribosomal_reference=reference$ribosomal,
      ribosomal_rna_definition="Cached gene_biotype rRNA, Mt_rRNA and rRNA_pseudogene, separately; candidate agreement required"),
    summary=data.frame(sample=record$sample_id,features=nrow(counts),cells=ncol(counts),
      matched=sum(mapped$annotation_status=="matched"),ambiguous=sum(mapped$annotation_status=="ambiguous"),
      partial=sum(mapped$annotation_status=="partial"),unmapped=sum(mapped$annotation_status=="unmapped"),
      mt_features=sum(mapped$mitochondrial),hb_features=sum(mapped$hemoglobin),
      mito_status=mt_status,hemoglobin_status=hb_status,preserved_conflicting_fields=paste(conflicts,collapse=",")),
    ribosomal_summary=data.frame(sample=record$sample_id,group=.qc_ribosomal_fields,
      annotated_features=vapply(.qc_ribosomal_fields,function(g) sum(mapped[[g]] %in% TRUE),integer(1)),
      unresolved_features=vapply(.qc_ribosomal_fields,function(g) sum(is.na(mapped[[g]])),integer(1))),
    legacy_field_conflicts=conflicts)
}

# Materialise the baseline once. Reuse imported SCE assays/reductions/metadata;
# for matrix inputs create the SCE here. Subsetting preserves all cells' metadata.
make_baseline_sce <- function(record, prepared) {
  stopifnot(identical(record$status,"ready"),identical(prepared$metrics$cells,colnames(record$counts)))
  .check_qc_namespace(record)
  object <- record$source_object
  if (is.null(object)) object <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=record$counts),
    rowData=S4Vectors::DataFrame(record$feature_data,row.names=rownames(record$counts)),
    colData=S4Vectors::DataFrame(record$cell_data,row.names=colnames(record$counts)))
  # Dimname changes only implement previously recorded missing/duplicate-ID
  # decisions. Counts are never aggregated, rounded or reannotated by renaming.
  dimnames(object) <- dimnames(record$counts)
  SummarizedExperiment::assay(object,"counts") <- record$counts
  for (field in .qc_metric_names) {
    object[[paste0("qc_",field)]] <- prepared$metrics[[field]]
    if (!field %in% names(SummarizedExperiment::colData(object))) object[[field]] <- prepared$metrics[[field]]
  }
  object$qc_cell_id <- record$cell_identity$internal_cell_id
  object$qc_generated_cell_id <- record$cell_identity$generated
  for (field in .qc_annotation_fields) SummarizedExperiment::rowData(object)[[paste0("qc_",field)]] <- prepared$feature_annotation[[field]]
  keep <- prepared$metrics$nUMI>=100
  md <- S4Vectors::metadata(object)
  md$scrna_explorer_qc <- list(producer="scRNA-Explorer",schema_version=1L,config=prepared$config,
    previous_qc=md[["scrna_explorer_qc"]],
    mapping_summary=prepared$summary,
    import_source=record$source,import_decisions=record$decisions,import_sample_names=record$sample_name_history,
    legacy_field_conflicts=prepared$legacy_field_conflicts,
    initial_filter=list(rule="nUMI >= 100",input_cells=ncol(object),retained_cells=sum(keep),
      excluded_cells=colnames(object)[!keep]))
  S4Vectors::metadata(object) <- md
  object[,keep,drop=FALSE]
}

# Read current metrics by name, never by the first nine colData positions. This
# also works after the legacy filter's SCE subset, without recomputing expression.
baseline_qc_metrics <- function(object) {
  fields <- lapply(.qc_metric_names,function(field) {
    key <- if (paste0("qc_",field) %in% names(SummarizedExperiment::colData(object))) paste0("qc_",field) else field
    value <- object[[key]]
    if (is.null(value)) stop("Missing QC metric: ",field)
    value
  })
  names(fields) <- .qc_metric_names
  as.data.frame(fields,stringsAsFactors=FALSE,row.names=colnames(object))
}
