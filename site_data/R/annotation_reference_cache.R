# Reference maintenance is separate from reading counts and preparing QC.
# Only explicit maintenance actions access providers. Downloads and supplements
# save separate annotation versions; neither action selects one for analysis.
annotation_cache_dir <- function(app_root=".") {
  getOption("scrnaExplorer.annotationCache",file.path(normalizePath(app_root,mustWork=TRUE),"..","analysis","local","annotation-cache"))
}

annotation_genome_build <- function(x) {
  aliases <- c(hg19="GRCh37",GRCh37="GRCh37",hg38="GRCh38",GRCh38="GRCh38",
    mm10="GRCm38",GRCm38="GRCm38",mm39="GRCm39",GRCm39="GRCm39")
  result <- unname(aliases[as.character(x)])
  # Unknown assemblies cannot be advertised as compatible.
  result[is.na(result)] <- NA_character_
  result
}

.annotation_species <- function(organism) {
  if (length(organism)!=1L || is.na(organism) || !organism %in% c("Homo sapiens","Mus musculus"))
    stop("Select human or mouse explicitly.")
  organism
}

.annotation_target <- function(organism,genome_build) {
  .annotation_species(organism)
  build <- annotation_genome_build(genome_build)
  allowed <- if (organism=="Homo sapiens") c("GRCh37","GRCh38") else c("GRCm38","GRCm39")
  if (length(build)!=1L || is.na(build) || !build %in% allowed)
    stop("Select a known, species-matched genome build before checking or downloading annotations.")
  build
}

# Older human assemblies are available as GTF-derived GRanges in AnnotationHub,
# not necessarily as EnsDb databases. Accept gene annotations, never arbitrary
# GRanges such as peaks. The downloaded object is checked again before export.
annotation_reference_kind <- function(rows) {
  kind <- rep(NA_character_,nrow(rows))
  if (!all(c("dataprovider","rdataclass","title") %in% names(rows))) return(kind)
  ensembl <- !is.na(rows$dataprovider) & rows$dataprovider=="Ensembl"
  kind[which(ensembl & rows$rdataclass=="EnsDb")] <- "EnsDb"
  gtf <- ensembl & rows$rdataclass=="GRanges" &
    grepl("\\.[0-9]+\\.gtf(\\.gz)?$",rows$title,ignore.case=TRUE)
  if ("sourcetype" %in% names(rows)) gtf <- gtf & !is.na(rows$sourcetype) & rows$sourcetype=="GTF"
  kind[which(gtf)] <- "GTF"
  kind
}

# Export unique gene-annotation combinations, not exon/transcript counts. Keep
# conflicting mappings for the existing ambiguity review rather than choosing
# the first gene name. GTFs lack Entrez links and descriptions: record these as
# unavailable; do not fabricate them from a newer, different annotation release.
annotation_gtf_table <- function(resource,candidate) {
  if (!inherits(resource,"GRanges")) stop("Downloaded GTF annotation is not a GRanges.")
  if (!identical(annotation_reference_kind(candidate),"GTF")) stop("Select an Ensembl GTF annotation.")
  columns <- S4Vectors::mcols(resource)
  if (!all(c("gene_id","gene_name","gene_biotype") %in% names(columns)))
    stop("GTF annotation lacks gene_id, gene_name or gene_biotype.")
  object_builds <- unique(as.character(GenomeInfoDb::genome(resource)))
  object_builds <- object_builds[!is.na(object_builds) & nzchar(object_builds)]
  build <- .annotation_target(candidate$species[[1]],candidate$genome[[1]])
  if (length(object_builds) && any(is.na(annotation_genome_build(object_builds)) |
      annotation_genome_build(object_builds)!=build)) stop("GTF sequence genome differs from the checked catalogue build.")
  table <- unique(data.frame(gene_id=as.character(columns$gene_id),gene_name=as.character(columns$gene_name),
    gene_biotype=as.character(columns$gene_biotype),seq_name=as.character(GenomicRanges::seqnames(resource)),
    stringsAsFactors=FALSE))
  table <- table[!is.na(table$gene_id) & nzchar(table$gene_id),,drop=FALSE]
  if (!nrow(table)) stop("GTF annotation contains no gene identifiers.")
  table <- table[order(table$gene_id,table$gene_name,table$gene_biotype,table$seq_name),,drop=FALSE]
  rownames(table) <- NULL
  table$description <- NA_character_
  table$entrezid <- NA_character_
  release <- sub("^.*\\.([0-9]+)\\.gtf(\\.gz)?$","\\1",candidate$title[[1]],ignore.case=TRUE)
  list(table=table,species=candidate$species[[1]],genome_build=build,annotation_release=release,
    source_metadata=list(resource_class="GRanges",catalogue_candidate=candidate,sequence_genomes=object_builds,
      build_evidence=if (length(object_builds)) "Catalogue and sequence genome agree" else "Catalogue build; sequence genome unavailable",
      unavailable_fields=c("description","entrezid"),export="Distinct gene-level annotation combinations; conflicts retained"))
}

bundled_annotation_provenance <- function(organism,app_root=".") {
  .annotation_species(organism)
  files <- c("Homo sapiens"="annotations_hsap.Rda","Mus musculus"="annotations_mmus.Rda")
  path <- normalizePath(file.path(app_root,files[[organism]]),mustWork=TRUE)
  md5 <- unname(tools::md5sum(path))
  result <- list(source="AnnotationHub export (origin confirmed by the user)",path=path,md5=md5,
    resource_id=NA_character_,snapshot_date=NA_character_,annotation_release=NA_character_,
    genome_build=NA_character_,provenance_status="Exact resource not recovered",reference_id="bundled")
  sidecar <- file.path(app_root,"annotation-provenance.csv")
  if (file.exists(sidecar)) {
    rows <- utils::read.csv(sidecar,stringsAsFactors=FALSE,na.strings=c("","NA"))
    hit <- rows[rows$species==organism & rows$md5==md5,,drop=FALSE]
    if (nrow(hit)==1L) for (field in intersect(names(result),names(hit))) result[[field]] <- hit[[field]][1]
    # A changed file cannot inherit provenance from the old file's sidecar.
    if (nrow(hit)!=1L) result$provenance_status <- "File differs from recorded export; resource details unverified"
    result$provenance_manifest_md5 <- unname(tools::md5sum(sidecar))
  }
  result
}

.annotation_read_entry <- function(directory) {
  manifest <- file.path(directory,"manifest.rds")
  if (!file.exists(manifest)) return(NULL)
  entry <- tryCatch(readRDS(manifest),error=function(e) NULL)
  if (!is.list(entry) || !identical(entry$schema_version,1L) ||
      !identical(entry$reference_id,basename(directory)) || !is.list(entry$provenance)) return(NULL)
  entry$path <- file.path(directory,"reference.rds")
  if (!file.exists(entry$path)) return(NULL)
  entry
}

list_annotation_references <- function(organism,app_root=".",cache_dir=annotation_cache_dir(app_root)) {
  provenance <- bundled_annotation_provenance(organism,app_root)
  result <- list(bundled=list(reference_id="bundled",label="Bundled AnnotationHub export (release/build unverified)",
    provenance=provenance,path=provenance$path))
  root <- file.path(cache_dir,"references")
  if (dir.exists(root)) for (directory in list.dirs(root,recursive=FALSE,full.names=TRUE)) {
    if (!grepl("^AH[0-9]+_",basename(directory))) next
    entry <- .annotation_read_entry(directory)
    if (!is.null(entry) && identical(entry$species,organism)) result[[entry$reference_id]] <- entry
  }
  result
}

# Validate hashes at load time, not just when listing choices. A changed cached
# reference must never acquire the old provenance or silently alter a result.
load_versioned_annotation <- function(reference_id,organism,app_root=".",cache_dir=annotation_cache_dir(app_root)) {
  entries <- list_annotation_references(organism,app_root,cache_dir)
  entry <- entries[[reference_id]]
  if (is.null(entry)) stop("Selected annotation reference is unavailable for this species.")
  if (!identical(unname(tools::md5sum(entry$path)),entry$provenance$md5)) stop("Annotation reference checksum mismatch.")
  if (reference_id=="bundled") {
    e <- new.env(parent=baseenv()); load(entry$path,envir=e)
    if (!exists("annotations",envir=e,inherits=FALSE)) stop("The cache has no annotations table.")
    table <- e$annotations
  } else {
    table <- readRDS(entry$path)
    if (!identical(digest::digest(table,algo="sha256"),entry$table_hash)) stop("Annotation table hash mismatch.")
  }
  list(table=table,species=organism,provenance=entry$provenance)
}

.annotation_check_path <- function(organism,build,cache_dir) {
  file.path(cache_dir,"checks",paste0(gsub(" ","_",organism),"_",build,".rds"))
}

annotation_update_state <- function(organism,genome_build,cache_dir,today=Sys.Date()) {
  build <- .annotation_target(organism,genome_build)
  path <- .annotation_check_path(organism,build,cache_dir)
  state <- if (file.exists(path)) readRDS(path) else list(last_success=as.Date(NA),last_attempt=as.Date(NA),candidates=NULL)
  # Monthly means one calendar month after the last successful check. Failures
  # do not advance this date or claim that the catalogue is up to date.
  next_check <- if (is.na(state$last_success)) as.Date(NA) else {
    first <- as.Date(format(state$last_success,"%Y-%m-01"))
    next_first <- seq(first,by="month",length.out=2)[2]
    after_next <- seq(next_first,by="month",length.out=2)[2]
    next_first+min(as.integer(format(state$last_success,"%d")),as.integer(after_next-next_first))-1
  }
  state$next_check <- next_check
  state$due <- is.na(next_check) || as.Date(today)>=next_check
  state
}

# This may refresh/download AnnotationHub's catalogue database into the
# project-only hub directory. No EnsDb resource is downloaded by this function.
# Results are scoped to the installed Bioconductor version and hub snapshot;
# they are not a claim to cover every provider's newest release.
annotationhub_catalog <- function(cache_dir) {
  if (!requireNamespace("AnnotationHub",quietly=TRUE)) stop("Existing AnnotationHub installation required.")
  hub <- withCallingHandlers(AnnotationHub::AnnotationHub(cache=file.path(cache_dir,"hub"),ask=FALSE,localHub=FALSE),
    message=function(m) {
      if (grepl("Cannot connect|using 'localHub=TRUE'",conditionMessage(m),ignore.case=TRUE))
        stop("Online catalogue unavailable; existing analysis references are unchanged.")
    })
  rows <- as.data.frame(S4Vectors::mcols(hub),stringsAsFactors=FALSE)
  rows$resource_id <- names(hub)
  list(rows=rows,snapshot_date=as.character(AnnotationHub::snapshotDate(hub)),
    annotationhub_version=as.character(utils::packageVersion("AnnotationHub")))
}

check_annotation_updates <- function(organism,genome_build,app_root=".",cache_dir=annotation_cache_dir(app_root),
    force=FALSE,today=Sys.Date(),catalog_fetch=annotationhub_catalog) {
  build <- .annotation_target(organism,genome_build)
  state <- annotation_update_state(organism,build,cache_dir,today)
  if (!force && !state$due) return(state)
  state$last_attempt <- as.Date(today)
  state <- tryCatch({
    catalog <- catalog_fetch(cache_dir)
    rows <- catalog$rows
    required <- c("resource_id","species","genome","rdataclass","dataprovider","title")
    if (!is.data.frame(rows) || !all(required %in% names(rows)) ||
        length(catalog$snapshot_date)!=1L || is.na(as.Date(catalog$snapshot_date))) stop("Incomplete AnnotationHub catalogue metadata.")
    compatible <- !is.na(rows$species) & rows$species==organism &
      !is.na(annotation_reference_kind(rows)) &
      !is.na(annotation_genome_build(rows$genome)) & annotation_genome_build(rows$genome)==build
    selected <- rows[which(compatible),intersect(c(required,"rdatadateadded","description","sourceurl","sourcesize","sourcetype","rdatapath"),names(rows)),drop=FALSE]
    if (anyNA(selected$resource_id) || any(!grepl("^AH[0-9]+$",selected$resource_id)) || anyDuplicated(selected$resource_id))
      stop("Invalid or duplicated AnnotationHub resource identifiers.")
    selected$genome <- rep(build,nrow(selected))
    if ("rdatadateadded" %in% names(selected)) selected <- selected[order(selected$rdatadateadded,decreasing=TRUE),,drop=FALSE]
    state$candidates <- selected
    state$snapshot_date <- catalog$snapshot_date
    state$annotationhub_version <- catalog$annotationhub_version
    state$last_success <- as.Date(today)
    state$error <- NULL
    state
  },error=function(e) {state$error<-conditionMessage(e);state})
  path <- .annotation_check_path(organism,build,cache_dir)
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  # Rename a complete small record; readers never see a half-written RDS.
  tmp <- tempfile("check-",tmpdir=dirname(path)); on.exit(unlink(tmp),add=TRUE)
  saveRDS(state,tmp)
  if (!file.rename(tmp,path)) stop("Cannot save annotation check history.")
  annotation_update_state(organism,build,cache_dir,today)
}

# Retrieve exactly the chosen resource in the catalogue snapshot. Selecting a
# new catalogue entry is not enough: resource metadata must independently agree.
annotationhub_fetch <- function(candidate,snapshot_date,cache_dir,hub=NULL) {
  needed <- c("AnnotationHub","ensembldb","AnnotationDbi","GenomicRanges","GenomeInfoDb")
  if (!all(vapply(needed,requireNamespace,logical(1),quietly=TRUE))) stop("Existing AnnotationHub/ensembldb/AnnotationDbi installations required; no packages installed.")
  # Ordinary reference loading stays offline. A caller may supply an offline
  # AnnotationHub to reuse an already cached resource without downloading it.
  if (is.null(hub)) hub <- AnnotationHub::AnnotationHub(cache=file.path(cache_dir,"hub"),ask=FALSE,localHub=FALSE)
  if (!identical(as.character(AnnotationHub::snapshotDate(hub)),as.character(as.Date(snapshot_date))))
    AnnotationHub::snapshotDate(hub) <- as.character(as.Date(snapshot_date))
  id <- candidate$resource_id
  if (!id %in% names(hub)) stop("Chosen resource no longer available in the checked snapshot.")
  catalog_row <- as.data.frame(S4Vectors::mcols(hub[id]),stringsAsFactors=FALSE)
  if (!identical(as.character(catalog_row$species),candidate$species) ||
      !identical(annotation_genome_build(catalog_row$genome),annotation_genome_build(candidate$genome)) ||
      !identical(annotation_reference_kind(catalog_row),annotation_reference_kind(candidate)) ||
      !identical(as.character(catalog_row$title),candidate$title)) stop("Resource catalogue metadata changed; check again before downloading.")
  resource <- hub[[id]]
  if (identical(annotation_reference_kind(candidate),"GTF")) return(annotation_gtf_table(resource,candidate))
  if (!inherits(resource,"EnsDb")) stop("Downloaded reference is not an EnsDb.")
  metadata <- AnnotationDbi::metadata(resource)
  on.exit(try(DBI::dbDisconnect(AnnotationDbi::dbconn(resource)),silent=TRUE),add=TRUE)
  value <- function(key) {x<-metadata$value[metadata$name==key];if(length(x)==1L) as.character(x) else NA_character_}
  table <- as.data.frame(ensembldb::genes(resource,return.type="data.frame"))
  list(table=table[,c("gene_id","gene_name","gene_biotype","seq_name","description","entrezid"),drop=FALSE],
    species=value("Organism"),genome_build=value("genome_build"),annotation_release=value("ensembl_version"),
    source_metadata=metadata)
}

download_annotation_reference <- function(candidate,snapshot_date,reason,app_root=".",cache_dir=annotation_cache_dir(app_root),
    resource_fetch=annotationhub_fetch) {
  if (!is.data.frame(candidate) || nrow(candidate)!=1L || !all(c("resource_id","species","genome","title","rdataclass","dataprovider") %in% names(candidate)))
    stop("Select exactly one checked catalogue candidate.")
  id <- candidate$resource_id[[1]]
  kind <- annotation_reference_kind(candidate)
  if (is.na(id) || !grepl("^AH[0-9]+$",id) || is.na(kind)) stop("Select an Ensembl EnsDb or GTF gene-annotation resource.")
  build <- .annotation_target(candidate$species[[1]],candidate$genome[[1]])
  if (length(snapshot_date)!=1L || is.na(as.Date(snapshot_date))) stop("Record the checked snapshot date.")
  if (length(reason)!=1L || is.na(reason) || !nzchar(trimws(reason))) stop("Provide a reason before downloading a reference.")
  payload <- resource_fetch(candidate,snapshot_date,cache_dir)
  if (!identical(payload$species,candidate$species[[1]]) ||
      !identical(annotation_genome_build(payload$genome_build),build)) stop("Downloaded reference species/build does not match the reviewed candidate.")
  required <- c("gene_id","gene_name","gene_biotype","seq_name","description","entrezid")
  if (!is.data.frame(payload$table) || !nrow(payload$table) || !all(required %in% names(payload$table))) stop("Downloaded annotation table lacks required fields.")
  if (length(payload$annotation_release)!=1L || is.na(payload$annotation_release) || !nzchar(payload$annotation_release))
    stop("Downloaded reference lacks an annotation release; nothing selected or overwritten.")
  table <- payload$table[,required,drop=FALSE]
  table_hash <- digest::digest(table,algo="sha256")
  reference_id <- paste(id,as.character(as.Date(snapshot_date)),substr(table_hash,1,16),sep="_")
  root <- file.path(cache_dir,"references"); dir.create(root,recursive=TRUE,showWarnings=FALSE)
  destination <- file.path(root,reference_id)
  if (dir.exists(destination)) {
    entry <- .annotation_read_entry(destination)
    if (is.null(entry) || !identical(entry$table_hash,table_hash) ||
        !identical(unname(tools::md5sum(entry$path)),entry$provenance$md5)) stop("Existing cached reference failed verification; nothing overwritten.")
    return(entry)
  }
  staging <- tempfile("pending-",tmpdir=root); dir.create(staging)
  on.exit(unlink(staging,recursive=TRUE),add=TRUE) # only this function's unfinished files
  path <- file.path(staging,"reference.rds"); saveRDS(table,path)
  provenance <- list(source=paste("AnnotationHub Ensembl",kind),reference_id=reference_id,resource_id=id,
    snapshot_date=as.character(as.Date(snapshot_date)),annotation_release=payload$annotation_release,
    genome_build=build,md5=unname(tools::md5sum(path)),path=file.path(destination,"reference.rds"),
    retrieved_at=format(Sys.time(),tz="UTC",usetz=TRUE),reason=trimws(reason),
    source_metadata=payload$source_metadata,catalogue_candidate=candidate)
  entry <- list(schema_version=1L,reference_id=reference_id,species=payload$species,table_hash=table_hash,
    label=paste(build,"— Ensembl",payload$annotation_release,paste0("(",id,", ",kind,")"),"snapshot",snapshot_date),provenance=provenance)
  saveRDS(entry,file.path(staging,"manifest.rds"))
  if (!file.rename(staging,destination)) stop("Cannot save versioned reference; existing references unchanged.")
  .annotation_read_entry(destination)
}
