# GRCh37 enrichment is explicit maintenance, never part of ordinary annotation
# loading. Gene models remain those of the pinned AnnotationHub reference.
# Only missing descriptions/Entrez links are filled; expression is never read.
.annotation_grch37_parent <- function(reference) {
  p <- reference$provenance
  if (!identical(reference$species,"Homo sapiens") ||
      !identical(p$genome_build,"GRCh37") || !identical(p$annotation_release,"75") ||
      !identical(p$source,"AnnotationHub Ensembl GTF") || !is.null(p$supplement))
    stop("Choose an unsupplemented human GRCh37 Ensembl 75 AnnotationHub GTF reference.")
  invisible(TRUE)
}

.annotation_validate_biomart <- function(payload) {
  m <- payload$metadata
  if (!is.list(m) || !identical(m$species,"Homo sapiens") ||
      !identical(m$genome_build,"GRCh37") || !identical(m$annotation_release,"75") ||
      !identical(m$dataset,"hsapiens_gene_ensembl") ||
      length(m$host)!=1L || is.na(m$host) ||
      !grepl("^https://grch37\\.ensembl\\.org(:443)?/biomart/martservice/?$",m$host) ||
      length(m$dataset_version)!=1L || is.na(m$dataset_version) ||
      !grepl("^GRCh37(\\.p[0-9]+)?$",m$dataset_version) ||
      length(m$mart_version)!=1L || is.na(m$mart_version) ||
      !grepl("^Ensembl Genes [0-9]+$",m$mart_version) ||
      !identical(m$release_evidence_url,"https://jun2026.archive.ensembl.org/info/website/tutorials/grch37.html"))
    stop("BioMart metadata must verify human GRCh37 and Ensembl release 75; no supplement saved.")
  d <- payload$rows
  if (!is.data.frame(d) || !nrow(d) || !all(c("gene_id","description","entrezid") %in% names(d)))
    stop("BioMart returned no usable description/Entrez table.")
  d <- d[,c("gene_id","description","entrezid"),drop=FALSE]
  d[] <- lapply(d,as.character)
  if (anyNA(d$gene_id) || any(!nzchar(d$gene_id))) stop("BioMart returned empty gene identifiers.")
  d$description[!is.na(d$description) & !nzchar(trimws(d$description))] <- NA_character_
  d$entrezid[!is.na(d$entrezid) & !nzchar(trimws(d$entrezid))] <- NA_character_
  if (any(!is.na(d$entrezid) & !grepl("^[0-9]+$",d$entrezid))) stop("BioMart returned invalid Entrez identifiers.")
  d <- unique(d)
  d <- d[order(d$gene_id,d$description,d$entrezid,na.last=TRUE),,drop=FALSE]
  rownames(d) <- NULL
  d
}

# Use BioMart's documented XML interface through existing R packages. Disable
# redirects so a service redirect cannot silently change assembly. No global
# HTTP/SSL options or package settings are modified. A completion stamp below
# distinguishes a complete query from an interrupted/partial HTTP 200 response.
.annotation_biomart_request <- function(parameters) {
  if (!requireNamespace("httr2",quietly=TRUE)) stop("Existing httr2 installation required; no installation attempted.")
  request <- httr2::request("https://grch37.ensembl.org/biomart/martservice")
  request <- do.call(httr2::req_url_query,c(list(request),parameters,list(redirect="no")))
  request <- httr2::req_options(request,followlocation=FALSE)
  request <- httr2::req_user_agent(request,"scRNA-Explorer annotation maintenance (R httr2)")
  request <- httr2::req_timeout(request,120)
  response <- httr2::req_perform(request)
  if (httr2::resp_status(response)!=200L) stop("GRCh37 BioMart did not return HTTP 200; no redirect followed.")
  httr2::resp_body_string(response)
}

grch37_biomart_fetch <- function(request=.annotation_biomart_request) {
  if (!requireNamespace("xml2",quietly=TRUE)) stop("Existing xml2 installation required; no installation attempted.")
  registry_text <- request(list(type="registry"))
  registry <- xml2::read_xml(registry_text)
  entry <- xml2::xml_find_all(registry,".//MartURLLocation[@name='ENSEMBL_MART_ENSEMBL']")
  if (length(entry)!=1L || !identical(xml2::xml_attr(entry,"host"),"grch37.ensembl.org"))
    stop("Cannot verify the GRCh37 BioMart registry host.")
  tsv <- function(text,fill=FALSE) {
    lines <- strsplit(text,"\n",fixed=TRUE)[[1]]
    lines <- lines[nzchar(trimws(lines))]
    utils::read.table(text=paste(lines,collapse="\n"),sep="\t",header=FALSE,quote="",comment.char="",
      fill=fill,colClasses="character",stringsAsFactors=FALSE,na.strings="")
  }
  dataset_text <- request(list(type="datasets",mart="ENSEMBL_MART_ENSEMBL"))
  datasets <- tsv(dataset_text,fill=TRUE)
  if (ncol(datasets)<5L) stop("Incomplete GRCh37 BioMart dataset metadata.")
  dataset <- datasets[which(datasets[[2]]=="hsapiens_gene_ensembl"),,drop=FALSE]
  if (nrow(dataset)!=1L) stop("Cannot verify the human BioMart dataset.")
  metadata <- list(source="Ensembl GRCh37 BioMart",species="Homo sapiens",genome_build="GRCh37",
    annotation_release="75",host="https://grch37.ensembl.org/biomart/martservice",dataset="hsapiens_gene_ensembl",
    dataset_version=dataset[[5]],mart_version=xml2::xml_attr(entry,"displayName"),
    release_evidence_url="https://jun2026.archive.ensembl.org/info/website/tutorials/grch37.html",
    release_basis="Ensembl documents GRCh37 gene annotation as frozen at release 75; BioMart and cross-reference releases are recorded separately, not assumed to be historical release-75 snapshots",
    registry_entry=as.character(entry),dataset_metadata=as.list(dataset[1,,drop=FALSE]),
    retrieved_at=format(Sys.time(),tz="UTC",usetz=TRUE),
    transport="BioMart XML API via R httr2",httr2_version=as.character(utils::packageVersion("httr2")))
  # Validate metadata before making the potentially larger gene query.
  .annotation_validate_biomart(list(metadata=metadata,
    rows=data.frame(gene_id="validation",description=NA_character_,entrezid=NA_character_)))
  attribute_text <- request(list(type="attributes",dataset="hsapiens_gene_ensembl"))
  attributes <- tsv(attribute_text,fill=TRUE)[[1]]
  entrez <- intersect(c("entrezgene_id","entrezgene"),attributes)
  if (!all(c("ensembl_gene_id","description") %in% attributes) || !length(entrez))
    stop("GRCh37 BioMart does not expose the required annotation fields.")
  selected <- c("ensembl_gene_id","description",entrez[1])
  query <- xml2::xml_new_root("Query",virtualSchemaName="default",formatter="TSV",header="0",
    uniqueRows="1",datasetConfigVersion="0.6",completionStamp="1")
  dataset_node <- xml2::xml_add_child(query,"Dataset",name="hsapiens_gene_ensembl",interface="default")
  for (field in selected) xml2::xml_add_child(dataset_node,"Attribute",name=field)
  xml <- gsub("[\r\n]","",as.character(query))
  result <- request(list(query=xml))
  lines <- strsplit(result,"\n",fixed=TRUE)[[1]]
  lines <- sub("\r$","",lines)
  while (length(lines) && !nzchar(utils::tail(lines,1L))) lines <- utils::head(lines,-1L)
  if (!length(lines) || !identical(utils::tail(lines,1L),"[success]")) stop("BioMart query incomplete or unsuccessful; no supplement saved.")
  rows <- tsv(paste(utils::head(lines,-1L),collapse="\n"))
  if (ncol(rows)!=3L) stop("Unexpected BioMart query columns.")
  names(rows) <- c("gene_id","description","entrezid")
  metadata$attributes <- selected
  metadata$metadata_response_hashes <- vapply(list(registry=registry_text,datasets=dataset_text,attributes=attribute_text),
    digest::digest,character(1),algo="sha256")
  metadata$query_xml <- xml
  metadata$query_scope <- "All human genes; join to the parent reference by exact Ensembl gene ID"
  result <- list(rows=rows,metadata=metadata)
  result$rows <- .annotation_validate_biomart(result)
  result
}

# Keep parent rows/order and all Entrez candidates. A contradictory description
# remains unresolved. Existing nonempty values are never overwritten, even when
# the supplement disagrees; the review records both candidate sets.
merge_grch37_biomart <- function(reference,payload) {
  .annotation_grch37_parent(reference)
  rows <- .annotation_validate_biomart(payload)
  table <- reference$table
  required <- c("gene_id","gene_name","gene_biotype","seq_name","description","entrezid")
  if (!is.data.frame(table) || !all(required %in% names(table))) stop("Parent reference lacks required annotation fields.")
  groups <- split(seq_len(nrow(rows)),rows$gene_id)
  values <- function(x) sort(unique(as.character(x[!is.na(x) & nzchar(trimws(as.character(x)))])))
  descriptions <- lapply(groups,function(i) values(rows$description[i]))
  entrez <- lapply(groups,function(i) values(rows$entrezid[i]))
  # List-column links retain one-to-many gene/Entrez relationships without
  # duplicating expression features or selecting the first mapping.
  links <- lapply(seq_len(nrow(table)),function(i) values(unlist(table$entrezid[i],use.names=FALSE)))
  original_description <- description <- as.character(table$description)
  original_entrezid <- vapply(links,paste,character(1),collapse=";")
  desc_status <- id_status <- desc_candidates <- id_candidates <- character(nrow(table))
  # Match once: repeated name-based list lookup and data-frame column updates
  # inside a 60,000-gene loop otherwise make explicit maintenance needlessly slow.
  matched <- match(as.character(table$gene_id),names(groups))
  for (i in seq_len(nrow(table))) {
    j <- matched[i]
    desc <- if (is.na(j)) character() else descriptions[[j]]
    ids <- if (is.na(j)) character() else entrez[[j]]
    old <- values(description[i])
    desc_candidates[i] <- paste(desc,collapse=" | ")
    id_candidates[i] <- paste(ids,collapse=";")
    if (length(old)) {
      desc_status[i] <- if (length(desc) && !setequal(old,desc)) "existing_conflict_preserved" else "existing_preserved"
    } else if (length(desc)==1L) {
      description[i] <- desc; desc_status[i] <- "filled"
    } else desc_status[i] <- if (length(desc)) "conflicting_unresolved" else "unmapped"
    if (length(links[[i]])) {
      id_status[i] <- if (length(ids) && !setequal(links[[i]],ids)) "existing_conflict_preserved" else "existing_preserved"
    } else if (length(ids)) {
      links[[i]] <- ids; id_status[i] <- if (length(ids)>1L) "filled_multiple" else "filled"
    } else id_status[i] <- "unmapped"
  }
  table$description <- description
  table$entrezid <- I(links)
  review <- data.frame(source_row=seq_len(nrow(table)),gene_id=as.character(table$gene_id),
    original_description=original_description,original_entrezid=original_entrezid,
    description_status=desc_status,entrezid_status=id_status,
    description_candidates=desc_candidates,entrezid_candidates=id_candidates,stringsAsFactors=FALSE)
  list(table=table,review=review,rows=rows,summary=list(
    parent_rows=nrow(table),descriptions_filled=sum(review$description_status=="filled"),
    descriptions_conflicting=sum(review$description_status=="conflicting_unresolved"),
    descriptions_unmapped=sum(review$description_status=="unmapped"),
    entrez_links_filled=sum(review$entrezid_status %in% c("filled","filled_multiple")),
    entrez_multiple=sum(review$entrezid_status=="filled_multiple"),
    entrez_unmapped=sum(review$entrezid_status=="unmapped")))
}

# Save a separate immutable derivative. The parent table and saved analyses are
# untouched. Ordinary loading of either reference remains completely offline.
supplement_grch37_reference <- function(reference_id,reason,app_root=".",cache_dir=annotation_cache_dir(app_root),
    supplement_fetch=grch37_biomart_fetch) {
  if (length(reason)!=1L || is.na(reason) || !nzchar(trimws(reason))) stop("Provide a reason for the BioMart supplement.")
  reference <- load_versioned_annotation(reference_id,"Homo sapiens",app_root,cache_dir)
  .annotation_grch37_parent(reference)
  payload <- supplement_fetch()
  merged <- merge_grch37_biomart(reference,payload)
  table_hash <- digest::digest(merged$table,algo="sha256")
  source <- payload$metadata
  source$response_hash <- digest::digest(merged$rows,algo="sha256")
  identity <- source[setdiff(names(source),"retrieved_at")]
  suffix <- substr(digest::digest(list(table_hash,identity),algo="sha256"),1,16)
  id <- paste(reference_id,"biomart",suffix,sep="_")
  root <- file.path(cache_dir,"references"); dir.create(root,recursive=TRUE,showWarnings=FALSE)
  destination <- file.path(root,id)
  if (dir.exists(destination)) {
    entry <- .annotation_read_entry(destination)
    load_versioned_annotation(id,"Homo sapiens",app_root,cache_dir)
    load_annotation_supplement_review(id,app_root,cache_dir)
    return(entry)
  }
  staging <- tempfile("supplement-",tmpdir=root);dir.create(staging)
  on.exit(unlink(staging,recursive=TRUE),add=TRUE)
  saveRDS(merged$table,file.path(staging,"reference.rds"))
  saveRDS(list(rows=merged$rows,metadata=source),file.path(staging,"supplement.rds"))
  utils::write.csv(merged$review,file.path(staging,"supplement-review.csv"),row.names=FALSE,na="")
  p <- reference$provenance
  p$parent_reference <- list(reference_id=reference_id,md5=p$md5,source=p$source,
    resource_id=p$resource_id,snapshot_date=p$snapshot_date)
  p$source <- "AnnotationHub Ensembl GTF + Ensembl GRCh37 BioMart"
  p$reference_id <- id; p$path <- file.path(destination,"reference.rds")
  p$md5 <- unname(tools::md5sum(file.path(staging,"reference.rds")))
  p$supplement <- c(source,list(summary=merged$summary,reason=trimws(reason),
    policy="Fill missing fields only; exact gene-ID joins; conflicting descriptions unresolved; all Entrez candidates retained"))
  entry <- list(schema_version=1L,reference_id=id,species="Homo sapiens",table_hash=table_hash,
    label=paste("GRCh37 — Ensembl 75",paste0("(",p$resource_id,", GTF + BioMart)")),provenance=p,
    review_md5=unname(tools::md5sum(file.path(staging,"supplement-review.csv"))),
    supplement_md5=unname(tools::md5sum(file.path(staging,"supplement.rds"))))
  saveRDS(entry,file.path(staging,"manifest.rds"))
  if (!file.rename(staging,destination)) stop("Cannot save supplemented reference; existing references unchanged.")
  .annotation_read_entry(destination)
}

load_annotation_supplement_review <- function(reference_id,app_root=".",cache_dir=annotation_cache_dir(app_root)) {
  entries <- list_annotation_references("Homo sapiens",app_root,cache_dir)
  entry <- entries[[reference_id]]
  if (is.null(entry) || is.null(entry$review_md5)) stop("Selected reference has no BioMart supplement review.")
  load_versioned_annotation(reference_id,"Homo sapiens",app_root,cache_dir)
  directory <- dirname(entry$path)
  path <- file.path(directory,"supplement-review.csv")
  if (!identical(unname(tools::md5sum(path)),entry$review_md5) ||
      !identical(unname(tools::md5sum(file.path(directory,"supplement.rds"))),entry$supplement_md5))
    stop("BioMart supplement audit checksum mismatch.")
  utils::read.csv(path,stringsAsFactors=FALSE,na.strings="")
}
