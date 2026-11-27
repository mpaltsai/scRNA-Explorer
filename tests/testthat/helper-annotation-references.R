# Shared synthetic catalogue/payload; reference tests perform no downloads.
reference_cache_fixture <- function() {
  root <- tempfile("reference-cache-"); dir.create(root)
  app_root <- normalizePath(file.path("..","..","site_data"))
  api <- .annotation_cache_api(app_root)
  rows <- data.frame(resource_id=c("AH101","AH102","AH103","AH104","AH105","AH106"),
    species=c(rep("Homo sapiens",4),"Mus musculus","Homo sapiens"),
    genome=c("GRCh37","hg19","GRCh38",NA,"GRCm38","GRCh37"),
    rdataclass=c(rep("EnsDb",5),"OrgDb"),dataprovider=rep("Ensembl",6),
    title=paste("Ensembl",101:106),rdatadateadded=as.character(as.Date("2025-01-01")+1:6))
  fetch <- function(cache_dir) list(rows=rows,snapshot_date="2026-04-23",annotationhub_version="test")
  payload <- function(candidate,snapshot_date,cache_dir) {
    a <- annotation_fixture()$table; a$entrezid <- I(lapply(seq_len(nrow(a)),as.integer))
    list(table=a,species="Homo sapiens",genome_build="GRCh37",annotation_release="105",
      source_metadata=data.frame(name=c("Organism","genome_build","ensembl_version"),value=c("Homo sapiens","GRCh37","105")))
  }
  list(root=root,app_root=app_root,api=api,rows=rows,fetch=fetch,payload=payload)
}
