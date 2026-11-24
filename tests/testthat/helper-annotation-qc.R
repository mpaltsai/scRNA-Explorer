source(file.path("..","..","site_data","R","gene_annotations.R"))
source(file.path("..","..","site_data","R","annotation_reference_server.R"))
source(file.path("..","..","site_data","R","annotation_qc_server.R"))

annotation_fixture <- function() {
  list(table=data.frame(gene_id=c("ENSG1","ENSG2","ENSG3","ENSG4","ENSG5"),
    gene_name=c("MT-A","HBA","GeneA","Shared","Shared"),
    gene_biotype=rep("protein_coding",5),seq_name=c("MT","1","1","2","3"),
    description=c("mitochondrial gene","hemoglobin test","ordinary","ordinary","ordinary")),
    species="Homo sapiens",provenance=list(source="Synthetic test reference"))
}
