cell_cycle_fixture_annotation <- function(reference) {
  genes <- unique(c(reference$s_genes,reference$g2m_genes))
  list(table=data.frame(gene_name=genes,gene_id=sprintf("ENSG%06d",seq_along(genes))),
    species=reference$species,feature_id_type="gene_name",provenance=list(source="Synthetic test annotation"))
}
