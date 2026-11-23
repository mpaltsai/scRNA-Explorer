registry_module <- file.path("..", "..", "analysis", "R", "gse176031_registry.R")
if (file.exists(registry_module)) source(registry_module)

registry_fixture <- function() {
  root <- tempfile("gse176031-fixture-")
  dir.create(root)
  data_root <- file.path(root, "dge")
  dir.create(data_root)
  write_dge <- function(name, cells = c("shared", "cell-b")) {
    con <- gzfile(file.path(data_root, name), "wt")
    on.exit(close(con))
    writeLines(c(paste(c("GENE", cells), collapse = "\t"),
                 paste(c("A1BG", rep(1, length(cells))), collapse = "\t")), con)
  }
  write_dge("GSM101_tissue_dge.txt.gz")
  write_dge("GSM102_organoid_dge.txt.gz")
  metadata <- file.path(root, "series.txt.gz")
  write_metadata <- function(extra = character()) {
    con <- gzfile(metadata, "wt")
    on.exit(close(con))
    quote_row <- function(key, values) {
      paste(c(key, paste0('"', values, '"')), collapse = "\t")
    }
    writeLines(c(
      quote_row("!Series_geo_accession", "GSE176031"),
      quote_row("!Series_overall_design", "Seq-Well S^3 protocol"),
      quote_row("!Sample_geo_accession", c("GSM102", "GSM101")),
      quote_row("!Sample_title", c("organoid", "tissue")),
      quote_row("!Sample_characteristics_ch1", c("patient: P1", "patient: P1")),
      quote_row("!Sample_characteristics_ch1", c("tissue: prostate", "tissue: prostate")),
      quote_row("!Sample_source_name_ch1", c("Organoid", "Radical Prostatectomy Paired Normal")),
      quote_row("!Sample_organism_ch1", rep("Homo sapiens", 2)),
      quote_row("!Sample_platform_id", rep("GPL11154", 2)),
      extra, "!series_matrix_table_begin", "!series_matrix_table_end"
    ), con)
  }
  write_metadata()
  list(root = root, data_root = data_root, metadata = metadata,
       write_dge = write_dge, write_metadata = write_metadata)
}
