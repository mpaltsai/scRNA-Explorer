# Demonstration step for the decision runner. No thresholds or exclusions.
qc_metric_step <- function(object, parameters) {
  if (inherits(object, "SingleCellExperiment")) {
    counts <- SummarizedExperiment::assay(object, "counts")
  } else if (inherits(object, "Seurat")) {
    layers <- SeuratObject::Layers(object[["RNA"]], search = "^counts($|\\.)")
    if (length(layers) != 1L) stop("The example requires one RNA counts layer.")
    counts <- SeuratObject::LayerData(object, assay = "RNA", layer = layers)
  } else stop("Supply an existing SCE or Seurat object.")
  umi <- Matrix::colSums(counts)
  genes <- Matrix::colSums(counts > 0)
  if (!length(umi)) stop("The example requires at least one cell.")
  object$review_nUMI <- as.numeric(umi)
  object$review_nGene <- as.numeric(genes)
  plot_path <- parameters$plot_path
  if (!is.character(plot_path) || length(plot_path) != 1L || is.na(plot_path) || file.exists(plot_path)) {
    stop("Supply a new plot_path; existing evidence is never overwritten.")
  }
  grDevices::pdf(plot_path, width = 12, height = 4)
  on.exit(grDevices::dev.off())
  graphics::par(mfrow = c(1, 3))
  graphics::hist(log10(umi + 1), main = "Library size", xlab = "log10(UMIs + 1)", col = "grey80")
  graphics::hist(log10(genes + 1), main = "Detected genes", xlab = "log10(genes + 1)", col = "grey80")
  graphics::plot(log10(umi + 1), log10(genes + 1), pch = 16, cex = 0.4,
                 xlab = "log10(UMIs + 1)", ylab = "log10(genes + 1)", main = "Complexity")
  list(object = object,
       diagnostics = list(cells = ncol(counts), features = nrow(counts),
                          total_umis = sum(umi), median_umis = stats::median(umi),
                          median_genes = stats::median(genes), excluded_cells = 0L),
       artifacts = plot_path)
}
