# Compose QC evidence into one review stage before choosing exclusions.
# Source analysis_runner.R, legacy_qc_adapter.R and qc_optional_checks.R first.
preprocessing_qc_step <- function(object,parameters) {
  legacy <- make_legacy_qc_steps(parameters$app_root)
  result <- legacy$metadata(object,parameters$legacy)
  checks <- list()
  ambient_method <- if (identical(parameters$ambient$method,"DecontX")) decontx_qc_step else ambient_rna_readiness_step
  methods <- list(cell_cycle=cell_cycle_qc_step,doublets=doublet_qc_step,ambient=ambient_method)
  for (name in names(methods)) {
    if (is.null(parameters[[name]])) {
      checks[[name]] <- list(status="disabled",excluded_cells=0L)
    } else {
      step <- methods[[name]](result$object,parameters[[name]])
      result$object <- step$object
      checks[[name]] <- step$diagnostics
    }
  }
  result$diagnostics$checks <- checks
  result$diagnostics$optional_checks <- do.call(rbind,lapply(names(checks),function(name) {
    x <- checks[[name]]
    data.frame(check=name,status=x$status,notes=if (is.null(x$reasons)) {
      if (identical(name,"cell_cycle") && identical(x$status,"scored")) {
        "S/G2M scores and predicted phases; no regression or exclusions"
      } else if (identical(name,"ambient") && identical(x$status,"estimated")) x$interpretation else ""
    } else paste(x$reasons,collapse=" | "))
  }))
  directory <- parameters$diagnostic_dir
  if (!is.null(directory)) {
    if (dir.exists(directory)) stop("Use a new diagnostic directory.")
    dir.create(directory,recursive=TRUE)
    directory <- normalizePath(directory)
    metrics <- as.data.frame(SummarizedExperiment::colData(result$object))
    columns <- intersect(c("cells","sample","nUMI","nGene","mitoRatio","rbcRatio","log10GenesPerUMI",
      "cell_cycle_S_score","cell_cycle_G2M_score","cell_cycle_phase","scDblFinder.score","scDblFinder.class",
      "ambient_decontX_fraction","ambient_decontX_preliminary_cluster"),names(metrics))
    metric_path <- file.path(directory,"qc-check-metrics.csv")
    status_path <- file.path(directory,"optional-check-status.csv")
    utils::write.csv(metrics[,columns,drop=FALSE],metric_path,row.names=FALSE)
    utils::write.csv(result$diagnostics$optional_checks,status_path,row.names=FALSE)
    result$artifacts <- c(result$artifacts,metric_path,status_path)
    if (identical(checks$cell_cycle$status,"scored")) {
      path <- file.path(directory,"cell-cycle-scores.png")
      plot <- ggplot2::ggplot(metrics,ggplot2::aes(x=cell_cycle_S_score,y=cell_cycle_G2M_score,color=cell_cycle_phase)) +
        ggplot2::geom_point(size=0.9,alpha=0.6) + ggplot2::theme_bw() +
        ggplot2::labs(title="Cell-cycle program scores · predicted phases",x="S score",y="G2M score",color="Predicted phase")
      ggplot2::ggsave(path,plot=plot,width=9,height=6,dpi=110)
      result$artifacts <- c(result$artifacts,path)
      reference_path <- file.path(directory,"cell-cycle-marker-coverage.csv")
      coverage <- checks$cell_cycle$marker_mapping
      coverage$in_matrix <- coverage$status=="matched"
      utils::write.csv(coverage,reference_path,row.names=FALSE)
      annotation_path <- file.path(directory,"cell-cycle-annotation-provenance.txt")
      writeLines(capture.output(dput(checks$cell_cycle$annotation)),annotation_path)
      result$artifacts <- c(result$artifacts,reference_path,annotation_path)
    }
    if (identical(checks$doublets$status,"scored")) {
      path <- file.path(directory,"doublet-scores.png")
      plot <- ggplot2::ggplot(metrics,ggplot2::aes(x=nUMI,y=scDblFinder.score,color=scDblFinder.class)) +
        ggplot2::geom_point(size=0.9,alpha=0.6) + ggplot2::scale_x_log10() + ggplot2::theme_bw() +
        ggplot2::labs(title="Doublet predictions · no cells excluded",x="UMIs",y="scDblFinder score",color="Prediction")
      ggplot2::ggsave(path,plot=plot,width=9,height=6,dpi=110)
      result$artifacts <- c(result$artifacts,path)
    }
    if (identical(checks$ambient$status,"estimated")) {
      path <- file.path(directory,"ambient-decontx-fraction.png")
      plot <- ggplot2::ggplot(metrics,ggplot2::aes(x=nUMI,y=ambient_decontX_fraction,color=ambient_decontX_preliminary_cluster)) +
        ggplot2::geom_point(size=0.9,alpha=0.6) + ggplot2::scale_x_log10() + ggplot2::theme_bw() +
        ggplot2::labs(title="DecontX estimates · original counts retained",x="UMIs",y="Estimated contamination fraction",color="Preliminary cluster")
      ggplot2::ggsave(path,plot=plot,width=9,height=6,dpi=110)
      result$artifacts <- c(result$artifacts,path)
    }
  }
  result
}
