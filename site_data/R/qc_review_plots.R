# Plot functions consume metric tables only. Facets keep baseline/retained
# scales identical. Percentages are display units; stored fractions stay 0..1.
qc_comparison_metrics <- function(metrics,selection=NULL,after_label="Preview retained") {
  before <- metrics
  before$stage <- "Baseline (UMIs >= 100)"
  before$selection <- if (is.null(selection)) "Not previewed" else ifelse(selection$keep,"Retained","Excluded")
  before$reasons <- if (is.null(selection)) "" else selection$per_cell$reasons
  if (is.null(selection)) return(before)
  after <- before[selection$keep,,drop=FALSE];after$stage<-rep(after_label,nrow(after))
  result <- rbind(before,after)
  result$stage <- factor(result$stage,levels=c("Baseline (UMIs >= 100)",after_label))
  result
}

.qc_long_metrics <- function(metrics,fields,labels=fields) {
  do.call(rbind,lapply(seq_along(fields),function(i) {
    x <- metrics[[fields[i]]]
    if (fields[i] %in% c("nUMI","nGene")) x <- log10(x)
    if (grepl("^(mitoRatio|rbcRatio|ribo_|annotation_)",fields[i])) x <- 100*x
    stages <- if("stage" %in% names(metrics)) as.character(metrics$stage) else rep("Baseline",nrow(metrics))
    stage_levels <- if(is.factor(metrics$stage)) levels(metrics$stage) else unique(stages)
    data.frame(sample=metrics$sample,stage=factor(stages,levels=stage_levels),
      metric=labels[i],value=x)
  }))
}

qc_plot_overview <- function(metrics,compare=FALSE) {
  d <- .qc_long_metrics(metrics,c("nUMI","nGene","mitoRatio","rbcRatio","log10GenesPerUMI"),
    c("log10 UMIs","log10 detected genes","Mitochondrial %","Hemoglobin %","log10(genes) / log10(UMIs)"))
  d$category <- if (compare) d$stage else d$sample
  ggplot2::ggplot(d,ggplot2::aes(x=category,y=value))+
    ggplot2::geom_boxplot(outlier.alpha=0.25,na.rm=TRUE)+
    ggplot2::scale_x_discrete(drop=FALSE)+
    ggplot2::facet_wrap(~metric,scales="free_y",ncol=3)+ggplot2::theme_bw()+
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=25,hjust=1))+
    ggplot2::labs(x=NULL,y=NULL,title=if(compare) "Selected sample: baseline and retained cells" else "Baseline QC across samples")
}

qc_joint_plot_data <- function(metrics,x="nGene") {
  if (!x %in% c("nUMI","nGene")) stop("Choose UMIs or detected genes.")
  metrics$qc_x <- log10(metrics[[x]]);metrics$qc_y <- 100*metrics$mitoRatio
  metrics
}

qc_plot_joint <- function(data,thresholds,x="nGene") {
  p <- ggplot2::ggplot(data,ggplot2::aes(x=qc_x,y=qc_y,color=selection))+
    ggplot2::geom_point(alpha=0.55,size=1,na.rm=TRUE)+ggplot2::facet_wrap(~stage,drop=FALSE)+
    ggplot2::scale_color_manual(values=c(Retained="#287d5e",Excluded="#c54b42",`Not previewed`="#426c9a"),drop=FALSE)+
    ggplot2::theme_bw()+ggplot2::labs(x=if(x=="nGene") "log10 detected genes" else "log10 UMIs",y="Mitochondrial UMIs (%)",color=NULL)
  if (!is.null(thresholds)) {
    xvals <- unlist(thresholds[if(x=="nGene") c("minGenesPerCell","maxGenesPerCell") else c("minUMIsPerCell","maxUMIsPerCell")])
    p <- p+ggplot2::geom_vline(xintercept=log10(xvals[xvals>0]),linetype=2)+
      ggplot2::geom_hline(yintercept=100*c(thresholds$minMitoPerCell,thresholds$maxMitoPerCell),linetype=2)
  }
  p
}

qc_plot_ribosomal <- function(metrics) {
  fields <- paste0("ribo_",.qc_ribosomal_fields)
  d <- .qc_long_metrics(metrics,fields,gsub("_"," ",.qc_ribosomal_fields))
  ggplot2::ggplot(d,ggplot2::aes(x=stage,y=value))+
    ggplot2::geom_boxplot(outlier.alpha=0.3,na.rm=TRUE)+
    ggplot2::scale_x_discrete(drop=FALSE)+
    ggplot2::facet_wrap(~metric,scales="free_y",ncol=3)+ggplot2::theme_bw()+
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=25,hjust=1))+
    ggplot2::labs(x=NULL,y="Annotated member UMIs / all UMIs (%)",
      caption="Annotated shares can be incomplete. Inspect unresolved UMI shares below. Protein parent/subgroups overlap; do not add them.")
}

qc_ribosomal_coverage <- function(metrics) {
  do.call(rbind,lapply(.qc_ribosomal_fields,function(group) {
    known <- metrics[[paste0("ribo_",group)]]; unknown <- metrics[[paste0("ribo_unknown_",group)]]
    data.frame(group=group,assessed_cells=sum(is.finite(known)),unavailable_cells=sum(!is.finite(known)),
      median_annotated_percent=if(any(is.finite(known))) stats::median(100*known,na.rm=TRUE) else NA_real_,
      median_unresolved_percent=if(any(is.finite(unknown))) stats::median(100*unknown,na.rm=TRUE) else NA_real_)
  }))
}

qc_plot_annotation_coverage <- function(metrics) {
  fields <- paste0("annotation_",c("matched","ambiguous","partial","unmapped"))
  d <- .qc_long_metrics(metrics,fields,sub("annotation_","",fields))
  ggplot2::ggplot(d,ggplot2::aes(x=stage,y=value))+
    ggplot2::geom_boxplot(outlier.alpha=0.25,na.rm=TRUE)+ggplot2::facet_wrap(~metric,ncol=4)+
    ggplot2::scale_x_discrete(drop=FALSE)+
    ggplot2::theme_bw()+ggplot2::labs(x=NULL,y="UMIs / all UMIs (%)",title="Annotation coverage weighted by observed counts")
}

qc_plot_rule_impact <- function(proposal) {
  d <- proposal$rule_counts
  ggplot2::ggplot(d,ggplot2::aes(x=stats::reorder(rule,failed),y=failed))+
    ggplot2::geom_col(fill="#426c9a")+ggplot2::coord_flip()+ggplot2::theme_bw()+
    ggplot2::labs(x=NULL,y="Cells failing rule",caption="Rules overlap. Their counts must not be added to obtain total exclusions.")
}

qc_plot_overlap <- function(proposal) {
  d <- utils::head(proposal$overlap,15L)
  ggplot2::ggplot(d,ggplot2::aes(x=stats::reorder(combination,cells),y=cells))+
    ggplot2::geom_col(fill="#c54b42")+ggplot2::coord_flip()+ggplot2::theme_bw()+
    ggplot2::labs(x=NULL,y="Cells with this exact rule combination",caption="Top 15 combinations; each excluded cell belongs to exactly one combination. Full reasons are in the cell table/download.")
}

qc_plot_sensitivity <- function(data) {
  ggplot2::ggplot(data,ggplot2::aes(x=threshold,y=100*retained/baseline_cells))+
    ggplot2::geom_line()+ggplot2::geom_vline(ggplot2::aes(xintercept=current),linetype=2,color="#c54b42")+
    ggplot2::facet_wrap(~rule,scales="free_x",ncol=3)+ggplot2::coord_cartesian(ylim=c(0,100))+
    ggplot2::theme_bw()+ggplot2::labs(x="Cutoff (fractions use 0–1)",y="Baseline cells retained (%)",
      caption="One cutoff varied at a time; all others fixed at the selected preview/applied values. Dashed lines mark those values.")
}
