# Join by stable cell identifier; never join scores by the displayed row number.
# Baseline fits are reused for retained subsets, with common scales per plot.
qc_diagnostic_comparison <- function(comparison,result) {
  if(!result$status %in% c("scored","estimated")) stop("This check has no current fitted evidence.")
  if(anyDuplicated(result$values$cell_id)) stop("Duplicate optional-QC cell identifiers.")
  i<-match(comparison$cell_id,result$values$cell_id)
  if(anyNA(i)) stop("Diagnostic results do not cover the current comparison cells.")
  cbind(comparison,result$values[i,setdiff(names(result$values),c("cell_id","cells")),drop=FALSE])
}

.qc_doublet_plot_groups <- function(d) {
  groups<-if("qc_diag_doublet_group" %in% names(d)) d$qc_diag_doublet_group else d$qc_diag_capture
  exploratory<-if("qc_diag_doublet_mode" %in% names(d)) d$qc_diag_doublet_mode=="exploratory" else rep(FALSE,nrow(d))
  d$qc_plot_group<-paste(ifelse(exploratory,"Provisional group:","Declared capture:"),groups)
  d
}

qc_plot_doublet_distribution <- function(d) ggplot2::ggplot(.qc_doublet_plot_groups(d),ggplot2::aes(x=qc_diag_doublet_score,fill=qc_diag_doublet_class))+
  ggplot2::geom_histogram(bins=30,position="identity",alpha=.6)+ggplot2::facet_grid(qc_plot_group~stage,drop=FALSE)+
  ggplot2::theme_bw()+ggplot2::labs(x="scDblFinder score",y="Cells",fill="Prediction",
    caption="Groups are scored separately. Exploratory groups are provisional, with unknown capture identity. Predictions do not remove cells.")

qc_plot_doublet_joint <- function(d,x="nUMI") {
  d<-.qc_doublet_plot_groups(d)
  d$diagnostic_x<-log10(d[[x]])
  ggplot2::ggplot(d,ggplot2::aes(x=diagnostic_x,y=qc_diag_doublet_score,color=qc_diag_doublet_class))+
    ggplot2::geom_point(alpha=.6,size=1)+ggplot2::facet_grid(qc_plot_group~stage,drop=FALSE)+
    ggplot2::theme_bw()+ggplot2::labs(x=paste("log10",if(x=="nUMI") "UMIs" else "detected genes"),y="scDblFinder score",color="Prediction")
}

qc_plot_doublet_counts <- function(d) {
  d<-.qc_doublet_plot_groups(d)
  counts<-as.data.frame(table(group=d$qc_plot_group,stage=d$stage,prediction=d$qc_diag_doublet_class))
  ggplot2::ggplot(counts,ggplot2::aes(x=group,y=Freq,fill=prediction))+
    ggplot2::geom_col()+ggplot2::facet_wrap(~stage,drop=FALSE)+ggplot2::theme_bw()+
    ggplot2::labs(x="Analysis group",y="Cells",fill="Prediction")
}

qc_plot_ambient_distribution <- function(d) ggplot2::ggplot(d,ggplot2::aes(x=100*qc_diag_ambient_fraction))+
  ggplot2::geom_histogram(bins=30,fill="#647fa3")+ggplot2::facet_grid(qc_diag_ambient_fit_group~stage,drop=FALSE)+
  ggplot2::theme_bw()+ggplot2::labs(x="Estimated contamination (%)",y="Cells",
    caption="Exploratory, population-dependent DecontX estimates. Original counts remain selected; no ambient correction is applied.")

qc_plot_ambient_joint <- function(d,x="nUMI") {
  d$diagnostic_x<-log10(d[[x]])
  ggplot2::ggplot(d,ggplot2::aes(x=diagnostic_x,y=100*qc_diag_ambient_fraction,color=qc_diag_ambient_population))+
    ggplot2::geom_point(alpha=.6,size=1)+ggplot2::facet_grid(qc_diag_ambient_fit_group~stage,drop=FALSE)+ggplot2::theme_bw()+
    ggplot2::labs(x=paste("log10",if(x=="nUMI") "UMIs" else "detected genes"),y="Estimated contamination (%)",color="Preliminary population")
}

qc_plot_ambient_populations <- function(d) ggplot2::ggplot(d,ggplot2::aes(x=qc_diag_ambient_population,y=100*qc_diag_ambient_fraction))+
  ggplot2::geom_boxplot(outlier.alpha=.3)+ggplot2::facet_grid(qc_diag_ambient_fit_group~stage,drop=FALSE)+ggplot2::theme_bw()+
  ggplot2::labs(x="Preliminary population (not final biological clusters)",y="Estimated contamination (%)")

qc_ambient_fit_table <- function(result) {
  fits<-result$diagnostics$fit_diagnostics
  if(!length(fits)) return(data.frame())
  do.call(rbind,lapply(seq_along(fits),function(i) {
    f<-fits[[i]]
    data.frame(fit=if(length(names(fits)) && nzchar(names(fits)[i])) names(fits)[i] else as.character(i),
      reported_iteration=f$reported_iteration,iteration_limit=f$iteration_limit,reached_iteration_limit=isTRUE(f$reached_iteration_limit))
  }))
}

qc_plot_ambient_fit <- function(result) {
  fits<-result$diagnostics$fit_diagnostics
  d<-do.call(rbind,lapply(seq_along(fits),function(i) {
    values<-as.numeric(fits[[i]]$log_likelihood)
    if(!length(values)) return(NULL)
    data.frame(fit=if(length(names(fits)) && nzchar(names(fits)[i])) names(fits)[i] else as.character(i),iteration=seq_along(values),log_likelihood=values)
  }))
  if(is.null(d) || !nrow(d)) stop("Likelihood history is unavailable in this DecontX result; inspect the iteration table and warnings.")
  ggplot2::ggplot(d,ggplot2::aes(x=iteration,y=log_likelihood))+ggplot2::geom_line()+
    ggplot2::facet_wrap(~fit,scales="free_y")+ggplot2::theme_bw()+
    ggplot2::labs(x="Iteration",y="Reported log likelihood",caption="A stopped fit is not independent validation of contamination estimates.")
}

qc_plot_cycle_scores <- function(d) ggplot2::ggplot(d,ggplot2::aes(x=qc_diag_S_score,y=qc_diag_G2M_score,color=qc_diag_phase))+
  ggplot2::geom_point(alpha=.6,size=1)+ggplot2::facet_wrap(~stage,drop=FALSE)+ggplot2::theme_bw()+
  ggplot2::labs(x="S program score",y="G2M program score",color="Predicted phase",
    caption="Scores and predicted phases are diagnostic. No regression or cell exclusions are applied.")

qc_plot_cycle_counts <- function(d) {
  counts<-as.data.frame(table(stage=d$stage,phase=d$qc_diag_phase))
  ggplot2::ggplot(counts,ggplot2::aes(x=phase,y=Freq,fill=phase))+ggplot2::geom_col()+
    ggplot2::facet_wrap(~stage,drop=FALSE)+ggplot2::theme_bw()+ggplot2::labs(x="Predicted phase",y="Cells",fill="Phase")
}

qc_plot_cycle_qc <- function(d) {
  specs<-c(nUMI="log10 UMIs",nGene="log10 genes",mitoRatio="Mitochondrial %",rbcRatio="Hemoglobin %")
  x<-do.call(rbind,lapply(names(specs),function(n) data.frame(phase=d$qc_diag_phase,stage=d$stage,metric=specs[[n]],
    value=if(n %in% c("nUMI","nGene")) log10(d[[n]]) else 100*d[[n]])))
  ggplot2::ggplot(x,ggplot2::aes(x=phase,y=value))+ggplot2::geom_boxplot(na.rm=TRUE,outlier.alpha=.25)+
    ggplot2::facet_grid(metric~stage,scales="free_y",drop=FALSE)+ggplot2::theme_bw()+
    ggplot2::labs(x="Predicted phase",y=NULL,caption="QC relationships can reflect biology. Phase alone is not an exclusion rule.")
}
