source(file.path("..","..","site_data","R","qc_diagnostics.R"))
source(file.path("..","..","site_data","R","qc_diagnostic_controls.R"))
source(file.path("..","..","site_data","R","qc_diagnostic_plots.R"))
source(file.path("..","..","site_data","R","qc_diagnostic_server.R"))

qc_diag_fixture <- function() {
  root<-tempfile("optional-qc-");dir.create(root)
  genes<-c("MT-A","HBA","PCNA","MCM5","TOP2A","MKI67",paste0("Background",1:16))
  table<-data.frame(gene_id=sprintf("ENSG%06d",seq_along(genes)),gene_name=genes,gene_biotype="protein_coding",
    seq_name=c("MT",rep("1",length(genes)-1)),description=c("mitochondrial","hemoglobin test",rep("ordinary",length(genes)-2)))
  path<-file.path(root,"reference.rds");saveRDS(table,path)
  reference<-list(table=table,species="Homo sapiens",provenance=list(source="Synthetic pinned annotation",
    path=path,md5=unname(tools::md5sum(path))))
  counts<-Matrix::Matrix(matrix(10*rep(1:11,length.out=length(genes)*8),nrow=length(genes),
    dimnames=list(genes,paste0("cell",1:8))),sparse=TRUE)
  object<-SingleCellExperiment::SingleCellExperiment(assays=list(counts=counts,extra=counts))
  object$donor<-rep("fixture donor",8)
  r<-new_count_import(counts,"A",cell_data=SummarizedExperiment::colData(object),object=object)
  baseline<-make_baseline_sce(r,prepare_sample_qc(r,reference,"gene_name"))
  list(root=root,baseline=baseline,reference=reference,
    app_root=normalizePath(file.path("..","..","site_data")))
}

qc_diag_config <- function(method,app_root) {
  cfg<-list(seed=1,library_path=qc_optional_library(app_root),parameters=list(),
    capture=list(mode="unresolved",id="",column="",complete=FALSE,evidence=""),technology="simulation")
  if(method=="cell_cycle") {
    cfg$minimum_markers<-1
    cfg$parameters<-list(reference=list(name="Synthetic programs",species="Homo sapiens",
      s_genes=c("PCNA","MCM5"),g2m_genes=c("TOP2A","MKI67")),nbin=2,ctrl=2,scale_factor=10000)
  } else if(method=="doublets") {
    cfg$rates_mode<-"shared";cfg$expected_rate<-.1;cfg$uncertainty<-.02;cfg$rate_source<-"Simulation parameter"
  } else {
    cfg$ambient_method<-"DecontX";cfg$exploratory_pooled<-FALSE;cfg$clusters_column<-""
    cfg$parameters$model_args<-list(maxIter=20,delta=c(10,10),estimateDelta=TRUE,convergence=.001,varGenes=20,dbscanEps=1)
  }
  cfg
}

qc_diag_fake_fit <- function(object,parameters) {
  if(!is.null(parameters$annotation)) {
    object$cell_cycle_S_score<-seq_len(ncol(object))/10
    object$cell_cycle_G2M_score<-rev(seq_len(ncol(object)))/10
    object$cell_cycle_phase<-rep(c("S","G2M"),length.out=ncol(object))
    diagnostics<-list(status="scored",regression_applied=FALSE,excluded_cells=0L)
  } else if(!is.null(parameters$rate_table) || !is.null(parameters$group_column) || !is.null(parameters$capture_column)) {
    object$scDblFinder.score<-seq(.05,.95,length.out=ncol(object))
    object$scDblFinder.class<-ifelse(object$scDblFinder.score>.5,"doublet","singlet")
    diagnostics<-list(status="scored",excluded_cells=0L)
  } else {
    object$ambient_decontX_fraction<-seq(.01,.1,length.out=ncol(object))
    object$ambient_decontX_preliminary_cluster<-rep(c("preliminary1","preliminary2"),length.out=ncol(object))
    diagnostics<-list(status="estimated",excluded_cells=0L,counts_corrected=FALSE,
      fit_diagnostics=list(pooled=list(reported_iteration=10,iteration_limit=20,reached_iteration_limit=FALSE,log_likelihood=1:10)))
  }
  list(object=object,diagnostics=diagnostics)
}

qc_diag_set_controls <- function(session,controller,method,overrides=list(),app_root=normalizePath(file.path("..","..","site_data"))) {
  p<-shiny::isolate(controller$prefix())
  values<-c(list(seed=1,use_library=TRUE,library_path=qc_optional_library(app_root),
    capture_mode="unresolved",capture_id="",capture_column="",capture_evidence="",capture_complete=FALSE,technology="simulation"),
    switch(method,doublets=list(doublet_mode="capture",provisional_group_mode="input",provisional_group_column="",rates_mode="shared",expected_rate=.1,uncertainty=.02,rate_source="Simulation parameter"),
      ambient=list(ambient_method="DecontX",pooled=FALSE,clusters_column="",maxIter=20,delta1=10,delta2=10,estimateDelta=TRUE,convergence=.001,varGenes=20,dbscanEps=1,droplet_compatible=FALSE),
      cell_cycle=list(marker_mode="seurat",marker_name="",marker_species="",nbin=2,ctrl_auto=FALSE,ctrl=2,scale_factor=10000,minimum_markers=1)))
  values[names(overrides)]<-overrides
  values$ready_common<-p;values[[paste0("ready_",method)]]<-p
  do.call(session$setInputs,stats::setNames(values,paste0(p,names(values))))
}
