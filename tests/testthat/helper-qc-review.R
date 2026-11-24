source(file.path("..","..","site_data","R","qc_review.R"))
source(file.path("..","..","site_data","R","qc_review_plots.R"))
source(file.path("..","..","site_data","R","qc_review_server.R"))
source(file.path("..","..","site_data","R","createQCPlots.R"))

qc_review_fixture <- function(sample="A") {
  counts <- Matrix::Matrix(matrix(c(10,10,80, 30,10,160, 150,0,150, 20,40,340),nrow=3,
    dimnames=list(c("MT-A","HBA","GeneA"),c("at100","at200","highMito","highHb"))),sparse=TRUE)
  object <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=counts,other=counts))
  object$donor <- rep("donor1",4);object$nUMI<-rep(999,4)
  ref <- annotation_fixture()
  ref$table$entrezid <- as.list(as.character(1:5))
  ref$ribosomal <- list(status="available",species="Homo sapiens",groups=list(
    ribosomal_protein="3",cytosolic_ribosomal_protein="3",mitochondrial_ribosomal_protein=character()))
  record <- new_count_import(counts,sample,cell_data=SummarizedExperiment::colData(object),object=object)
  make_baseline_sce(record,prepare_sample_qc(record,ref,"gene_name"))
}

qc_test_set_draft <- function(session,controller,thresholds) {
  ids <- shiny::isolate(controller$control_ids())[names(thresholds)]
  do.call(session$setInputs,stats::setNames(thresholds,unname(ids)))
}
