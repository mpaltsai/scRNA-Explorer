testthat::test_that("shared optional methods source from outside the repo and preserve duplicate-feature uncertainty", {
  path<-normalizePath(file.path("..","..","analysis","R","qc_optional_checks.R"))
  withr::local_dir(tempdir());e<-new.env(parent=baseenv());source(path,local=e)
  a<-list(table=data.frame(gene_id=c("ENSG1","ENSG2"),gene_name=c("Sgene","Ggene")),
    species="Homo sapiens",feature_id_type="gene_name",provenance=list(source="fixture"),
    feature_map=data.frame(feature=c("Sgene","Sgene_1","Ggene"),lookup_id=c("Sgene","Sgene","Ggene")))
  m<-e$.qc_cell_cycle_mapping(c("Sgene","Sgene_1","Ggene"),list(S="Sgene",G2M="Ggene"),a,"Homo sapiens")
  testthat::expect_identical(m$mapping$status,c("ambiguous_features","matched"))
  testthat::expect_length(m$matched$S,0)
  testthat::expect_identical(m$matched$G2M,"Ggene")
})

testthat::test_that("runtime restores RNG and library paths on both success and failure", {
  libraries<-.libPaths();set.seed(923);seed<-.Random.seed
  qc_optional_runtime("",1,function() {runif(4);invisible(NULL)})
  testthat::expect_identical(.Random.seed,seed);testthat::expect_identical(.libPaths(),libraries)
  testthat::expect_error(qc_optional_runtime(tempfile("absent-library"),1,function() NULL),"unavailable")
  testthat::expect_identical(.Random.seed,seed);testthat::expect_identical(.libPaths(),libraries)
})

testthat::test_that("diagnostic CSV uploads preserve literal capture IDs and reject invalid rates", {
  root<-tempfile("qc-tables-");dir.create(root);on.exit(unlink(root,recursive=TRUE))
  p<-file.path(root,"upload")
  writeLines(c("capture,expected_rate,uncertainty,source","001,0.1,0.02,simulation"),p)
  u<-data.frame(name="rates.csv",datapath=p)
  r<-qc_read_diagnostic_table(u)
  testthat::expect_identical(r$table$capture,"001")
  testthat::expect_equal(r$table$expected_rate,.1)
  writeLines(c("capture,expected_rate,uncertainty,source","001,unknown,0.02,simulation"),p)
  testthat::expect_error(qc_read_diagnostic_table(u),"must be numeric")
})

testthat::test_that("capture declarations preserve unknowns and reject shared incomplete inputs", {
  f<-qc_diag_fixture();on.exit(unlink(f$root,recursive=TRUE))
  cfg<-qc_diag_config("doublets",f$app_root)
  testthat::expect_error(qc_capture_assignment(f$baseline,cfg$capture),"unresolved")
  cfg$capture<-list(mode="single",id="capture1",complete=TRUE,evidence="Known simulation")
  testthat::expect_identical(qc_capture_assignment(f$baseline,cfg$capture),rep("capture1",8))
  testthat::expect_error(qc_capture_assignment(f$baseline,cfg$capture,list(list(declared_captures="capture1"))),"multiple")
  b<-f$baseline;b$capture<-rep("capture1",8)
  cfg$capture$mode<-"metadata";cfg$capture$column<-"capture"
  testthat::expect_error(qc_capture_assignment(b,cfg$capture,list(b)),"spans multiple")
  b$capture<-rep(NA_character_,8)
  testthat::expect_error(qc_capture_assignment(b,cfg$capture,list(b)),"unresolved")
})

testthat::test_that("real cell-cycle scoring reuses pinned annotation without changing baseline or original fields", {
  f<-qc_diag_fixture();on.exit(unlink(f$root,recursive=TRUE))
  b<-f$baseline;b$cell_cycle_S_score<-rep(99,ncol(b));before<-digest::digest(b);hash<-qc_baseline_signature(b)
  cfg<-qc_diag_config("cell_cycle",f$app_root)
  mapped<-qc_cycle_marker_review(b,cfg,f$app_root)
  testthat::expect_equal(mapped$coverage$matched,c(2,2))
  set.seed(811);seed<-.Random.seed
  r<-qc_run_optional(b,"cell_cycle",cfg,hash,f$app_root)
  testthat::expect_identical(r$status,"scored")
  testthat::expect_true(all(is.finite(r$values$qc_diag_S_score)))
  testthat::expect_identical(.Random.seed,seed)
  testthat::expect_identical(digest::digest(b),before)
  testthat::expect_identical(b$cell_cycle_S_score,rep(99,ncol(b)))
  testthat::expect_false("object" %in% names(r))
  cfg$parameters$reference$species<-"Mus musculus"
  wrong<-qc_run_optional(b,"cell_cycle",cfg,hash,f$app_root)
  testthat::expect_identical(wrong$status,"not_assessed")
  testthat::expect_true(all(is.na(wrong$values$qc_diag_S_score)))
  saveRDS(data.frame(changed=TRUE),f$reference$provenance$path)
  testthat::expect_error(qc_pinned_cycle_annotation(b),"missing or changed")
})

testthat::test_that("doublet and ambient evidence abstain before fitting when prerequisites are unknown", {
  f<-qc_diag_fixture();on.exit(unlink(f$root,recursive=TRUE))
  calls<-0L;fit<-function(...) {calls<<-calls+1L;qc_diag_fake_fit(...)}
  for(method in c("doublets","ambient")) {
    cfg<-qc_diag_config(method,f$app_root)
    r<-qc_run_optional(f$baseline,method,cfg,qc_baseline_signature(f$baseline),f$app_root,execute=fit)
    testthat::expect_identical(r$status,"not_assessed")
    testthat::expect_true(all(is.na(r$values[[.qc_optional_columns[[method]][1]]])))
  }
  testthat::expect_equal(calls,0)
  cfg<-qc_diag_config("ambient",f$app_root);cfg$ambient_method<-"SoupX"
  r<-qc_run_optional(f$baseline,"ambient",cfg,qc_baseline_signature(f$baseline),f$app_root)
  testthat::expect_identical(r$status,"not_assessed")
  testthat::expect_match(paste(r$reasons,collapse=" "),"background")
})

testthat::test_that("the adapter runs the existing doublet method on a known synthetic capture", {
  testthat::skip_if_not_installed("scDblFinder")
  f<-qc_diag_fixture();on.exit(unlink(f$root,recursive=TRUE))
  set.seed(14)
  raw<-scDblFinder::mockDoubletSCE(ncells=c(60,60),ngenes=300,dbl.rate=.1)
  rownames(raw)[1:2]<-c("MT-A","HBA")
  record<-new_count_import(SummarizedExperiment::assay(raw,"counts"),"simulation",object=raw)
  b<-make_baseline_sce(record,prepare_sample_qc(record,annotation_fixture(),"gene_name"))
  before<-digest::digest(b)
  cfg<-qc_diag_config("doublets",f$app_root);cfg$capture<-list(mode="single",id="simulation",complete=TRUE,evidence="Known synthetic capture")
  r<-qc_run_optional(b,"doublets",cfg,qc_baseline_signature(b),f$app_root)
  testthat::expect_identical(r$status,"scored")
  testthat::expect_true(all(is.finite(r$values$qc_diag_doublet_score)))
  testthat::expect_identical(r$diagnostics$rate_table$capture,"simulation")
  testthat::expect_equal(nrow(r$values),ncol(b))
  testthat::expect_identical(digest::digest(b),before)
})

testthat::test_that("compact optional evidence joins by cell ID and survives exports without replacing counts", {
  f<-qc_diag_fixture();on.exit(unlink(f$root,recursive=TRUE))
  b<-f$baseline;hash<-qc_baseline_signature(b);cfg<-qc_diag_config("ambient",f$app_root);cfg$exploratory_pooled<-TRUE
  r<-qc_run_optional(b,"ambient",cfg,hash,f$app_root,execute=qc_diag_fake_fit)
  testthat::expect_identical(r$status,"estimated")
  shuffled<-r;shuffled$values<-r$values[8:1,,drop=FALSE]
  m<-qc_review_metrics(b);p<-qc_filter_proposal(m,qc_default_thresholds(m),hash)
  d<-qc_diagnostic_comparison(qc_comparison_metrics(m,p),shuffled)
  testthat::expect_equal(d$qc_diag_ambient_fraction[1:8],r$values$qc_diag_ambient_fraction)
  subset<-b[,c(8,2,4)]
  attached<-qc_attach_optional(subset,list(ambient=shuffled),hash)
  testthat::expect_equal(attached$qc_diag_ambient_fraction,r$values$qc_diag_ambient_fraction[c(8,2,4)])
  testthat::expect_identical(SummarizedExperiment::assay(attached,"counts"),SummarizedExperiment::assay(subset,"counts"))
  testthat::expect_identical(SummarizedExperiment::assay(attached,"extra"),SummarizedExperiment::assay(subset,"extra"))
  testthat::expect_identical(attached$donor,subset$donor)
  pth<-file.path(f$root,"export.rds");saveRDS(attached,pth)
  testthat::expect_identical(readRDS(pth)$qc_diag_ambient_fraction,attached$qc_diag_ambient_fraction)
  b$qc_diag_ambient_fraction<-rep(99,8)
  testthat::expect_error(qc_attach_optional(b,list(ambient=r),hash),"none were overwritten")
  S4Vectors::metadata(b)$scrna_explorer_qc$optional_diagnostics<-list(fields="qc_diag_ambient_fraction")
  testthat::expect_error(qc_attach_optional(b,list(ambient=r),hash),"ownership/version")
})

testthat::test_that("a single DecontX population without background is unavailable evidence, not zero contamination", {
  f<-qc_diag_fixture();on.exit(unlink(f$root,recursive=TRUE))
  cfg<-qc_diag_config("ambient",f$app_root);cfg$exploratory_pooled<-TRUE
  fit<-function(...) {r<-qc_diag_fake_fit(...);r$object$ambient_decontX_preliminary_cluster<-rep("only population",ncol(r$object));r}
  r<-qc_run_optional(f$baseline,"ambient",cfg,qc_baseline_signature(f$baseline),f$app_root,execute=fit)
  testthat::expect_identical(r$status,"not_assessed")
  testthat::expect_true(all(is.na(r$values$qc_diag_ambient_fraction)))
  testthat::expect_match(r$reasons,"Insufficient preliminary population")
})

testthat::test_that("all diagnostic plot types build and share baseline/retained axes", {
  f<-qc_diag_fixture();on.exit(unlink(f$root,recursive=TRUE))
  b<-f$baseline;hash<-qc_baseline_signature(b);m<-qc_review_metrics(b)
  p<-qc_filter_proposal(m,qc_default_thresholds(m),hash);p$keep[1:2]<-FALSE
  comparison<-qc_comparison_metrics(m,p)
  for(method in c("doublets","ambient","cell_cycle")) {
    cfg<-qc_diag_config(method,f$app_root)
    if(method=="doublets") cfg$capture<-list(mode="single",id="simulation",complete=TRUE,evidence="Known simulation")
    if(method=="ambient") cfg$exploratory_pooled<-TRUE
    r<-qc_run_optional(b,method,cfg,hash,f$app_root,execute=qc_diag_fake_fit)
    testthat::expect_true(r$status %in% c("scored","estimated"))
    d<-qc_diagnostic_comparison(comparison,r)
    plots<-switch(method,doublets=list(qc_plot_doublet_distribution(d),qc_plot_doublet_joint(d),qc_plot_doublet_counts(d)),
      ambient=list(qc_plot_ambient_distribution(d),qc_plot_ambient_joint(d),qc_plot_ambient_populations(d),qc_plot_ambient_fit(r)),
      cell_cycle=list(qc_plot_cycle_scores(d),qc_plot_cycle_counts(d),qc_plot_cycle_qc(d)))
    for(plot in plots) testthat::expect_s3_class(ggplot2::ggplot_build(plot),"ggplot_built")
    if(method=="doublets") {
      built<-ggplot2::ggplot_build(qc_plot_doublet_joint(d))
      testthat::expect_length(built$layout$panel_scales_x,1)
      testthat::expect_length(built$layout$panel_scales_y,1)
    }
  }
})
