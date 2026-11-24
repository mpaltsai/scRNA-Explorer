testthat::test_that("review uses canonical metrics and preserves partial ribosomal/annotation evidence", {
  b<-qc_review_fixture();counts<-SummarizedExperiment::assay(b,"counts")
  before<-digest::digest(b);m<-qc_review_metrics(b)
  testthat::expect_equal(m$nUMI,c(100,200,300,400))
  testthat::expect_equal(m$ribo_cytosolic_ribosomal_protein,c(.8,.8,.5,.85))
  testthat::expect_equal(m$ribo_mitochondrial_ribosomal_protein,rep(0,4))
  testthat::expect_equal(unname(rowSums(m[paste0("annotation_",c("matched","ambiguous","partial","unmapped"))])),rep(1,4))
  testthat::expect_identical(digest::digest(b),before)
  SummarizedExperiment::rowData(b)$qc_cytosolic_ribosomal_protein[1]<-NA
  SummarizedExperiment::rowData(b)$qc_annotation_status[3]<-"unmapped"
  q<-qc_review_metrics(b)
  testthat::expect_equal(q$ribo_cytosolic_ribosomal_protein,m$ribo_cytosolic_ribosomal_protein)
  testthat::expect_equal(q$ribo_unknown_cytosolic_ribosomal_protein,c(.1,.15,.5,.05))
  testthat::expect_equal(q$annotation_unmapped,c(.8,.8,.5,.85))
  SummarizedExperiment::rowData(b)$qc_ribosomal_protein<-rep(NA,3)
  testthat::expect_true(all(is.na(qc_review_metrics(b)$ribo_ribosomal_protein)))
  testthat::expect_identical(SummarizedExperiment::assay(b,"counts"),counts)
})

testthat::test_that("inclusive previews match the unchanged legacy filter and count overlaps once", {
  b<-qc_review_fixture();m<-qc_review_metrics(b);t<-qc_default_thresholds(m)
  t$maxUMIsPerCell<-300;t$maxMitoPerCell<-.15;t$maxRbcPerCell<-.1
  p<-qc_filter_proposal(m,t,qc_baseline_signature(b))
  testthat::expect_identical(p$keep,c(TRUE,TRUE,FALSE,FALSE))
  testthat::expect_identical(p$per_cell$cell_id,as.character(b$qc_cell_id))
  e<-new.env(parent=baseenv());e$`%>%`<-dplyr::`%>%`;e$pull<-dplyr::pull
  sys.source(file.path("..","..","site_data","R","filter_cells.R"),envir=e)
  legacy<-do.call(e$filter_cells,c(list(raw_data=b,metrics=m),t))$se
  result<-qc_apply_filter(b,p,"Keep the inclusive boundary cells")
  testthat::expect_identical(SummarizedExperiment::assay(result$object,"counts"),SummarizedExperiment::assay(legacy,"counts"))
  testthat::expect_identical(result$object$donor,rep("donor1",2))
  testthat::expect_identical(result$object$nUMI,rep(999,2))
  testthat::expect_identical(SummarizedExperiment::assay(result$object,"other"),SummarizedExperiment::assay(b,"other")[,1:2])
  testthat::expect_identical(SummarizedExperiment::rowData(result$object),SummarizedExperiment::rowData(b))
  t$maxRbcPerCell<-.09
  overlap<-qc_filter_proposal(m,t,p$baseline_hash)
  testthat::expect_equal(overlap$summary$excluded,3)
  testthat::expect_equal(sum(overlap$rule_counts$failed),4)
  testthat::expect_equal(sum(overlap$overlap$cells),3)
  testthat::expect_equal(overlap$per_cell$failed_rules, c(1,0,1,2))
  single<-m[1,,drop=FALSE]
  testthat::expect_identical(qc_filter_proposal(single,qc_default_thresholds(single),"single")$keep,TRUE)
})

testthat::test_that("invalid, stale or tampered decisions cannot alter an SCE", {
  b<-qc_review_fixture();m<-qc_review_metrics(b);t<-qc_default_thresholds(m);hash<-qc_baseline_signature(b)
  p<-qc_filter_proposal(m,t,hash);before<-digest::digest(b)
  testthat::expect_identical(names(p$overlap),c("combination","cells"))
  testthat::expect_s3_class(ggplot2::ggplot_build(qc_plot_overlap(p)),"ggplot_built")
  testthat::expect_error(qc_apply_filter(b,p," "),"Provide a reason")
  bad<-t;bad$minUMIsPerCell<-99
  testthat::expect_error(qc_filter_proposal(m,bad,hash),"cannot be lowered")
  bad<-t;bad$maxUMIsPerCell<-99
  testthat::expect_error(qc_filter_proposal(m,bad,hash),"minimum")
  bad<-t;bad$minUMIsPerCell<-100.5
  testthat::expect_error(qc_filter_proposal(m,bad,hash),"integers")
  testthat::expect_error(qc_filter_proposal(m,c(t,t[1]),hash),"nine named")
  missing<-m;missing$mitoRatio<-rep(NA_real_,4)
  testthat::expect_error(qc_filter_proposal(missing,t,hash),"mitoRatio")
  changed<-b;SummarizedExperiment::assay(changed,"counts")[1,1]<-11
  testthat::expect_error(qc_apply_filter(changed,p,"Changed input"),"not for this baseline")
  p$keep[1]<-FALSE
  testthat::expect_error(qc_apply_filter(b,p,"Tampered mask"),"differs")
  testthat::expect_identical(digest::digest(b),before)
})

testthat::test_that("decisions retain lineage without expression copies and empty filters block handoff", {
  b<-qc_review_fixture();m<-qc_review_metrics(b);t<-qc_default_thresholds(m);hash<-qc_baseline_signature(b)
  first<-qc_apply_filter(b,qc_filter_proposal(m,t,hash),"Baseline-only acceptance")
  t$maxUMIsPerCell<-100
  second<-qc_apply_filter(b,qc_filter_proposal(m,t,hash),"Inspect one boundary cell",first$history)
  testthat::expect_identical(second$decision$previous_revision,first$decision$revision)
  testthat::expect_length(second$history,2)
  testthat::expect_false(any(vapply(second$history,function(d) any(vapply(d,function(x) inherits(x,"SingleCellExperiment"),logical(1))),logical(1))))
  path<-tempfile(fileext=".rds");on.exit(unlink(path));saveRDS(second$object,path)
  restored<-S4Vectors::metadata(readRDS(path))$scrna_explorer_qc$filtering
  testthat::expect_identical(restored$history,second$history)
  testthat::expect_identical(colnames(first$object),colnames(b))
  testthat::expect_error(qc_require_applied(list(A=b,B=b),list(A=second$object,B=NULL)),"Pending/empty: B")
  t$minGenesPerCell<-4;t$maxGenesPerCell<-4
  empty<-qc_apply_filter(b,qc_filter_proposal(m,t,hash),"Inspect empty-selection handling")
  testthat::expect_equal(ncol(empty$object),0)
  testthat::expect_error(qc_require_applied(list(A=b),list(A=empty$object)),"Pending/empty: A")
})

testthat::test_that("sensitivity and plot data agree with the preview, including an empty retained panel", {
  b<-qc_review_fixture();m<-qc_review_metrics(b);t<-qc_default_thresholds(m)
  t$maxUMIsPerCell<-300
  p<-qc_filter_proposal(m,t,qc_baseline_signature(b))
  s<-qc_threshold_sensitivity(m,t)
  testthat::expect_true(all(s$retained[s$threshold==s$current]==sum(p$keep)))
  d<-qc_comparison_metrics(m,p)
  testthat::expect_equal(nrow(d),nrow(m)+sum(p$keep))
  testthat::expect_identical(d$cell_id[as.character(d$stage)=="Preview retained"],p$per_cell$cell_id[p$keep])
  joint<-qc_joint_plot_data(d)
  build<-ggplot2::ggplot_build(qc_plot_joint(joint,t))
  testthat::expect_length(build$layout$panel_scales_x,1)
  testthat::expect_length(build$layout$panel_scales_y,1)
  testthat::expect_equal(nrow(build$layout$layout),2)
  testthat::expect_equal(joint$qc_y,100*d$mitoRatio)
  for(plot in list(qc_plot_overview(m),qc_plot_overview(d,TRUE),qc_plot_ribosomal(d),
    qc_plot_annotation_coverage(d),qc_plot_rule_impact(p),qc_plot_overlap(p),qc_plot_sensitivity(s)))
    testthat::expect_s3_class(ggplot2::ggplot_build(plot),"ggplot_built")
  p$keep[]<-FALSE
  empty<-qc_comparison_metrics(m,p)
  testthat::expect_equal(nrow(ggplot2::ggplot_build(qc_plot_joint(qc_joint_plot_data(empty),t))$layout$layout),2)
})
