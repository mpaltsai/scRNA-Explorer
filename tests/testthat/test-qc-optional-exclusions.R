qc_exclusion_evidence <- function(f) {
  b<-f$baseline;hash<-qc_baseline_signature(b)
  d<-qc_diag_config("doublets",f$app_root)
  d$analysis_mode<-"exploratory";d$group_mode<-"input";d$rates_mode<-"unknown"
  a<-qc_diag_config("ambient",f$app_root);a$exploratory_pooled<-TRUE
  r<-list(doublets=qc_run_optional(b,"doublets",d,hash,f$app_root,execute=qc_diag_fake_fit),
    ambient=qc_run_optional(b,"ambient",a,hash,f$app_root,execute=qc_diag_fake_fit))
  r$doublets$values$qc_diag_doublet_score<-c(.1,.6,.4,.8,.1,.7,.2,.9)
  r$doublets$values$qc_diag_doublet_class<-ifelse(r$doublets$values$qc_diag_doublet_score>.5,"doublet","singlet")
  r$ambient$values$qc_diag_ambient_fraction<-c(.05,.1,.5,.4,.05,.1,.5,.4)
  r
}

testthat::test_that("optional exclusions combine by cell ID with overlap reasons and exact boundaries", {
  f<-qc_diag_fixture();on.exit(unlink(f$root,recursive=TRUE))
  b<-f$baseline;before<-digest::digest(b);hash<-qc_baseline_signature(b)
  m<-qc_review_metrics(b);t<-qc_default_thresholds(m);r<-qc_exclusion_evidence(f)
  for(method in names(r)) r[[method]]$values<-r[[method]]$values[8:1,,drop=FALSE]
  rules<-qc_optional_filter_defaults();rules$doublets<-"prediction";rules$ambient<-TRUE;rules$ambient_fraction<-.3
  p<-qc_filter_proposal(m,t,hash,rules,r)
  testthat::expect_identical(p$keep,c(TRUE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE))
  testthat::expect_equal(p$summary$excluded,6)
  testthat::expect_equal(sum(p$rule_counts$failed),8)
  testthat::expect_match(p$per_cell$reasons[4],"Predicted doublet; High estimated ambient fraction",fixed=TRUE)
  testthat::expect_identical(p$optional_evidence$doublets$analysis_context$capture_status,"unknown")
  applied<-qc_apply_filter(b,p,"Use reviewed exploratory findings for this synthetic test",optional_results=r)
  testthat::expect_identical(colnames(applied$object),colnames(b)[c(1,5)])
  testthat::expect_identical(SummarizedExperiment::assay(applied$object,"counts"),SummarizedExperiment::assay(b,"counts")[,c(1,5)])
  testthat::expect_identical(applied$decision$optional_evidence_used,p$optional_evidence)
  testthat::expect_identical(applied$decision$per_cell$qc_diag_doublet_class_used,
    r$doublets$values$qc_diag_doublet_class[match(m$cell_id,r$doublets$values$cell_id)])
  s<-qc_threshold_sensitivity(m,t,optional_keep=p$optional_keep)
  testthat::expect_true(all(s$retained[s$threshold==s$current]==sum(p$keep)))
  rules$doublets<-"score";rules$doublet_score<-.6
  scored<-qc_filter_proposal(m,t,hash,rules,r)
  testthat::expect_true(scored$keep[2]) # Exactly at the maximum stays.
  testthat::expect_equal(sum(scored$keep),3)
  rules$doublets<-"none";rules$ambient_fraction<-.4
  ambient<-qc_filter_proposal(m,t,hash,rules,r)
  testthat::expect_identical(ambient$keep,c(TRUE,TRUE,FALSE,TRUE,TRUE,TRUE,FALSE,TRUE))
  testthat::expect_identical(digest::digest(b),before)
})

testthat::test_that("enabled exclusions reject unavailable, stale, partial and tampered evidence", {
  f<-qc_diag_fixture();on.exit(unlink(f$root,recursive=TRUE))
  b<-f$baseline;hash<-qc_baseline_signature(b);m<-qc_review_metrics(b);t<-qc_default_thresholds(m)
  r<-qc_exclusion_evidence(f);rules<-qc_optional_filter_defaults();rules$doublets<-"score";rules$doublet_score<-.6
  testthat::expect_true(all(qc_filter_proposal(m,t,hash)$keep))
  testthat::expect_error(qc_filter_proposal(m,t,hash,rules),"current doublets")
  for(change in c("stale_settings","not_assessed","failed")) {
    bad<-r;bad$doublets$status<-change
    testthat::expect_error(qc_filter_proposal(m,t,hash,rules,bad),"current doublets")
  }
  bad<-r;bad$doublets$baseline_hash<-"other"
  testthat::expect_error(qc_filter_proposal(m,t,hash,rules,bad),"this baseline")
  bad<-r;bad$doublets$values$cell_id[1]<-bad$doublets$values$cell_id[2]
  testthat::expect_error(qc_filter_proposal(m,t,hash,rules,bad),"cell identifiers")
  bad<-r;bad$doublets$values$qc_diag_doublet_score[1]<-NA_real_
  testthat::expect_error(qc_filter_proposal(m,t,hash,rules,bad),"missing/invalid")
  badrules<-rules;badrules$doublet_score<-NA_real_
  testthat::expect_error(qc_filter_proposal(m,t,hash,badrules,r),"no cutoff is chosen")
  p<-qc_filter_proposal(m,t,hash,rules,r)
  bad<-r;bad$doublets$values$qc_diag_doublet_score[1]<-.95 # Same run ID, changed values.
  testthat::expect_error(qc_apply_filter(b,p,"Changed evidence",optional_results=bad),"differs")
  bad<-r;bad$doublets$revision<-"rerun"
  testthat::expect_error(qc_apply_filter(b,p,"New run",optional_results=bad),"differs")
  p$keep[1]<-!p$keep[1]
  testthat::expect_error(qc_apply_filter(b,p,"Tampered mask",optional_results=r),"differs")
  testthat::expect_true(all(qc_filter_proposal(m,t,hash,optional_results=bad)$keep))
})

testthat::test_that("optional filtering stays per sample and preserves the applied result when evidence changes", {
  f<-qc_diag_fixture();on.exit(unlink(f$root,recursive=TRUE))
  b2<-f$baseline;b2$qc_cell_id<-paste0("B::",colnames(b2));b2$qc_sample<-rep("B",ncol(b2))
  evidence<-qc_exclusion_evidence(f)
  server<-function(input,output,session) {
    r<-shiny::reactiveVal(evidence)
    controller<-qc_review_server(input,output,session,function() list(A=f$baseline,B=b2),
      optional_results=function(n) if(n=="A") r() else list())
  }
  shiny::testServer(server,{
    action<-function(n) structure(as.integer(n),class="shinyActionButtonValue")
    session$setInputs(selectSample="A",qcPreview=action(0),QCFiltering=action(0),qcView="preview")
    qc_test_set_draft(session,controller,qc_default_thresholds(controller$metrics()))
    ids<-controller$control_ids()
    do.call(session$setInputs,stats::setNames(list("prediction",TRUE,.3),ids[c("doubletRule","ambientExclude","ambientMaxFraction")]))
    session$setInputs(qcPreview=action(1))
    testthat::expect_equal(sum(controller$preview()$keep),2)
    testthat::expect_null(controller$filtered()$A)
    session$setInputs(QCFiltering=action(1),qcApplyReason="Explicitly use the reviewed model findings")
    testthat::expect_equal(ncol(controller$filtered()$A),2)
    testthat::expect_true(all(is.na(controller$filtered()$A$qc_diag_capture)))
    applied<-digest::digest(controller$filtered()$A)
    do.call(session$setInputs,stats::setNames(list(NULL),ids[["ambientMaxFraction"]]))
    testthat::expect_null(controller$preview())
    testthat::expect_match(output$qcOptionalFilterNotice,"Enter a maximum retained ambient fraction")
    testthat::expect_identical(digest::digest(controller$filtered()$A),applied)
    do.call(session$setInputs,stats::setNames(list(.3),ids[["ambientMaxFraction"]]))
    changed<-r();changed$ambient$revision<-"new-ambient-run";r(changed);session$flushReact()
    testthat::expect_null(controller$preview())
    session$setInputs(QCFiltering=action(2))
    testthat::expect_match(output$QCDataFilt,"Preview the current thresholds")
    testthat::expect_identical(digest::digest(controller$filtered()$A),applied)
    session$setInputs(qcView="applied")
    testthat::expect_equal(sum(controller$selection()$keep),2)
    testthat::expect_true(all(controller$sensitivity()$retained[controller$sensitivity()$threshold==controller$sensitivity()$current]==2))
    session$setInputs(selectSample="B",qcView="preview")
    testthat::expect_identical(controller$optional_draft(),qc_optional_filter_defaults())
    qc_test_set_draft(session,controller,qc_default_thresholds(controller$metrics()))
    session$setInputs(qcPreview=action(2))
    testthat::expect_equal(sum(controller$preview()$keep),8)
    session$setInputs(selectSample="A")
    testthat::expect_identical(controller$optional_draft()$doublets,"prediction")
    changed<-r();changed$ambient$status<-"stale_settings";r(changed);session$flushReact()
    session$setInputs(qcPreview=action(3))
    testthat::expect_match(output$QCDataFilt,"current ambient")
    testthat::expect_null(controller$preview())
    testthat::expect_identical(digest::digest(controller$filtered()$A),applied)
    ids<-controller$control_ids()
    do.call(session$setInputs,stats::setNames(list("none",FALSE),ids[c("doubletRule","ambientExclude")]))
    session$setInputs(qcPreview=action(4))
    testthat::expect_equal(sum(controller$preview()$keep),8)
  })
})

testthat::test_that("ambient UMI-gene cluster colours remain aligned across baseline and retained panels", {
  f<-qc_diag_fixture();on.exit(unlink(f$root,recursive=TRUE))
  b<-f$baseline;hash<-qc_baseline_signature(b);m<-qc_review_metrics(b);r<-qc_exclusion_evidence(f)$ambient
  r$values<-r$values[8:1,,drop=FALSE]
  t<-qc_default_thresholds(m);t$maxUMIsPerCell<-floor(stats::median(m$nUMI))
  d<-qc_diagnostic_comparison(qc_comparison_metrics(m,qc_filter_proposal(m,t,hash)),r)
  plot<-qc_plot_ambient_umi_genes_clusters(d);built<-ggplot2::ggplot_build(plot)
  testthat::expect_identical(as.character(plot$data$ambient_plot_cluster),r$values$qc_diag_ambient_population[match(d$cell_id,r$values$cell_id)])
  testthat::expect_equal(as.numeric(built$data[[1]]$x),log10(d$nUMI))
  testthat::expect_equal(as.numeric(built$data[[1]]$y),log10(d$nGene))
  testthat::expect_length(built$layout$panel_scales_x,1L)
  testthat::expect_length(built$layout$panel_scales_y,1L)
  testthat::expect_identical(built$data[[1]]$colour,built$plot$scales$get_scales("colour")$map(plot$data$ambient_plot_cluster))
})
