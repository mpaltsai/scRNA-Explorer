qc_exploratory_config <- function(app_root) {
  x<-qc_diag_config("doublets",app_root)
  x$analysis_mode<-"exploratory";x$group_mode<-"input";x$group_column<-"";x$rates_mode<-"unknown"
  x$expected_rate<-NULL;x$uncertainty<-NULL;x$rate_source<-NULL
  x
}

testthat::test_that("exploratory groups do not fabricate capture identity or an expected rate", {
  f<-qc_diag_fixture();on.exit(unlink(f$root,recursive=TRUE))
  b<-f$baseline;before<-digest::digest(b);cfg<-qc_exploratory_config(f$app_root)
  seen<-NULL
  r<-qc_run_optional(b,"doublets",cfg,qc_baseline_signature(b),f$app_root,
    execute=function(object,parameters) {seen<<-parameters;qc_diag_fake_fit(object,parameters)})
  testthat::expect_identical(r$status,"scored")
  testthat::expect_identical(seen$analysis_mode,"exploratory")
  testthat::expect_identical(seen$rate_mode,"unknown")
  testthat::expect_null(seen$rate_table)
  testthat::expect_null(seen$capture_column)
  testthat::expect_true(all(is.na(r$values$qc_diag_capture)))
  testthat::expect_true(all(r$values$qc_diag_capture_status=="unknown"))
  testthat::expect_identical(unique(r$values$qc_diag_doublet_group),"A")
  testthat::expect_identical(r$analysis_context$group_basis,"Selected input matrix: A")
  attached<-qc_attach_optional(b,list(doublets=r),qc_baseline_signature(b))
  testthat::expect_true(all(is.na(attached$qc_diag_capture)))
  testthat::expect_identical(SummarizedExperiment::assay(attached,"counts"),SummarizedExperiment::assay(b,"counts"))
  testthat::expect_identical(digest::digest(b),before)
  comparison<-qc_comparison_metrics(qc_review_metrics(b))
  d<-qc_diagnostic_comparison(comparison,r)
  testthat::expect_true(all(grepl("Provisional group:",.qc_doublet_plot_groups(d)$qc_plot_group,fixed=TRUE)))
  for(plot in list(qc_plot_doublet_distribution(d),qc_plot_doublet_joint(d),qc_plot_doublet_counts(d)))
    testthat::expect_s3_class(ggplot2::ggplot_build(plot),"ggplot_built")
})

testthat::test_that("metadata grouping needs complete labels but not capture declarations", {
  f<-qc_diag_fixture();on.exit(unlink(f$root,recursive=TRUE))
  b<-f$baseline;b$geo_accession<-rep(c("GSM1","GSM2"),each=4)
  cfg<-qc_exploratory_config(f$app_root);cfg$group_mode<-"metadata";cfg$group_column<-"geo_accession"
  r<-qc_run_optional(b,"doublets",cfg,qc_baseline_signature(b),f$app_root,execute=qc_diag_fake_fit)
  testthat::expect_identical(r$status,"scored")
  testthat::expect_identical(r$values$qc_diag_doublet_group,as.character(b$geo_accession))
  testthat::expect_true(all(is.na(r$values$qc_diag_capture)))
  testthat::expect_match(r$analysis_context$group_basis,"geo_accession")
  b$geo_accession[1]<-NA_character_
  failed<-qc_run_optional(b,"doublets",cfg,qc_baseline_signature(b),f$app_root,
    execute=function(...) stop("Missing labels must not reach the model"))
  testthat::expect_identical(failed$status,"not_assessed")
  testthat::expect_match(paste(failed$reasons,collapse=" "),"missing labels")
  testthat::expect_true(all(is.na(failed$values$qc_diag_doublet_score)))
})

testthat::test_that("explicit group rates preserve old capture tables and reject incomplete assumptions", {
  f<-qc_diag_fixture();on.exit(unlink(f$root,recursive=TRUE))
  cfg<-qc_exploratory_config(f$app_root);cfg$rates_mode<-"shared"
  r<-qc_run_optional(f$baseline,"doublets",cfg,qc_baseline_signature(f$baseline),f$app_root,execute=qc_diag_fake_fit)
  testthat::expect_identical(r$status,"not_assessed")
  testthat::expect_match(paste(r$reasons,collapse=" "),"rate table")
  cfg$rates_mode<-"csv";cfg$rate_table<-data.frame(capture="A",expected_rate=.05,uncertainty=.02,source="Exploratory sensitivity assumption")
  r<-qc_run_optional(f$baseline,"doublets",cfg,qc_baseline_signature(f$baseline),f$app_root,execute=qc_diag_fake_fit)
  testthat::expect_identical(r$status,"scored")
  cfg$analysis_mode<-"capture";cfg$rates_mode<-"unknown"
  strict<-qc_run_optional(f$baseline,"doublets",cfg,qc_baseline_signature(f$baseline),f$app_root,execute=qc_diag_fake_fit)
  testthat::expect_identical(strict$status,"not_assessed")
  testthat::expect_match(paste(strict$reasons,collapse=" "),"Capture identity is unresolved")
})

testthat::test_that("unknown-rate controls ignore unused numbers and capture declarations", {
  p<-"form_"
  values<-list(ready_common=p,ready_doublets=p,seed=1,use_library=FALSE,doublet_mode="exploratory",
    provisional_group_mode="input",provisional_group_column="",rates_mode="unknown",expected_rate=.05,
    uncertainty=.02,rate_source="",capture_mode="single",capture_id="unverified",capture_complete=TRUE,capture_evidence="")
  cfg<-qc_diagnostic_config(stats::setNames(values,paste0(p,names(values))),p,"doublets",list())
  testthat::expect_null(cfg$invalid_reason)
  testthat::expect_null(cfg$capture)
  testthat::expect_null(cfg$expected_rate)
  testthat::expect_null(cfg$uncertainty)
  testthat::expect_identical(cfg$group_mode,"input")
  testthat::expect_identical(cfg$rates_mode,"unknown")
})

testthat::test_that("the shared method scores an unknown-rate provisional group without changing counts", {
  testthat::skip_if_not_installed("scDblFinder")
  f<-qc_diag_fixture();on.exit(unlink(f$root,recursive=TRUE))
  set.seed(14);raw<-scDblFinder::mockDoubletSCE(ncells=c(60,60),ngenes=300,dbl.rate=.1)
  record<-new_count_import(SummarizedExperiment::assay(raw,"counts"),"Public matrix",object=raw)
  b<-make_baseline_sce(record,prepare_sample_qc(record,annotation_fixture(),"gene_name"));before<-digest::digest(b)
  r<-qc_run_optional(b,"doublets",qc_exploratory_config(f$app_root),qc_baseline_signature(b),f$app_root)
  testthat::expect_identical(r$status,"scored")
  testthat::expect_true(all(is.finite(r$values$qc_diag_doublet_score)))
  testthat::expect_true(all(is.na(r$values$qc_diag_capture)))
  testthat::expect_identical(r$diagnostics$rate_mode,"unknown")
  testthat::expect_null(r$diagnostics$rate_table)
  testthat::expect_null(r$diagnostics$per_capture)
  testthat::expect_match(r$diagnostics$threshold_policy,"dbr.sd=1",fixed=TRUE)
  testthat::expect_true(is.numeric(r$diagnostics$internal_dbr_per1k_default))
  testthat::expect_identical(r$diagnostics$per_group$group,"Public matrix")
  testthat::expect_identical(digest::digest(b),before)
})

testthat::test_that("the exploratory Shiny mode reuses groups and invalidates changed scope", {
  f<-qc_diag_fixture();on.exit(unlink(f$root,recursive=TRUE))
  b<-f$baseline;b$geo_accession<-rep(c("GSM1","GSM2"),each=4)
  server<-function(input,output,session) {
    calls<-shiny::reactiveVal(0L)
    engine<-function(object,method,config,baseline_hash,app_root,cohort) qc_run_optional(object,method,config,baseline_hash,app_root,cohort,
      execute=function(...) {calls(calls()+1L);qc_diag_fake_fit(...)})
    review<-qc_review_server(input,output,session,function() list(A=b),optional_results=function(n) diagnostics$results(n))
    diagnostics<-qc_diagnostics_server(input,output,session,function() list(A=b),review,f$app_root,services=list(run=engine))
  }
  shiny::testServer(server,{
    action<-function(n) structure(as.integer(n),class="shinyActionButtonValue")
    session$setInputs(selectSample="A",qcReviewTab="optional",qcOptionalTab="doublets",qcRun_doublets=action(0))
    cfg<-list(doublet_mode="exploratory",rates_mode="unknown",expected_rate=NA_real_,uncertainty=NA_real_,rate_source="")
    qc_diag_set_controls(session,diagnostics,"doublets",cfg)
    session$setInputs(qcRun_doublets=action(1));session$flushReact()
    testthat::expect_identical(diagnostics$results("A")$doublets$status,"scored")
    testthat::expect_true(all(is.na(diagnostics$cell_table()$qc_diag_capture)))
    for(id in c("qcDoubletDistribution","qcDoubletJoint","qcDoubletCounts")) testthat::expect_true(is.list(output[[id]]))
    testthat::expect_equal(calls(),1L)
    cfg$provisional_group_mode<-"metadata";cfg$provisional_group_column<-"geo_accession"
    qc_diag_set_controls(session,diagnostics,"doublets",cfg)
    testthat::expect_identical(diagnostics$results("A")$doublets$status,"stale_settings")
    testthat::expect_true(all(is.na(diagnostics$cell_table()$qc_diag_doublet_score)))
    session$setInputs(qcRun_doublets=action(2));session$flushReact()
    testthat::expect_setequal(diagnostics$cell_table()$qc_diag_doublet_group,c("GSM1","GSM2"))
    cfg$doublet_mode<-"capture";qc_diag_set_controls(session,diagnostics,"doublets",cfg)
    testthat::expect_identical(diagnostics$results("A")$doublets$status,"stale_settings")
    session$setInputs(qcRun_doublets=action(3));session$flushReact()
    testthat::expect_identical(diagnostics$results("A")$doublets$status,"not_assessed")
    testthat::expect_equal(calls(),2L)
    testthat::expect_null(review$filtered()$A)
  })
})
