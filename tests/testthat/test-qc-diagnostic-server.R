qc_diag_action<-function(n) structure(as.integer(n),class="shinyActionButtonValue")

testthat::test_that("unstarted optional plots keep Shiny validation with Bioconductor packages attached", {
  # The app attaches SingleCellExperiment. Its S4 generics can wrap a validation
  # condition raised while forcing a plotting argument into an ordinary error.
  # Namespace-only fixtures miss that difference in the app's search path.
  suppressPackageStartupMessages(withr::local_package("SingleCellExperiment"))
  f<-qc_diag_fixture();on.exit(unlink(f$root,recursive=TRUE))
  before<-digest::digest(f$baseline)
  server<-function(input,output,session) {
    calls<-shiny::reactiveVal(0L)
    engine<-function(...) {calls(calls()+1L);stop("Opening a panel must not run a model.")}
    review<-qc_review_server(input,output,session,function() list(A=f$baseline),
      optional_results=function(n) diagnostics$results(n))
    diagnostics<-qc_diagnostics_server(input,output,session,function() list(A=f$baseline),
      review,f$app_root,services=list(run=engine))
  }
  shiny::testServer(server,{
    panels<-list(doublets=c("qcDoubletDistribution","qcDoubletJoint","qcDoubletCounts","qcDoubletUmiGenes"),
      ambient=c("qcAmbientDistribution","qcAmbientJoint","qcAmbientPopulations","qcAmbientFit"),
      cell_cycle=c("qcCycleScores","qcCycleCounts","qcCycleQC"))
    for(method in names(panels)) {
      session$setInputs(selectSample="A",qcReviewTab="optional",qcOptionalTab=method)
      for(id in panels[[method]]) {
        condition<-tryCatch({output[[id]];NULL},error=identity)
        testthat::expect_true(inherits(condition,"shiny.silent.error"),info=id)
        testthat::expect_match(conditionMessage(condition),"Review settings and run this optional check first.",fixed=TRUE,info=id)
      }
    }
    testthat::expect_equal(calls(),0L)
    testthat::expect_length(diagnostics$results("A"),0L)
    testthat::expect_null(review$filtered()$A)
    testthat::expect_identical(digest::digest(f$baseline),before)
  })
})

testthat::test_that("optional doublet review preserves filtering state and invalidates changed or shared-capture evidence", {
  f<-qc_diag_fixture();on.exit(unlink(f$root,recursive=TRUE))
  b2<-f$baseline;b2$qc_sample<-rep("B",ncol(b2));b2$qc_cell_id<-paste0("B::",colnames(b2))
  source_hash<-digest::digest(f$baseline)
  server<-function(input,output,session) {
    b<-shiny::reactiveVal(list(A=f$baseline,B=b2))
    fits<-shiny::reactiveVal(0L)
    engine<-function(object,method,config,baseline_hash,app_root,cohort) qc_run_optional(object,method,config,baseline_hash,app_root,cohort,
      execute=function(...) {fits(fits()+1L);qc_diag_fake_fit(...)})
    review<-qc_review_server(input,output,session,function() b(),optional_results=function(n) diagnostics$results(n))
    diagnostics<-qc_diagnostics_server(input,output,session,function() b(),review,f$app_root,services=list(run=engine))
  }
  shiny::testServer(server,{
    session$setInputs(selectSample="A",qcView="preview",qcReviewTab="overview",qcRun_doublets=qc_diag_action(0),qcPreview=qc_diag_action(0),QCFiltering=qc_diag_action(0))
    testthat::expect_equal(fits(),0)
    qc_diag_set_controls(session,diagnostics,"doublets")
    session$setInputs(qcRun_doublets=qc_diag_action(1));session$flushReact()
    r<-diagnostics$results("A")$doublets
    testthat::expect_identical(r$status,"not_assessed")
    testthat::expect_match(paste(r$reasons,collapse=" "),"unresolved")
    testthat::expect_true(all(is.na(r$values$qc_diag_doublet_score)))
    testthat::expect_equal(fits(),0)
    qc_diag_set_controls(session,diagnostics,"doublets",list(capture_mode="single",capture_id="capA",capture_complete=TRUE,capture_evidence="Synthetic complete capture"))
    session$setInputs(qcRun_doublets=qc_diag_action(2));session$flushReact()
    r<-diagnostics$results("A")$doublets
    testthat::expect_identical(r$status,"scored")
    testthat::expect_equal(fits(),1)
    testthat::expect_null(review$filtered()$A)
    t<-qc_default_thresholds(review$metrics());qc_test_set_draft(session,review,t)
    session$setInputs(qcPreview=qc_diag_action(1));session$setInputs(QCFiltering=qc_diag_action(1),qcApplyReason="Keep baseline despite optional flags")
    testthat::expect_equal(ncol(review$filtered()$A),8)
    testthat::expect_true(any(review$filtered()$A$qc_diag_doublet_class=="doublet"))
    testthat::expect_identical(review$history()$A[[1]]$optional_evidence$doublets$revision,r$revision)
    applied_hash<-digest::digest(review$filtered()$A)
    session$setInputs(qcReviewTab="optional",qcOptionalTab="doublets")
    testthat::expect_true(is.list(output$qcDoubletDistribution))
    testthat::expect_true(is.list(output$qcDoubletJoint))
    testthat::expect_true(is.list(output$qcDoubletCounts))
    testthat::expect_true(is.list(output$qcDoubletUmiGenes))
    p<-diagnostics$prefix();do.call(session$setInputs,stats::setNames(list(2),paste0(p,"seed")))
    stale<-diagnostics$results("A")$doublets
    testthat::expect_identical(stale$status,"stale_settings")
    testthat::expect_true(all(is.na(stale$values$qc_diag_doublet_score)))
    testthat::expect_identical(digest::digest(review$filtered()$A),applied_hash)
    exported<-diagnostics$decorate(review$filtered()$A,"A")
    testthat::expect_true(all(is.na(exported$qc_diag_doublet_score)))
    do.call(session$setInputs,stats::setNames(list(1),paste0(p,"seed")))
    testthat::expect_identical(diagnostics$results("A")$doublets$status,"scored")
    session$setInputs(selectSample="B",qcReviewTab="overview")
    testthat::expect_false(identical(p,diagnostics$prefix()))
    session$setInputs(qcRun_doublets=qc_diag_action(3));session$flushReact()
    testthat::expect_null(diagnostics$results("B")$doublets)
    qc_diag_set_controls(session,diagnostics,"doublets",list(capture_mode="single",capture_id="capA",capture_complete=TRUE,capture_evidence="Conflicting synthetic declaration"))
    testthat::expect_identical(diagnostics$results("A")$doublets$status,"stale_settings")
    session$setInputs(qcRun_doublets=qc_diag_action(4));session$flushReact()
    testthat::expect_identical(diagnostics$results("B")$doublets$status,"not_assessed")
    testthat::expect_match(paste(diagnostics$results("B")$doublets$reasons,collapse=" "),"multiple")
    testthat::expect_equal(fits(),1)
    testthat::expect_identical(digest::digest(f$baseline),source_hash)
    # Even a previously visited sample needs freshly bound controls. An old
    # sample's queued Run click must not start work on the new selection.
    session$setInputs(selectSample="A",qcRun_doublets=qc_diag_action(5));session$flushReact()
    testthat::expect_equal(fits(),1)
    testthat::expect_match(output$qcOptionalNotice,"Wait for the selected sample")
    b(NULL);session$flushReact()
    testthat::expect_length(diagnostics$results("A"),0)
  })
})

testthat::test_that("ambient and cycle panels keep fit provenance and compare subsets without refitting", {
  f<-qc_diag_fixture();on.exit(unlink(f$root,recursive=TRUE))
  server<-function(input,output,session) {
    calls<-shiny::reactiveVal(0L)
    engine<-function(object,method,config,baseline_hash,app_root,cohort) qc_run_optional(object,method,config,baseline_hash,app_root,cohort,
      execute=function(...) {calls(calls()+1L);qc_diag_fake_fit(...)})
    review<-qc_review_server(input,output,session,function() list(A=f$baseline),optional_results=function(n) diagnostics$results(n))
    diagnostics<-qc_diagnostics_server(input,output,session,function() list(A=f$baseline),review,f$app_root,services=list(run=engine))
  }
  shiny::testServer(server,{
    session$setInputs(selectSample="A",qcReviewTab="overview",qcRun_ambient=qc_diag_action(0),qcRun_cell_cycle=qc_diag_action(0),
      qcReviewCycleMarkers=qc_diag_action(0),qcPreview=qc_diag_action(0),QCFiltering=qc_diag_action(0))
    qc_diag_set_controls(session,diagnostics,"ambient")
    session$setInputs(qcRun_ambient=qc_diag_action(1));session$flushReact()
    testthat::expect_identical(diagnostics$results("A")$ambient$status,"not_assessed")
    testthat::expect_equal(calls(),0)
    qc_diag_set_controls(session,diagnostics,"ambient",list(pooled=TRUE))
    session$setInputs(qcRun_ambient=qc_diag_action(2));session$flushReact()
    testthat::expect_identical(diagnostics$results("A")$ambient$status,"estimated")
    testthat::expect_equal(calls(),1)
    qc_diag_set_controls(session,diagnostics,"cell_cycle")
    session$setInputs(qcReviewCycleMarkers=qc_diag_action(1));session$flushReact()
    testthat::expect_true(all(diagnostics$marker_review()$coverage$matched>0))
    testthat::expect_equal(calls(),1)
    session$setInputs(qcRun_cell_cycle=qc_diag_action(1));session$flushReact()
    testthat::expect_identical(diagnostics$results("A")$cell_cycle$status,"scored")
    testthat::expect_equal(calls(),2)
    t<-qc_default_thresholds(review$metrics());t$maxUMIsPerCell<-stats::median(review$metrics()$nUMI)
    t$maxUMIsPerCell<-floor(t$maxUMIsPerCell);qc_test_set_draft(session,review,t)
    session$setInputs(qcPreview=qc_diag_action(1),qcView="preview")
    rows<-diagnostics$compared("cell_cycle")
    testthat::expect_equal(nrow(rows),8+sum(review$preview()$keep))
    testthat::expect_equal(calls(),2)
    session$setInputs(QCFiltering=qc_diag_action(1),qcApplyReason="Apply the reviewed simple QC limits",qcView="applied")
    testthat::expect_true(all(c("qc_diag_S_score","qc_diag_ambient_fraction") %in% names(SummarizedExperiment::colData(review$filtered()$A))))
    session$setInputs(qcReviewTab="optional",qcOptionalTab="ambient")
    testthat::expect_true(is.list(output$qcAmbientDistribution))
    testthat::expect_true(is.list(output$qcAmbientJoint))
    testthat::expect_true(is.list(output$qcAmbientPopulations))
    testthat::expect_true(is.list(output$qcAmbientFit))
    session$setInputs(qcOptionalTab="cell_cycle")
    testthat::expect_true(is.list(output$qcCycleScores))
    testthat::expect_true(is.list(output$qcCycleCounts))
    testthat::expect_true(is.list(output$qcCycleQC))
    testthat::expect_equal(calls(),2)
    testthat::expect_true(all(is.finite(diagnostics$cell_table()$qc_diag_ambient_fraction)))
    testthat::expect_length(diagnostics$history()$A,3)
  })
})
