testthat::test_that("deferred optional checks clean up outside reactive consumers and allow correction", {
  suppressPackageStartupMessages(withr::local_package("SingleCellExperiment"))
  f<-qc_diag_fixture();on.exit(unlink(f$root,recursive=TRUE))
  before<-digest::digest(f$baseline)
  control<-new.env(parent=emptyenv());control$fail<-FALSE;control$fits<-0L
  # MockShinySession normally invokes onFlushed callbacks inside isolate(),
  # unlike the real session. Queue the real controller callback and invoke it
  # below with a session domain but deliberately no reactive consumer.
  session_class<-R6::R6Class(inherit=shiny::MockShinySession,portable=FALSE,lock_objects=FALSE,
    parent_env=asNamespace("shiny"),public=list(qcPendingCallbacks=list(),
    onFlushed=function(fun,once=TRUE) {
      self$qcPendingCallbacks<-c(self$qcPendingCallbacks,list(fun))
      invisible(function() NULL)
    }))
  session<-session_class$new();on.exit(session$close(),add=TRUE)
  shiny::withReactiveDomain(session,{
    review<-qc_review_server(session$input,session$output,session,function() list(A=f$baseline),
      optional_results=function(n) diagnostics$results(n))
    engine<-function(object,method,config,baseline_hash,app_root,cohort) {
      if(control$fail) stop("Simulated diagnostic failure")
      qc_run_optional(object,method,config,baseline_hash,app_root,cohort,
        execute=function(...) {control$fits<-control$fits+1L;qc_diag_fake_fit(...)})
    }
    diagnostics<-qc_diagnostics_server(session$input,session$output,session,
      function() list(A=f$baseline),review,f$app_root,services=list(run=engine))
  })
  complete_run<-function() {
    testthat::expect_length(session$qcPendingCallbacks,1L)
    callback<-session$qcPendingCallbacks[[1L]];session$qcPendingCallbacks<-session$qcPendingCallbacks[-1L]
    error<-tryCatch({shiny::withReactiveDomain(session,callback());NULL},error=identity)
    testthat::expect_true(is.null(error),info=if(inherits(error,"error")) conditionMessage(error) else NULL)
    if(!is.null(error)) return(FALSE)
    session$flushReact()
    testthat::expect_length(shiny::isolate(diagnostics$busy()),0L)
    testthat::expect_false(session$isClosed())
    TRUE
  }
  action<-function(n) structure(as.integer(n),class="shinyActionButtonValue")
  session$setInputs(selectSample="A",qcRun_doublets=action(0),qcReviewTab="optional",qcOptionalTab="doublets")
  for(stage in c("missing_capture","missing_rates","runner_error","corrected")) {
    settings<-list(expected_rate=NA_real_,uncertainty=NA_real_,rate_source="")
    if(stage!="missing_capture") settings<-c(settings,list(capture_mode="single",capture_id="simulation",
      capture_complete=TRUE,capture_evidence="Known synthetic capture"))
    if(stage %in% c("runner_error","corrected")) {
      settings$expected_rate<-.1;settings$uncertainty<-.02;settings$rate_source<-"Simulation settings"
    }
    control$fail<-stage=="runner_error"
    qc_diag_set_controls(session,diagnostics,"doublets",settings,app_root=f$app_root)
    session$setInputs(qcRun_doublets=action(match(stage,c("missing_capture","missing_rates","runner_error","corrected"))))
    if(!complete_run()) break
    notice<-session$getOutput("qcOptionalNotice")
    if(stage %in% c("missing_capture","missing_rates")) {
      result<-shiny::isolate(diagnostics$results("A")$doublets)
      testthat::expect_identical(result$status,"not_assessed")
      testthat::expect_true(all(is.na(result$values$qc_diag_doublet_score)))
      testthat::expect_match(notice,if(stage=="missing_capture") "Capture identity is unresolved" else "rate table")
      testthat::expect_equal(control$fits,0L)
    } else if(stage=="runner_error") {
      testthat::expect_match(notice,"Optional check stopped: Simulated diagnostic failure",fixed=TRUE)
      testthat::expect_equal(control$fits,0L)
    } else {
      testthat::expect_identical(shiny::isolate(diagnostics$results("A")$doublets$status),"scored")
      testthat::expect_equal(control$fits,1L)
    }
  }
  if(control$fits==1L) {
    # Marker review returns early without fitting a model, but its callback
    # must still clear the busy state safely outside a reactive consumer.
    session$setInputs(qcReviewCycleMarkers=action(0))
    qc_diag_set_controls(session,diagnostics,"cell_cycle",app_root=f$app_root)
    session$setInputs(qcReviewCycleMarkers=action(1))
    if(complete_run()) {
      testthat::expect_match(session$getOutput("qcOptionalNotice"),"Marker coverage reviewed",fixed=TRUE)
      testthat::expect_equal(control$fits,1L)
    }
  }
  testthat::expect_null(shiny::isolate(review$filtered()$A))
  testthat::expect_identical(digest::digest(f$baseline),before)
})
