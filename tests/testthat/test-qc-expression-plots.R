# Exercise the real controller and unchanged scater plotting helpers. Resolve
# their dependencies privately instead of attaching the app's entire package set.
qc_expression_fixture <- function(sample="A") {
  # Enough distinct observations for the legacy plots' linear-fit intervals;
  # the regression does not depend on degenerate two-point fits.
  x<-outer(seq_len(24),seq_len(12),function(g,c) ifelse(g<=c+6,10*(1+(g*c)%%7),0))
  dimnames(x)<-list(c("MT-A","HBA","GeneA",paste0("Other",4:24)),paste0("cell",1:12))
  r<-new_count_import(Matrix::Matrix(x,sparse=TRUE),sample)
  make_baseline_sce(r,prepare_sample_qc(r,annotation_fixture(),"gene_name"))
}

qc_expression_server <- function(baselines) {
  plot_env<-new.env(parent=asNamespace("ggplot2"))
  plot_env$`%>%`<-dplyr::`%>%`;plot_env$lm<-stats::lm
  plot_env$plotScater<-scater::plotScater
  plot_env$plotHighestExprs<-scater::plotHighestExprs
  sys.source(file.path("..","..","site_data","R","createQCPlots.R"),envir=plot_env)
  e<-new.env(parent=environment(qc_review_server))
  for(name in ls(plot_env,pattern="^create")) e[[name]]<-plot_env[[name]]
  controller_function<-qc_review_server;environment(controller_function)<-e
  function(input,output,session) {
    controller<-controller_function(input,output,session,function() baselines)
  }
}

qc_action <- function(n) structure(as.integer(n),class="shinyActionButtonValue")
qc_has_image <- function(output,id) tryCatch({
  z<-output[[id]];is.list(z) && length(z$src)==1L && nzchar(z$src)
},error=function(e) FALSE)

# The installed scater already warns about this deprecated internal call.
# Keep other warnings visible so this helper cannot mask rendering regressions.
qc_scater_warnings <- function(expr) withCallingHandlers(expr,warning=function(w) {
  if(grepl("perCellQCMetrics.*deprecated",conditionMessage(w))) invokeRestart("muffleWarning")
})

testthat::test_that("both applied images render when plot clicks and Apply share an input batch", {
  testthat::skip_if_not_installed("scater")
  b<-list(A=qc_expression_fixture("A"));before<-digest::digest(b)
  qc_scater_warnings(shiny::testServer(qc_expression_server(b),{
    session$setInputs(selectSample="A",qcReviewTab="overview",qcPreview=qc_action(0),QCFiltering=qc_action(0),
      makeLibraryFiltered=qc_action(0),makeHighExprFiltered=qc_action(0))
    t<-qc_default_thresholds(controller$metrics());t$minGenesPerCell<-1
    qc_test_set_draft(session,controller,t);session$setInputs(qcPreview=qc_action(1))
    # Real action-button values and the failing delivery order: requests can be
    # received before Apply in the same browser input batch.
    session$setInputs(makeLibraryFiltered=qc_action(1),makeHighExprFiltered=qc_action(1),
      QCFiltering=qc_action(1),qcApplyReason="Accept fixture baseline")
    session$setInputs(qcReviewTab="existing")
    testthat::expect_equal(ncol(controller$filtered()$A),12)
    testthat::expect_true(qc_has_image(output,"libraryPlotFiltered"))
    testthat::expect_true(qc_has_image(output,"highExprPlotFiltered"))
    revision<-controller$revisions()[["A"]]
    t$maxUMIsPerCell<-sort(controller$metrics()$nUMI)[6];qc_test_set_draft(session,controller,t)
    session$setInputs(qcPreview=qc_action(2))
    session$setInputs(makeLibraryFiltered=qc_action(2),makeHighExprFiltered=qc_action(2),
      QCFiltering=qc_action(2),qcApplyReason="Accept fixture upper cutoff")
    testthat::expect_equal(ncol(controller$filtered()$A),6)
    testthat::expect_false(identical(controller$revisions()[["A"]],revision))
    testthat::expect_true(qc_has_image(output,"libraryPlotFiltered"))
    testthat::expect_true(qc_has_image(output,"highExprPlotFiltered"))
    session$setInputs(qcView="applied")
    testthat::expect_true(qc_has_image(output,"libraryPlotFiltered"))
    testthat::expect_true(qc_has_image(output,"highExprPlotFiltered"))
  }))
  testthat::expect_identical(digest::digest(b),before)
})

testthat::test_that("Preview, Apply and plotting in one batch follow the decision order", {
  testthat::skip_if_not_installed("scater")
  qc_scater_warnings(shiny::testServer(qc_expression_server(list(A=qc_expression_fixture("A"))),{
    session$setInputs(selectSample="A",qcReviewTab="overview",qcPreview=qc_action(0),QCFiltering=qc_action(0),
      makeLibraryFiltered=qc_action(0),makeHighExprFiltered=qc_action(0))
    t<-qc_default_thresholds(controller$metrics());t$minGenesPerCell<-1
    qc_test_set_draft(session,controller,t)
    session$setInputs(makeLibraryFiltered=qc_action(1),makeHighExprFiltered=qc_action(1),
      QCFiltering=qc_action(1),qcPreview=qc_action(1),qcApplyReason="Batched explicit decision")
    session$setInputs(qcReviewTab="existing")
    testthat::expect_equal(ncol(controller$filtered()$A),12)
    testthat::expect_length(controller$history()$A,1)
    testthat::expect_true(qc_has_image(output,"libraryPlotFiltered"))
    testthat::expect_true(qc_has_image(output,"highExprPlotFiltered"))
  }))
})

testthat::test_that("missing, empty and stale applied plot requests explain the next action", {
  testthat::skip_if_not_installed("scater")
  qc_scater_warnings(shiny::testServer(qc_expression_server(list(A=qc_expression_fixture("A"),B=qc_expression_fixture("B"))),{
    session$setInputs(selectSample="A",qcReviewTab="overview",qcPreview=qc_action(0),QCFiltering=qc_action(0),
      makeLibrary=qc_action(0),makeHighExpr=qc_action(0),makeLibraryFiltered=qc_action(0),makeHighExprFiltered=qc_action(0))
    t<-qc_default_thresholds(controller$metrics());t$minGenesPerCell<-1
    qc_test_set_draft(session,controller,t)
    session$setInputs(qcReviewTab="existing")
    testthat::expect_error(output$libraryPlotFiltered,"Apply filtering")
    testthat::expect_error(output$highExprPlotFiltered,"Apply filtering")
    session$setInputs(qcPreview=qc_action(1));session$setInputs(QCFiltering=qc_action(1),qcApplyReason="Keep A")
    testthat::expect_error(output$libraryPlotFiltered,"Click the plot button")
    session$setInputs(makeLibrary=qc_action(1),makeHighExpr=qc_action(1),makeLibraryFiltered=qc_action(1),makeHighExprFiltered=qc_action(1))
    testthat::expect_true(qc_has_image(output,"libraryPlot"))
    testthat::expect_true(qc_has_image(output,"highExprPlot"))
    testthat::expect_true(qc_has_image(output,"libraryPlotFiltered"))
    testthat::expect_true(qc_has_image(output,"highExprPlotFiltered"))
    t$maxUMIsPerCell<-sort(controller$metrics()$nUMI)[6];qc_test_set_draft(session,controller,t)
    session$setInputs(qcPreview=qc_action(2));session$setInputs(QCFiltering=qc_action(2),qcApplyReason="Revise A")
    testthat::expect_error(output$libraryPlotFiltered,"changed.*Click")
    testthat::expect_error(output$highExprPlotFiltered,"changed.*Click")
    # A changed applied revision leaves valid baseline plot requests intact.
    testthat::expect_true(qc_has_image(output,"libraryPlot"))
    testthat::expect_true(qc_has_image(output,"highExprPlot"))
    session$setInputs(makeLibraryFiltered=qc_action(2),makeHighExprFiltered=qc_action(2))
    testthat::expect_true(qc_has_image(output,"libraryPlotFiltered"))
    testthat::expect_true(qc_has_image(output,"highExprPlotFiltered"))
    session$setInputs(selectSample="B",qcReviewTab="overview")
    tb<-qc_default_thresholds(controller$metrics());tb$minGenesPerCell<-1;qc_test_set_draft(session,controller,tb)
    session$setInputs(qcReviewTab="existing")
    testthat::expect_error(output$libraryPlotFiltered,"Apply filtering")
    testthat::expect_error(output$highExprPlotFiltered,"Apply filtering")
    session$setInputs(selectSample="A",qcReviewTab="overview")
    t$minGenesPerCell<-4;t$maxGenesPerCell<-4;qc_test_set_draft(session,controller,t)
    session$setInputs(qcPreview=qc_action(3));session$setInputs(QCFiltering=qc_action(3),qcApplyReason="Test empty selection")
    session$setInputs(qcReviewTab="existing")
    session$setInputs(makeLibraryFiltered=qc_action(3),makeHighExprFiltered=qc_action(3))
    testthat::expect_equal(ncol(controller$filtered()$A),0)
    testthat::expect_error(output$libraryPlotFiltered,"No cells remain")
    testthat::expect_error(output$highExprPlotFiltered,"No cells remain")
  }))
})
