testthat::test_that("standard QC plots handle missing and invalidated baselines with Bioconductor attached", {
  # Match the app's attached packages, not just loaded namespaces: S4 dispatch
  # must not turn an ordinary Shiny readiness condition into a console error.
  suppressPackageStartupMessages(withr::local_package("SingleCellExperiment"))
  baselines<-list(A=qc_review_fixture("A"));before<-digest::digest(baselines)
  server<-function(input,output,session) {
    b<-shiny::reactiveVal(NULL)
    controller<-qc_review_server(input,output,session,function() b())
  }
  shiny::testServer(server,{
    panels<-c(overview="qcSelectedOverview",ribosomal="qcRibosomalPlot",annotation="qcAnnotationCoverage")
    for(state in c("startup","ready","removed","restored")) {
      ready<-state %in% c("ready","restored")
      b(if(ready) baselines else NULL);session$flushReact()
      session$setInputs(selectSample="A")
      for(tab in names(panels)) {
        session$setInputs(qcReviewTab=tab)
        value<-tryCatch(output[[panels[[tab]]]],error=identity)
        if(ready) {
          testthat::expect_true(is.list(value) && length(value$src)==1L && nzchar(value$src),
            info=paste(state,tab))
        } else {
          testthat::expect_true(inherits(value,"shiny.silent.error"),info=paste(state,tab))
        }
      }
    }
    testthat::expect_null(controller$filtered()$A)
    testthat::expect_length(controller$history()$A,0L)
    testthat::expect_identical(digest::digest(baselines),before)
  })
})
