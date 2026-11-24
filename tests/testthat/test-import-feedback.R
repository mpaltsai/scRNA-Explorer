testthat::test_that("changing the input selection clears prior completion and readiness", {
  f <- input_fixture(); on.exit(unlink(f$root,recursive=TRUE))
  server <- function(input,output,session) {
    controller <- input_import_server(input,output,session,function(directory,selected_input) {
      list(sample=new_count_import(f$counts))
    })
  }
  shiny::testServer(server,{
    session$setInputs(testData="2",typeData="1",startQC=0)
    session$setInputs(startQC=1)
    session$flushReact()
    testthat::expect_match(output$readData,"Reading data completed")
    session$setInputs(typeData="2")
    testthat::expect_match(output$readData,"Not yet read")
    testthat::expect_null(controller$ready())
    testthat::expect_null(controller$records())
    session$setInputs(startQC=2)
    session$flushReact()
    testthat::expect_match(output$readData,"Reading data completed")
    session$setInputs(mtxInput=data.frame(name="new.mtx",datapath="/new-upload"))
    testthat::expect_match(output$readData,"Not yet read")
    testthat::expect_null(controller$ready())
  })
})

testthat::test_that("reading status is flushed before the reader and uses a selection snapshot", {
  f <- input_fixture(); on.exit(unlink(f$root,recursive=TRUE))
  server <- function(input,output,session) {
    probe <- new.env(parent=emptyenv())
    probe$status <- character(); probe$sent <- character(); probe$selected <- list(); probe$busy <- logical()
    controller <- input_import_server(input,output,session,function(directory,selected_input) {
      probe$status <- c(probe$status,controller$status())
      probe$sent <- c(probe$sent,output$readData)
      probe$busy <- c(probe$busy,controller$busy())
      probe$selected[[length(probe$selected)+1L]] <- selected_input
      list(sample=new_count_import(f$counts))
    })
  }
  shiny::testServer(server,{
    session$setInputs(testData="2",typeData="1",startQC=0)
    session$setInputs(startQC=1)
    session$flushReact()
    session$setInputs(startQC=2)
    session$flushReact()
    testthat::expect_length(probe$status,2)
    testthat::expect_true(all(grepl("Reading data",probe$status)))
    testthat::expect_false(any(grepl("completed",probe$status)))
    testthat::expect_true(all(grepl("Reading data",probe$sent)))
    testthat::expect_false(any(grepl("completed",probe$sent)))
    testthat::expect_identical(probe$selected[[2]]$typeData,"1")
    testthat::expect_match(output$readData,"Reading data completed")
    testthat::expect_true(all(probe$busy))
    testthat::expect_false(controller$busy())
  })
})

testthat::test_that("a selection superseded during a read cannot publish its old result", {
  f <- input_fixture(); on.exit(unlink(f$root,recursive=TRUE))
  server <- function(input,output,session) {
    controller <- input_import_server(input,output,session,function(directory,selected_input) {
      if (identical(selected_input$typeData,"1")) session$setInputs(typeData="2")
      list(sample=new_count_import(f$counts))
    })
  }
  shiny::testServer(server,{
    session$setInputs(testData="2",typeData="1",startQC=0)
    session$flushReact()
    session$setInputs(startQC=1)
    session$flushReact()
    testthat::expect_match(output$readData,"Not yet read")
    testthat::expect_null(controller$ready())
    testthat::expect_false(controller$busy())
    session$setInputs(startQC=2)
    session$flushReact()
    testthat::expect_match(output$readData,"Reading data completed")
    # Editing QC parameters or an inactive upload does not discard a valid read.
    session$setInputs(typeOrganism="Homo sapiens",csvFile=data.frame(name="hidden.csv",datapath="/hidden"))
    testthat::expect_false(is.null(controller$ready()))
    testthat::expect_match(output$readData,"Reading data completed")
  })
})
