testthat::test_that("Shiny import pauses duplicates, records reasons and clears stale results on failure", {
  testthat::skip_if_not_installed("shiny")
  f <- input_fixture(); on.exit(unlink(f$root,recursive=TRUE))
  p <- file.path(f$root,"counts.tsv")
  utils::write.table(data.frame(feature=c("gene","gene","other"),as.matrix(f$counts),check.names=FALSE),
    p,sep="\t",row.names=FALSE,quote=FALSE)
  server <- function(input,output,session) {
    controller <- input_import_server(input,output,session,function(directory,selected_input) read_app_input(selected_input,directory))
    # This is the same ready-record gate used by the app's annotation input.
    available <- shiny::reactive({x<-controller$ready();if(is.null(x)) NULL else import_counts(x)})
  }
  owned <- NULL
  shiny::testServer(server,{
    owned <<- controller$directory
    session$setInputs(testData="2",typeData="1",csvFile=data.frame(name="counts.tsv",datapath=p))
    session$setInputs(startQC=1)
    session$flushReact()
    testthat::expect_null(controller$ready())
    testthat::expect_null(available())
    testthat::expect_match(controller$status(),"paused")
    testthat::expect_match(output$importDuplicateTable,"gene",fixed=TRUE)
    session$setInputs(importDecisionSample="sample",importDecisionAction="keep_separate",importDecisionReason="")
    session$setInputs(resolveImport=1)
    testthat::expect_null(controller$ready())
    testthat::expect_match(controller$status(),"reason")
    session$setInputs(importDecisionReason="Retain source rows separately")
    session$setInputs(resolveImport=2)
    testthat::expect_false(is.null(available()))
    testthat::expect_equal(unname(as.matrix(available()$sample)),unname(as.matrix(f$counts)))
    testthat::expect_identical(controller$decisions()$reason,"Retain source rows separately")
    # Repeating a decision cannot change a ready import or invent another record.
    session$setInputs(resolveImport=3)
    testthat::expect_equal(nrow(controller$decisions()),1)
    testthat::expect_match(controller$status(),"No unresolved")
    session$setInputs(csvFile=data.frame(name="bad.tsv",datapath=file.path(f$root,"absent")))
    session$setInputs(startQC=2)
    session$flushReact()
    testthat::expect_null(controller$ready())
    testthat::expect_null(available())
    testthat::expect_match(controller$status(),"Import stopped")
    testthat::expect_false(controller$busy())
  })
  testthat::expect_false(dir.exists(owned))
})

testthat::test_that("all samples must resolve and cancellation cannot leak to annotation", {
  testthat::skip_if_not_installed("shiny")
  f <- input_fixture(); on.exit(unlink(f$root,recursive=TRUE))
  x <- f$counts; rownames(x) <- c("a","a","b")
  server <- function(input,output,session) {
    controller <- input_import_server(input,output,session,function(directory,selected_input)
      list(A=new_count_import(x,"A"),B=new_count_import(x,"B")))
  }
  shiny::testServer(server,{
    session$setInputs(startQC=0)
    session$flushReact()
    session$setInputs(startQC=1)
    session$flushReact()
    session$setInputs(importDecisionSample="A",importDecisionAction="keep_separate",importDecisionReason="Keep rows")
    session$setInputs(resolveImport=1)
    testthat::expect_null(controller$ready())
    testthat::expect_identical(controller$records()$A$status,"ready")
    session$setInputs(importDecisionSample="B",importDecisionAction="cancel",importDecisionReason="Correct B")
    session$setInputs(resolveImport=2)
    testthat::expect_null(controller$ready())
    testthat::expect_identical(controller$records()$B$status,"cancelled")
    testthat::expect_equal(nrow(controller$decisions()),2)
  })
})

testthat::test_that("Shiny sessions own distinct working directories", {
  testthat::skip_if_not_installed("shiny")
  directories <- character()
  server <- function(input,output,session) {
    controller <- input_import_server(input,output,session,function(directory,selected_input) NULL)
  }
  for (i in 1:2) shiny::testServer(server,{
    directories <<- c(directories,controller$directory)
    testthat::expect_true(dir.exists(controller$directory))
  })
  testthat::expect_identical(anyDuplicated(directories),0L)
  testthat::expect_false(any(dir.exists(directories)))
})

testthat::test_that("the edited app server and input UI parse successfully", {
  testthat::expect_silent(parse(file.path("..","..","site_data","server.R")))
  testthat::expect_silent(parse(file.path("..","..","site_data","tab_1.R")))
})

testthat::test_that("the actual app server retains imported metadata and invalidates data after a failed read", {
  testthat::skip_if_not_installed("shiny")
  testthat::skip_if_not_installed("shinyjs")
  f <- input_fixture(); on.exit(unlink(f$root,recursive=TRUE))
  object <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=f$counts))
  object$geo_accession <- rep("GSM101",ncol(object))
  S4Vectors::metadata(object)$geo <- list(source="fixture")
  p <- file.path(f$root,"sce.rds"); saveRDS(object,p)
  env <- new.env(parent=environment())
  env$`%>%` <- dplyr::`%>%`
  # Shiny runs the application with its own directory as the working directory.
  withr::local_dir(file.path("..","..","site_data"))
  app_server <- source("server.R",local=env)$value
  shiny::testServer(app_server,{
    session$setInputs(testData="2",typeData="3",startQC=0,applySampleNames=0)
    session$setInputs(rdsFile=data.frame(name="sce.rds",datapath=p))
    session$setInputs(startQC=1)
    session$flushReact()
    testthat::expect_identical(my_react_vals$d_to_metadata$sample,f$counts)
    testthat::expect_identical(my_react_vals$import_records$sample$source_object,object)
    stable_ids <- my_react_vals$import_records$sample$cell_identity
    my_react_vals$filt_samples_list <- list(sample="old result")
    my_react_vals$loaded.dataSO.combined.no.cluster <- "old Seurat result"
    session$setInputs(importSampleName_1="Patient 1",applySampleNames=1)
    testthat::expect_identical(names(my_react_vals$d_to_metadata),"Patient 1")
    testthat::expect_identical(my_react_vals$import_records[["Patient 1"]]$cell_identity,stable_ids)
    testthat::expect_identical(my_react_vals$import_records[["Patient 1"]]$source_object,object)
    testthat::expect_identical(my_react_vals$import_records[["Patient 1"]]$source$upload_names$rdsFile,"sce.rds")
    testthat::expect_null(my_react_vals$filt_samples_list)
    testthat::expect_null(my_react_vals$loaded.dataSO.combined.no.cluster)
    my_react_vals$filt_samples_list <- list("Patient 1"="old result")
    my_react_vals$loaded.dataSO.combined.no.cluster <- "old Seurat result"
    session$setInputs(rdsFile=data.frame(name="missing.rds",datapath=file.path(f$root,"missing")))
    session$setInputs(startQC=2)
    session$flushReact()
    testthat::expect_null(my_react_vals$d_to_metadata)
    testthat::expect_null(my_react_vals$import_records)
    testthat::expect_null(my_react_vals$filt_samples_list)
    testthat::expect_null(my_react_vals$loaded.dataSO.combined.no.cluster)
    testthat::expect_match(output$readData,"Import stopped")
  })
})
