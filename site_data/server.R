#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyjs)
#library(tools)
#library(dplyr)
#library(utils)


options(shiny.maxRequestSize=15000*1024^2)  #max upload file size 15GB

function(input, output, session){
  
  #reactive values to store reactive objects
  my_react_vals <- reactiveValues()
  
  
  ##### Upload, setting parameters and read data routine #####
  
  #Validate uploaded .csv or .tsv file
  output$uploadCsvValidation <- reactive({
           if(tools::file_ext(input$csvFile$name)!= "csv" & tools::file_ext(input$csvFile$name)!= "tsv"){
           validate("Invalid file; Please upload a .csv or .tsv file")}
           
    }) %>% bindEvent(input$csvFile, ignoreInit = TRUE, ignoreNULL=TRUE)
  
  #Validate uploaded mtx (10X matrix)
  output$uploadmtxValidation <- reactive({
    if(tools::file_ext(input$mtxInput$name) != "mtx" & 
       substr(input$mtxInput$name, nchar(input$mtxInput$name)-6, nchar(input$mtxInput$name)-0) != ".mtx.gz"){
      validate("Invalid file; Please upload a .mtx or .mtx.gz file")}
  }) %>% bindEvent(input$mtxInput, ignoreInit = TRUE, ignoreNULL=TRUE)

  
  #Validate uploaded genes 10X
  output$uploadgenes10XValidation <- reactive({
    if(substr(input$genes10XInput$name, nchar(input$genes10XInput$name)-3, nchar(input$genes10XInput$name)-0) != ".tsv" & 
       substr(input$genes10XInput$name, nchar(input$genes10XInput$name)-6, nchar(input$genes10XInput$name)-0)!= ".tsv.gz"){
      validate("Invalid file; Please upload a .tsv or .tsv.gz file")}
  }) %>% bindEvent(input$genes10XInput, ignoreInit = TRUE, ignoreNULL=TRUE)

  
  #Validate uploaded barcodes 10X
  output$uploadbarcodes10XValidation <- reactive({
    if(substr(input$barcodes10XInput$name, nchar(input$barcodes10XInput$name)-3, nchar(input$barcodes10XInput$name)-0) != ".tsv" & 
       substr(input$barcodes10XInput$name, nchar(input$barcodes10XInput$name)-6, nchar(input$barcodes10XInput$name)-0)!= ".tsv.gz"){
      validate("Invalid file; Please upload a .tsv or .tsv.gz file")}
  }) %>% bindEvent(input$barcodes10XInput, ignoreInit = TRUE, ignoreNULL=TRUE)
  
  
  #Validate uploaded .rds file
  output$uploadRdsValidation <- reactive({
    validate(need(identical(tools::file_ext(input$rdsFile$name),"rds"),"Invalid extension (not a .rds file)"))
  }) %>% bindEvent(input$rdsFile, ignoreInit = TRUE, ignoreNULL=TRUE)
  
  #Validate uploaded .zip file
  output$uploadZipValidation <- reactive({
    validate(need(identical(tools::file_ext(input$zipFile$name),"zip"),"Invalid extension (not a .zip file)"))
  }) %>% bindEvent(input$zipFile, ignoreInit = TRUE, ignoreNULL=TRUE)
  
  #Check if there ism't any file uploaded or state which kind of input the user has provided
    user_input <- reactive({
      #if(length(input$testData)==0){
       # validate("You haven't chosen your data type yet")
      validate(need(input$testData, message = FALSE))
      if(input$testData == 1){
        count.matrix <- TRUE
        seurat <- TRUE
        inputFile <- "raw_se_S190.rds"
        inputDir <- FALSE
        my_samples_list <- c(my_list = readData(count.matrix, seurat, inputFile, inputDir), ############## <<- #############################
                    mes = "Reading data completed. The origin organism is Mus musculus. Please select the appropriate button below.")
        
      }else{
        validate(need(input$typeData, message = FALSE))
        if(input$typeData == 1){
          #req(input$csvFile$datapath)
          validate(need(input$csvFile$datapath, message = "Please upload a file"))
          count.matrix <- TRUE
          seurat <- FALSE
          inputFile <- input$csvFile$datapath
          inputDir <- FALSE
          my_samples_list <- c(my_list = readData(count.matrix, seurat, inputFile, inputDir),
                               mes = "Reading data completed ")
          
        }else if(input$typeData==3){
          req(input$rdsFile$datapath)
          count.matrix <- TRUE
          seurat <- TRUE
          inputFile <- input$rdsFile$datapath
          inputDir <- FALSE
          my_samples_list <- c(my_list = readData(count.matrix, seurat, inputFile, inputDir),
                               mes = "Reading data completed ")
          
        }else if(input$typeData==2){
          req(input$mtxInput$datapath)
          req(input$genes10XInput$datapath)
          my_samples_list <- c(my_list = my_read10X(mtxPath = input$mtxInput$datapath, 
                                                    genesPath = input$genes10XInput$datapath, 
                                                    barcodesPath = input$barcodes10XInput$datapath),
                               mes = "Reading data completed ")
          
          ###i put it inside my_read10X at the end 
          #counts <- as(input_counts, "dgCMatrix")
          
        }else if(input$typeData==4){
          req(input$zipFile$datapath)
          if(isTRUE(file.exists("./userFile"))){
            unlink("./userFile", recursive = TRUE)
          }
          u <- unzip(input$zipFile$datapath, exdir = "userFile")
          user_list_dir <- list.dirs("./userFile")
          my_samples_list <- c(my_list = list(readZip(user_list_dir)),
                               mes = "Reading data completed ")
        }
      } 
    })%>% bindEvent(input$startQC , ignoreInit = TRUE, ignoreNULL = TRUE)

   # output$readData <- renderText({
  #   need(input$typeData, message=FALSE)
      #req(input$typeData)
   #   my_datatest <<- user_input()$my_list
      #dd <- req(input$typeData)
      #if(dd==4){
    #    cond_satisf = c()
     #   for(x in user_input()$my_list){
      #    my_condition = is.null(base::rownames(x))
       #   cond_satisf <- c(cond_satisf, my_condition)
        #}
    #    if(any(cond_satisf)==TRUE){
     #     paste0("Genes names are absent for counts data in subdirectory ", which(cond_satisf==TRUE), ". Please upload a file with gene names as row names if a count matrix or check your data in any other case")
    #    }else{
    #      user_input()$mes
     #   }

    # })
    
    #observeEvent(input$typeData,{
      #hide("readData")
    #  output$readData <- renderText({NULL})
    #  output$paramsAddMetadata <- renderText({NULL})
   #   output$initDataFilt <- renderText({NULL})
   # })
    
    #empty outputs and reactive values when changing input data
    observeEvent(input$startQC,{
      my_react_vals$d_to_filter <- NULL
      my_react_vals$d_to_metadata <- NULL
    #hide("readData")
      if(input$startQC>0){
        output$readData <- renderText({NULL})
        output$paramsAddMetadata <- renderText({NULL})
        output$initDataFilt <- renderText({NULL})
        #shinyjs::reset("reset_inputs")
        
        my_react_vals$filt_samples_list <- NULL
        my_react_vals$loaded.dataSO.combined.no.cluster <- NULL
        my_react_vals$loaded.dataSO.combined.markers <- NULL
        my_react_vals$loaded.dataSO.combined.markerstop <- NULL
        my_react_vals$gen_cluster_res <- NULL
        my_react_vals$loaded.dataSO.combined <- NULL
        my_react_vals$iRNA_cor_tables <- NULL
        my_react_vals$iRNA_enrich_tables <- NULL
        #metrics_clean() <- NULL
        #se_c() <- NULL
        

      }
      
    })
    
    #to prevent withSpinner running when starting the app
    output$readData <- renderText({NULL})
    output$paramsAddMetadata <- renderText({NULL})
    output$initDataFilt <- renderText({NULL})
    
    
    #output readData and check if rownames (genes) are present
    observeEvent(input$startQC,{
        #validate(need(input$typeData, message=FALSE))
      validate(need(user_input(), message = FALSE))
      validate(need(input$testData, message=FALSE))
      
      my_react_vals$d_to_metadata <- user_input()$my_list
      
      if(input$testData==2){
          
          output$readData <- renderText({
            validate(need(input$typeData, message=FALSE))
            
            #my_datatest <<- user_input()$my_list ############## <<- #############################
            
            cond_satisf = c()
            for(x in user_input()$my_list){
              my_condition = is.null(base::rownames(x))
              cond_satisf <- c(cond_satisf, my_condition)
            }
            if(any(cond_satisf)==TRUE){
              paste0("Genes names are absent for counts data in subdirectory ", which(cond_satisf==TRUE), ". Please upload a file with gene names as row names if a count matrix or check your data in any other case")
            }else{
              user_input()$mes
            }
            
          })
      }else if(input$testData == 1){
        output$readData <- renderText({
          user_input()$mes
        })
        }
      #show("readData")
      }, ignoreNULL = FALSE)
    #})
  ##### Preprocessing #####
    
    #add metadata and update counts, sometimes we have to subset the count.matrix because we don't have matches between gene names and ensembl IDs
    calc_metadata <- reactive({
      #validate(need(try(user_input()), message = FALSE))
      validate(need(try(my_react_vals$d_to_metadata), message = "Please read data first"))
      if(input$testData == 1){#test dataset
        if(length(input$typeOrganism) == 0){
          #output$paramsAddMetadata <- renderText(
          my_metadata_list <- list(mes = "You have to select the origin organism")
          #)
        }else{
          #d_to_metadata <- user_input()$my_list
          my_metadata_list <- c(my_list = list(createMetadata(list_samples= my_react_vals$d_to_metadata,"Mus musculus", FALSE)), ############## <<- #############################
                                 mes = "Metadata added")
          
        }

      }else if(input$typeData == 1 | input$typeData == 3){#count.matrix or rds
        if(length(input$typeOrganism) == 0){
          #output$paramsAddMetadata <- renderText(
          my_metadata_list <- list(mes = "You have to select the origin organism")
          #)
        }else{
          #d_to_metadata <- user_input()$my_list
          my_metadata_list <- c(my_list = list(createMetadata(list_samples = my_react_vals$d_to_metadata,input$typeOrganism, FALSE)),
                        mes = "Metadata added")
          
        }
      }else if(input$typeData == 2){#10X
        if(length(input$typeOrganism) == 0){
          
          my_metadata_list <- list(mes = "You have to select the origin organism")
        }else{
          #d_to_metadata <- user_input()$my_list
          my_metadata_list <- c(my_list = list(createMetadata(list_samples = my_react_vals$d_to_metadata,input$typeOrganism, FALSE)),
            mes = "Metadata added")
          
        }
         
        }else if(input$typeData == 4){#multiple samples
          
          if(length(input$typeOrganism) == 0){
            
            my_metadata_list <- list(mes = "You have to select the origin organism")
          }else{
            #d_to_metadata <- user_input()$my_list
            my_metadata_list <- c(my_list = list(createMetadata(list_samples = my_react_vals$d_to_metadata,input$typeOrganism, TRUE)),
              mes = "Metadata added")
          }
        }
        
    })%>% bindEvent(input$addMetadata , ignoreInit = TRUE, ignoreNULL = TRUE)
    
    observeEvent(input$addMetadata, {
      
      my_react_vals$d_to_filter <- calc_metadata()$my_list
      
      output$paramsAddMetadata <- renderText({
        
          calc_metadata()$mes

      })
    }, ignoreNULL = FALSE)
    
    #minimal filtering
    calc_se <- reactive({

      validate(need(try(my_react_vals$d_to_filter), message = "Please add metadata first")) ########################
      #validate(need(calc_metadata()$my_list, message = FALSE))
      #metaData <<- calc_metadata()$my_list
      if(input$testData == 1){#test dataset
        #d_to_filter <- calc_metadata()$my_list ########################
        my_filter_list <- c(my_list = list(initial_filt(list_samples = my_react_vals$d_to_filter, FALSE)), ############## <<- #############################
                            mes = "Minimal filtering completed")
        
        
      }else if(input$testData == 2){
        if(input$typeData == 1 | input$typeData == 3){#count.matrix or rds
          if(length(input$typeOrganism) == 0){
            #output$paramsAddMetadata <- renderText(
            my_filter_list <- list(mes = "You have to select the origin organism") ############## <<- #############################
            #)
          }else{
            #d_to_filter <- calc_metadata()$my_list
            my_filter_list <- c(my_list = list(initial_filt(list_samples = my_react_vals$d_to_filter, FALSE)),
                                mes = "Minimal filtering completed")
            
          }
        }else if(input$typeData == 2){#10X
          if(length(input$typeOrganism) == 0){
            
            my_filter_list <- list(mes = "You have to select the origin organism")
          }else{
            #d_to_filter <- calc_metadata()$my_list
            my_filter_list <- c(my_list = list(initial_filt(list_samples = my_react_vals$d_to_filter, FALSE)),
                                mes = "Minimal filtering completed")
            
          }
          
        }else if(input$typeData == 4){#multiple samples
          
          if(length(input$typeOrganism) == 0){
            
            my_filter_list <- list(mes = "You have to select the origin organism")
          }else{
            #d_to_filter <- calc_metadata()$my_list
            my_filter_list <- c(my_list = list(initial_filt(list_samples = my_react_vals$d_to_filter, TRUE)),
                                mes = "Minimal filtering completed")
          }
        }
      }
    #if(input$typeData == 4){
    #    d_to_filter <- calc_metadata()$my_list
     #   my_filter_list <- c(my_list = list(initial_filt(list_samples = d_to_filter, TRUE)),
    #                           mes = "Minimal filtering completed")
    #    
    #  }else{
    #    d_to_filter <- calc_metadata()$my_list
    #    my_filter_list <- c(my_list = list(initial_filt(list_samples = d_to_filter, FALSE)),
    #                         mes = "Minimal filtering completed")
     # }

      })%>% bindEvent(input$initFiltering, ignoreInit = TRUE, ignoreNULL = TRUE)
    
    observeEvent(input$initFiltering, {
      #my_calc_se <<- calc_se() ############## <<- #############################
      output$initDataFilt<- renderText({

        calc_se()$mes
      
      }) 
    }, ignoreNULL = FALSE)
    
    #save and download as rds
    output$initDownload <- 
      downloadHandler(
        filename = function() {
          paste0("init_filt_data_",  Sys.Date(), ".rds")
        },
        content = function(file) {
          saveRDS(calc_se()$my_list, file)
        }
      )
  
    
  ###### QC routine #####
  
  #update selectInput for samples 
  observeEvent(input$initFiltering,{
    
    my_choices <-  names(calc_se()$my_list)
    updateSelectInput(session, "selectSample",
                      label = "Choose a sample",
                      choices =  my_choices,
                      selected = NULL)
  }, ignoreInit = TRUE, ignoreNULL = FALSE)#%>% bindEvent(input$showSamples , ignoreInit = TRUE, ignoreNULL = TRUE)
    
    #make selectSample input reactive
    my_sample_name <- reactive({
      req(input$selectSample)
      my_name <- input$selectSample
    })
  

    #create metrics
    calc_metrics <- reactive({
      validate(need(calc_se()$my_list, message="Please complete the initial filtering before proceeding")) 
      validate(need(input$selectSample, message=FALSE))
     #my_se<<- calc_se()$my_list ############################ ATTENTION #########################
      
      my_metrics <- c(as.data.frame(colData(calc_se()$my_list[[my_sample_name()]]$se)), sample = my_sample_name(),
        mes = paste0("Metrics created for ", my_sample_name()))
      
    })%>% bindEvent(input$makeMetrics , ignoreInit = TRUE, ignoreNULL = TRUE)
    
    output$createdMetrics <- renderText({
      validate(need(calc_se()$my_list, message="Please complete the initial filtering before proceeding")) 
        calc_metrics()$mes
      }) 

    output$dimsMetrics <- renderText({
      validate(need(calc_metrics(), message= FALSE))
      paste0("Dimensions of ", my_sample_name(), ":",
             "\nFeatures: ", dim(calc_se()$my_list[[my_sample_name()]]$se)[1],
             "\nCells: ", dim(calc_se()$my_list[[my_sample_name()]]$se)[2])
    })
    
  #Update slider min & max values when metrics object is created
    observeEvent(input$drawPlots,{
      validate(need(input$makeMetrics, message=FALSE))
      validate(need(calc_metrics(), message= FALSE))
      validate(need(input$selectSample, message="You have to pick a sample to work with"))
      #update slider min and max values of UMIs per cell
      updateSliderInput(session,"nUMIs", "Number of UMIs:",
                        min = 0, max = max(calc_metrics()$nUMI), value = c(min(calc_metrics()$nUMI),max(calc_metrics()$nUMI)), step = 10)
      #nGenes histogram
      updateSliderInput(session, "nGenes", "Number of genes:",
                        min = 0, max = max(calc_metrics()$nGene), value = c(min(calc_metrics()$nGene),max(calc_metrics()$nGene)), step = 10)
      #UMIs vs nGenes with mitoRatio
      updateSliderInput(session,"nUMIsVSnGenes", "Genes detected:",
                        min = 0, max = max(calc_metrics()$nGene), value = c(min(calc_metrics()$nGene),max(calc_metrics()$nGene)), step = 10)
      #UMIs vs nGenes with hemoRatio
      updateSliderInput(session,"nUMIsVSnGenesRBC", "Genes detected:",
                        min = 0, max = max(calc_metrics()$nGene), value = c(min(calc_metrics()$nGene),max(calc_metrics()$nGene)), step = 10)
      #mitoRatio
      updateSliderInput(session,"mitoRatio", "Mitochondrial ratio:",
                        min = 0, max = round(max(calc_metrics()$mitoRatio), digits = 5), value = c(round(min(calc_metrics()$mitoRatio),digits = 5),round(max(calc_metrics()$mitoRatio), digits = 5)), step = 0.0001)
      #hemoRatio
      updateSliderInput(session,"rbcRatio", "Hemoglobins ratio:",
                        min = 0, max = round(max(calc_metrics()$rbcRatio), digits = 5), value = c(round(min(calc_metrics()$rbcRatio),digits = 5),round(max(calc_metrics()$rbcRatio), digits = 5)), step = 0.0001)
      #Genes vs UMIs (novelty)
      updateSliderInput(session,"novelty", "Genes detected per UMI:",
                        min = 0, max = round(max(calc_metrics()$log10GenesPerUMI), digits = 2), value = round(min(calc_metrics()$log10GenesPerUMI),digits = 2), step = 0.01)
      
    })
    
  #number of UMIs per cell
  #update slider min and max values
   #observe({
  #  updateSliderInput(session,"nUMIs", "Number of UMIs:",
   #                   min = 0, max = max(calc_metrics()$nUMI), value = c(min(calc_metrics()$nUMI),max(calc_metrics()$nUMI)), step = 10)
  #}) %>% bindEvent(input$makeNBUMIs , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  
  output$UMIsPlot <- renderPlot({
      validate(need(calc_metrics(), message= "You have to compute metrics for your sample first"))

      createUMIsPlot(calc_metrics(), input$nUMIs[1], input$nUMIs[2])
      
    })%>% bindEvent(c(input$drawPlots, input$nUMIs) , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #})%>% bindEvent(c(input$makeNBUMIs, input$nUMIs) , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #nGenes histogram
  #observe({
    #updateSliderInput(session, "nGenes", "Number of genes:",
    #                  min = 0, max = max(calc_metrics()$nGene), value = c(min(calc_metrics()$nGene),max(calc_metrics()$nGene)), step = 10)
  #}) %>% bindEvent(input$makeNBgenesHist , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  output$GenesHistPlot<- renderPlot({
    validate(need(calc_metrics(), message= "You have to compute metrics for your sample first"))
    createGenesHistPlot(calc_metrics(), input$nGenes[1], input$nGenes[2])
  })%>% bindEvent(c(input$drawPlots, input$nGenes) , ignoreInit = TRUE, ignoreNULL = TRUE)
  #})%>% bindEvent(c(input$makeNBgenesHist, input$nGenes) , ignoreInit = TRUE, ignoreNULL = TRUE)

  #nGenes boxplot
  output$GenesBoxPlot <- renderPlot({
    validate(need(calc_metrics(), message= "You have to compute metrics for your sample first"))
    createGenesBoxPlot(calc_metrics())
  })%>% bindEvent(c(input$drawPlots) , ignoreInit = TRUE, ignoreNULL = TRUE)  
  #})%>% bindEvent(c(input$makeNBgenesBoxPlot) , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #UMIs vs nGenes with mitoRatio
  #observe({
    #updateSliderInput(session,"nUMIsVSnGenes", "Genes detected:",
   #                   min = 0, max = max(calc_metrics()$nGene), value = c(min(calc_metrics()$nGene),max(calc_metrics()$nGene)), step = 10)
  #})%>% bindEvent(input$makeUMIsGenes , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  output$UMIsGenesPlot <- renderPlot({
    validate(need(calc_metrics(), message= "You have to compute metrics for your sample first"))
    createUMIsGenesPlot(calc_metrics(), input$nUMIsVSnGenes[1], input$nUMIsVSnGenes[2])
  })%>% bindEvent(c(input$drawPlots, input$nUMIsVSnGenes) , ignoreInit = TRUE, ignoreNULL = TRUE)  
  #})%>% bindEvent(c(input$makeUMIsGenes, input$nUMIsVSnGenes) , ignoreInit = TRUE, ignoreNULL = TRUE)
 
  
  #UMIs vs nGenes with hemoRatio
  #observe({
  #  updateSliderInput(session,"nUMIsVSnGenesRBC", "Genes detected:",
   #                   min = 0, max = max(calc_metrics()$nGene), value = c(min(calc_metrics()$nGene),max(calc_metrics()$nGene)), step = 10)
  #})%>% bindEvent(input$makeUMIsGenesRBC , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  output$UMIsGenesPlotRBC <- renderPlot({
    validate(need(calc_metrics(), message= "You have to compute metrics for your sample first"))
    createUMIsGenesPlotRBC(calc_metrics(), input$nUMIsVSnGenesRBC[1], input$nUMIsVSnGenesRBC[2])
  })%>% bindEvent(c(input$drawPlots, input$nUMIsVSnGenesRBC) , ignoreInit = TRUE, ignoreNULL = TRUE)    
  #})%>% bindEvent(c(input$makeUMIsGenesRBC, input$nUMIsVSnGenesRBC) , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #plot mitoRatio
  #observe({
  #  updateSliderInput(session,"mitoRatio", "Mitochondrial ratio:",
   #                   min = 0, max = round(max(calc_metrics()$mitoRatio), digits = 5), value = c(round(min(calc_metrics()$mitoRatio),digits = 5),round(max(calc_metrics()$mitoRatio), digits = 5)), step = 0.0001)
 # })%>% bindEvent(input$makemitoRatio , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  output$mitoRatioPlot <- renderPlot({
    validate(need(calc_metrics(), message= "You have to compute metrics for your sample first"))
    createmitoRatioPlot(calc_metrics(), input$mitoRatio[1], input$mitoRatio[2])
  })%>% bindEvent(c(input$drawPlots, input$mitoRatio) , ignoreInit = TRUE, ignoreNULL = TRUE)  
  #})%>% bindEvent(c(input$makemitoRatio, input$mitoRatio) , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #plot hemoRatio
  #observe({
   # updateSliderInput(session,"rbcRatio", "Hemoglobins ratio:",
  #                    min = 0, max = round(max(calc_metrics()$rbcRatio), digits = 5), value = c(round(min(calc_metrics()$rbcRatio),digits = 5),round(max(calc_metrics()$rbcRatio), digits = 5)), step = 0.0001)
  #})%>% bindEvent(input$makerbcRatio , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  output$rbcRatioPlot <- renderPlot({
    validate(need(calc_metrics(), message= "You have to compute metrics for your sample first"))
    createrbcRatioPlot(calc_metrics(), input$rbcRatio[1], input$rbcRatio[2])
  })%>% bindEvent(c(input$drawPlots, input$rbcRatio) , ignoreInit = TRUE, ignoreNULL = TRUE)  
  #})%>% bindEvent(c(input$makerbcRatio, input$rbcRatio) , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #Genes vs UMIs
  #observe({
   # updateSliderInput(session,"novelty", "Genes detected per UMI:",
   #                   min = 0, max = round(max(calc_metrics()$log10GenesPerUMI), digits = 2), value = round(min(calc_metrics()$log10GenesPerUMI),digits = 2), step = 0.01)
  #})%>% bindEvent(input$makeNovelty , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  output$noveltyPlot <- renderPlot({
    validate(need(calc_metrics(), message= "You have to compute metrics for your sample first"))
    createnoveltyPlot(calc_metrics(), input$novelty[1], input$novelty[2])
  })%>% bindEvent(input$drawPlots, input$novelty, ignoreInit = TRUE, ignoreNULL = TRUE)  
  #})%>% bindEvent(input$makeNovelty, input$novelty, ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #plot scaterPlot
  scater_plot <- reactive({
    validate(need(calc_se()$my_list, message= FALSE))
    createlibraryPlot(calc_se()$my_list[[my_sample_name()]]$se)
  #})%>%bindEvent(input$drawPlots, ignoreInit = TRUE, ignoreNULL = TRUE)  
  })%>%bindEvent(input$makeLibrary, ignoreInit = TRUE, ignoreNULL = TRUE)
  
  
  output$libraryPlot <- renderPlot({
    validate(need(calc_se()$my_list, message= "You have to complete initial filtering first"))
    scater_plot()
  })
  
  #plot highest expressed features
  highest_expr_plot <- reactive({
    validate(need(calc_se()$my_list, message= FALSE))
    createhighExprPlot(calc_se()$my_list[[my_sample_name()]]$se)
  #})%>%bindEvent(input$drawPlots, ignoreInit = TRUE, ignoreNULL = TRUE)  
  })%>%bindEvent(input$makeHighExpr, ignoreInit = TRUE, ignoreNULL = TRUE)
  
  output$highExprPlot <- renderPlot({
    validate(need(calc_se()$my_list, message= "You have to complete initial filtering first"))
    highest_expr_plot()
  })
  
  #zoom on plot
  observeEvent(input$zoomhighExprPlot, {
    validate(need(calc_se()$my_list, message= FALSE))
    showModal(modalDialog(
      renderPlot({
        highest_expr_plot()
        },
        width = 600,
        height = 1200,
        res = 300),
      footer = NULL,
      size = "l",
      easyClose = TRUE
    ))
  })
  
  ###### Filtering thresholds #####
  
  #update filtering thresholds when metrics are computed to avoid wrong initial values as input to plots
  observe({
    validate(need(calc_se(), message=FALSE))
    validate(need(calc_metrics(), message=FALSE))
    
    updateNumericInput(session, "minUMIsPerCell", value = min(calc_metrics()$nUMI))
    updateNumericInput(session, "maxUMIsPerCell", value = max(calc_metrics()$nUMI))
    updateNumericInput(session, "minGenesPerCell", value = min(calc_metrics()$nGene))
    updateNumericInput(session, "maxGenesPerCell", value = max(calc_metrics()$nGene))
    updateNumericInput(session, "minMitoPerCell", value = round(min(calc_metrics()$mitoRatio),digits = 5))
    updateNumericInput(session, "maxMitoPerCell", value = round(max(calc_metrics()$mitoRatio),digits = 5))
    updateNumericInput(session, "minRbcPerCell", value = round(min(calc_metrics()$rbcRatio),digits = 5))
    updateNumericInput(session, "maxRbcPerCell", value = round(max(calc_metrics()$rbcRatio),digits = 5))
    updateNumericInput(session, "noveltyPerCell", value = round(min(calc_metrics()$log10GenesPerUMI),digits = 5))
    
  }) %>% bindEvent(input$makeMetrics, ignoreInit = TRUE, ignoreNULL = TRUE)
  
  
  #update filtering thresholds from sliders upon user request via action button "updateThresholds"
  observe({
    updateNumericInput(session, "minUMIsPerCell", value = input$nUMIs[1])
    updateNumericInput(session, "maxUMIsPerCell", value = input$nUMIs[2])
    updateNumericInput(session, "minGenesPerCell", value = input$nGenes[1])
    updateNumericInput(session, "maxGenesPerCell", value = input$nGenes[2])
    updateNumericInput(session, "minMitoPerCell", value = input$mitoRatio[1])
    updateNumericInput(session, "maxMitoPerCell", value = input$mitoRatio[2])
    updateNumericInput(session, "minRbcPerCell", value = input$rbcRatio[1])
    updateNumericInput(session, "maxRbcPerCell", value = input$rbcRatio[2])
    updateNumericInput(session, "noveltyPerCell", value = input$novelty)
    
  }) %>% bindEvent(input$updateThresholds, ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #User defined thresholds from numericInput boxes
  output$minUMIsPerCellThreshold <- renderText({ paste0("We'll keep cells with number of UMIs detected per cell above: ", input$minUMIsPerCell) })
  output$maxUMIsPerCellThreshold <- renderText({ paste0("We'll keep cells with number of UMIs detected per cell below: ", input$maxUMIsPerCell) })
  output$minGenesPerCellThreshold <- renderText({ paste0("We'll keep cells with number of genes detected per cell above: ", input$minGenesPerCell) })
  output$maxGenesPerCellThreshold <- renderText({ paste0("We'll keep cells with number of genes detected per cell below: ", input$maxGenesPerCell) })
  output$minMitoPerCellThreshold <- renderText({ paste0("We'll keep cells with mitochondrial ratio above: ", input$minMitoPerCell) })
  output$maxMitoPerCellThreshold <- renderText({ paste0("We'll keep cells with mitochondrial ratio below: ", input$maxMitoPerCell) })
  output$minRbcPerCellThreshold <- renderText({ paste0("We'll keep cells with hemoglobins ratio above: ", input$minRbcPerCell) })
  output$maxRbcPerCellThreshold <- renderText({ paste0("We'll keep cells with hemoglobins ratio below: ", input$maxRbcPerCell) })
  output$noveltyPerCellThreshold <- renderText({ paste0("We'll keep cells with number of genes detected per UMI per cell (novelty) above: ", input$noveltyPerCell) })
  
  #QC Filtering
  se_c <- reactive({
    if(is.null(calc_metrics())){
      validate(FALSE, message=paste0("Doesn't seem you've created the metrics for this sample yet (input= ", my_sample_name(), ")"))
    }
    validate(need(calc_metrics(), message=FALSE))
    metrics <- as.data.frame(calc_metrics())
    raw_data <- calc_se()$my_list[[my_sample_name()]]$se
    #filter_cells(raw_data, metrics, input$minUMIsPerCell, input$maxUMIsPerCell, input$minGenesPerCell, input$maxGenesPerCell, 
     # input$noveltyPerCell, input$minMitoPerCell, input$maxMitoPerCell, input$minRbcPerCell, input$maxRbcPerCell)
    tryCatch( 
      {
        filter_cells(raw_data, metrics, input$minUMIsPerCell, input$maxUMIsPerCell, input$minGenesPerCell, input$maxGenesPerCell, 
                     input$noveltyPerCell, input$minMitoPerCell, input$maxMitoPerCell, input$minRbcPerCell, input$maxRbcPerCell)
      },
      error = function(e) {
        validate(FALSE,paste0("Doesn't seem you've created the metrics for this sample yet (input= ", my_sample_name(), ")"))
      }
    )

  })%>% bindEvent(input$QCFiltering , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  output$QCDataFilt <- reactive({
      se_c()$mes
  })
  
  #create a list to keep filtered data for each sample
  #create the names of the list
  filt_samples_list_names <- reactive({
    validate(need(calc_se()$my_list, message=FALSE))
    names(calc_se()$my_list)
  })
  #reactive values to store the samples list
  #my_react_vals <- reactiveValues()#(filt_samples_list = NULL, loaded.dataSO.combined.no.cluster = NULL, loaded.dataSO.combined.markers = NULL, loaded.dataSO.combined.markerstop = NULL)
  
  #initialize a samples' list for multiple samples
  observeEvent(input$initFiltering,{
    validate(need(input$typeData, message = FALSE))
    if(input$typeData==4){
      my_react_vals$filt_samples_list <- vector(mode="list", length=length(filt_samples_list_names()))
      names(my_react_vals$filt_samples_list) <- filt_samples_list_names()
    }#else{
      #my_react_vals$filt_samples_list <- list(se_c()$se)
      #names(my_react_vals$filt_samples_list)[1] <- "sample"
    #}
  })
   # initialize a samples' list for testData and update it with filtered testData
  observeEvent(input$QCFiltering,{
    validate(need(input$testData, message = FALSE))
    if(input$testData==1){
    my_react_vals$filt_samples_list <- list(se_c()$se)
    names(my_react_vals$filt_samples_list)[1] <- "sample"
    }
  })
  
  #update the filtered samples list
  observeEvent(input$QCFiltering,{
    validate(need(se_c(), message=FALSE))
    validate(need(calc_se()$my_list, message=FALSE))
    validate(need(input$typeData, message = FALSE))
    if(input$typeData==4){
      #my_react_vals$filt_samples_list <- vector(mode="list", length=length(filt_samples_list_names()))
      my_react_vals$filt_samples_list[[my_sample_name()]] <- se_c()$se
    }else{
      my_react_vals$filt_samples_list <- list(se_c()$se)
      names(my_react_vals$filt_samples_list)[1] <- "sample"
    }
  })

  # Save subset to new metrics variable
  metrics_clean <- reactive({
    validate(need(input$selectSample, message=FALSE))
    validate(need(input$initFiltering, message="Please complete the initial filtering before proceeding")) 
    validate(need(se_c()$se, message = FALSE))######################## ATTENTION ################################################
    #test_se_c <<- se_c()$se ######################## ATTENTION ################################################
    #my_filt_samples_list <<- my_react_vals$filt_samples_list
    my_metrics_clean <- c(as.data.frame(colData(se_c()$se)), sample = my_sample_name(),
    mes = paste0("Metrics recomputed for ", my_sample_name()))
    
  })%>% bindEvent(input$filt_metrics, ignoreInit = TRUE, ignoreNULL = TRUE)

   
  output$metricsFilt<- renderText({

    metrics_clean()$mes
 
  })
  
  output$dimsMetricsClean <- renderText({
    validate(need(metrics_clean(), message= FALSE))
    paste0("Dimensions of filtered ", my_sample_name(), ":",
           "\nFeatures: ", dim(my_react_vals$filt_samples_list[[my_sample_name()]])[1],
           "\nCells: ", dim(my_react_vals$filt_samples_list[[my_sample_name()]])[2])
  })
  #save and download as rds
  output$QCDownload <- 
    #validate(need(se_c()$se, message = "You have to filter out cells first"))
    downloadHandler(
      filename = function() {
        paste0("filtered_se_", my_sample_name(), "_", Sys.Date(), ".rds")
      },
      content = function(file) {
       # se= calc_se()$se
       # keep_index = keep_index()$keep_index
        saveRDS(se_c()$se, file)
      }
    )

  
  ##### Redraw Plots with filtered data #####
  #UMIs Plot
  output$UMIsPlotFiltered <- renderPlot({
    validate(need(metrics_clean(), message = "You have to compute metrics for filtered cells first"))
    
    createUMIsPlot(metrics_clean(), min(metrics_clean()$nUMI),max(metrics_clean()$nUMI))
    
  })%>% bindEvent(input$redrawPlots , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #nGenes histogram
  output$GenesHistPlotFiltered <- renderPlot({
    validate(need(metrics_clean(), message = "You have to compute metrics for filtered cells first"))
    
    createGenesHistPlot(metrics_clean(), min(metrics_clean()$nGene),max(metrics_clean()$nGene))
    
  })%>% bindEvent(input$redrawPlots , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #nGenes Boxplot
  output$GenesBoxPlotFiltered <- renderPlot({
    validate(need(metrics_clean(), message = "You have to compute metrics for filtered cells first"))
    
    createGenesBoxPlot(metrics_clean())
    
  })%>% bindEvent(input$redrawPlots , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #UMIs vs genes with mitoRatio
  output$UMIsGenesPlotFiltered <- renderPlot({
    validate(need(metrics_clean(), message = "You have to compute metrics for filtered cells first"))
    
    createUMIsGenesPlot(metrics_clean(), min(metrics_clean()$nGene), max(metrics_clean()$nGene))
    
  })%>% bindEvent(input$redrawPlots , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #UMIs vs genes with rbcRatio
  output$UMIsGenesPlotRBCFiltered <- renderPlot({
    validate(need(metrics_clean(), message = "You have to compute metrics for filtered cells first"))
    
    createUMIsGenesPlotRBC(metrics_clean(), min(metrics_clean()$nGene), max(metrics_clean()$nGene))
    
  })%>% bindEvent(input$redrawPlots , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #mitoRatio
  output$mitoRatioPlotFiltered <- renderPlot({
    validate(need(metrics_clean(), message = "You have to compute metrics for filtered cells first"))
    
    createmitoRatioPlot(metrics_clean(), min(metrics_clean()$mitoRatio), max(metrics_clean()$mitoRatio))
    
  })%>% bindEvent(input$redrawPlots , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #rbcRatio
  output$rbcRatioPlotFiltered <- renderPlot({
    validate(need(metrics_clean(), message = "You have to compute metrics for filtered cells first"))
    
    createrbcRatioPlot(metrics_clean(), min(metrics_clean()$rbcRatio), max(metrics_clean()$rbcRatio))
    
  })%>% bindEvent(input$redrawPlots , ignoreInit = TRUE, ignoreNULL = TRUE)
  
   
  #novelty
  output$noveltyPlotFiltered <- renderPlot({
    validate(need(metrics_clean(), message = "You have to compute metrics for filtered cells first"))
    
    createnoveltyPlot(metrics_clean(), min(metrics_clean()$log10GenesPerUMI), max(metrics_clean()$log10GenesPerUMI))
    
  })%>% bindEvent(input$redrawPlots , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #library plot
  output$libraryPlotFiltered <- renderPlot({
    validate(need(se_c()$se, message = "You have to filter out cells first"))
    createlibraryPlot(se_c()$se)
    
  })%>% bindEvent(input$makeLibraryFiltered , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #plot highest expressed features
  highest_expr_plot_filtered <- reactive({
    validate(need(se_c()$se, message = FALSE))
    createhighExprPlot(se_c()$se)
  })%>%bindEvent(input$makeHighExprFiltered , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  output$highExprPlotFiltered <- renderPlot({
    validate(need(se_c()$se, message = "You have to filter out cells first"))
    highest_expr_plot_filtered()
  })
  
  #zoom on plot
  observeEvent(input$zoomhighExprPlotFiltered, {
    #validate(need(se_c()$se, message = "You have to filter out cells first"))
    showModal(modalDialog(
      renderPlot({
        highest_expr_plot_filtered()
      },
      width = 600,
      height = 1200,
      res = 300),
      footer = NULL,
      size = "l",
      easyClose = TRUE
    ))
  })
  
  ##### Cell clustering #####
  
  #I have to initialize a variable that corresponds to input$typeData == 4, else I have a problem when working with testData, because input$typeData is NULL in this case
  observeEvent(input$startQC,{
    validate(need(input$testData, message = FALSE))
    if(input$testData == 1){
      my_react_vals$multiple_case <- FALSE
    }else if(input$testData == 2){
      validate(need(input$typeData, message = FALSE))
      if(input$typeData== 4){
        my_react_vals$multiple_case <- TRUE
      }else{
        my_react_vals$multiple_case <- FALSE
      }
    }
  }, ignoreNULL = TRUE)
  
  ##Normalize data and store the clustering Seurat object as a reactiveVal object so I can update it at the next steps, avoiding the creation of new objects
  observeEvent(input$normData,{
    validate(need(se_c(), message=FALSE))
    validate(need(calc_se()$my_list, message=FALSE))
    #validate(need(my_react_vals$filt_samples_list, message = FALSE))
    #validate(need(any(mapply(is.null, my_react_vals$filt_samples_list )) == FALSE, message = "It seems you forgot to filter out cells in some/all of your samples"))
    #filt_samples_list_pre_norm <<- my_react_vals$filt_samples_list ############## <<- #############################
    #if(input$typeData==4){
    if(my_react_vals$multiple_case){
      validate(need(my_react_vals$filt_samples_list, message = FALSE))
      validate(need(any(mapply(is.null, my_react_vals$filt_samples_list )) == FALSE, message = "It seems you forgot to filter out cells in some/all of your samples"))
      
      spsComps::shinyCatch({
        
        my_react_vals$loaded.dataSO.combined.no.cluster <- norm_data_and_pca_multiple(my_react_vals$filt_samples_list)
        
      },
      # blocking recommended
      blocking_level = "error",
      prefix = "My-project" #change console prefix if you don't want "SPS"
      )
      
    }else{
      validate(need(my_react_vals$filt_samples_list, message = "Please filter out cells first before proceeding to normalization"))
      validate(need(any(mapply(is.null, my_react_vals$filt_samples_list )) == FALSE, message = "It seems you forgot to filter out cells in your sample"))
      
      spsComps::shinyCatch({
        my_react_vals$loaded.dataSO.combined.no.cluster <- norm_data_and_pca(my_react_vals$filt_samples_list)
      },
      # blocking recommended
      blocking_level = "error",
      prefix = "My-project" #change console prefix if you don't want "SPS"
      )
      
    }
  })
  
  #Normalize data
 # norm_pca_data <- reactive({
  #  if(input$typeData == 4){
  #    validate(need(my_react_vals$filt_samples_list, message = FALSE))
  #    validate(need(any(mapply(is.null, my_react_vals$filt_samples_list )) == FALSE, message = "It seems you forgot to filter out cells in some/all of your samples"))

   #   norm_data_and_pca_multiple(my_react_vals$filt_samples_list)
  #  }else{
   #   validate(need(my_react_vals$filt_samples_list, message = "Please filter out cells first before proceeding to normalization"))
   #   validate(need(any(mapply(is.null, my_react_vals$filt_samples_list )) == FALSE, message = "It seems you forgot to filter out cells in your sample"))
  #    norm_data_and_pca(my_react_vals$filt_samples_list)
  #  }
 #}) %>% bindEvent(input$normData , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  output$normalizeData <- renderText({
    #my_filt_samples_list <<- my_react_vals$filt_samples_list
    #loaded.dataSO.combined.no.cluster<<- my_react_vals$loaded.dataSO.combined.no.cluster
    validate(need(my_react_vals$filt_samples_list, message = "Please filter out cells first before proceeding to normalization"))
    validate(need(any(mapply(is.null, my_react_vals$filt_samples_list )) == FALSE, message = "It seems you forgot to filter out cells in some/all of your samples"))
    if(!is.null(my_react_vals$loaded.dataSO.combined.no.cluster )){
      "Normalization and PCA completed successfully"
    }else if(is.null(my_react_vals$loaded.dataSO.combined.no.cluster )){
      #null_samples <- which(is.null(my_react_vals$loaded.dataSO.combined.no.cluster ))
      #validate(FALSE, paste0("You have to complete the filtering step first"))
      validate(FALSE,"You have to complete the filtering step first")
    }
  })%>% bindEvent(input$normData , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #elbow plot
  output$elbowPlot <- renderPlot({ 
    validate(need(my_react_vals$loaded.dataSO.combined.no.cluster, message = FALSE))
    tryCatch({
      ElbowPlot(my_react_vals$loaded.dataSO.combined.no.cluster, ndims = 30, reduction = "pca")
    },
    error = function(e) {
      validate(FALSE,"Please run normalization step first")
    })
  })
  
  #Find Neighbors
  observeEvent(input$findNeigh,{
    validate(need(my_react_vals$loaded.dataSO.combined.no.cluster, message=FALSE ))
    req(input$nPCs)
    #if(input$typeData==4){
    if(my_react_vals$multiple_case){
      validate(need(my_react_vals$filt_samples_list, message = "Please filter out cells first before proceeding to cell type clustering"))
      validate(need(my_react_vals$loaded.dataSO.combined.no.cluster@reductions$pca, message = "It seems you didn't normalize data first"))

      obj <- IntegrateLayers(
        object = my_react_vals$loaded.dataSO.combined.no.cluster, method = CCAIntegration,
        orig.reduction = "pca", new.reduction = "integrated.cca",
        verbose = FALSE)

      spsComps::shinyCatch({

        my_react_vals$loaded.dataSO.combined.no.cluster <- FindNeighbors(obj, reduction = "integrated.cca", dims = 1:input$nPCs)

      },
      # blocking recommended
      blocking_level = "error",
      prefix = "My-project" #change console prefix if you don't want "SPS"
      )

      #my_merged_seurat_clust <<- my_react_vals$loaded.dataSO.combined.no.cluster
      
    }else{
      validate(need(my_react_vals$filt_samples_list, message = "Please filter out cells first before proceeding to cell type clustering"))
      validate(need(my_react_vals$loaded.dataSO.combined.no.cluster@reductions$pca, message = "It seems you didn't normalize data first"))
      
      spsComps::shinyCatch({
        
        my_react_vals$loaded.dataSO.combined.no.cluster <- FindNeighbors(my_react_vals$loaded.dataSO.combined.no.cluster, dims = 1:input$nPCs, reduction = "pca")
        
        },
      # blocking recommended
      blocking_level = "error",
      prefix = "My-project" #change console prefix if you don't want "SPS"
      )
      
      
    }
  })
  
  #find_neigh_data <- reactive({
  #  
  #  validate(need(my_react_vals$loaded.dataSO.combined.no.cluster, message="Please normalize your data first" ))
  #  req(input$nPCs)
  #  
  #  if(input$typeData==4){
  #    
  #    obj <- IntegrateLayers(
  #      object = my_react_vals$loaded.dataSO.combined.no.cluster, method = CCAIntegration,
  #      orig.reduction = "pca", new.reduction = "integrated.cca",
  #      verbose = FALSE
  #    )
  #    FindNeighbors(obj, reduction = "integrated.cca", dims = 1:input$nPCs)
  #    
  #  }else{
  #    FindNeighbors(my_react_vals$loaded.dataSO.combined.no.cluster, dims = 1:input$nPCs, reduction = "pca")
  #  }
    
  #})%>% bindEvent(input$findNeigh , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  output$findNeighRes <- renderText({
    #loaded.dataSO.combined.no.cluster<<- my_react_vals$loaded.dataSO.combined.no.cluster
    validate(need(my_react_vals$loaded.dataSO.combined.no.cluster, message = "You have to complete the normalization and reduction step first"))
    validate(need(grep("RNA_nn", names(my_react_vals$loaded.dataSO.combined.no.cluster))!=0 | grep("RNA_snn", names(my_react_vals$loaded.dataSO.combined.no.cluster))!=0, message = FALSE))
 
    if(is.null(my_react_vals$loaded.dataSO.combined.no.cluster)){
      validate(FALSE,"You have to complete the normalization and reduction step first")
    }
    if(grep("RNA_nn", names(my_react_vals$loaded.dataSO.combined.no.cluster))!=0 | grep("RNA_snn", names(my_react_vals$loaded.dataSO.combined.no.cluster))!=0){
      "Nearest neighbors computed"
    }else{
      validate(FALSE,"You have to complete the normalization and reduction step first")
    }
  })%>% bindEvent(input$findNeigh , ignoreInit = TRUE, ignoreNULL = TRUE)

  #Find Clusters
  observeEvent(input$findClusters,{
    validate(need(my_react_vals$loaded.dataSO.combined.no.cluster, message=FALSE))
    #validate(need(!is.null(my_react_vals$loaded.dataSO.combined.no.cluster@commands$FindVariableFeatures.RNA), message = "Please compute neighbors first"))
    req(input$resolutionValue)
    #if(input$typeData==4){
    if(my_react_vals$multiple_case){
      validate(need(any(mapply(is.null, my_react_vals$filt_samples_list )) == FALSE, message = "It seems you forgot to filter out cells in some/all of your samples"))
      #validate(need(my_react_vals$filt_samples_list, message = "Please filter out cells first before proceeding to cell type clustering"))
      
      if(length(grep("RNA_nn", names(my_react_vals$loaded.dataSO.combined.no.cluster)))==0 | length(grep("RNA_snn", names(my_react_vals$loaded.dataSO.combined.no.cluster)))==0){
        validate(FALSE, message = "Please compute neighbors first")
      }
      
      spsComps::shinyCatch({

        my_react_vals$loaded.dataSO.combined.no.cluster <- FindClusters(my_react_vals$loaded.dataSO.combined.no.cluster, resolution = input$resolutionValue)#, cluster.name = "cca_clusters")

        },
      # blocking recommended
      blocking_level = "error",
      prefix = "My-project" #change console prefix if you don't want "SPS"
      )
      }else{
      validate(need(any(mapply(is.null, my_react_vals$filt_samples_list )) == FALSE, message = "It seems you forgot to filter out cells in some/all of your samples"))
      #validate(need(my_react_vals$filt_samples_list, message = "Please filter out cells first before proceeding to cell type clustering"))
      validate(need(my_react_vals$loaded.dataSO.combined.no.cluster@reductions$pca, message = "It seems you didn't normalize data first"))
      if(length(grep("RNA_nn", names(my_react_vals$loaded.dataSO.combined.no.cluster)))==0 | length(grep("RNA_snn", names(my_react_vals$loaded.dataSO.combined.no.cluster)))==0){
        validate(FALSE, message = "Please compute neighbors first")}
      
      spsComps::shinyCatch({
        my_react_vals$loaded.dataSO.combined.no.cluster <- FindClusters(my_react_vals$loaded.dataSO.combined.no.cluster, resolution = input$resolutionValue)#, cluster.name = "pca_clusters")
      },
      # blocking recommended
      blocking_level = "error",
      prefix = "My-project" #change console prefix if you don't want "SPS"
      )
      
    }
  })
  
  
 # find_clusters_data <- reactive({
    
  #  validate(need(find_neigh_data(), message="Please find neighbors first" ))
  #  req(input$resolutionValue)
    
  #  FindClusters(find_neigh_data(), resolution  = input$resolutionValue)
    # t-SNE and Clustering

 # })%>% bindEvent(input$findNeigh , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  output$findClustersRes <- renderText({
    validate(need(my_react_vals$loaded.dataSO.combined.no.cluster, message = "You have to complete the normalization and reduction step first"))
    validate(need((length(grep("RNA_nn", names(my_react_vals$loaded.dataSO.combined.no.cluster)))==0 | length(grep("RNA_snn", names(my_react_vals$loaded.dataSO.combined.no.cluster)))==0)== FALSE, 
                  "You have to compute nearest neighbors first"))
    validate(need(my_react_vals$loaded.dataSO.combined.no.cluster@meta.data$seurat_clusters, message = FALSE))
    #my_merged_seurat_clust <<- my_react_vals$loaded.dataSO.combined.no.cluster
    if(!is.null(my_react_vals$loaded.dataSO.combined.no.cluster@meta.data$seurat_clusters)){
      "Clusters identified"
    }else{
      validate(FALSE, "You have to compute nearest neighbors first")
    }
  })%>% bindEvent(input$findClusters , ignoreInit = TRUE, ignoreNULL = TRUE)

  
  observeEvent(input$cellMarkers, {

      validate(need(my_react_vals$loaded.dataSO.combined.no.cluster, message= FALSE))
      if(is.null(my_react_vals$loaded.dataSO.combined.no.cluster@meta.data$seurat_clusters)){
        message("You have to identify clusters at the previous step first")
      }else if(length(grep("RNA_nn", names(my_react_vals$loaded.dataSO.combined.no.cluster)))==0 | length(grep("RNA_snn", names(my_react_vals$loaded.dataSO.combined.no.cluster)))==0){
        message("You have to compute nearest neighbors first")
      }else if(is.null(my_react_vals$loaded.dataSO.combined.no.cluster)){
        "You have to complete the normalization and reduction step first"
      }else{
        #if(input$typeData==4){
        if(my_react_vals$multiple_case){
          my_react_vals$loaded.dataSO.combined.no.cluster <- JoinLayers(my_react_vals$loaded.dataSO.combined.no.cluster)

        }
        ############ 
        spsComps::shinyCatch({

          data_markers <- FindAllMarkers(my_react_vals$loaded.dataSO.combined.no.cluster,  only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

        },
        # blocking recommended
        blocking_level = "error",
        prefix = "My-project" #change console prefix if you don't want "SPS"
        )
          
       my_react_vals$loaded.dataSO.combined.markers <- data_markers
      }

  })
  
  #Define cell type markers
  #####################################I have to JoinLayers in integrated option
  #data.markers <- reactive({
  #  if(input$typeData == 4){
  #    RunUMAP(my_react_vals$loaded.dataSO.combined.no.cluster, reduction = "integrated.cca", dims = 1:30)
  #    FindAllMarkers(my_react_vals$loaded.dataSO.combined.no.cluster, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  #  }else{
      #DefaultAssay(loaded.dataSO.combined.no.cluster) <<- "RNA"
  #    FindAllMarkers(loaded.dataSO.combined.no.cluster, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  #  }
    #DefaultAssay(loaded.dataSO.combined.no.cluster) <<- "RNA"
  #  FindAllMarkers(loaded.dataSO.combined.no.cluster, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  #})
  
  output$cellTypeMarkers <- reactive({
    validate(need(my_react_vals$loaded.dataSO.combined.no.cluster, message = "You have to complete the normalization and reduction step first"))
    validate(need((length(grep("RNA_nn", names(my_react_vals$loaded.dataSO.combined.no.cluster)))==0 | length(grep("RNA_snn", names(my_react_vals$loaded.dataSO.combined.no.cluster)))==0)== FALSE, 
                  "You have to compute nearest neighbors first"))
    validate(need(my_react_vals$loaded.dataSO.combined.no.cluster@meta.data$seurat_clusters, message = "You have to identify clusters at the previous step first"))

    
    if(length(dim(my_react_vals$loaded.dataSO.combined.markers)[1]==0)==0){
      "You have to complete previous steps in Cell type clustering analysis"
    }else if(dim(my_react_vals$loaded.dataSO.combined.markers)[1]==0){
      "No clusters to show. Either DE genes not detected and/or couldn't find cells with identity class. You have to change parameters at previous steps, before proceeding"
    }else{
    "Marker definition completed"}
  })%>% bindEvent(input$cellMarkers , ignoreInit = TRUE, ignoreNULL = TRUE)
  
    #define top markers
  top_mark_def <- reactive({
    req(my_react_vals$loaded.dataSO.combined.markers)
    my_react_vals$loaded.dataSO.combined.markers %>%
      group_by(cluster) %>%
      slice_max(n = input$nTopMarkers, order_by = avg_log2FC) 
   })
  
  observeEvent(input$topMarkers ,{
    top_mark_def()
    my_react_vals$loaded.dataSO.combined.markerstop <- top_mark_def() ############## <<- #############################
  })
  
  
  output$userTopMarkers <- reactive({
    validate(
      need(my_react_vals$loaded.dataSO.combined.no.cluster, "It seems you have skipped some previous steps!")
    )
    validate(
      need(try(my_react_vals$loaded.dataSO.combined.markers), "You have to define cell type markers first")
    )
    #my_react_vals$loaded.dataSO.combined.markerstop<- top_mark_def()
    #message printed
    paste0("We'll use top ", input$nTopMarkers, " marker genes for cell type annotation")
    
  }) %>% bindEvent(input$topMarkers , ignoreInit = TRUE, ignoreNULL = TRUE)
  
    #Plot heatmap of top 20 markers 
  output$HeatmapMarkers <- renderPlot({
    validate(
      need(my_react_vals$loaded.dataSO.combined.no.cluster, "It seems you have skipped some previous steps!")
    )
    validate(
      need(try(my_react_vals$loaded.dataSO.combined.markers), "You have to define cell type markers first and then choose the number of marker genes before plotting the heatmap")
    )
    validate(
      need(try(my_react_vals$loaded.dataSO.combined.markerstop), "You have to choose the number of marker genes before plotting the heatmap")
    )
    #loaded.dataSO.combined.markerstop<<- top_mark_def()
    #plot heatmap
    DoHeatmap(my_react_vals$loaded.dataSO.combined.no.cluster, features = my_react_vals$loaded.dataSO.combined.markerstop$gene[1:20]) + NoLegend()
    
  })%>% bindEvent(input$plotHeatmap , ignoreInit = TRUE, ignoreNULL = TRUE)
 ##############################################################################################################################################
  #annotate gene clusters
  observeEvent(input$annotClusters,{
    validate(need(my_react_vals$loaded.dataSO.combined.markerstop, message = FALSE))
    spsComps::shinyCatch({
      my_react_vals$gen_cluster_res <- annotateMyClusters(input$annotLibrary, input$typeOrganism, my_react_vals$loaded.dataSO.combined.markerstop, my_react_vals$loaded.dataSO.combined.no.cluster, my_react_vals$loaded.dataSO.combined.markers)
    },
    # blocking recommended
    blocking_level = "error",
    prefix = "My-project" #change console prefix if you don't want "SPS"
    )
  })
  
  #gen_cluster_res <- reactive({
   # validate(need(my_react_vals$loaded.dataSO.combined.markerstop, message = FALSE))
   # spsComps::shinyCatch({
   #   annotateMyClusters(input$annotLibrary, input$typeOrganism, my_react_vals$loaded.dataSO.combined.markerstop, my_react_vals$loaded.dataSO.combined.no.cluster, my_react_vals$loaded.dataSO.combined.markers)
   # },
    # blocking recommended
   # blocking_level = "error",
   # prefix = "My-project" #change console prefix if you don't want "SPS"
   # )
    
 # }) #%>% bindEvent(input$annotClusters , ignoreInit = TRUE, ignoreNULL = TRUE)
  
   output$anot_completed <- renderText({
     validate(need(my_react_vals$loaded.dataSO.combined.no.cluster, message = "You have to complete the normalization and reduction step first"))
     validate(
       need(try(my_react_vals$loaded.dataSO.combined.markers), "You have to define cell type markers first and then choose the number of marker genes before proceeding")
     )
     validate(need(my_react_vals$loaded.dataSO.combined.markerstop, message = "You have to choose how many markers to use first"))
     validate(
       need(try(my_react_vals$gen_cluster_res), message = FALSE)
     )
    #gen_cluster_res <- gen_cluster_res()
   "Annotation completed"
  }) %>% bindEvent(input$annotClusters , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  
  #update selectInput for gene clusters 
  observe({
    validate(need(my_react_vals$gen_cluster_res, message = FALSE))
    updateSelectInput(session, "selectCluster",
                      label = "Choose a gene cluster",
                      choices =  levels(my_react_vals$gen_cluster_res[1][[1]]$celltype))
  })%>% bindEvent(input$showClusters , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  
  #Plot annotation results 
  output$anotResults <- renderPlot({
    validate(need(my_react_vals$loaded.dataSO.combined.no.cluster, message = "You have to complete the normalization and reduction step first"))
    validate(
      need(try(my_react_vals$loaded.dataSO.combined.markers), "You have to define cell type markers first and then choose the number of marker genes before proceeding")
    )
    validate(need(my_react_vals$loaded.dataSO.combined.markerstop, message = "You have to choose how many markers to use first"))
    validate(need(my_react_vals$gen_cluster_res, message = FALSE))
    
    enriched <- my_react_vals$gen_cluster_res[4][[1]]
    if(input$typeOrganism == "Mus musculus"){
      if (getOption("enrichR.live")){
        #for mouse
       plot(plotEnrich(enriched[[paste("PanglaoDB_Augmented_2021_",input$selectCluster,sep="")]],  numChar = 40, y = "Count", orderBy = "P.value",title = paste("PanglaoDB"))+plotEnrich(enriched[[paste("CellMarker_Augmented_2021_",input$selectCluster,sep="")]], numChar = 40, y = "Count", orderBy = "P.value",title = paste("CellMarker"))+plotEnrich(enriched[[paste("Tabula_Muris_",input$selectCluster,sep="")]], numChar = 40, y = "Count", orderBy = "P.value",title = paste("Tabula_Muris")))
      }else{message("Cannot resolve host: maayanlab.cloud. Check your internet connection and try again")}
    }else{
      if (getOption("enrichR.live")){
        # for human
        plot(plotEnrich(enriched[[paste("PanglaoDB_Augmented_2021_",input$selectCluster,sep="")]],  numChar = 40, y = "Count", orderBy = "P.value",title = paste("PanglaoDB"))+plotEnrich(enriched[[paste("CellMarker_Augmented_2021_",input$selectCluster,sep="")]], numChar = 40, y = "Count", orderBy = "P.value",title = paste("CellMarker"))+plotEnrich(enriched[[paste("Tabula_Sapiens_",input$selectCluster,sep="")]], numChar = 40, y = "Count", orderBy = "P.value",title = paste("Tabula_Sapiens")))
      }else{message("Cannot resolve host: maayanlab.cloud. Check your internet connection and try again")}
    }
  })%>% bindEvent(input$plotAnnot , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #UMAP plot
  #compute UMAP first
  #loaded.dataSO.combined <- reactive({
    observeEvent(input$plotUMAP, {
    req(my_react_vals$gen_cluster_res)
    req(input$nPCs)
    #if(input$typeData==4){
    if(my_react_vals$multiple_case){
      spsComps::shinyCatch({

        my_react_vals$loaded.dataSO.combined <- RunUMAP(my_react_vals$gen_cluster_res[1][[1]], reduction = "integrated.cca", dims = 1:input$nPCs, reduction.name = "umap")

        },
      # blocking recommended
      blocking_level = "error",
      prefix = "My-project" #change console prefix if you don't want "SPS"
      )
    }else{
      spsComps::shinyCatch({
        my_react_vals$loaded.dataSO.combined <- RunUMAP(my_react_vals$gen_cluster_res[1][[1]], reduction = "pca", dims = 1:input$nPCs, reduction.name = "umap")
      },
      # blocking recommended
      blocking_level = "error",
      prefix = "My-project" #change console prefix if you don't want "SPS"
      )
    }
  })#%>% bindEvent(input$plotUMAP , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #prevent withSpinner running at output$vizClusters when app starts
  output$vizClusters <- renderText({NULL})
  
  observeEvent(input$plotUMAP, {
    #req(loaded.dataSO.combined())
  output$vizClusters <- renderPlot({
    
    validate(need(my_react_vals$loaded.dataSO.combined.no.cluster, message = "You have to complete the normalization and reduction step first"))
    validate(
      need(try(my_react_vals$loaded.dataSO.combined.markers), "You have to define cell type markers first and then choose the number of marker genes before proceeding")
    )
    validate(need(my_react_vals$loaded.dataSO.combined.markerstop, message = "You have to choose how many markers to use first"))
    validate(need(my_react_vals$gen_cluster_res, message = "You have to complete the cluster annotation step first"))
    validate(need(my_react_vals$loaded.dataSO.combined, message = FALSE))
    
    DimPlot(my_react_vals$loaded.dataSO.combined, reduction = "umap", group.by = c("orig.ident", "celltype") ,combine = TRUE, label.size = 2)
    })#%>% bindEvent(input$plotUMAP , ignoreInit = TRUE, ignoreNULL = TRUE)
  },ignoreInit = TRUE, ignoreNULL = TRUE) ################################################ check if ignoreInit works
  #})
  
  
  #review expression of top markers
  output$vizTopMarkersExpr <- renderPlot({
    validate(need(my_react_vals$loaded.dataSO.combined.markerstop, message = "You have to complete previous steps first"))
    validate(need(my_react_vals$gen_cluster_res, message = "You have to complete the cluster annotation step first"))
    validate(need(my_react_vals$loaded.dataSO.combined, message = "You have to compute UMAP at the previous step"))
    #loaded.dataSO.combined.markerstop <<- gen_cluster_res()[3][[1]]
    markers.to.plot <- unique(my_react_vals$loaded.dataSO.combined.markerstop$gene) 
    DotPlot(my_react_vals$loaded.dataSO.combined, features = markers.to.plot[1:10], cols = c("blue", "red"), dot.scale = 8, group.by = "celltype") +
      RotatedAxis()
  })%>% bindEvent(input$topMarkersExpr , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  
  #Vizualize markers
  #update dropdown menu of selectizeInput for available marker genes
  observe({
    validate(need(my_react_vals$loaded.dataSO.combined.markerstop, message = FALSE))
    #loaded.dataSO.combined.markerstop <<- gen_cluster_res()[3][[1]]
    updateSelectizeInput(session, "selectMarkerGene",
                      label = "Choose a gene",
                      choices =  unique(my_react_vals$loaded.dataSO.combined.markerstop$gene), server =TRUE )
  })%>% bindEvent(input$showMarkerGene , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #Violin plot 
  output$vizGeneViolin <- renderPlot({
    validate(
      need(try(my_react_vals$loaded.dataSO.combined), "You have to complete previous steps first")
    )
    validate(
      need(try(input$selectMarkerGene), "Please choose a marker from the dropdown menu")
    )
    #loaded.dataSO.combined <<- gen_cluster_res()[1][[1]]
    VlnPlot(my_react_vals$loaded.dataSO.combined, features = input$selectMarkerGene)
  })%>% bindEvent(input$plotViolin , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #Feature plot 
  output$vizGeneFeature <- renderPlot({
    validate(
      need(try(my_react_vals$loaded.dataSO.combined), "You have to complete previous steps first")
    )
    validate(
      need(try(input$selectMarkerGene), "Please choose a marker from the dropdown menu")
    )
    #loaded.dataSO.combined <<- gen_cluster_res()[1][[1]]
    FeaturePlot(my_react_vals$loaded.dataSO.combined, features = input$selectMarkerGene,label = T)& theme(legend.position = c(0.1,0.2))
  })%>% bindEvent(input$plotFeature , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #Proportion of cells
  output$vizPropCells <- renderPlot({
    validate(
      need(try(my_react_vals$loaded.dataSO.combined), "You have to complete previous steps first")
    )
    #loaded.dataSO.combined <<- gen_cluster_res()[1][[1]]
    proportions<-as.data.frame(prop.table(table(my_react_vals$loaded.dataSO.combined$celltype, my_react_vals$loaded.dataSO.combined$orig.ident), margin = 2))
    colnames(proportions) <- c("CellIDs", "Condition", "Proportion")
    #pie chart with proportions
    #ggplot(proportions, aes(x="", y=Proportion, fill=CellIDs)) +
    #  geom_bar(stat="identity", width=1, color="white") +
    #  coord_polar("y", start=0) +
    #  theme_void() # remove background, grid, numeric labels
    
    #barplot works better with multiple samples
    ggplot(proportions, aes(Condition, Proportion, fill = CellIDs))+
      geom_col() +
      #geom_text(aes(label = scales::percent(Proportion)),
      #          position = position_stack(vjust = 0.5)) +
      scale_y_continuous(labels = scales::percent)
    
  })%>% bindEvent(input$plotPropCells , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  
  ##### Correlation #####
  #observeEvent(input$corData,{
   observe({ 
    #request minimal filtered data
    validate(need(calc_se()$my_list, FALSE))
     validate(need(input$corData, message = FALSE))
    #update list of sample and list of genes
    #if(input$typeData==4){
      
      if(input$corData=='filtered'){
        #validate(need(try(se_c()), message="You have to complete filtering steps first before running a correlation analysis"))
        validate(need(try(my_react_vals$filt_samples_list), message = "You have to complete filtering steps first before running a correlation analysis"))
        validate(need(any(mapply(is.null, my_react_vals$filt_samples_list )) == FALSE, message = "It seems you forgot to filter out cells in some/all of your samples"))
        #update samples choices  
        updateSelectInput(inputId="targetGenesSampleFiltered",
                    label = "Choose a sample",
                    choices =  names(my_react_vals$filt_samples_list),
                    selected = NULL)
        
        #update gene choices
        #updateSelectizeInput(session, "targetGenes",
        #                     label = "'Bait' gene to run correlation analysis with",
       #                      choices =  array(my_react_vals$filt_samples_list[[input$targetGenesSample]]$se@assays@data@listData$counts@Dimnames[1][[1]]), server = TRUE)
        #my_react_vals$cor_counts <- my_react_vals$filt_samples_list[[input$targetGenesSample]]$se@assays@data@listData$counts
        
      }else if(input$corData=='unfiltered'){

        validate(need(try(calc_se()$my_list), "You have to complete initial filtering at the preprocessing stage before running a correlation analysis"))
        #cor_unfiltered1 <<- calc_se()$my_list
        #update samples choices  
        updateSelectInput(inputId="targetGenesSampleUnfiltered",
                          label = "Choose a sample",
                          choices =  names(calc_se()$my_list),
                          selected = NULL)
        
        
        #update gene choices
        #validate(need(input$targetGenesSample, message = FALSE))
        #updateSelectizeInput(session, "targetGenes",
        #                     label = "'Bait' gene to run correlation analysis with",
        #                     choices =  array(calc_se()$my_list[[cor_sample]]$se@assays@data@listData$counts@Dimnames[1][[1]]), server = TRUE)
        #my_react_vals$cor_counts <- calc_se()$my_list[[cor_sample]]$se@assays@data@listData$counts
        
      }else if(input$corData=='clustered'){
        validate(need(try(my_react_vals$gen_cluster_res[1][[1]]), "You have to complete the annotation of clusters at the cell type clustering stage before running a correlation analysis"))
        validate(need(try(my_react_vals$filt_samples_list), message = "You have to complete filtering steps first before running a correlation analysis"))
        validate(need(any(mapply(is.null, my_react_vals$filt_samples_list )) == FALSE, message = "It seems you forgot to filter out cells in some/all of your samples"))
        validate(need(unique(my_react_vals$loaded.dataSO.combined.no.cluster$orig.ident), message = FALSE))
        #update samples choices  
        updateSelectInput(inputId="targetGenesSampleClustered",
                          label = "Choose a sample",
                          choices =  unique(my_react_vals$loaded.dataSO.combined.no.cluster$orig.ident),
                          selected = NULL)
        
        #update cell cluster choices
        #updateSelectInput(inputId="targetGenesCluster",
        #            label = "Choose a cell cluster",
        #            choices =  levels(subset(gen_cluster_res()[1][[1]], subset = orig.ident == input$targetGenesSample)$celltype))
        #update gene choices
        #updateSelectizeInput(session, "targetGenes",
        #                     label = "'Bait' gene to run correlation analysis with",
        #                     choices =  array(subset(gen_cluster_res()[1][[1]],  orig.ident == input$targetGenesSample & celltype == input$targetGenesCluster)@assays$RNA$counts@Dimnames[1][[1]]), server = TRUE)
        #my_react_vals$cor_counts = subset(gen_cluster_res()[1][[1]],  orig.ident == input$targetGenesSample & celltype == input$targetGenesCluster)@assays$RNA$counts
      }
      
   # }else{
      
   # }
    
    #updateSelectizeInput(session, "targetGenes",
    #                     label = "'Bait' gene to run correlation analysis with",
    #                     choices =  array(calc_se()@assays$RNA$counts@Dimnames[1]), server = TRUE)
    
    
   # validate(need(try(levels(my_react_vals$loaded.dataSO.combined$celltype)), "You have to complete gene cluster annotation steps first!"))
    #validate(need(try(loaded.dataSO.combined.no.cluster), "You have to complete gene cluster annotation steps first!"))
    
    #if(input$corData=="clustered"){
    #  updateSelectInput(session, "targetGenesCluster",
     #                   label = "Choose a gene cluster",
     #                   choices =  levels(my_react_vals$loaded.dataSO.combined$celltype.label))
      
    #}
     #})#, ignoreInit = TRUE, ignoreNULL = TRUE)
  })%>% bindEvent(input$corData , ignoreInit = TRUE, ignoreNULL = TRUE)
  
    #observe({
    observeEvent(input$targetGenesSampleFiltered,{
      validate(need(my_react_vals$filt_samples_list, message = FALSE))
      validate(need(input$targetGenesSampleFiltered, message = FALSE))
      updateSelectizeInput(session, "targetGenesFiltered",
                           label = "'Bait' gene to run correlation analysis with",
                           choices =  array(my_react_vals$filt_samples_list[[input$targetGenesSampleFiltered]]@assays@data@listData$counts@Dimnames[1][[1]]), server = TRUE)
      #my_react_vals$cor_counts <- my_react_vals$filt_samples_list[[input$targetGenesSampleFiltered]]$se@assays@data@listData$counts
    })
    observeEvent(input$targetGenesSampleUnfiltered,{
    #observe({
      validate(need(calc_se()$my_list, message= FALSE))
      validate(need(input$targetGenesSampleUnfiltered, message = FALSE))
      #cor_unfiltered2 <<- calc_se()$my_list
      updateSelectizeInput(session, "targetGenesUnfiltered",
                           label = "'Bait' gene to run correlation analysis with",
                           choices =  array(calc_se()$my_list[[input$targetGenesSampleUnfiltered]]$se@assays@data@listData$counts@Dimnames[1][[1]]), server = TRUE)
      #my_react_vals$cor_counts <- calc_se()$my_list[[input$targetGenesSampleUnfiltered]]$se@assays@data@listData$counts
    })
    observeEvent(input$targetGenesSampleClustered,{
   #observe({
      validate(need(my_react_vals$gen_cluster_res[1][[1]], message = FALSE))
      updateSelectInput(inputId="targetGenesCluster",
                        label = "Choose a cell cluster",
                        choices =  levels(subset(my_react_vals$gen_cluster_res[1][[1]], subset = orig.ident == input$targetGenesSampleClustered)$celltype),
                        selected = NULL)  
      
    })#%>% bindEvent(input$targetGenesSampleClustered, ignoreInit = TRUE, ignoreNULL = TRUE)
      
    observeEvent(input$targetGenesCluster,{
     #observe({ 
      validate(need(my_react_vals$gen_cluster_res[1][[1]], message = FALSE))
      #update gene choices
      updateSelectizeInput(session, "targetGenes",
                           label = "'Bait' gene to run correlation analysis with",
                           choices =  array(subset(my_react_vals$gen_cluster_res[1][[1]],  orig.ident == input$targetGenesSampleClustered & celltype == input$targetGenesCluster)@assays$RNA$counts@Dimnames[1][[1]]), server = TRUE)
      #my_react_vals$cor_counts = subset(gen_cluster_res()[1][[1]],  orig.ident == input$targetGenesSampleClustered & celltype == input$targetGenesCluster)@assays$RNA$counts
    })#%>% bindEvent(input$targetGenesCluster, ignoreInit = TRUE, ignoreNULL = TRUE)#, ignoreInit = TRUE, ignoreNULL = TRUE)
    
    
    #observeEvent(my_react_vals$cor_sample,{
       
        #validate(need(my_react_vals$filt_samples_list, message = FALSE))
        #validate(need(calc_se()$my_list, message = FALSE))
        #validate(need(gen_cluster_res()[1][[1]], message = FALSE))
        
        #if(input$corData=='filtered'){
          #update gene choices
         # updateSelectizeInput(session, "targetGenes",
         #                      label = "'Bait' gene to run correlation analysis with",
         #                      choices =  array(my_react_vals$filt_samples_list[[input$targetGenesSample]]$se@assays@data@listData$counts@Dimnames[1][[1]]), server = TRUE)
        #  my_react_vals$cor_counts <- my_react_vals$filt_samples_list[[my_react_vals$cor_sample]]$se@assays@data@listData$counts
        #}else if(input$corData=='unfiltered'){
         # cor_unfiltered2 <<- calc_se()$my_list
         # updateSelectizeInput(session, "targetGenes",
          #                     label = "'Bait' gene to run correlation analysis with",
         #                      choices =  array(calc_se()$my_list[[my_react_vals$cor_sample]]$se@assays@data@listData$counts@Dimnames[1][[1]]), server = TRUE)
         # my_react_vals$cor_counts <- calc_se()$my_list[[my_react_vals$cor_sample]]$se@assays@data@listData$counts
        #}else if(input$corData=='clustered'){
          #update cell cluster choices
         # updateSelectInput(inputId="targetGenesCluster",
         #                   label = "Choose a cell cluster",
         #                   choices =  levels(subset(gen_cluster_res()[1][[1]], subset = orig.ident == my_react_vals$cor_sample)$celltype))  
         # my_react_vals$cor_cluster <- isolate(my_react_vals$cor_sample)
       # }
      #}, ignoreInit = TRUE, ignoreNULL = TRUE)
      
      #observeEvent(my_react_vals$cor_cluster,{
     #   cor_cluster <<- my_react_vals$cor_cluster
        #validate(need(gen_cluster_res()[1][[1]], message = FALSE))
        #update gene choices
      #  updateSelectizeInput(session, "targetGenes",
      #                       label = "'Bait' gene to run correlation analysis with",
     #                        choices =  array(subset(gen_cluster_res()[1][[1]],  orig.ident == my_react_vals$cor_cluster & celltype == my_react_vals$cor_cluster)@assays$RNA$counts@Dimnames[1][[1]]), server = TRUE)
     #   my_react_vals$cor_counts = subset(gen_cluster_res()[1][[1]],  orig.ident == my_react_vals$cor_cluster & celltype == my_react_vals$cor_cluster)@assays$RNA$counts
    #  }, ignoreInit = TRUE, ignoreNULL = TRUE)

  ########### how to subset layered seurat:
  # subset(loaded.dataSO.combined, subset = orig.ident == "S190_rds")
  
  observeEvent(input$runCor, {
      withCallingHandlers({
        shinyjs::html("cor_completed", "")
    if(length(input$corData) == 0 | length(input$method) == 0){
      message("You have to define which method and/or which data you want to run the analysis with")
    }else{
    if(input$corData=="filtered"){
      #cor_dat_filtered <<-my_react_vals$filt_samples_list[[input$targetGenesSampleFiltered]]@assays@data@listData$counts
      my_react_vals$iRNA_cor_tables <- tryCatch( ############## <<- #############################
        {
          iRNA(my_react_vals$filt_samples_list[[input$targetGenesSampleFiltered]]@assays@data@listData$counts, input$method, input$alpha, input$geneSparsity, input$targetGenesFiltered, input$correlation_threshold, input$corData)
        },
        error = function(e) {
         message("No filtered data available. Did you forget to filter out cells after plotting metrics?")
        }
      )
      
      
    }else if(input$corData=="unfiltered"){
      #"Running analysis with raw data"
      #cor_dat_unfiltered <<- calc_se()$my_list[[input$targetGenesSampleUnfiltered]]$se@assays@data@listData$counts
      my_react_vals$iRNA_cor_tables <- tryCatch({ ############## <<- #############################
        iRNA(calc_se()$my_list[[input$targetGenesSampleUnfiltered]]$se@assays@data@listData$counts, input$method, input$alpha, input$geneSparsity, input$targetGenesUnfiltered, input$correlation_threshold, input$corData)
      },
      error = function(e) {
        message("Seems initial filtering at the preprocessing stage is not complete. Please go back and follow all the steps of preprocessing")
      })
      
    }else{
      cat(paste0("Running analysis with cell type cluster:", input$targetGenesCluster))
      #cor_dat_clustered <<-subset(gen_cluster_res()[1][[1]],  orig.ident == input$targetGenesSampleClustered & celltype == input$targetGenesCluster)@assays$RNA$counts
      my_react_vals$iRNA_cor_tables <- tryCatch({ ############## <<- #############################
        iRNA_cluster(subset(my_react_vals$gen_cluster_res[1][[1]],  orig.ident == input$targetGenesSampleClustered & celltype == input$targetGenesCluster)@assays$RNA$counts, input$method, input$alpha, input$geneSparsity, input$targetGenes, input$correlation_threshold, input$targetGenesCluster)
        #iRNA_cluster(loaded.dataSO.combined, input$method, input$alpha, input$geneSparsity, input$targetGenes, input$correlation_threshold, input$selectCluster)
        
        },
      error = function(e) {
        message("Cell type clusters not computed yet. Please go back and run cell clustering steps")
        })
    }}
  },
  message = function(m) {
    shinyjs::html(id = "cor_completed", html = paste0(m$message, '<br>'), add = TRUE)
      })
    })
  
  #empty outputs when changing input data
  observeEvent(input$runCor,{
    #hide("readData")
    output$dt1 <- DT::renderDataTable({NULL})
    output$dt2 <- DT::renderDataTable({NULL})
    output$dt3 <- DT::renderDataTable({NULL})
  })
  
  #Show correlation results in tables
  #All correlated genes
  
  observeEvent(input$runCor, {
      validate(need(my_react_vals$iRNA_cor_tables, message = FALSE))
      output$dt1 <- DT::renderDataTable({
        if (identical(my_react_vals$iRNA_cor_tables, list())){
          
          return(NULL)
          #iRNA_cor_tables[[1]]
        }else if(is.null(my_react_vals$iRNA_cor_tables[[1]])){
          return(NULL)
         }else{ DT::datatable(my_react_vals$iRNA_cor_tables[[1]],
                        filter = 'top',
                        extensions = 'Buttons',
                        options = list(dom = 'Blfrtip',
                                       buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                                       lengthMenu = list(c(10, 20, 50, 100, -1), c('10', '20', '50', '100', 'All')),
                                       pageLength = 15),
                        #escape = c(3,4,5,6)
          )%>%
             formatSignif(columns = c("p_r", "p_w"))%>%
             formatRound(columns=c("r", "percentage in Ea", "percentage in Eb"), digits=3)}
        })#%>%bindEvent(input$tableCor, ignoreInit = TRUE, ignoreNULL = TRUE)
      
      #positive correlated genes
      output$dt2 <- DT::renderDataTable({
        if (identical(my_react_vals$iRNA_cor_tables, list())){
          return(NULL)
        }else if(is.null(my_react_vals$iRNA_cor_tables[[2]])){
          return(NULL)
        }else{ DT::datatable(my_react_vals$iRNA_cor_tables[[2]],
                             filter = 'top',
                             extensions = 'Buttons',
                             options = list(dom = 'Blfrtip',
                                            buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                                            lengthMenu = list(c(10, 20, 50, 100, -1), c('10', '20', '50', '100', 'All')),
                                            pageLength = 15),
                             #escape = c(3,4,5,6)
        )%>%
            formatSignif(columns = c("p_r", "p_w"))%>%
            formatRound(columns=c("r", "percentage in Ea", "percentage in Eb"), digits=3)}# DT::datatable
      })#%>%bindEvent(input$tableCor, ignoreInit = TRUE, ignoreNULL = TRUE)
      
      #negative correlated genes
      output$dt3 <- DT::renderDataTable({
        if (identical(my_react_vals$iRNA_cor_tables, list())){
          return(NULL)
        }else if(is.null(my_react_vals$iRNA_cor_tables[[3]])){
          return(NULL)
        }else{ DT::datatable(my_react_vals$iRNA_cor_tables[[3]],
                             filter = 'top',
                             extensions = 'Buttons',
                             options = list(dom = 'Blfrtip',
                                            buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                                            lengthMenu = list(c(10, 20, 50, 100, -1), c('10', '20', '50', '100', 'All')),
                                            pageLength = 15),
                             #escape = c(3,4,5,6)
        )%>%
            formatSignif(columns = c("p_r", "p_w"))%>%
            formatRound(columns=c("r", "percentage in Ea", "percentage in Eb"), digits=3)}
      })#%>%bindEvent(input$tableCor, ignoreInit = TRUE, ignoreNULL = TRUE)

  }, ignoreNULL = FALSE)
  
  
  ##### Gene enrichment #####
  observeEvent(input$runEnrich, {
    validate(need(my_react_vals$iRNA_cor_tables, message = FALSE))
    withCallingHandlers({
      shinyjs::html("enrich_completed", "")
      
      my_react_vals$iRNA_enrich_tables <- run_gene_enrich(my_react_vals$iRNA_cor_tables, input$typeOrganism) ############## <<- #############################
    },
    message = function(m) {
      shinyjs::html(id = "enrich_completed", html = paste0(m$message, '<br>'), add = TRUE)
    })
  })
  
  #Show enrichment results in tables
  observeEvent(input$runEnrich, {
    validate(need(my_react_vals$iRNA_enrich_tables, message = FALSE))
    #iRNA_enrich_tables  <<- my_react_vals$iRNA_enrich_tables ############## <<- #############################
  #All correlated genes
  output$dt1_enrich <- DT::renderDataTable({
    if (identical(my_react_vals$iRNA_enrich_tables, list())){
      
      return(NULL)
      #iRNA_cor_tables[[1]]
    }else if(is.null(dim(my_react_vals$iRNA_enrich_tables[[1]]))){
      return(NULL)
    }else{ DT::datatable(as.data.frame(my_react_vals$iRNA_enrich_tables[[1]]),
                         filter = 'top',
                         extensions = 'Buttons',
                         options = list(dom = 'Blfrtip',
                                        buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                                        lengthMenu = list(c(10, 20, 50, 100, -1), c('10', '20', '50', '100', 'All')),
                                        pageLength = 15),
                         #escape = c(3,4,5,6)
    )}# DT::datatable
  })#%>%bindEvent(input$tableEnrich, ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #positive correlated genes
  output$dt2_enrich <- DT::renderDataTable({
    if (identical(my_react_vals$iRNA_enrich_tables, list())){
      return(NULL)
    }else if(is.null(dim(my_react_vals$iRNA_enrich_tables[[2]]))){
      return(NULL)
    }else{ DT::datatable(as.data.frame(my_react_vals$iRNA_enrich_tables[[2]]),
                         filter = 'top',
                         extensions = 'Buttons',
                         options = list(dom = 'Blfrtip',
                                        buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                                        lengthMenu = list(c(10, 20, 50, 100, -1), c('10', '20', '50', '100', 'All')),
                                        pageLength = 15),
                         #escape = c(3,4,5,6)
    )}# DT::datatable
  })#%>%bindEvent(input$tableEnrich, ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #negative correlated genes
  output$dt3_enrich <- DT::renderDataTable({
    if (identical(my_react_vals$iRNA_enrich_tables, list())){
      return(NULL)
    }else if(is.null(dim(my_react_vals$iRNA_enrich_tables[[3]]))){
      return(NULL)
    }else{ DT::datatable(as.data.frame(my_react_vals$iRNA_enrich_tables[[3]]),
                         filter = 'top',
                         extensions = 'Buttons',
                         options = list(dom = 'Blfrtip',
                                        buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                                        lengthMenu = list(c(10, 20, 50, 100, -1), c('10', '20', '50', '100', 'All')),
                                        pageLength = 15),
                         #escape = c(3,4,5,6)
    )}# DT::datatable
  })#%>%bindEvent(input$tableEnrich, ignoreInit = TRUE, ignoreNULL = TRUE)
  }, ignoreNULL = FALSE)
  
  
  #delete all files
  session$onSessionEnded(function() { unlink(input, recursive = TRUE) })
  }


  #shinyApp(ui = ui, server = server)