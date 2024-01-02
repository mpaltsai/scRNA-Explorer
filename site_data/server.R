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
library(tools)
library(dplyr)
library(utils)


options(shiny.maxRequestSize=5000*1024^2)  #max upload file size 5GB


#load functions
#source("scRNA-Explorer/scripts/readData.R", local = TRUE)

function(input, output, session){
  
  ##### Upload, setting parameters and read data routine #####
  
  #Validate uploaded .csv file
  output$uploadCsvValidation <- reactive({
    validate(need(identical(tools::file_ext(input$csvFile$datapath),"csv"),"Invalid extension (not a .csv file)"))
    }) %>% bindEvent(input$csvFile, ignoreInit = TRUE, ignoreNULL=TRUE)
  
  #Validate uploaded .rds file
  output$uploadRdsValidation <- reactive({
    validate(need(identical(tools::file_ext(input$rdsFile$datapath),"rds"),"Invalid extension (not a .rds file)"))
  }) %>% bindEvent(input$rdsFile, ignoreInit = TRUE, ignoreNULL=TRUE)
  
  #Check if there ism't any file uploaded or state which kind of input the user has provided
  output$uploadSummary <- reactive({
      #need(input$dir, 'Provide a directory'),
      if(input$testData==1){
        count.matrix <<- TRUE
        seurat <<- TRUE
        inputFile <<- "raw_se_S190.rds"
        inputDir <<- FALSE
        paste0("Your input to the following QC pipeline is the test dataset.","\n",
               "You'll have to define origin organism and features annotation type below.", "\n",
               "This dataset comes from Mus mussculus and genes are annotated as gene names.")
      }else if(input$dir== "Enter directory path..." & all(is.null(c(input$csvFile, input$rdsFile)))){
        "Please provide an input (file or directory)"
      }else if(input$typeData == 1){
        count.matrix <<- TRUE
        seurat <<- FALSE
        inputFile <<- input$csvFile$datapath
        inputDir <<- FALSE
        paste0("Your input to the following QC pipeline is: ", input$csvFile$name)
      }else if(input$typeData ==2){
        count.matrix <<- FALSE
        seurat <<- FALSE
        inputFile <<- FALSE
        inputDir <<- TRUE
        paste0("The directory to search for additional files is: ", input$dir)
      } else if(input$typeData==3){
        count.matrix <<- TRUE
        seurat <<- TRUE
        inputFile <<- input$rdsFile$datapath
        inputDir <<- FALSE
        paste0("Your input to the following QC pipeline is: ", input$rdsFile$name)
      }
   #validate(need(identical(all(is.null(c(input$csvFile, input$dir, input$rdsFile))), NULL), "Please upload a file"))
      #combn(c(input$csvFile, input$dir, input$rdsFile), 2, simplify= FALSE)
    }) %>% bindEvent(input$checkInput, ignoreInit = TRUE, ignoreNULL = FALSE)
  
  ### Text box to check inputs and backend parameters
  #output$inputParameters <- reactive({
  #  paste0("count.matrix: ", count.matrix, " seurat: ", seurat, " input dir: " ,inputDir, " input file: ", inputFile)
  #  }) %>% bindEvent(input$checkInput, ignoreInit = TRUE, ignoreNULL = TRUE)
  
  ##### Preprocessing #####
  
  #read data
  output$readData <- reactive({
    #source("scripts/load_libraries.R", local=TRUE)
    countsData <<- readData()
   "Reading data completed"
  }) %>% bindEvent(input$startQC, ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #add metadata
  output$paramsAddMetadata <- reactive({
    if(length(input$typeOrganism) == 0 | length(input$typeGeneId) == 0){
      ("You have to select the organism and gene IDs type")
    }else{
      metaData <<- createMetadata(countsData,input$typeOrganism, input$typeGeneId)
      "Metadata added."
    }
  }) %>% bindEvent(input$addMetadata, ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #minimal filtering and save as an Rds object
  output$initDataFilt<- reactive({
    # Keep cells with nUMI greater than 100
    idx <- which(metaData$nUMI > 100)
    
    # Extract the counts for those cells
    counts_c <- countsData[, idx]
    
    # Extract the metadata for those cells
    metadata_c <- metaData[idx,]
    
    # Save data to single cell experiment variable
    se <<- SingleCellExperiment(assays=list(counts=counts_c), 
                               colData = metadata_c)
    #remove countsData and metaData since they are stored in se object
    rm(list=c("countsData", "metaData"), envir = .GlobalEnv)
    
    #and delete initial csv or rds file
   
    if (count.matrix) {
      if(seurat){
        if(input$testData==1){
          NULL
        }else{
        file.remove(input$rdsFile$datapath)
        rm(list="inputFile", envir = .GlobalEnv)
      }}else{file.remove(input$csvFile$datapath)
        rm(list="inputFile", envir = .GlobalEnv)}
    }
    
    "Data filtered"
  }) %>% bindEvent(input$initFiltering, ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #save and download as rds
  output$dataNameValidation <- reactive({
    if(input$dataName=="Enter a name..."){
      "Please provide a name"
    }else if (input$dataName==""){
      "Please provide a name"
    }else{
      # Create .RData object to load at any time
      saveRDS(se, paste0("raw_se_", input$dataName, ".rds"))
      "Saved"
    }
  }) %>% bindEvent(input$initDownload , ignoreInit = TRUE, ignoreNULL = TRUE)
    
  ###### QC routine #####
  
  #create metrics
  output$createdMetrics <- reactive({
    metrics <<- colData(se) %>%
      as.data.frame
      "Metrics created"
    
  }) %>% bindEvent(input$makeMetrics , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #Update slider max values when metrics object is created
  #metrics_created <- reactive({get(metrics)})
  
  
  #Make plots
  #number of cells in the sample
  #output$NBcellsPlot <- renderPlot({
  #  metrics %>% 
  #    ggplot(aes(x=sample, fill=sample)) + 
  #    geom_bar(width=.6) + 
  #    ggtitle("Number of Cells")
  #  output$NBinSample <-({
  #    renderText({
  #      paste0("Number of features (genes): ", dim(se)[1], " Number of cells: ", dim(se)[2])
  #    })
  #})%>% bindEvent(input$makeNBcells , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #number of UMIs per cell
  #update slider min and max values

  observe({
    updateSliderInput(session,"nUMIs", "Number of UMIs:",
                      min = 0, max = max(metrics$nUMI), value = c(min(metrics$nUMI),max(metrics$nUMI)), step = 10)
  }) %>% bindEvent(input$makeNBUMIs , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  output$UMIsPlot <- renderPlot({
    
      createUMIsPlot(metrics, input$nUMIs[1], input$nUMIs[2])

  })%>% bindEvent(c(input$makeNBUMIs, input$nUMIs) , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #nGenes histogram
  observe({
    updateSliderInput(session, "nGenes", "Number of genes:",
                      min = 0, max = max(metrics$nGene), value = c(min(metrics$nGene),max(metrics$nGene)), step = 10)
  }) %>% bindEvent(input$makeNBgenesHist , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  output$GenesHistPlot<- renderPlot({
    createGenesHistPlot(metrics, input$nGenes[1], input$nGenes[2])

  })%>% bindEvent(c(input$makeNBgenesHist, input$nGenes) , ignoreInit = TRUE, ignoreNULL = TRUE)

  #nGenes boxplot
  output$GenesBoxPlot <- renderPlot({
    
    createGenesBoxPlot(metrics)
    
  })%>% bindEvent(c(input$makeNBgenesBoxPlot) , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #UMIs vs nGenes with mitoRatio
  observe({
    updateSliderInput(session,"nUMIsVSnGenes", "Genes detected:",
                      min = 0, max = max(metrics$nGene), value = c(min(metrics$nGene),max(metrics$nGene)), step = 10)
  })%>% bindEvent(input$makeUMIsGenes , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  output$UMIsGenesPlot <- renderPlot({
    
    createUMIsGenesPlot(metrics, input$nUMIsVSnGenes[1], input$nUMIsVSnGenes[2])
    
  })%>% bindEvent(c(input$makeUMIsGenes, input$nUMIsVSnGenes) , ignoreInit = TRUE, ignoreNULL = TRUE)
 
  
  #UMIs vs nGenes with hemoRatio
  observe({
    updateSliderInput(session,"nUMIsVSnGenesRBC", "Genes detected:",
                      min = 0, max = max(metrics$nGene), value = c(min(metrics$nGene),max(metrics$nGene)), step = 10)
  })%>% bindEvent(input$makeUMIsGenesRBC , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  output$UMIsGenesPlotRBC <- renderPlot({
    
    createUMIsGenesPlotRBC(metrics, input$nUMIsVSnGenesRBC[1], input$nUMIsVSnGenesRBC[2])
      
  })%>% bindEvent(c(input$makeUMIsGenesRBC, input$nUMIsVSnGenesRBC) , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #plot mitoRatio
  observe({
    updateSliderInput(session,"mitoRatio", "Mitochondrial ratio:",
                      min = 0, max = round(max(metrics$mitoRatio), digits = 2), value = c(round(min(metrics$mitoRatio),digits = 2),round(max(metrics$mitoRatio), digits = 2)), step = 0.01)
  })%>% bindEvent(input$makemitoRatio , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  output$mitoRatioPlot <- renderPlot({
    
    createmitoRatioPlot(metrics, input$mitoRatio[1], input$mitoRatio[2])
    
  })%>% bindEvent(c(input$makemitoRatio, input$mitoRatio) , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #plot hemoRatio
  observe({
    updateSliderInput(session,"rbcRatio", "Hemoglobins ratio:",
                      min = 0, max = round(max(metrics$rbcRatio), digits = 2), value = c(round(min(metrics$rbcRatio),digits = 2),round(max(metrics$rbcRatio), digits = 2)), step = 0.01)
  })%>% bindEvent(input$makerbcRatio , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  output$rbcRatioPlot <- renderPlot({
    
    createrbcRatioPlot(metrics, input$rbcRatio[1], input$rbcRatio[2])
    
  })%>% bindEvent(c(input$makerbcRatio, input$rbcRatio) , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #Genes vs UMIs
  observe({
    updateSliderInput(session,"novelty", "Genes detected per UMI:",
                      min = 0, max = round(max(metrics$log10GenesPerUMI), digits = 2), value = round(min(metrics$log10GenesPerUMI),digits = 2), step = 0.01)
  })%>% bindEvent(input$makeNovelty , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  output$noveltyPlot <- renderPlot({
    createnoveltyPlot(metrics, input$novelty[1], input$novelty[2])
  })%>% bindEvent(input$makeNovelty, input$novelty, ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #plot scaterPlot
  output$libraryPlot <- renderPlot({
    createlibraryPlot(se)
  })%>% bindEvent(input$makeLibrary, ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #plot highest expressed features
  output$highExprPlot <- renderPlot({
    createhighExprPlot(se)
  })%>% bindEvent(input$makeHighExpr, ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #zoom on plot
  observeEvent(input$zoomhighExprPlot, {
    showModal(modalDialog(
      renderPlot({
        createhighExprPlot(se)
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
  output$QCDataFilt<- reactive({
    keep_index <<- metrics %>%
      dplyr::filter(nUMI > input$minUMIsPerCell , 
                    nUMI < input$maxUMIsPerCell,
                    nGene > input$minGenesPerCell,
                    nGene < input$maxGenesPerCell,
                    log10GenesPerUMI > input$noveltyPerCell,
                    mitoRatio > input$minMitoPerCell,
                    mitoRatio < input$maxMitoPerCell,
                    rbcRatio > input$minRbcPerCell,
                    rbcRatio < input$maxRbcPerCell
      ) %>% 
      pull(cells)
    
    # Subset the cells to only include those that meet the thresholds specified and save subset to new metrics variable
   metrics_clean <<- colData(se[ ,keep_index]) %>%
      as.data.frame()
    "Data filtered"
  })%>% bindEvent(input$QCFiltering, ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #save and download as rds
  output$dataNameValidationQC <- reactive({
    if(input$dataNameQC=="Enter a name..."){
      "Please provide a name"
    }else if (input$dataNameQC==""){
      "Please provide a name"
    }else{
      # Create .RData object to load at any time
      saveRDS(se[ ,keep_index], paste0("se_filtered_", input$dataNameQC, ".rds"))
      "Saved"
    }
  }) %>% bindEvent(input$QCDownload , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  ##### Redraw Plots with filtered data #####
  #UMIs Plot
  output$UMIsPlotFiltered <- renderPlot({
    
    createUMIsPlot(metrics_clean, min(metrics_clean$nUMI),max(metrics_clean$nUMI))
    
  })%>% bindEvent(input$redrawPlots , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #nGenes histogram
  output$GenesHistPlotFiltered <- renderPlot({
    
    createGenesHistPlot(metrics_clean, min(metrics_clean$nGene),max(metrics_clean$nGene))
    
  })%>% bindEvent(input$redrawPlots , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #nGenes Boxplot
  output$GenesBoxPlotFiltered <- renderPlot({
    
    createGenesBoxPlot(metrics_clean)
    
  })%>% bindEvent(input$redrawPlots , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #UMIs vs genes with mitoRatio
  output$UMIsGenesPlotFiltered <- renderPlot({
    
    createUMIsGenesPlot(metrics_clean, min(metrics_clean$nGene), max(metrics_clean$nGene))
    
  })%>% bindEvent(input$redrawPlots , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #UMIs vs genes with rbcRatio
  output$UMIsGenesPlotRBCFiltered <- renderPlot({
    
    createUMIsGenesPlotRBC(metrics_clean, min(metrics_clean$nGene), max(metrics_clean$nGene))
    
  })%>% bindEvent(input$redrawPlots , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #mitoRatio
  output$mitoRatioPlotFiltered <- renderPlot({
    
    createmitoRatioPlot(metrics_clean, min(metrics_clean$mitoRatio), max(metrics_clean$mitoRatio))
    
  })%>% bindEvent(input$redrawPlots , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #rbcRatio
  output$rbcRatioPlotFiltered <- renderPlot({
    
    createrbcRatioPlot(metrics_clean, min(metrics_clean$rbcRatio), max(metrics_clean$rbcRatio))
    
  })%>% bindEvent(input$redrawPlots , ignoreInit = TRUE, ignoreNULL = TRUE)
  
   
  #novelty
  output$noveltyPlotFiltered <- renderPlot({
    
    createnoveltyPlot(metrics_clean, min(metrics_clean$log10GenesPerUMI), max(metrics_clean$log10GenesPerUMI))
    
  })%>% bindEvent(input$redrawPlots , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #library plot
  output$libraryPlotFiltered <- renderPlot({
    
    createlibraryPlot(se[ ,keep_index])
    
  })%>% bindEvent(input$redrawPlots , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #plot highest expressed genes
  output$highExprPlotFiltered <- renderPlot({
    
    createhighExprPlot(se[ ,keep_index])
    
  })%>% bindEvent(input$redrawPlots , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #zoom on plot
  observeEvent(input$zoomhighExprPlotFiltered, {
    showModal(modalDialog(
      renderPlot({
        createhighExprPlot(se[ ,keep_index])
      },
      width = 600,
      height = 1200,
      res = 300),
      footer = NULL,
      size = "l",
      easyClose = TRUE
    ))
  })
  
  ##### Correlation #####
  observeEvent(input$runCor, {
      withCallingHandlers({
        shinyjs::html("cor_completed", "")
    if(length(input$corData) == 0 | length(input$method) == 0){
      message("You have to define which method and/or which data you want to run the analysis with")
    }else{
    if(input$corData=="filtered"){
      "Running analysis with filtered data"
      iRNA_cor_tables <<- iRNA(se[ ,keep_index], input$method, input$alpha, input$geneSparsity, input$targetGenes, input$correlation_threshold)
      
    }else{
      "Running analysis with raw data"
      iRNA_cor_tables <<- iRNA(se, input$method, input$alpha, input$geneSparsity, input$targetGenes, input$correlation_threshold)
      #output$dt1 <- renderDataTable(DT::datatable(
        #iRNA_cor_tables[[1]]))
     # output$dt2 <- renderDataTable(DT::datatable(
      #  iRNA_cor_tables[[2]]))
      #output$dt3 <- renderDataTable(DT::datatable(
     #   iRNA_cor_tables[[3]]))
      
    }}
  },
  message = function(m) {
    shinyjs::html(id = "cor_completed", html = paste0(m$message, '<br>'), add = TRUE)
      })
    })
  
  #Show correlation results in tables
  #All correlated genes
  output$dt1 <- DT::renderDataTable({
    if (identical(iRNA_cor_tables, list())){
      
      return(NULL)
      #iRNA_cor_tables[[1]]
    }else if(is.null(iRNA_cor_tables[[1]])){
      return(NULL)
     }else{ DT::datatable(iRNA_cor_tables[[1]],
                    filter = 'top',
                    extensions = 'Buttons',
                    options = list(dom = 'Blfrtip',
                                   buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                                   lengthMenu = list(c(10, 20, 50, 100, -1), c('10', '20', '50', '100', 'All')),
                                   pageLength = 15),
                    #escape = c(3,4,5,6)
      )}# DT::datatable
    })%>%bindEvent(input$tableCor, ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #positive correlated genes
  output$dt2 <- DT::renderDataTable({
    if (identical(iRNA_cor_tables, list())){
      return(NULL)
    }else if(is.null(iRNA_cor_tables[[2]])){
      return(NULL)
    }else{ DT::datatable(iRNA_cor_tables[[2]],
                         filter = 'top',
                         extensions = 'Buttons',
                         options = list(dom = 'Blfrtip',
                                        buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                                        lengthMenu = list(c(10, 20, 50, 100, -1), c('10', '20', '50', '100', 'All')),
                                        pageLength = 15),
                         #escape = c(3,4,5,6)
    )}# DT::datatable
  })%>%bindEvent(input$tableCor, ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #negative correlated genes
  output$dt3 <- DT::renderDataTable({
    if (identical(iRNA_cor_tables, list())){
      return(NULL)
    }else if(is.null(iRNA_cor_tables[[3]])){
      return(NULL)
    }else{ DT::datatable(iRNA_cor_tables[[3]],
                         filter = 'top',
                         extensions = 'Buttons',
                         options = list(dom = 'Blfrtip',
                                        buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                                        lengthMenu = list(c(10, 20, 50, 100, -1), c('10', '20', '50', '100', 'All')),
                                        pageLength = 15),
                         #escape = c(3,4,5,6)
    )}# DT::datatable
  })%>%bindEvent(input$tableCor, ignoreInit = TRUE, ignoreNULL = TRUE)

  
  ##### Gene enrichment #####
  observeEvent(input$runEnrich, {
    withCallingHandlers({
      shinyjs::html("enrich_completed", "")
      
      iRNA_enrich_tables <<- run_gene_enrich(iRNA_cor_tables, input$typeOrganism)
    },
    message = function(m) {
      shinyjs::html(id = "enrich_completed", html = paste0(m$message, '<br>'), add = TRUE)
    })
  })
  
  #Show enrichment results in tables
  #All correlated genes
  output$dt1_enrich <- DT::renderDataTable({
    if (identical(iRNA_enrich_tables, list())){
      
      return(NULL)
      #iRNA_cor_tables[[1]]
    }else if(is.null(iRNA_enrich_tables[[1]])){
      return(NULL)
    }else{ DT::datatable(as.data.frame(iRNA_enrich_tables[[1]]),
                         filter = 'top',
                         extensions = 'Buttons',
                         options = list(dom = 'Blfrtip',
                                        buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                                        lengthMenu = list(c(10, 20, 50, 100, -1), c('10', '20', '50', '100', 'All')),
                                        pageLength = 15),
                         #escape = c(3,4,5,6)
    )}# DT::datatable
  })%>%bindEvent(input$tableEnrich, ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #positive correlated genes
  output$dt2_enrich <- DT::renderDataTable({
    if (identical(iRNA_enrich_tables, list())){
      return(NULL)
    }else if(is.null(iRNA_enrich_tables[[2]])){
      return(NULL)
    }else{ DT::datatable(as.data.frame(iRNA_enrich_tables[[2]]),
                         filter = 'top',
                         extensions = 'Buttons',
                         options = list(dom = 'Blfrtip',
                                        buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                                        lengthMenu = list(c(10, 20, 50, 100, -1), c('10', '20', '50', '100', 'All')),
                                        pageLength = 15),
                         #escape = c(3,4,5,6)
    )}# DT::datatable
  })%>%bindEvent(input$tableEnrich, ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #negative correlated genes
  output$dt3_enrich <- DT::renderDataTable({
    if (identical(iRNA_enrich_tables, list())){
      return(NULL)
    }else if(is.null(iRNA_enrich_tables[[3]])){
      return(NULL)
    }else{ DT::datatable(as.data.frame(iRNA_enrich_tables[[3]]),
                         filter = 'top',
                         extensions = 'Buttons',
                         options = list(dom = 'Blfrtip',
                                        buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                                        lengthMenu = list(c(10, 20, 50, 100, -1), c('10', '20', '50', '100', 'All')),
                                        pageLength = 15),
                         #escape = c(3,4,5,6)
    )}# DT::datatable
  })%>%bindEvent(input$tableEnrich, ignoreInit = TRUE, ignoreNULL = TRUE)
  #insertUI(
   # selector = "#add",
    #where = "afterEnd",
   # ui = textInput(paste0("txt", input$add),
                 #  "Insert some text")
  #)
  #delete all files
  session$onSessionEnded(function() { unlink(input, recursive = TRUE) })
  }


  #shinyApp(ui = ui, server = server)