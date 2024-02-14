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


options(shiny.maxRequestSize=10000*1024^2)  #max upload file size 5GB

function(input, output, session){
  
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
  
  #Check if there ism't any file uploaded or state which kind of input the user has provided
    user_input <<- reactive({
      #if(length(input$testData)==0){
       # validate("You haven't chosen your data type yet")
      if(input$testData == 1){
        count.matrix <- TRUE
        seurat <- TRUE
        inputFile <- "raw_se_S190.rds"
        inputDir <- FALSE
        my_counts <- c(counts = readData(count.matrix, seurat, inputFile, inputDir),
                    mes = "Reading data completed. Features are in Ensembl IDs and come from mouse. These parameters are pre-selected, no need to define them below.")
        
      }else if(input$typeData == 1){
        count.matrix <- TRUE
        seurat <- FALSE
        inputFile <- input$csvFile$datapath
        inputDir <- FALSE
        my_counts <- c(counts = readData(count.matrix, seurat, inputFile, inputDir),
                   mes = "Reading data completed ")
       
      }else if(input$typeData==3){
        count.matrix <- TRUE
        seurat <- TRUE
        inputFile <- input$rdsFile$datapath
        inputDir <- FALSE
        my_counts <- c(counts = readData(count.matrix, seurat, inputFile, inputDir),
                    mes = "Reading data completed ")

      }else if(input$typeData==2){
        my_counts <- c(counts = my_read10X(mtxPath = input$mtxInput$datapath, 
                           genesPath = input$genes10XInput$datapath, 
                           barcodesPath = input$barcodes10XInput$datapath),
                    mes = "Reading data completed ")
          
        ###i put it inside my_read10X at the end 
        #counts <- as(input_counts, "dgCMatrix")
        
      }
    })%>% bindEvent(input$startQC , ignoreInit = TRUE, ignoreNULL = TRUE)

    output$readData <- renderText({
      user_input()$mes
    })
    
  ##### Preprocessing #####
    
    #add metadata and update counts, sometimes we have to subset the count.matrix because we don't have matches between gene names and ensembl IDs
    calc_metadata <- reactive({
      if(input$testData == 1){#test dataset
        
        c(createMetadata(countsData= user_input()$counts,"Mus musculus", TRUE),
          mes = "Metadata added")
        
      }else if(input$typeData == 1 | input$typeData == 3){#count.matrix or rds
        if(length(input$typeOrganism) == 0 | length(input$typeGeneId) == 0){
          #output$paramsAddMetadata <- renderText(
          message("You have to select the organism and/or gene IDs type")
          #)
        }else{
           c(createMetadata(countsData = user_input()$counts,input$typeOrganism, input$typeGeneId),
                        mes = "Metadata added")
        }
      }else{#10X
        if(length(input$typeOrganism) == 0){
          
          message("You have to select the organism")
         
        }else{
          
          c(createMetadata(countsData = user_input()$counts,input$typeOrganism, TRUE),
            mes = "Metadata added")
          
        }
        
      }
    })%>% bindEvent(input$addMetadata , ignoreInit = TRUE, ignoreNULL = TRUE)
    
    output$paramsAddMetadata <- renderText({
        calc_metadata()$mes
      

    })
    
    #minimal filtering
    calc_se <- reactive({
      metaData <- calc_metadata()
      # Keep cells with nUMI greater than 100
      idx <- which(metaData$nUMI > 100) 
      
      # Extract the counts for those cells
      counts_c <- metaData$counts[, idx] 
      
      # Extract the metadata for those cells
      metadata_c <- lapply(metaData[1:10], "[", idx)  
      
      # Save data to single cell experiment variable
      c(se = SingleCellExperiment::SingleCellExperiment(assays=list(counts=counts_c), 
                                        colData = metadata_c),  
              mes = "Data filtered")
    })%>% bindEvent(input$initFiltering, ignoreInit = TRUE, ignoreNULL = TRUE)
    
    output$initDataFilt<- renderText({
      calc_se()$mes
      
      }) 
    
    #save and download as rds
    output$initDownload <- 
      downloadHandler(
        filename = function() {
          paste0("raw_se_",  Sys.Date(), ".rds")
        },
        content = function(file) {
          saveRDS(calc_se()$se, file)
        }
      )
  
    
  ###### QC routine #####
  
  #create metrics
    calc_metrics <- reactive({
      #need(try(calc_se()$se), "Please make the initial filtering before proceeding")
      my_se<- calc_se()$se
      c(as.data.frame(colData(calc_se()$se)),
        mes = "Metrics created")
    })%>% bindEvent(input$makeMetrics , ignoreInit = TRUE, ignoreNULL = TRUE)
    
  output$createdMetrics <- renderText({
    calc_metrics()$mes

  }) 
  
  #Update slider max values when metrics object is created

  #number of UMIs per cell
  #update slider min and max values
   observe({
    updateSliderInput(session,"nUMIs", "Number of UMIs:",
                      min = 0, max = max(calc_metrics()$nUMI), value = c(min(calc_metrics()$nUMI),max(calc_metrics()$nUMI)), step = 10)
  }) %>% bindEvent(input$makeNBUMIs , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  
  
  output$UMIsPlot <- renderPlot({
      createUMIsPlot(calc_metrics()[1:10], input$nUMIs[1], input$nUMIs[2])

  })%>% bindEvent(c(input$makeNBUMIs, input$nUMIs) , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #nGenes histogram
  observe({
    updateSliderInput(session, "nGenes", "Number of genes:",
                      min = 0, max = max(calc_metrics()$nGene), value = c(min(calc_metrics()$nGene),max(calc_metrics()$nGene)), step = 10)
  }) %>% bindEvent(input$makeNBgenesHist , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  output$GenesHistPlot<- renderPlot({
    createGenesHistPlot(calc_metrics()[1:10], input$nGenes[1], input$nGenes[2])

  })%>% bindEvent(c(input$makeNBgenesHist, input$nGenes) , ignoreInit = TRUE, ignoreNULL = TRUE)

  #nGenes boxplot
  output$GenesBoxPlot <- renderPlot({

    createGenesBoxPlot(calc_metrics()[1:10])
    
  })%>% bindEvent(c(input$makeNBgenesBoxPlot) , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #UMIs vs nGenes with mitoRatio
  observe({

    updateSliderInput(session,"nUMIsVSnGenes", "Genes detected:",
                      min = 0, max = max(calc_metrics()$nGene), value = c(min(calc_metrics()$nGene),max(calc_metrics()$nGene)), step = 10)
  })%>% bindEvent(input$makeUMIsGenes , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  output$UMIsGenesPlot <- renderPlot({
    
    createUMIsGenesPlot(calc_metrics()[1:10], input$nUMIsVSnGenes[1], input$nUMIsVSnGenes[2])
    
  })%>% bindEvent(c(input$makeUMIsGenes, input$nUMIsVSnGenes) , ignoreInit = TRUE, ignoreNULL = TRUE)
 
  
  #UMIs vs nGenes with hemoRatio
  observe({
    
    updateSliderInput(session,"nUMIsVSnGenesRBC", "Genes detected:",
                      min = 0, max = max(calc_metrics()$nGene), value = c(min(calc_metrics()$nGene),max(calc_metrics()$nGene)), step = 10)
  })%>% bindEvent(input$makeUMIsGenesRBC , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  output$UMIsGenesPlotRBC <- renderPlot({
    
    createUMIsGenesPlotRBC(calc_metrics()[1:10], input$nUMIsVSnGenesRBC[1], input$nUMIsVSnGenesRBC[2])
      
  })%>% bindEvent(c(input$makeUMIsGenesRBC, input$nUMIsVSnGenesRBC) , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #plot mitoRatio
  observe({
    
    updateSliderInput(session,"mitoRatio", "Mitochondrial ratio:",
                      min = 0, max = round(max(calc_metrics()$mitoRatio), digits = 2), value = c(round(min(calc_metrics()$mitoRatio),digits = 2),round(max(calc_metrics()$mitoRatio), digits = 2)), step = 0.01)
  })%>% bindEvent(input$makemitoRatio , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  output$mitoRatioPlot <- renderPlot({
    
    createmitoRatioPlot(calc_metrics()[1:10], input$mitoRatio[1], input$mitoRatio[2])
    
  })%>% bindEvent(c(input$makemitoRatio, input$mitoRatio) , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #plot hemoRatio
  observe({
    
    updateSliderInput(session,"rbcRatio", "Hemoglobins ratio:",
                      min = 0, max = round(max(calc_metrics()$rbcRatio), digits = 2), value = c(round(min(calc_metrics()$rbcRatio),digits = 2),round(max(calc_metrics()$rbcRatio), digits = 2)), step = 0.01)
  })%>% bindEvent(input$makerbcRatio , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  output$rbcRatioPlot <- renderPlot({
    
    createrbcRatioPlot(calc_metrics()[1:10], input$rbcRatio[1], input$rbcRatio[2])
    
  })%>% bindEvent(c(input$makerbcRatio, input$rbcRatio) , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #Genes vs UMIs
  observe({
    
    updateSliderInput(session,"novelty", "Genes detected per UMI:",
                      min = 0, max = round(max(calc_metrics()$log10GenesPerUMI), digits = 2), value = round(min(calc_metrics()$log10GenesPerUMI),digits = 2), step = 0.01)
  })%>% bindEvent(input$makeNovelty , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  output$noveltyPlot <- renderPlot({
    
    createnoveltyPlot(calc_metrics()[1:10], input$novelty[1], input$novelty[2])
  })%>% bindEvent(input$makeNovelty, input$novelty, ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #plot scaterPlot
  scater_plot <- reactive({
    createlibraryPlot(calc_se()$se)
  })%>%bindEvent(input$makeLibrary, ignoreInit = TRUE, ignoreNULL = TRUE)
  
  
  output$libraryPlot <- renderPlot({
    scater_plot()
  })
  
  #plot highest expressed features
  highest_expr_plot <- reactive({
    createhighExprPlot(calc_se()$se)
  })%>%bindEvent(input$makeHighExpr, ignoreInit = TRUE, ignoreNULL = TRUE)
  
  output$highExprPlot <- renderPlot({
    highest_expr_plot()
  })
  
  #zoom on plot
  observeEvent(input$zoomhighExprPlot, {
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
    metrics <- as.data.frame(calc_metrics()[1:10])
    filt_cells <-    metrics %>%
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
    my_se_c <- calc_se()$se
    #return filtered data
   c(se = my_se_c[ ,filt_cells],
     mes= "Data filtered")
    
  })%>% bindEvent(input$QCFiltering , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  output$QCDataFilt <- reactive({
    se_c()$mes
  })
  # Save subset to new metrics variable
  metrics_clean <- reactive({
    c(as.data.frame(colData(se_c()$se)),
    mes = "Metrics recomputed")
  })%>% bindEvent(input$filt_metrics, ignoreInit = TRUE, ignoreNULL = TRUE)
  
  output$metricsFilt<- renderText({
    metrics_clean()$mes
  })
  
  #save and download as rds
  output$QCDownload <- 
    downloadHandler(
      filename = function() {
        paste0("filtered_se_",  Sys.Date(), ".rds")
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
    
    createUMIsPlot(metrics_clean()[1:10], min(metrics_clean()$nUMI),max(metrics_clean()$nUMI))
    
  })%>% bindEvent(input$redrawPlots , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #nGenes histogram
  output$GenesHistPlotFiltered <- renderPlot({
    
    createGenesHistPlot(metrics_clean()[1:10], min(metrics_clean()[1:10]$nGene),max(metrics_clean()[1:10]$nGene))
    
  })%>% bindEvent(input$redrawPlots , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #nGenes Boxplot
  output$GenesBoxPlotFiltered <- renderPlot({
    
    createGenesBoxPlot(metrics_clean()[1:10])
    
  })%>% bindEvent(input$redrawPlots , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #UMIs vs genes with mitoRatio
  output$UMIsGenesPlotFiltered <- renderPlot({
    
    createUMIsGenesPlot(metrics_clean()[1:10], min(metrics_clean()[1:10]$nGene), max(metrics_clean()[1:10]$nGene))
    
  })%>% bindEvent(input$redrawPlots , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #UMIs vs genes with rbcRatio
  output$UMIsGenesPlotRBCFiltered <- renderPlot({
    
    createUMIsGenesPlotRBC(metrics_clean()[1:10], min(metrics_clean()[1:10]$nGene), max(metrics_clean()[1:10]$nGene))
    
  })%>% bindEvent(input$redrawPlots , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #mitoRatio
  output$mitoRatioPlotFiltered <- renderPlot({
    
    createmitoRatioPlot(metrics_clean()[1:10], min(metrics_clean()[1:10]$mitoRatio), max(metrics_clean()[1:10]$mitoRatio))
    
  })%>% bindEvent(input$redrawPlots , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #rbcRatio
  output$rbcRatioPlotFiltered <- renderPlot({
    
    createrbcRatioPlot(metrics_clean()[1:10], min(metrics_clean()[1:10]$rbcRatio), max(metrics_clean()[1:10]$rbcRatio))
    
  })%>% bindEvent(input$redrawPlots , ignoreInit = TRUE, ignoreNULL = TRUE)
  
   
  #novelty
  output$noveltyPlotFiltered <- renderPlot({
    
    createnoveltyPlot(metrics_clean()[1:10], min(metrics_clean()[1:10]$log10GenesPerUMI), max(metrics_clean()[1:10]$log10GenesPerUMI))
    
  })%>% bindEvent(input$redrawPlots , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #library plot
  output$libraryPlotFiltered <- renderPlot({

    createlibraryPlot(se_c()$se)
    
  })%>% bindEvent(input$redrawPlots , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #plot highest expressed features
  highest_expr_plot_filtered <- reactive({

    createhighExprPlot(se_c()$se)
  })%>%bindEvent(input$redrawPlots , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  output$highExprPlotFiltered <- renderPlot({
    highest_expr_plot_filtered()
  })
  
  #zoom on plot
  observeEvent(input$zoomhighExprPlotFiltered, {
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
  
  #Normalize data
  output$normalizeData <- reactive({
    
   
    loaded.dataSO <<- tryCatch( 
      {
       CreateSeuratObject(counts = as.matrix(se_c()$se@assays@data@listData$counts)[which(!duplicated(row.names(as.matrix(se_c()$se@assays@data@listData$counts)))),], project = "sample")
      },
      error = function(e) {
        calc_se()$se
        CreateSeuratObject(counts = as.matrix(calc_se()$se@assays@data@listData$counts)[which(!duplicated(row.names(as.matrix(calc_se()$se@assays@data@listData$counts)))),], project = "sample")

      }
    )

    loaded.dataSO <- NormalizeData(loaded.dataSO)
    #Identification of highly variable features (feature selection) *
    loaded.dataSO <- FindVariableFeatures(loaded.dataSO, selection.method = "vst", nfeatures = 2000)
    loaded.dataSO.combined.no.cluster <<-loaded.dataSO
    "Data normalized"
    
  }) %>% bindEvent(input$normData , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #Dimensionality reduction and clustering
  dim_red_data <- reactive({
    loaded.dataSO.combined.no.cluster <<- ScaleData(loaded.dataSO.combined.no.cluster, verbose = FALSE)
    loaded.dataSO.combined.no.cluster <<- RunPCA(loaded.dataSO.combined.no.cluster, npcs = 30, verbose = FALSE)
    # t-SNE and Clustering
    loaded.dataSO.combined.no.cluster <<- RunUMAP(loaded.dataSO.combined.no.cluster, reduction = "pca", dims = 1:20)
    loaded.dataSO.combined.no.cluster <<- FindNeighbors(loaded.dataSO.combined.no.cluster, reduction = "pca", dims = 1:20)
    #loaded.dataSO.combined.no.cluster <<- FindClusters(loaded.dataSO.combined.no.cluster, resolution  = 0.3)#was at 0.3
    FindClusters(loaded.dataSO.combined.no.cluster, resolution  = 0.3)#was at 0.3
  })
  output$DimRedCluster <- reactive({

    loaded.dataSO.combined.no.cluster <<- dim_red_data()
    
    "Clustering completed"
    
  }) %>% bindEvent(input$dimRedClust , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #Define cell type markers
  data.markers <- reactive({
    loaded.dataSO.combined.no.cluster <<- dim_red_data()
    DefaultAssay(loaded.dataSO.combined.no.cluster) <<- "RNA"
    FindAllMarkers(loaded.dataSO.combined.no.cluster, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  })
  
  output$cellTypeMarkers <- reactive({
    loaded.dataSO.combined.markers <<- data.markers()
    if(dim(loaded.dataSO.combined.markers)[1]==0){
      "No clusters to show. Either DE genes not detected and/or couldn't find cells with identity class"
    }else{
    "Marker definition completed"}
  })%>% bindEvent(input$cellMarkers , ignoreInit = TRUE, ignoreNULL = TRUE)
  
    #define top markers
  top_mark_def <- reactive({
    loaded.dataSO.combined.markers <- data.markers()
    loaded.dataSO.combined.markers %>%
      group_by(cluster) %>%
      slice_max(n = input$nTopMarkers, order_by = avg_log2FC) 
   })
  
  output$userTopMarkers <- reactive({
    validate(
      need(try(loaded.dataSO.combined.markers <<- data.markers()), "You have to define cell type markers first")
    )
    loaded.dataSO.combined.markerstop<<- top_mark_def()
    #message printed
    paste0("We'll use top ", input$nTopMarkers, " marker genes for cell type annotation")
    
  }) %>% bindEvent(input$topMarkers , ignoreInit = TRUE, ignoreNULL = TRUE)
  
    #Plot heatmap of top 20 markers 
  output$HeatmapMarkers <- renderPlot({
    validate(
      need(try(loaded.dataSO.combined.markers <<- data.markers()), "You have to define cell type markers first and then choose the number of marker genes before plotting the heatmap")
    )
    validate(
      need(try(loaded.dataSO.combined.no.cluster), "It seems you have skipped some previous steps!")
    )
    loaded.dataSO.combined.markerstop<<- top_mark_def()
    #plot heatmap
    DoHeatmap(loaded.dataSO.combined.no.cluster, features = loaded.dataSO.combined.markerstop$gene[1:20]) + NoLegend()
    
  })%>% bindEvent(input$plotHeatmap , ignoreInit = TRUE, ignoreNULL = TRUE)
 
  #annotate gene clusters
  gen_cluster_res <- reactive({
    loaded.dataSO.combined.markerstop<<- top_mark_def()
    annotateMyClusters(input$annotLibrary, input$typeOrganism, loaded.dataSO.combined.markerstop, loaded.dataSO.combined.no.cluster)
  }) %>% bindEvent(input$annotClusters , ignoreInit = TRUE, ignoreNULL = TRUE)
  
   output$anot_completed <- reactive({
    "Annotation started"
    gen_cluster_res <- gen_cluster_res()
   "Annotation completed"
  }) 
  
  
  #update selectInput for gene clusters 
  observe({
    loaded.dataSO.combined <<- gen_cluster_res()[1][[1]]
    updateSelectInput(session, "selectCluster",
                      label = "Choose a gene cluster",
                      choices =  levels(loaded.dataSO.combined$celltype))
  })%>% bindEvent(input$showClusters , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  
  #Plot annotation results 
  output$anotResults <- renderPlot({
    enriched <<- gen_cluster_res()[4][[1]]
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
  output$vizClusters <- renderPlot({
    loaded.dataSO.combined <<- gen_cluster_res()[1][[1]]
    DimPlot(loaded.dataSO.combined, reduction = "umap", split.by = "label",label = TRUE)
  })%>% bindEvent(input$plotUMAP , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #review expression of top markers
  output$vizTopMarkersExpr <- renderPlot({
    loaded.dataSO.combined.markerstop <<- gen_cluster_res()[3][[1]]
    markers.to.plot <- unique(loaded.dataSO.combined.markerstop$gene) 
    DotPlot(loaded.dataSO.combined, features = markers.to.plot[1:10], cols = c("blue", "red"), dot.scale = 8, split.by = "label") +
      RotatedAxis()
  })%>% bindEvent(input$topMarkersExpr , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  
  #Vizualize markers
  #update dropdown menu of selectizeInput for available marker genes
  observe({
    loaded.dataSO.combined.markerstop <<- gen_cluster_res()[3][[1]]
    updateSelectizeInput(session, "selectMarkerGene",
                      label = "Choose a gene",
                      choices =  unique(loaded.dataSO.combined.markerstop$gene), server =TRUE )
  })%>% bindEvent(input$showMarkerGene , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #Violin plot 
  output$vizGeneViolin <- renderPlot({
    loaded.dataSO.combined <<- gen_cluster_res()[1][[1]]
    VlnPlot(loaded.dataSO.combined, features = input$selectMarkerGene)
  })%>% bindEvent(input$plotViolin , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #Feature plot 
  output$vizGeneFeature <- renderPlot({
    loaded.dataSO.combined <<- gen_cluster_res()[1][[1]]
    FeaturePlot(loaded.dataSO.combined, features = input$selectMarkerGene,label = T)& theme(legend.position = c(0.1,0.2))
  })%>% bindEvent(input$plotFeature , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  #Proportion of cells
  output$vizPropCells <- renderPlot({
    loaded.dataSO.combined <<- gen_cluster_res()[1][[1]]
    proportions<-as.data.frame(prop.table(table(loaded.dataSO.combined$celltype.label, loaded.dataSO.combined$orig.ident), margin = 2))
    colnames(proportions) <- c("CellIDs", "Condition", "Proportion")
    #pie chart with proportions
    ggplot(proportions, aes(x="", y=Proportion, fill=CellIDs)) +
      geom_bar(stat="identity", width=1, color="white") +
      coord_polar("y", start=0) +
      theme_void() # remove background, grid, numeric labels
    
  })%>% bindEvent(input$plotPropCells , ignoreInit = TRUE, ignoreNULL = TRUE)
  
  
  ##### Correlation #####
  
  observe({
    #validate(need(try(keep_index), "It seems you haven't filtered your data yet."))
    #validate(need(try(myReacValues$se), "You'll have to complete the preprocessing steps before running a correlation analysis."))
    #validate(need(try(loaded.dataSO.combined <<- gen_cluster_res()[1][[1]]), "You have to complete gene cluster annotation steps first!"))
    #validate(need(try(loaded.dataSO.combined.no.cluster), "You have to complete gene cluster annotation steps first!"))
    
    if(input$corData=="clustered"){
      updateSelectInput(session, "targetGenesCluster",
                        label = "Choose a gene cluster",
                        choices =  levels(loaded.dataSO.combined$celltype))
      
    }
    updateSelectizeInput(session, "targetGenes",
                         label = "'Bait' gene to run correlation analysis with",
                         choices =  array(calc_se()$se@assays@data@listData$counts@Dimnames[1][[1]]), server = TRUE)
  })%>% bindEvent(input$corData , ignoreInit = TRUE, ignoreNULL = TRUE)
  

  
  
  observeEvent(input$runCor, {
      withCallingHandlers({
        shinyjs::html("cor_completed", "")
    if(length(input$corData) == 0 | length(input$method) == 0){
      message("You have to define which method and/or which data you want to run the analysis with")
    }else{
    if(input$corData=="filtered"){
      iRNA_cor_tables <<- tryCatch( 
        {
          iRNA(se_c()$se, input$method, input$alpha, input$geneSparsity, input$targetGenes, input$correlation_threshold)
        },
        error = function(e) {
         message("No filtered data available. Will run analysis with raw data")
          iRNA(calc_se()$se, input$method, input$alpha, input$geneSparsity, input$targetGenes, input$correlation_threshold)
        }
      )
      
      
    }else if(input$corData=="unfiltered"){
      #"Running analysis with raw data"
      iRNA_cor_tables <<- tryCatch({
        iRNA(calc_se()$se, input$method, input$alpha, input$geneSparsity, input$targetGenes, input$correlation_threshold)
      },
      error = function(e) {
        message("Seems initial filtering at the preprocessing stage is not complete. Please go back and follow all the steps of preprocessing")
      })
      
    }else{
      cat(paste0("Running analysis with cell type cluster:", input$selectCluster))
      iRNA_cor_tables <<- tryCatch({
        iRNA_cluster(gen_cluster_res()[1][[1]], input$method, input$alpha, input$geneSparsity, input$targetGenes, input$correlation_threshold, input$selectCluster)
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

  #delete all files
  session$onSessionEnded(function() { unlink(input, recursive = TRUE) })
  }


  #shinyApp(ui = ui, server = server)