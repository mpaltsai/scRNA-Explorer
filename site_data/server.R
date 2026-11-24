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


options(shiny.maxRequestSize=10000*1024^2)  #max upload file size 10GB

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
  
  # Input readers retain structured records (including original SCE metadata).
  # Only fully reviewed records can supply the existing annotation workflow.
  importer <- input_import_server(input, output, session, function(directory, selected_input) {
    read_app_input(selected_input, directory)
  })
  user_input <- reactive({
    records <- importer$ready()
    req(records)
    list(my_list = import_counts(records), mes = importer$status())
  })

  # Starting/reviewing an import invalidates dependent data, even on failure.
  # State follows the import result rather than a second independent Read event.
  observeEvent(importer$records(), {
    fields <- c("d_to_metadata", "d_to_filter", "filt_samples_list",
      "loaded.dataSO.combined.no.cluster", "loaded.dataSO.combined.markers",
      "loaded.dataSO.combined.markerstop", "gen_cluster_res",
      "loaded.dataSO.combined", "iRNA_cor_tables", "iRNA_enrich_tables")
    for (field in fields) my_react_vals[[field]] <- NULL
    my_react_vals$import_records <- importer$records()
    records <- importer$ready()
    if (!is.null(records)) my_react_vals$d_to_metadata <- import_counts(records)
  }, ignoreNULL = FALSE)

  ##### Gene annotation, QC metrics and per-sample baseline SCEs #####
  annotator <- annotation_qc_server(input, output, session, importer$ready)
  calc_metadata <- reactive({
    p <- annotator$prepared(); req(p)
    list(my_list = p, mes = annotator$status())
  })
  calc_se <- reactive({
    b <- annotator$baseline(); req(b)
    list(my_list = lapply(b, function(object) list(se = object)),
      mes = annotator$baseline_status())
  })

  # A changed reference/baseline invalidates active downstream inputs. The
  # baseline SCE remains separate from later user-filtered SCEs.
  observeEvent(annotator$baseline(), {
    b <- annotator$baseline()
    fields <- c("loaded.dataSO.combined.no.cluster",
      "loaded.dataSO.combined.markers", "loaded.dataSO.combined.markerstop",
      "gen_cluster_res", "loaded.dataSO.combined", "iRNA_cor_tables", "iRNA_enrich_tables")
    for (field in fields) my_react_vals[[field]] <- NULL
    my_react_vals$multiple_case <- !is.null(b) && length(b) > 1L
  }, ignoreNULL = FALSE)

    #save and download as rds
    output$initDownload <- 
      
      downloadHandler(
        filename = function() {
          paste0("init_filt_data_",  Sys.Date(), ".rds")
        },
        content = function(filename) {
          if(length(calc_se()$my_list)>1L){
            objects<-calc_se()$my_list
            for(n in names(objects)) objects[[n]]$se<-optional_qc$decorate(objects[[n]]$se,n)
            saveRDS(objects, filename)
          }else{
            n<-names(calc_se()$my_list)[1]
            saveRDS(optional_qc$decorate(calc_se()$my_list[[1]]$se,n), filename)
          }
          #saveRDS(calc_se()$my_list, file)
        }
      )
  
    
  ##### Per-sample QC review: immutable baseline -> Preview -> Apply #####
  qc_review <- qc_review_server(input, output, session, annotator$baseline,
    optional_results=function(sample) optional_qc$results(sample))
  optional_qc <- qc_diagnostics_server(input,output,session,annotator$baseline,qc_review)
  my_sample_name <- qc_review$selected
  # Compatibility readers now follow the selected sample, rather than caching
  # metrics behind separate button events. These contain tables only.
  calc_metrics <- reactive({
    m <- qc_review$metrics(); req(m)
    c(as.list(m[.qc_metric_names]),list(mes=paste("Metrics ready for",my_sample_name())))
  })
  se_c <- reactive({
    n <- my_sample_name(); req(n)
    object <- qc_review$filtered()[[n]]; req(object)
    list(se=object,mes="Applied QC-filtered data")
  })
  metrics_clean <- reactive({
    object <- se_c()$se
    c(as.list(baseline_qc_metrics(object)),list(mes=paste("Applied metrics for",my_sample_name())))
  })
  observeEvent(qc_review$revisions(), {
    # Preview/draft edits leave revisions unchanged. An Apply replaces only
    # its sample and invalidates every dependent combined result.
    revisions <- qc_review$revisions()
    if (!identical(revisions,isolate(my_react_vals$qc_applied_revisions))) {
      my_react_vals$filt_samples_list <- qc_review$filtered()
      my_react_vals$qc_applied_revisions <- revisions
      for(field in c("loaded.dataSO.combined.no.cluster","loaded.dataSO.combined.markers",
        "loaded.dataSO.combined.markerstop","gen_cluster_res","loaded.dataSO.combined",
        "iRNA_cor_tables","iRNA_enrich_tables")) my_react_vals[[field]] <- NULL
    }
  }, ignoreNULL=FALSE, priority=90)

  ##### Cell clustering #####
  
  # All samples need an applied, nonempty QC revision before normalization.
  observeEvent(input$normData,{
    inputs <- tryCatch(qc_review$require_ready(),error=function(e) {
      validate(need(FALSE,conditionMessage(e)))
    })
    revisions <- qc_review$revisions()
    # Refresh optional evidence without refitting on a retained population.
    for(n in names(inputs)) inputs[[n]]<-optional_qc$decorate(inputs[[n]],n)
    spsComps::shinyCatch({
      object <- if(length(inputs)>1L) norm_data_and_pca_multiple(inputs) else norm_data_and_pca(inputs)
      object@misc$scrna_explorer_qc <- list(filtered_sample_revisions=revisions)
      my_react_vals$loaded.dataSO.combined.no.cluster <- object
    }, blocking_level="error", prefix="scRNA-Explorer")
  })

  output$normalizeData <- renderText({
    ready <- tryCatch({qc_review$require_ready();TRUE},error=function(e) FALSE)
    if(!ready) return("Apply nonempty QC filtering to every sample before normalization.")
    object <- my_react_vals$loaded.dataSO.combined.no.cluster
    if(is.null(object)) return("Current QC revisions have not been normalized. Run normalization/PCA.")
    "Normalization and PCA completed for the recorded applied QC revisions."
  })

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
    validate(need(tryCatch({qc_review$require_ready();TRUE},error=function(e) FALSE), message="Apply nonempty QC filtering to every sample first."))
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
    validate(need(tryCatch({qc_review$require_ready();TRUE},error=function(e) FALSE), message="Apply nonempty QC filtering to every sample first."))
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
      paste0("Clusters identified n=", length(levels(my_react_vals$loaded.dataSO.combined.no.cluster@meta.data$seurat_clusters)))
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
    validate(need(tryCatch({qc_review$require_ready();TRUE},error=function(e) FALSE), message="Apply nonempty QC filtering to every sample first."))
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
    validate(need(tryCatch({qc_review$require_ready();TRUE},error=function(e) FALSE), message="Apply nonempty QC filtering to every sample first."))
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
  })%>% bindEvent(input$corData, ignoreInit = TRUE, ignoreNULL = TRUE)
  
    #observe({
    #observeEvent(input$targetGenesSampleFiltered,{
      observeEvent(
        list(
          input$targetGenesSampleFiltered,
          my_react_vals$filt_samples_list
        ), {
      validate(need(my_react_vals$filt_samples_list, message = FALSE))
      validate(need(input$targetGenesSampleFiltered, message = FALSE))
      updateSelectizeInput(session, "targetGenesFiltered",
                           label = "'Bait' gene to run correlation analysis with",
                           choices =  array(my_react_vals$filt_samples_list[[input$targetGenesSampleFiltered]]@assays@data@listData$counts@Dimnames[1][[1]]), server = TRUE)
      #my_react_vals$cor_counts <- my_react_vals$filt_samples_list[[input$targetGenesSampleFiltered]]$se@assays@data@listData$counts
    })
    observeEvent(list(
      input$targetGenesSampleUnfiltered,
      calc_se()
    ), {
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
