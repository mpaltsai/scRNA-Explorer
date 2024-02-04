tab2_cell_clusters <-tabPanel("Cell type clustering",
                #for zoom in on images from https://stackoverflow.com/questions/62611894/adjust-image-that-is-shown-by-showmodal-on-shiny
                tags$head(tags$style("#shiny-modal img { max-width: 100%; }")),
                
                #Normalize Data
                fluidRow(
                  column(4,
                         wellPanel(
                           
                           helpText(strong("Now we have filtered out low quality cells we can proceed to cell type clustering.")),
                           helpText(strong("Firstly, we'll normalize expression values and look for outlier features")),
                           actionButton("normData", 'Normalize data'),
                           span(textOutput("normalizeData") %>% withSpinner(), style = "color:dodgerblue"),
                         )       
                  ),
                  
                  column(4,
                         #helpText(""),
                         #textOutput("")
                  ),
                  column(4,
                         #helpText("")
                  )
                ),
                
                #Dimensionality reduction
                fluidRow(
                  column(4,
                         wellPanel(
                           
                           helpText(strong("At this step we'll perform a dimensionality reduction procedure to define cell clusters")),
                           actionButton("dimRedClust", 'Dimensionality reduction and clustering'),
                           span(textOutput("DimRedCluster") %>% withSpinner(), style = "color:dodgerblue")
                         )       
                  ),
                  
                  column(4,
                         #helpText("Pre-filtering"),
                         #textOutput("NBinSample")
                  ),
                  column(4,
                         #helpText("Post-filtering")
                  )
                ),
                
                #Identify conserved cell type markers
                fluidRow(
                  column(4,
                         wellPanel(
                           
                           helpText(strong("Since we have defined cell type clusters we'll identify cell type markers")),
                           actionButton("cellMarkers", 'Define cell type markers'),
                           span(textOutput("cellTypeMarkers")%>% withSpinner(), style = "color:dodgerblue") 
                         )       
                  ),
                  
                  column(4,
                         #helpText("Pre-filtering"),
                         #textOutput("NBinSample")
                  ),
                  column(4,
                         #helpText("Post-filtering")
                  )
                ),
                
                #Use of number of markers for annotation
                fluidRow(
                  column(4,
                         wellPanel(
                           
                           helpText(strong("We have to define how many markers we are going to use for cell type annotation at the next step. You can choose from 1 to 200 top marker genes (default to 100). You can return and change your preference in case you don't find the annotation suitable.")),
                           numericInput(inputId= "nTopMarkers", "Define number of genes (up to 200)", value = 100, min = 1, max = 200, step = 1),
                           actionButton("topMarkers", 'Define number of markers'),
                           textOutput("userTopMarkers") 
                         )       
                  ),
                  
                  column(4,
                         #helpText("Pre-filtering"),
                         #textOutput("NBinSample")
                  ),
                  column(4,
                         #helpText("Post-filtering")
                  )
                ),
                
                
                #Plot a heatmap of top 20 genes in each cluster
                fluidRow(
                  column(4,
                         wellPanel(
                           
                           helpText(strong("Plot expression of top 20 genes in a heatmap (scaled data)")),
                           actionButton("plotHeatmap", 'Plot heatmap of 20 top markers')
                         )       
                  ),
                  
                  column(8,
                         wellPanel(
                           plotOutput("HeatmapMarkers")%>% withSpinner()

                  )
                  )
                ),
                
                #Cluster annotation
                fluidRow(
                  column(4,
                         wellPanel(
                           
                           helpText(strong("Now we can annotate our gene clusters. We are going to use the following libraries 
                           (see https://maayanlab.cloud/Enrichr/#libraries for further details). By default, the Tabula Sapiens (human genes) or Tabula Muris (mouse genes) are used 
                           to annotate the clusters ONLY for visualization purposes (i.e. name of clusters on plots). In any case, you can change this behavior by selecting the preferred 
                           library below.")),
                           radioButtons(inputId= "annotLibrary", label = "Which annotation to use on plots?",
                                        choices = c("PanglaoDB_Augmented_2021"= "1", "CellMarker_Augmented_2021" = "2", "Tabula_Sapiens or Tabula_Muris" = "3"), selected = "3", inline=FALSE),
                           actionButton("annotClusters", 'Annotate clusters'),
                           textOutput("anot_completed") %>% withSpinner()

                         )       
                  ),
                  
                  column(4,
                         #helpText("Pre-filtering"),
                         #textOutput("NBinSample")
                  ),
                  column(4,
                         #helpText("Post-filtering")
                  )
                ),
                
                #Plot annotation results
                fluidRow(
                  column(4,
                         wellPanel(
                           
                           helpText(strong("At this stage we can review annotation results")),
                           actionButton("showClusters", 'Press to provide clusters in the list below'),
                           selectInput("selectCluster", "Choose a gene cluster", ""),#choices = levels(loaded.dataSO.combined$celltype)),
                           actionButton("plotAnnot", 'Review annotation')
                           
                         )       
                  ),
                  
                  column(8,
                         wellPanel(
                         plotOutput("anotResults")%>% withSpinner()
                  ),
                  )
                ),
                
                #UMAP plot
                fluidRow(
                  column(4,
                         wellPanel(
                           
                           helpText(strong("We can visualize the clusters with UMAP")),
                           actionButton("plotUMAP", 'UMAP plot')
                           
                         )       
                  ),
                  
                  column(8,
                         wellPanel(
                           
                           plotOutput("vizClusters")%>% withSpinner()
                         ),
                  )
                ),
                
                #review expression of top markers
                fluidRow(
                  column(4,
                         wellPanel(
                           
                           helpText(strong("We can review the expression of the top 10 marker genes across all clusters")),
                           actionButton("topMarkersExpr", 'Dot plot of marker genes expression')
                         )       
                  ),
                  
                  column(8,
                         plotOutput("vizTopMarkersExpr")%>% withSpinner()
                         #textOutput("NBinSample")
                  )
                ),
                
                #Vizualize markers
                fluidRow(
                  column(4,
                         wellPanel(
                           
                           helpText(strong("We can visualize some individual markers")),
                           actionButton("showMarkerGene", 'Press to select a gene from the list below'),
                           selectInput("selectMarkerGene", "Choose a gene", ""),
                           actionButton("plotViolin", 'Violin plot')
                           
                         )       
                  ),
                  
                  column(8,
                         wellPanel(
                           
                           plotOutput("vizGeneViolin")%>% withSpinner()
                         ),
                  )
                ),
                
               #Vizualize markers feature Plot
                fluidRow(
                  column(4,
                         wellPanel(

                           actionButton("plotFeature", 'Feature plot')
                           
                         )       
                  ),
                  
                  column(8,
                         wellPanel(
                           
                           plotOutput("vizGeneFeature")%>% withSpinner()
                         ),
                  )
                ),
                
                
                #What proportion of cells are in each cluster by condition?
                fluidRow(
                  column(4,
                         wellPanel(
                           
                           helpText(strong("Additionally we can see what proportion of cells are in each cluster")),
                           actionButton("plotPropCells", 'Proportion of cells')
                           
                         )       
                  ),
                  
                  column(8,
                         wellPanel(
                           
                           plotOutput("vizPropCells")%>% withSpinner()
                         ),
                  )
                ),
                
                
)#end of tabPanel