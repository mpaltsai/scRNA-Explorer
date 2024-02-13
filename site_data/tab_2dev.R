tab2 <-tabPanel("Quality Control Plots",
  
  #pull metrics from the rds object
  fluidRow(
    column(4,
           wellPanel(
             #QC plots
             #create metrics
             helpText(strong("First of all we have to create the metrics for visualizations")),
             actionButton("makeMetrics", 'Create metrics'),
             span(textOutput("createdMetrics") %>% withSpinner(), style = "color:dodgerblue"),
           )       
    ),
    
    column(4,
           wellPanel(
           helpText(strong("Pre-filtering"))
           #textOutput("NBinSample")
    )),
    column(4,
           wellPanel(
             helpText(strong("Post-filtering"))
           )
  )),

  #UMIs per cell
  fluidRow(
    column(4,
           wellPanel(
             #make number of UMIs/transcripts plot
             helpText("Visualize UMI counts (transcripts) per cell"),
             actionButton("makeNBUMIs", 'Plot UMIs/transcripts per cell'),
             conditionalPanel(
               condition="input.makeNBUMIs>0",
               helpText("You can slide the vertical lines on the plot by changing the values below"),
               #uiOutput("myControl")
               sliderInput("nUMIs", "Number of UMIs:",
                          min = 10, max = 1000, value = c(10,1000), step = 10)
           )
      
    )),
    column(4,
           wellPanel(
           #number of UMIs/transcripts per cell
           plotOutput("UMIsPlot")
           )
      
    ),
    column(4,
           wellPanel(
             plotOutput("UMIsPlotFiltered")
           )
      
    )
  
  
  ),
  
  #genes per cell histogram
  fluidRow(
    column(4,
           wellPanel(
             #make number of genes HistPlot
             helpText("Visualize the distribution of genes detected per cell via histogram"),
             actionButton("makeNBgenesHist", 'Plot genes detected per cell in a histogram'),
             conditionalPanel(
               condition = "input.makeNBgenesHist>0",
               helpText("You can slide the vertical lines on the plot by changing the values below"),
               sliderInput("nGenes", "Number of genes:",
                           min = 10, max = 1000, value = c(10,1000), step = 10)
             )
             
           )),
    column(4,
           wellPanel(
             #number of genes per cell histogram
             plotOutput("GenesHistPlot")
           )
        ),
    column(4,
           wellPanel(
           plotOutput("GenesHistPlotFiltered")))
      ),
  
  #genes per cell boxplot
  fluidRow(
    column(4,
           wellPanel(
             #make number of genes BoxPlot
             helpText("Visualize the distribution of genes detected per cell via boxplot"),
             actionButton("makeNBgenesBoxPlot", 'Plot genes detected per cell in a boxplot')
             )
             
           ),
    column(4,
           wellPanel(
             #number of genes per cell boxplot
             plotOutput("GenesBoxPlot")
           )
    ),
    column(4,
           wellPanel(
             plotOutput("GenesBoxPlotFiltered")
           ))
  ),

#UMIs vs Genes with mitoRatio
  fluidRow(
    column(4,
           wellPanel(
             helpText("Visualize the correlation between genes detected and number of UMIs colored with mitochondrial ratio for each cell"),
             actionButton("makeUMIsGenes", 'Plot UMIs vs. genes detected'),
             conditionalPanel(
               condition="input.makeUMIsGenes>0",
               helpText("You can slide the horizontal lines on the plot by changing the values below"),
               sliderInput("nUMIsVSnGenes", "Number of genes:",
                           min = 10, max = 1000, value = c(10,1000), step = 10)
           )
           
    )),
  column(4,
         wellPanel(
           #number of UMIs vs Genes per cell 
           plotOutput("UMIsGenesPlot")
         )
  ),
  column(4,
         wellPanel(
           plotOutput("UMIsGenesPlotFiltered")
         ))
  ),

  #UMIs vs Genes with rbcRatio
  fluidRow(
    column(4,
           wellPanel(
             helpText("Visualize the correlation between genes detected and number of UMIs colored with hemoglobin ratio for each cell"),
             actionButton("makeUMIsGenesRBC", 'Plot UMIs vs. genes detected'),
             conditionalPanel(
               condition="input.makeUMIsGenesRBC>0",
               helpText("You can slide the horizontal lines on the plot by changing the values below"),
               sliderInput("nUMIsVSnGenesRBC", "Number of genes:",
                           min = 10, max = 1000, value = c(10,1000), step = 10)
             )
             
           )),
    column(4,
           wellPanel(
             #number of UMIs vs Genes per cell with rbcRatio
             plotOutput("UMIsGenesPlotRBC")
           )
    ),
    column(4,
           wellPanel(
             plotOutput("UMIsGenesPlotRBCFiltered")
           ))
  ),
  
  #Mitochondrial counts ratio
  fluidRow(
    column(4,
           wellPanel(
             helpText("Visualize the distribution of mitochondrial gene expression detected per cell"),
             actionButton("makemitoRatio", 'Plot mitochondrial ratio per cell'),
             conditionalPanel(
               condition="input.makemitoRatio>0",
               helpText("You can slide the vertical lines on the plot by changing the values below"),
               sliderInput("mitoRatio", "Mitochondrial ratio:",
                           min = 0, max = 1, value = c(0,1), step = 0.0001)
             )
             
           )),
    column(4,
           wellPanel(
             #number of UMIs vs Genes per cell with rbcRatio
             plotOutput("mitoRatioPlot")
           )
    ),
    column(4,
           wellPanel(
             plotOutput("mitoRatioPlotFiltered")
           ))
  ),

  #Hemoglobin counts ratio
  fluidRow(
    column(4,
           wellPanel(
             helpText("Visualize the distribution of hemoglobins expression detected per cell"),
             actionButton("makerbcRatio", 'Plot hemoglobins ratio per cell'),
             conditionalPanel(
               condition="input.makerbcRatio>0",
               helpText("You can slide the vertical lines on the plot by changing the values below"),
               sliderInput("rbcRatio", "Hemoglobins ratio:",
                           min = 0, max = 1, value = c(0,1), step = 0.0001)
             )
             
           )),
    column(4,
           wellPanel(
             #number of UMIs vs Genes per cell with rbcRatio
             plotOutput("rbcRatioPlot")
           )
    ),
    column(4,
           wellPanel(
             plotOutput("rbcRatioPlotFiltered")
           ))
  ),

  #Novelty
  fluidRow(
    column(4,
           wellPanel(
             helpText("Visualize the overall novelty of the gene expression by visualizing the genes detected per UMI"),
             actionButton("makeNovelty", 'Plot genes detected per UMI'),
             conditionalPanel(
               condition="input.makeNovelty>0",
               helpText("You can slide the vertical lines on the plot by changing the value below"),
               sliderInput("novelty", "Genes detected per UMI:",
                           min = 0, max = 1, value = 0, step = 0.01)
             )
             
           )),
    column(4,
           wellPanel(
             #number of UMIs vs Genes per cell with rbcRatio
             plotOutput("noveltyPlot")
           )
    ),
    column(4,
          wellPanel(
            plotOutput("noveltyPlotFiltered")
          ))
  ),
  #Plot an overview of expression for each cell
  fluidRow(
    column(4,
           wellPanel(
           helpText("Plot the relative proportion of the library size that is accounted for by the 300 most highly expressed features for each cell"),
           actionButton("makeLibrary", "Plot overview of expression per cell")
            )
           ),
    column(4,
           wellPanel(
             plotOutput("libraryPlot")%>% withSpinner()
             )
           
           ),
    column(4,
           wellPanel(
             plotOutput("libraryPlotFiltered")%>% withSpinner()
           ))
  ),

  #Plot highest expressed genes
  fluidRow(
    column(4,
           wellPanel(
             helpText("Plot highest expressed genes"),
             actionButton("makeHighExpr", "Plot highest expressed genes")
             
           )),
    column(4,
           wellPanel(
             actionButton("zoomhighExprPlot", "Zoom in"),
             plotOutput("highExprPlot")%>% withSpinner()
             
           )),
    column(4,
           wellPanel(
             actionButton("zoomhighExprPlotFiltered", "Zoom in"),
             plotOutput("highExprPlotFiltered")%>% withSpinner()
           ))
  ),

  #Filtering thresholds

  fluidRow(
    column(4,
           wellPanel(
             helpText(strong("Now we'll set thresholds to filter out low quality cells")),
             helpText("In case you want to use thresholds defined by the sliders of each plot above please press the button below and values will be updated"),
             helpText("Note: Min and max values for UMIs and genes are taken from the first and second plot respectively. 
                      Mitochondrial ratio min and max values from the mitochondrial gene expression plot, whereas hemoglobins ratio from the hemoglobins gene expression plot.
                      log10GenesPerUMI minimum value is taken from the genes detected per UMI plot."),
             actionButton("updateThresholds", "Set thresholds from sliders"),
             helpText("Else, provide the desired thresholds from the boxes below:"),
             numericInput("minUMIsPerCell", "Min number of UMIs detected per cell:", 1000, min = 0, max = as.integer(.Machine$integer.max)),
             numericInput("maxUMIsPerCell", "Max number of UMIs detected per cell:", 1000, min = 0, max = as.integer(.Machine$integer.max)),
             numericInput("minGenesPerCell", "Min number of genes detected per cell:", 1000, min = 0, max = as.integer(.Machine$integer.max)),
             numericInput("maxGenesPerCell", "Max number of genes detected per cell:", 1000, min = 0, max = as.integer(.Machine$integer.max)),
             numericInput("minMitoPerCell", "Min mitochondrial ratio detected per cell:", 1000, min = 0, max = as.integer(.Machine$integer.max), step = 0.0001),
             numericInput("maxMitoPerCell", "Max mitochondrial ratio detected per cell:", 1000, min = 0, max = as.integer(.Machine$integer.max), step = 0.0001),
             numericInput("minRbcPerCell", "Min hemoglobin ratio detected per cell:", 1000, min = 0, max = as.integer(.Machine$integer.max), step = 0.0001),
             numericInput("maxRbcPerCell", "Max hemoglobin ratio detected per cell:", 1000, min = 0, max = as.integer(.Machine$integer.max), step = 0.0001),
             numericInput("noveltyPerCell", "Min log10GenesPerUMI per cell:", 1000, min = 0, max = as.integer(.Machine$integer.max), step = 0.01)
             
           )),
    column(4,
           wellPanel(
             helpText(strong("Filtering thresholds")),
             verbatimTextOutput("minUMIsPerCellThreshold"),
             verbatimTextOutput("maxUMIsPerCellThreshold"),
             verbatimTextOutput("minGenesPerCellThreshold"),
             verbatimTextOutput("maxGenesPerCellThreshold"),
             verbatimTextOutput("minMitoPerCellThreshold"),
             verbatimTextOutput("maxMitoPerCellThreshold"),
             verbatimTextOutput("minRbcPerCellThreshold"),
             verbatimTextOutput("maxRbcPerCellThreshold"),
             verbatimTextOutput("noveltyPerCellThreshold")
                ),
             
           wellPanel(
             helpText(strong("Now we can filter out low quality reads using selected thresholds")),
             actionButton("QCFiltering", "Filter cells"),
             span(textOutput("QCDataFilt") %>% withSpinner(), style = "color:dodgerblue"),
             
             #save data as rds object
             helpText(strong("You can save filtered data with metadata in a Rds object")),
             #textInput(inputId="dataNameQC", "Please provide a name for your dataset", 
                    #   value = "Enter a name..."),
             helpText(strong("You can download the Rds object")),
             downloadButton("QCDownload", 'Save Rds object')

             
           ),
           wellPanel(
             helpText(strong("Now we have filtered out cells we can redraw plots to assess our thresholds")),
             helpText(strong("But first we have to re-compute our metrics")),
             actionButton("filt_metrics", "Metrics of filtered data"),
             span(textOutput("metricsFilt") %>% withSpinner(), style = "color:dodgerblue"),
             actionButton("redrawPlots", "Plots with filtered cells")
           ))
           
  )



)#end of tabPanel