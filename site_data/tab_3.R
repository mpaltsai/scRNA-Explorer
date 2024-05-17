tab3 <- tabPanel("Gene correlation",
                 sidebarLayout(
                   sidebarPanel(
                     ### Correlation method ###
                     radioButtons(inputId= "method", label = "Which correlation method do you want to apply?",
                                  choices = c("Pearson"= "pearson", "Kendall" = "kendall", "Spearman" = "spearman"), selected = character(0), inline=TRUE),
                     numericInput(inputId="correlation_threshold", "At which correlation level (correlation coefficient threshold)?", 
                                  value = 0.3, min = 0.01, max = 1, step = 0.01),
                     
                     numericInput(inputId="alpha", "At which level of significance (p value threshold)?", value = 0.05, min = 0.001, max = 1, step = 0.001),
                     radioButtons(inputId= "corData", label = "With which data?",
                                  choices = c("Filtered (across all cells)"= "filtered", "Unfiltered (across all cells)" = "unfiltered", "Cell type cluster" = "clustered"), selected = character(0)),


                     conditionalPanel(
                      condition= "input.corData=='filtered'",
                      helpText("Please select one sample from the dropdown menu"),
                      selectInput(inputId="targetGenesSampleFiltered",
                                  label = "Choose a sample",
                                  choices =  "",
                                  selected = NULL),
                      selectizeInput(inputId="targetGenesFiltered", label="'Bait' gene to run correlation analysis with", choices=NULL),
                      helpText("Note that some genes will be excluded from the subsequent correlation analysis due to limited expression among single cells"),
                      
                    ),
                    
                    conditionalPanel(
                      condition= "input.corData=='unfiltered'",
                      helpText("Please select one sample from the dropdown menu"),
                      selectInput(inputId="targetGenesSampleUnfiltered",
                                  label = "Choose a sample",
                                  choices =  "",
                                  selected = NULL),
                      selectizeInput(inputId="targetGenesUnfiltered", label="'Bait' gene to run correlation analysis with", choices=NULL),
                      helpText("Note that some genes will be excluded from the subsequent correlation analysis due to limited expression among single cells"),
                      
                      
                    ),

                    
                    conditionalPanel(
                      condition= "input.corData=='clustered'",
                      helpText("Please select one sample from the dropdown menu"),
                      selectInput(inputId="targetGenesSampleClustered",
                                  label = "Choose a sample",
                                  choices =  "",
                                  selected = NULL),
                      helpText("Please select one cluster from the dropdown menu"),
                      selectInput(inputId="targetGenesCluster",
                                  label = "Choose a cell cluster",
                                  choices =  "",
                                  selected = NULL),
                      selectizeInput(inputId="targetGenes", label="'Bait' gene to run correlation analysis with", choices=NULL),
                      helpText("Note that some genes will be excluded from the subsequent correlation analysis due to limited expression among single cells"),
                      
                    ),
                     
                     
                     numericInput(inputId= "geneSparsity", "Define the percentage of cells a gene is expressed (at least 2%)", value = 0.02, min = 0.019, max = 1, step = 0.01),
                     actionButton("runCor", "Run correlation analysis"),
                     span(textOutput("cor_completed"), style = "color:dodgerblue"),
                     
                    helpText(strong("Note")),
                    helpText("Table columns represent the following values:"),
                    helpText("r: correlation coefficient"),
                    helpText("p_r: p value of correlation coefficient"),
                    helpText("p_w: p value of Wilcoxon test"),
                    helpText("Ea: cells where the 'bait' gene is NOT expressed"),
                    helpText("percentage in Ea: percentage of cells in Ea where each correlated gene is expressed"),
                    helpText("Eb: cells where the 'bait' gene IS expressed"),
                    helpText("percentage in Eb: percentage of cells in Eb where each correlated gene is expressed")
                             
                   ),#end of sidebarPanel
                   
                   mainPanel(
                     
                     helpText(strong("All correlated genes")),
                     DT::dataTableOutput("dt1"),
                     HTML("<br/>"),
                     helpText(strong("Only Positive correlated genes")),
                     DT::dataTableOutput("dt2"),
                     HTML("<br/>"),
                     helpText(strong("Only Negative correlated genes")),
                     DT::dataTableOutput("dt3")
                   )#end of mainPanel
                   
                 )#end of sidebarLayout
                 
                 
                 
)