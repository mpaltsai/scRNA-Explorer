tab3 <- tabPanel("Gene correlation across all cells",
                 sidebarLayout(
                   sidebarPanel(
                     ### Correlation method ###
                     radioButtons(inputId= "method", label = "Which correlation method do you want to apply?",
                                  choices = c("Pearson"= "pearson", "Kendall" = "kendall", "Spearman" = "spearman"), selected = character(0), inline=TRUE),
                     numericInput(inputId="correlation_threshold", "At which correlation level (correlation coefficient threshold)?", 
                                  value = 0.3, min = 0.01, max = 1, step = 0.01),
                     
                     numericInput(inputId="alpha", "At which level of significance (p value threshold)?", value = 0.05, min = 0.001, max = 1, step = 0.001),
                     radioButtons(inputId= "corData", label = "With which data?",
                                  choices = c("Filtered"= "filtered", "Unfiltered" = "unfiltered"), selected = character(0), inline=TRUE),
                     numericInput(inputId= "geneSparsity", "Define the percentage of cells a gene is expressed (at least 2%)", value = 0.02, min = 0.019, max = 1, step = 0.01),
                     textInput(inputId="targetGenes", label="'Bait' gene to run correlation analysis with", value="Erf"),
                     helpText("Please provide one gene with a gene name ID"),
                     actionButton("runCor", "Run correlation analysis"),
                     textOutput("cor_completed"),
                     
                    helpText(strong("Note")),
                    helpText("Table columns represent the following values:"),
                    helpText("r: correlation coefficient"),
                    helpText("p_r: p value of correlation coefficient"),
                    helpText("p_w: p value of Wilcoxon test"),
                    helpText("percentage in Ea: percentage of cells where the selected gene is NOT expressed"),
                    helpText("percentage in Eb: percentage of cells where the selected gene is expressed"),
                    
                    actionButton("tableCor", "Show results of correlation analysis")
                             
                   ),#end of sidebarPanel
                   
                   mainPanel(
                     
                     helpText(strong("All correlated genes")),
                     DT::dataTableOutput("dt1"),
                     HTML("<br/>"),
                     helpText(strong("Positive correlated genes")),
                     DT::dataTableOutput("dt2"),
                     HTML("<br/>"),
                     helpText(strong("Negative correlated genes")),
                     DT::dataTableOutput("dt3")
                   )#end of mainPanel
                   
                 )#end of sidebarLayout
                 
                 
                 
)