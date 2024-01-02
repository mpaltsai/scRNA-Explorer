tab4<- tabPanel("Gene enrichment analysis",
                   sidebarLayout(
                     sidebarPanel(
                       actionButton("runEnrich", "Run enrichment analysis"),
                       textOutput("enrich_completed"),
                       actionButton("tableEnrich", "Show results of gene enrichment analysis")
                     ),
                     
                     mainPanel(
                       
                       
                       
                       helpText(strong("Enrichment results of all correlated genes")),
                       DT::dataTableOutput("dt1_enrich"),
                       HTML("<br/>"),
                       helpText(strong("Enrichment results of positive correlated genes")),
                       DT::dataTableOutput("dt2_enrich"),
                       HTML("<br/>"),
                       helpText(strong("Enrichment results of negative correlated genes")),
                       DT::dataTableOutput("dt3_enrich")
                     )#end of mainPanel
                     
                   )#end of sidebarLayout
  
)