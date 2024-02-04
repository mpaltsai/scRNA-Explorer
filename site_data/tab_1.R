tab1<-     tabPanel("Data input/preprocessing", 
                    sidebarLayout(
                      sidebarPanel(
                        
                        #Test dataset or file input
                        radioButtons(inputId= "testData", label = "Do you want to work with a test dataset?",
                                     choices = c("Yes"= "1", "No, I want to upload a dataset" = 0), selected = character(0), inline=TRUE),
                        conditionalPanel(
                          condition = "input.testData== 0",
                          ### Data upload ###
                          radioButtons(inputId= "typeData", label = "What type is your data?",
                                       choices = c("count matrix"= "1", "10X" = 2, "Rds object" = 3), selected = character(0), inline=TRUE),
                          helpText("Data in a count matrix format: a table with cell-barcodes as columns and features as rows."),
                          helpText("Data in 10X format: a directory with 3 separate files for counts, features and barcodes.")
                        ),
                        
                        conditionalPanel(
                          condition= "input.testData == 1",
                          #load test data
                          helpText(strong("We'll work with test data")) 
                          ),
                        ### Data upload ###
                        #radioButtons(inputId= "typeData", label = "What type is your data?",
                        #             choices = c("count matrix"= "1", "10X" = 2, "Rds object" = 3), selected = character(0), inline=TRUE),
                        
                        conditionalPanel(
                          condition = "input.typeData == 1",
                          fileInput(inputId="csvFile", NULL, buttonLabel = "Upload .csv",
                                    multiple = FALSE, accept = ".csv")),
                        textOutput("uploadCsvValidation"),
                        
                        conditionalPanel(
                          condition = "input.typeData == 2",
                          textInput(inputId="dir", "Please provide the absolute path of the directory containing the counts, features and barcodes files", 
                                    value = "Enter directory path...")),
                        
                        
                        conditionalPanel(
                          condition = "input.typeData == 3",
                          fileInput(inputId="rdsFile", NULL, buttonLabel = "Upload .rds",
                                    multiple = FALSE, accept = ".rds")),
                        textOutput("uploadRdsValidation"),
                        
                        
                        #helpText("Data in a count matrix format: a table with cell-barcodes as columns and features as rows."),
                        #helpText("Data in 10X format: a directory with 3 separate files for counts, features and barcodes."),
                        
                        #check user input
                        actionButton("checkInput", "Check your input"),
                        
                        span(textOutput("uploadSummary"), style="color:dodgerblue"),
                        
                        #set organism and gene annotation parameteres
                        radioButtons(inputId= "typeOrganism", label = "Which organism comes your data from?",
                                     choices = c("Homo sapiens", "Mus musculus" ), selected = character(0), inline=TRUE),
                        
                        radioButtons(inputId= "typeGeneId", label = "How are features (genes) annotated?",
                                     choices = c("Gene names" = TRUE, "Ensembl IDs" = FALSE ), selected = character(0), inline=TRUE),
                        
                        ###### QC analysis and  plots #####
                        
                        #textOutput("inputParameters"),
                        
                        #read data
                        actionButton("startQC", "Read data"),
                        
                        span(textOutput("readData") %>% withSpinner(), style = "color:dodgerblue"),
                        
                        #add metadata
                        helpText("We have to add some metadata regarding the genes present in your dataset"),
                        actionButton("addMetadata", 'Add gene annotations'),
                        
                        span(textOutput("paramsAddMetadata") %>% withSpinner(), style = "color:dodgerblue"),
                        
                        #inital minimal filtering
                        helpText("At this stage we'll filter out cells with less than 100 UMIs"),
                        actionButton("initFiltering", 'Minimal filtering'),
                        span(textOutput("initDataFilt") %>% withSpinner(), style = "color:dodgerblue"),
                        
                        #save data as rds object
                        helpText("Now we'll save data with metadata in a Rds object"),
                        #textInput(inputId="dataName", "Please provide a name for your dataset", 
                        #          value = "Enter a name..."),
                        helpText("You can download the Rds object (optional but helpful)"),
                        #actionButton("initDownload", 'Save Rds object'),
                        downloadButton("initDownload", 'Save Rds object')
                        #span(textOutput("dataNameValidation")%>% withSpinner(), style = "color:dodgerblue"),
                        
                        
                        
                      ),
                      
                      
                      mainPanel(
                        
                        
                      )
                    )
)