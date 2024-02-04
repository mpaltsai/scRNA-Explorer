#test 
#from server.R line 460
#Define cell type markers
observeEvent(input$cellMarkers, {
  withCallingHandlers({
    shinyjs::html("cellTypeMarkers", "")
    
    loaded.dataSO.combined.markers <<- FindAllMarkers(loaded.dataSO.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    
    
  },
  message = function(m) {
    shinyjs::html(id = "cellTypeMarkers", html = paste0(m$message, '<br>'), add = TRUE)
  })
  #shinyjs::html(id = "cellTypeMarkers", html = paste0("Marker definition completed", '<br>'), add = TRUE)
  
})

#from tab_2_cells.R line 64
#Identify conserved cell type markers
fluidRow(
  column(4,
         wellPanel(
           
           helpText(strong("Since we have defined cell type clusters we'll identify cell type markers")),
           actionButton("cellMarkers", 'Define cell type markers'),
           textOutput("cellTypeMarkers") 
         )       
  ),
  
  column(4,
         #helpText("Pre-filtering"),
         #textOutput("NBinSample")
  ),
  column(4,
         #helpText("Post-filtering")
  )
)