library(shiny)
library(shinyjs)
library(shinycssloaders)



source("tab_0.R", local=TRUE)
source("tab_1.R", local=TRUE)
source("tab_2dev.R", local=TRUE)
source("tab_2_cells.R", local=TRUE)
source("tab_3.R", local=TRUE)
source("tab_4.R", local=TRUE)
source("helpPage.R", local=TRUE)

ui <- fluidPage(
  
  titlePanel("scRNA-Explorer"),
  
  tags$head(
    tags$style(HTML("
      .shiny-output-error-validation {
        color: red;
      }
    "))
  ),
  shinyjs::useShinyjs(),
  tabsetPanel(
    tab0,
    tab1,
    tab2,
    tab2_cell_clusters,
    tab3,
    tab4,
    tab5

    )
)  
  
  