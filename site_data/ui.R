#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyjs)
library(shinycssloaders)

# Define UI for application that draws a histogram
#fluidPage(

    # Application title
#    titlePanel("Old Faithful Geyser Data"),

    # Sidebar with a slider input for number of bins
#    sidebarLayout(
#        sidebarPanel(
#            sliderInput("bins",
#                        "Number of bins:",
#                        min = 1,
#                       max = 50,
#                        value = 30)
#        ),

        # Show a plot of the generated distribution
#        mainPanel(
#            plotOutput("distPlot")
#        )
#    )
#)

source("tab_0.R", local=TRUE)
source("tab_1.R", local=TRUE)
source("tab_2dev.R", local=TRUE)
source("tab_2_cells.R", local=TRUE)
source("tab_3.R", local=TRUE)
source("tab_4.R", local=TRUE)

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
    #tabPanel("Gene correlation across cell types"),
    #tabPanel("Gene enrichment analysis"),
    #tabPanel("Gene set comparison")
    )
)  
  
  