library(markdown)
tab5 <- tabPanel("Help",

  fluidRow(

    column(6,
           includeMarkdown("helpPage.Rmd")
    )
  )
)




