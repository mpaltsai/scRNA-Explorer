#observeEvent(input$startQC, {
#  validate(need(input$testData, message=FALSE))
#  validate(need(user_input()$my_list, message = FALSE))
#  if(input$testData==2){
#    
#    output$readData <- renderText({
#      validate(need(input$typeData, message=FALSE))
#      
#      my_datatest <<- user_input()$my_list
#      
#      cond_satisf = c()
#      for(x in user_input()$my_list){
#        my_condition = is.null(base::rownames(x))
#        cond_satisf <- c(cond_satisf, my_condition)
#      }
#      if(any(cond_satisf)==TRUE){
#        paste0("Genes names are absent for counts data in subdirectory ", which(cond_satisf==TRUE), ". Please upload a file with gene names as row names if a count matrix or check your data in any other case")
#      }else{
#        user_input()$mes
#      }
#      
#    })
#  }else if(input$testData==1){
#    output$readData <- renderText({
#      user_input()$mes
#    })
#  }
  #show("readData")
#}, ignoreNULL = FALSE)


#if(input$testData==1){
#  validate(need(user_input()$mes, message = FALSE))
#  output$readData <- renderText({
#   user_input()$mes
 # })
#}else 