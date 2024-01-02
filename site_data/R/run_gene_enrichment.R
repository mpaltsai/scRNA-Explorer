run_gene_enrich <- function(cor.res, organismName){
  #Not for shiny:message("\n Reading data...\n")
  #Not for shiny:cor.res = readRDS(enrich.only)
  #Not for shiny:output.path = gsub("(.*)[/].*", "\\1", enrich.only)


library(gprofiler2)

message("\n Running enrichment analysis\n")

#set parameters
if (organismName=="Homo sapiens"){
  organism_name_gost = "hsapiens"
}else{
  organism_name_gost = "mmusculus"
}
#Not for shiny:#Read cor.res by 3 to follow the routine of enrichment analysis using the selectedGene variable
#Not for shiny:for (i in seq(1, length(cor.res))){
  
  selectedGene = names(cor.res[1])
  message("\n------Working on ", selectedGene, "------\n")
  
  g.sets = lapply(c(cor.res[selectedGene],
                    cor.res[paste0(selectedGene, "_pos.cor")],
                    cor.res[paste0(selectedGene, "_neg.cor")]), rownames)
  
  #gost_S198 = readRDS("gost_iRNA_filtered_S198.rds")
  
  gost_S198 = list()
  
  for (x in 1:length(g.sets)){
    
    if(length(g.sets[x][1][[1]])==0){
      message("\n ", names(g.sets[x]), " is empty")
      gost_S198[names(g.sets[x])] <- matrix(rep(NA,5),ncol=5, 
                                            dimnames = list(NULL, c("p_value","term_id", "source", "term_name", "intersection")))
      next
    }
    message("\n Awaiting a response from gProfiler server...")
    #a helper function to handle errors when calling gost
    gost_ja <- function(x){
      tryCatch(
        expr = {
          gost.res <<-gost(query= g.sets[x], organism = organism_name_gost, domain_scope = "annotated", significant = T, evcodes = TRUE,
                           sources = c("GO", "KEGG", "REAC", "WP", "MIRNA", "HPA", "CORUM", "HP"))
          message("\n Successfully executed server call.")
          
        },
        error = function(e){
          message('\n Error:')
          gost.response <<- TRUE
          
          message(e$message)
        },
        warning = function(w){
          message('\n Warning:')
          message(w)
        },
        finally = {
          message('\n')
        }
        
      ) 
      
    }
    gost.response <- FALSE
    gost_ja(x)
    b=0
    while(gost.response==TRUE && b<5) { 
      
      gost_ja(x)
      
      Sys.sleep(2)
      b= b+1
      
      if(b==5){
        message("\n Enrichment analysis failed. \n")
        gost_S198[names(g.sets[x])] <- matrix(rep(NA,5),ncol=5, 
                                              dimnames = list(NULL, c("p_value","term_id", "source", "term_name", "intersection")))
        stop("\n Exiting...",call. = FALSE)
        
      }
    }
    
    if(is.null(gost.res)){
      message("\n gene set for ", names(g.sets[x]), " didn't return results from enrichment analysis\n")
      gost_S198[names(g.sets[x])] <- matrix(rep(NA,5),ncol=5, 
                                      dimnames = list(NULL, c("p_value","term_id", "source", "term_name", "intersection")))
      next
    }
    
    message("\n enrichment analysis for ", names(g.sets[x])," completed successfully\n")
    gost_S198[names(g.sets[x])] <-list(as.matrix(gost.res$result[,c("p_value","term_id", "source", "term_name", "intersection")]))
    
    #Not for shiny:outputPath          = paste0(output.path,'/', selectedGene, '/')
    #Not for shiny:f.path = paste0(outputPath, paste0(names(g.sets[x]),"_", mitoRatio, "_", dataset.name, '_GO.csv'))
    #append parameters of analysis to file
    #Not for shiny: cat(paste0("paramCorrMethod: ",paramCorrMethod,"\n"), 
    #Not for shiny:   paste0("correlation level: ", c.r, "\n"),
    #Not for shiny:    paste0("level of significance: ", alpha.level, "\n"),
    #Not for shiny:    paste0("paramCorrPAdjust: ",paramCorrPAdjust,"\n"), 
    #Not for shiny:   paste0("paramGeneSparcityThr: ", paramGeneSparcityThr, "\n"), 
    #Not for shiny:   paste0("dataset: ", dataset.name, "\n"), 
        #paste0("percentage of cells where the selected gene is expressed: ", sel.gene.percent, "\n"),
    #Not for shiny:   paste0("mitoRatio: ", mitoRatio), file=f.path)
    #write enrichment analysis results
    #write.csv(gost_S198[genes.incl[x]], 
    #paste0(file.path(getwd(), outputPath), paste0(genes.incl[x],'_GO.csv')), row.names=FALSE)
    #Not for shiny: message("\n Writing results to file\n")
    
    #Not for shiny:table_enrich = gost_S198[names(g.sets[x])][[1]][order(as.numeric(gost_S198[names(g.sets[x])][[1]][,1])),  ,drop = FALSE]
    #gost_S198[names(g.sets[x])] <- table_enrich
    #Not for shiny:suppressWarnings(write.table(table_enrich, f.path, sep=",", append=TRUE, col.names=NA))
    
    #wait a second before calling again gost
    Sys.sleep(1)
    } 
  
  
#} 
return(gost_S198)
stop("\n Enrichment analysis completed")
}