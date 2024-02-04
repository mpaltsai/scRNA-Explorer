annotateMyClusters <- function(x, a, b, c){
  library(enrichR)
  input_library = x
  organismName = a
  loaded.dataSO.combined.markerstop = b
  loaded.dataSO.combined = c

  if(getOption("enrichR.live")){
    setEnrichrSite("Enrichr") # Human genes
    annotatedclusters<<-c()
    listofresults<<-c()
    #listOfClusterAnnotationImages <<-list()
    if(organismName == "Mus musculus"){
      #for mouse:
      dbstouse <<- c("PanglaoDB_Augmented_2021","CellMarker_Augmented_2021","Tabula_Muris")  
    }else{
      #for human:
      dbstouse <<- c("PanglaoDB_Augmented_2021","CellMarker_Augmented_2021","Tabula_Sapiens")
    }
    #m=2
    #input_library = "3"
    for(m in as.numeric(as.character(unique(loaded.dataSO.combined.markerstop$cluster)))){
    #for(m in c(1:length(as.numeric(as.character(unique(loaded.dataSO.combined.markerstop$cluster)))))){
      message(paste0("Working on cluster ", m+1, "/", length(unique(loaded.dataSO.combined.markerstop$cluster)), "\n"))
      markers_to_search<-loaded.dataSO.combined.markerstop$gene[loaded.dataSO.combined.markerstop$cluster==m]
      if (getOption("enrichR.live")) {
        #wait 1 sec
        Sys.sleep(1)
        enriched <- enrichr(markers_to_search, dbstouse)
        Sys.sleep(1)
      }else{message("Cannot resolve host: maayanlab.cloud. Check your internet connection and try again")}
      #Ismini: provide options for annotation DB e.g. celltype[1]== PanglaoDB_Augmented_2021, ccelltype[2]==CellMarker_Augmented_2021, celltype[3]==Tabula Sapiens or Muris
      celltype<-c(enriched[[1]]$Term[1],enriched[[2]]$Term[1],enriched[[3]]$Term[1])
      #print(paste("Cluster",m,celltype,sep = " "))
      if(input_library == "1"){
        annot_var = 1
      }else if(input_library == "2"){
        annot_var = 2
      }else{
        annot_var = 3
      }
      annotatedclusters<-c(annotatedclusters,toString(make.names(celltype[annot_var])))
      #annotatedclusters<-c(annotatedclusters,toString(celltype[annot_var]))#using to Tabula Sapiens only to annotate
      names(enriched)<-paste(names(enriched),toString(make.names(celltype[annot_var])),sep = "_")
      #names(enriched)<-paste(names(enriched),toString(celltype[annot_var]),sep = "_")
      listofresults<-c(listofresults,enriched)
      #jpeg(file=paste("Annotation",celltype[annot_var],".jpg",sep=""),
       #    width=1200, height=800)

    }
  }else{message("Cannot resolve host: maayanlab.cloud. Check your internet connection and try again")}
  message("Finished annotation")
  
  #Rename indents 
  new.cluster.ids <- make.names(annotatedclusters,unique = T)
  #We have to rename the names(listofresults) the same way so we can plot them because in case there are two clusters named the same way i.e. "Fibroblast.Heart.CL.0000057" the previous
  #will change the name of the clusters to "Fibroblast.Heart.CL.0000057" and "Fibroblast.Heart.CL.0000057.1" in loaded.dataSO.combined$celltype but they will stay identical i.e. "Fibroblast.Heart.CL.0000057" for both in listofresults
  names(listofresults) <- make.names(names(listofresults), unique = TRUE)
  
  names(new.cluster.ids) <- levels(loaded.dataSO.combined)
  loaded.dataSO.combined <- RenameIdents(loaded.dataSO.combined, new.cluster.ids)
  #Change the Labels
  loaded.dataSO.combined$celltype.label <- paste(Idents(loaded.dataSO.combined), gsub("[0-9]+","",loaded.dataSO.combined$orig.ident), sep = "_")
  loaded.dataSO.combined$label <- gsub("[0-9]+","",loaded.dataSO.combined$orig.ident)
  
  loaded.dataSO.combined$celltype <- Idents(loaded.dataSO.combined)
  
  message("Renamed clusters. Cluster annotation completed.")
  
  return(list(loaded.dataSO.combined, loaded.dataSO.combined.markers, loaded.dataSO.combined.markerstop, listofresults))
}