initial_filt <- function(list_samples, a){
  multiple.samples = a
  
  if(multiple.samples){
    my_keys<-c()
    for(x in names(list_samples)){
      my_key <- x
      my_keys <- c(my_keys, my_key)
    }
    samples.list <- vector(mode="list", length=length(my_keys))
    names(samples.list) <- my_keys
  }else{
    samples.list <- list()
  }
    
  i=1
  for(metaData in list_samples){
    
    # Keep cells with nUMI greater than 100
    idx <- which(metaData$nUMI > 100)
    # Extract the counts for those cells
    counts_c <- metaData$counts[, idx] 
    # Extract the metadata for those cells
    metadata_c <- lapply(metaData[1:9], "[", idx) #the list now has 9 elements not 10
    #metadata_c <- lapply(metaData[1:10], "[", idx)
    # Save data to single cell experiment variable
    
    if(multiple.samples){
      samples.list[[names(list_samples)[i]]] = c(se = SingleCellExperiment::SingleCellExperiment(assays=list(counts=counts_c), 
                                                                                                 colData = metadata_c))
      i = i+1
      
    }else{
      samples.list <- list(c(se = SingleCellExperiment::SingleCellExperiment(assays=list(counts=counts_c), 
                                                                             colData = metadata_c)))
      names(samples.list)[1] <- "sample"

    }
    
  }
  
  return(samples.list)
}
