my_read10X <- function (mtxPath, genesPath, barcodesPath = NULL, DGEList = FALSE) {
      
      
      
      feature.names <- read.delim(genesPath, 
                                header = FALSE,
                                stringsAsFactors = FALSE)
  

      my_counts <- readMM(mtxPath)


        rownames(my_counts) = feature.names$V2
         genes <- rownames(my_counts)

        colnames(my_counts) = c(1:dim(my_counts)[2]) 
         cell_ids <- colnames(my_counts)
         
         samples.list <- list(c(my_counts))
         names(samples.list)[1] <- "sample"
         return(samples.list)
}