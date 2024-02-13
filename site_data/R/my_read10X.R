my_read10X <- function (mtxPath, genesPath, barcodesPath = NULL, DGEList = FALSE) {

      feature.names <- read.delim(genesPath, 
                                header = FALSE,
                                stringsAsFactors = FALSE)
  

        counts <- readMM(mtxPath)


        rownames(counts) = feature.names$V2
         genes <- rownames(counts)

        colnames(counts) = c(1:dim(counts)[2]) 
         cell_ids <- colnames(counts)
         return(counts)
}