createMetadata <- function(countsData, a, b){
  ####
  # Check this:Cannot connect to AnnotationHub server, using 'localHub=TRUE' instead
  #Using 'localHub=TRUE'
  #If offline, please also see BiocManager vignette section on offline use
  #snapshotDate(): 2023-08-28
  #Warning: call dbDisconnect() when finished working with a connection
  ####
  
  organismName = a
  gene.name = b

  if(organismName == "Mus musculus"){
    load("annotations_mmus.Rda")
    if(identical(grep("ENSMUSG",rownames(countsData)), integer(0)) ){
      gene.name = TRUE
    }else{
      gene.name=FALSE
    }
  }else{
    load("annotations_hsap.Rda")
    if(identical(grep("ENSG",rownames(countsData)), integer(0)) ){
      gene.name = TRUE
    }else{
      gene.name=FALSE
    }
  }
  # Connect to AnnotationHub
  #ah <- AnnotationHub(localHub=FALSE)
  # Access the Ensembl database for organism
 
  #ahDb <- query(ah, 
   #             pattern = c(organismName, "EnsDb"), 
    #            ignore.case = TRUE)
  # Acquire the latest annotation files
  #id <- ahDb %>%
  #  mcols() %>%
  #  rownames() %>%
   # tail(n = 1)
 # edb <- ah[[id]]
  
  # Extract gene-level information from database
  #annotations <- genes(edb, 
   #                    return.type = "data.frame")
  ## Select annotations of interest
  #annotations <- annotations %>%
   # dplyr::select(gene_id, gene_name, gene_biotype, seq_name, description, entrezid)
  
  #check for duplicates
  countsData <- countsData[which(!duplicated(rownames(countsData))),]
 
  #check for NAs and "" in feature names and subset matrix
  if(length(which(is.na(rownames(countsData))))!=0 | length(which(rownames(countsData)==""))!=0){
    countsData <- countsData[which(!is.na(rownames(countsData))),]
    countsData <- countsData[which(rownames(countsData)!=""),]
    warning("NAs or empty names present as feature names")
  }
   
  # Extract IDs for mitochondrial genes and hemoglobins
    name_init <- "MT"
 
  if(gene.name){
    mt <- annotations %>% 
      dplyr::filter(seq_name == name_init) %>%
      dplyr::pull(gene_name)
  }else{
    rowNames <-  as.matrix(rownames(countsData) )
    rowNames2 <- annotations %>%
      dplyr::filter(gene_id %in% rowNames) %>%
      dplyr::select(gene_name, gene_id, seq_name)
    colnames(rowNames) = "gene_id"
    
    #xxx = merge(as.data.frame(rowNames), rowNames2, by.x = "gene_id", by.y = "gene_id", sort = FALSE)
    xxx <- merge(rowNames, rowNames2, by = "gene_id", all=TRUE, sort = FALSE)
    subset_counts <- FALSE
    rownames(countsData) <- tryCatch( 
      {
        xxx$gene_name
      },
      error = function(e) {
        subset_counts <- TRUE
      }
    )
    if(subset_counts){
      countsData <- countsData[which(!is.na(xxx$gene_name)),]
      rownames(countsData) <- xxx$gene_name[which(!is.na(xxx$gene_name))]
    }

    
    mt <- annotations %>% 
      dplyr::filter(seq_name == name_init) %>%
      dplyr::pull(gene_name)
  }
  rbc = grep("^hemoglobin", annotations$description, value = TRUE, ignore.case=TRUE)
  #rbc = rbc[-(grep("like", rbc, ignore.case = TRUE))]
  rbc = annotations[which(annotations$description %in% rbc), match("gene_name", colnames(annotations))]
  
  metadata <- data.frame(row.names = colnames(countsData), cells = colnames(countsData), stringsAsFactors = F)
  ## Add number of UMIs for each gene per cell to metadata
  metadata$nUMI <- Matrix::colSums(countsData)
  ## Add number of genes detected per cell to metadata
  metadata$nGene <- Matrix::colSums(countsData > 0)
  ## Add number of UMIs per gene for each cell to metadata
  metadata$log10GenesPerUMI <- log10(metadata$nGene) / log10(metadata$nUMI)
  metadata$sample <- "sample"
  # Number of UMIs assigned to mitochondrial genes
  metadata$mtUMI <- Matrix::colSums(countsData[which(rownames(countsData) %in% mt),], na.rm = T)
  #counts@Dimnames[1]
  # Ensure all NAs receive zero counts
  metadata$mtUMI[is.na(metadata$mtUMI)] <- 0
  
  # Calculate of mitoRatio per cell
  metadata$mitoRatio <- metadata$mtUMI/metadata$nUMI
  #Define Hemoglobins
  metadata$rbcUMI <- Matrix::colSums(countsData[which(rownames(countsData) %in% rbc),], na.rm = TRUE)
  metadata$rbcUMI[is.na(metadata$rbcUMI)] <- 0
  metadata$rbcRatio <- metadata$rbcUMI/metadata$nUMI
  
  return(c(metadata, counts = countsData))
}