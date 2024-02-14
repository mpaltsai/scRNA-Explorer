### read data

readData <- function(count.matrix, seurat, inputFile, inputDir){
  input = inputFile
  
  #provide the directory to access data from 10X platform
  input_dir = inputDir
  count.matrix = count.matrix
  seurat = seurat
  if(count.matrix){
    # Import count matrix
    #sm <- read.csv(input, header = TRUE, row.names = 1, sep = "\t")
    
    
    if(seurat){
      
      txi <- readRDS(input)
      #the following doesn't apply anymore (9/2023)
      #txi = txi@assays$RNA@counts
      txi <-txi@assays@data$counts
      counts <- as(as.matrix(txi[,]), "dgCMatrix")
      # Import barcodes
      colnames(counts) <- colnames(txi)
      # Import genes
      rownames(counts) <- rownames(txi)

      return(counts)
      
      
      #txi <- readRDS(input)
      #the following doesn't apply anymore (9/2023)
      #txi = txi@assays$RNA@counts
      #txi <-txi@assays@data$counts
      
      #return(txi)
      
    }else{
      #input = "input_data/scWCP10_S190_matrix.count.csv"
      txi <- as.matrix(data.table::fread(input),rownames=1)
      
      #txi <- read.csv(input, header = TRUE, row.names = 1, sep = "\t")
      # Convert to a sparse matrix for more efficient computation
      counts <- as(txi[,], "dgCMatrix")
      #counts <- mltools::sparsify(txi)
      # Import barcodes
      #cell_ids <- colnames(txi)
      # Import genes
      
      #rownames(counts) <-txi$V1

        return(counts)
      
    }
  }else{
    
    no_barcodes <- FALSE
    
    check.it <- tryCatch( Read10X(data.dir = input_dir), error = function(e) { 
      
      e
      #genes.xcl <<- c(genes.xcl, selectedGene)
      no_barcodes <<- TRUE 
      
      e$message <- message(e, "\nWill assign numeric indices from 1 to the number of cells as barcode names")
    })
    
    if(no_barcodes) { 
      
      #list files in directory
      files <- list.files(path="input_dir", full.names=TRUE, recursive=FALSE)
      for(x in files){
        if(grepl("feature", x)){
          feature.names <- read.delim(x, 
                                      header = FALSE,
                                      stringsAsFactors = FALSE)
        }
        if(grepl("count", x)){
          counts <- readMM(x)
        }
      }
      rownames(counts) = feature.names$V2
      genes <- rownames(counts)
      
      colnames(counts) = c(1:dim(counts)[2]) 
      cell_ids <- colnames(counts)
      return(counts)
    }else{
      
      # Read 10X data from the directory provided in input_dir
      counts <- Read10X(data.dir = input_dir)
      
      # Import barcodes
      cell_ids <- colnames(counts)
      # Import genes
      genes <- rownames(counts)
      return(counts)
    }
  }
}

