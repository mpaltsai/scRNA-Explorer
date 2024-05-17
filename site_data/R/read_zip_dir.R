readZip <- function(dirs){
  
  #initialize a list to store samples
  my_keys<-c()
  for(i in 3:length(dirs)){
    my_key <- basename(dirs[i])
    my_keys <- c(my_keys, my_key)
  }

  samples.list <- vector(mode="list", length=length(my_keys))
  names(samples.list) <- my_keys
  
  for(i in 3:length(dirs)){
    data_dir <- dirs[i]
    #check if dir contains seurat raw data
    filesindir<-list.files(data_dir)
    foundrds<-grep("rds",filesindir)
    #or 10X
    foundfalse<-which(c("barcodes.tsv","genes.tsv","matrix.mtx") %in% filesindir==FALSE)
    #or other types of data format
    #foundH5<-grep("h5",filesindir)

    if(length(foundfalse)==0){ #10X
      # Load the dataset
      #loaded.data <- Read10X(data.dir = data_dir)
      project_name<-basename(data_dir)
      print(project_name)
      # Read 10X data from the directory provided in input_dir
      my_counts <- Read10X(data.dir = data_dir)
      
      # Import barcodes
      cell_ids <- colnames(my_counts)
      # Import genes
      genes <- rownames(my_counts)
      # Initialize the Seurat object with the raw (non-normalized data).
      #loaded.dataSO <- CreateSeuratObject(counts = loaded.data, project = project_name, min.cells = 3, min.features = 200)
      #append to list
      #loaded.dataSO.list<-c(loaded.dataSO.list,loaded.dataSO)
      samples.list[[project_name]] = my_counts
    }else if(length(foundrds) !=0){ #rds
     
      txi<- readRDS(paste(data_dir,"/",filesindir[foundrds],sep=""))
      #the following doesn't apply anymore (9/2023)
      #txi = txi@assays$RNA@counts
      txi <-txi@assays@data$counts
      my_counts <- as(as.matrix(txi[,]), "dgCMatrix")
      # Import barcodes
      colnames(my_counts) <- colnames(txi)
      # Import genes
      rownames(my_counts) <- rownames(txi)
      
      # Initialize the Seurat object with the raw (non-normalized data).
      #loaded.dataSO<- readRDS(paste(data_dir,"/",filesindir[foundrds],sep=""))
      project_name<-basename(data_dir)
      #print(project_name)
      #DefaultAssay(loaded.dataSO)<-"RNA"
      #append to list
      #loaded.dataSO.list<-c(loaded.dataSO.list,loaded.dataSO)
      samples.list[[project_name]] = my_counts
    }else{ #count matrix in csv or tsv
      #check if it is a tsv or csv file
      if(tools::file_ext(filesindir)!= "csv" & tools::file_ext(filesindir)!= "tsv"){
        validate(paste0("Invalid file: ", filesindir,"; Please upload a .csv or .tsv file"))}
      
      txi <- as.matrix(data.table::fread(paste0(data_dir,"/",filesindir,sep="")),rownames=1)
      
      #txi <- read.csv(input, header = TRUE, row.names = 1, sep = "\t")
      # Convert to a sparse matrix for more efficient computation
      my_counts <- as(txi[,], "dgCMatrix")
      #append to list
      #loaded.dataSO.list<-c(loaded.dataSO.list,loaded.dataSO)
      project_name<-basename(data_dir)
      samples.list[[project_name]] = my_counts
    }
  }
  return(samples.list)
}
