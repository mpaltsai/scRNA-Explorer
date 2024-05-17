norm_data_and_pca_multiple <- function(my_filt_samples_list){

  #filtered data
  #my_react_vals$filt_samples_list
  list_samples <- my_filt_samples_list
  
    my_keys<-c()
    for(x in names(list_samples)){
      my_key <- x
      my_keys <- c(my_keys, my_key)
    }
    samples.list <- vector(mode="list", length=length(my_keys))
    names(samples.list) <- my_keys
    


    first_sample_name = names(samples.list)[1]
    
    first_sample <- CreateSeuratObject(counts = as.matrix(my_filt_samples_list[[first_sample_name]]@assays@data@listData$counts)[which(!duplicated(row.names(my_filt_samples_list[[first_sample_name]]@assays@data@listData$counts))),], project = first_sample_name)
    
    second_sample_name <- names(samples.list)[2]
    second_sample <- CreateSeuratObject(counts = as.matrix(my_filt_samples_list[[second_sample_name]]@assays@data@listData$counts)[which(!duplicated(row.names(my_filt_samples_list[[second_sample_name]]@assays@data@listData$counts))),], project = second_sample_name)

  merged_data <- merge(x = first_sample, y = second_sample)
  
  if(length(samples.list)>2){
    
    i=3
    for(x in samples.list[3:length(samples.list)]){
      sample_name = names(samples.list)[i]
      my_sample <- CreateSeuratObject(counts = as.matrix(my_filt_samples_list[[sample_name]]@assays@data@listData$counts)[which(!duplicated(row.names(my_filt_samples_list[[sample_name]]@assays@data@listData$counts))),], project = sample_name)
      merged_data <- merge(x = my_sample, y = merged_data)
      i = i+1
    }
    
  }else{
    next
  }
  
  merged_data <- NormalizeData(merged_data)
  merged_data <- FindVariableFeatures(merged_data)
  merged_data <- ScaleData(merged_data)
  merged_data <- RunPCA(merged_data, pcs = 30, verbose = FALSE)
  #merged_data <- RunUMAP(merged_data, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap")
  #ElbowPlot(merged_data, ndims = 30)
  return(merged_data)
}

norm_data_and_pca <- function(my_filt_samples_list){
  sample_name <- names(my_filt_samples_list)[1]
  my_sample <- CreateSeuratObject(counts = as.matrix(my_filt_samples_list[[1]]@assays@data@listData$counts)[which(!duplicated(row.names(my_filt_samples_list[[1]]@assays@data@listData$counts))),], project = sample_name)
  my_sample <- NormalizeData(my_sample)
  my_sample <- FindVariableFeatures(my_sample)
  my_sample <- ScaleData(my_sample)
  my_sample <- RunPCA(my_sample, pcs = 30, verbose = FALSE)
  #my_sample <- RunUMAP(my_react_vals$loaded.dataSO.combined.no.cluster, reduction = "pca", dims = 1:30, reduction.name = "umap")
  #ElbowPlot(merged_data, ndims = 30)
  return(my_sample)
}
