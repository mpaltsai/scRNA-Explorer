iRNA_cluster <- function(sce, method, alpha, geneSparsity, targetGenes, correlation_threshold, cluster){
  message('\nloading libraries, please wait')
  suppressPackageStartupMessages(library(stats))
  suppressPackageStartupMessages(library(psych))
  suppressPackageStartupMessages(library(plotly))
  suppressPackageStartupMessages(library(BioQC))
  
  message("libraries loaded\n")
  #set parameters
  #dataset.name        = opt$dataset.name 
  #mitoRatio           = opt$mitoRatio 
  paramCorrMethod     <- method 
  #paramCorrPAdjust    = opt$pAdjustMethod   
  #paramSeuMinFeatures = 200
  paramGeneSparcityThr <-  geneSparsity 
  alpha.level         <- alpha 
  #enrich              = opt$enrichment
  #gen.set.sign        = opt$signGeneSets
  #enrich.only         = opt$only_enrichment
  c.r                 <- correlation_threshold 
  cluster.type <- cluster

  
  #read data
  #loaded.dataSO.combined@assays$RNA@counts[,loaded.dataSO.combined$celltype==input$selectCluster]
  #seurat = CreateSeuratObject(counts = loaded.dataSO.combined@assays$RNA@counts[,loaded.dataSO.combined$celltype=="Fibroblast.Heart.CL.0000057"], project = "sample", min.cells = paramSeuMinCells)
  #targetGenes=c("Alpl")
  #paramCorrMethod = "pearson"
  paramSeuMinCells    = paramGeneSparcityThr*dim(sce)[2]
  
  #create a seurat object to work with but check first if it is already a seurat object
  seurat = CreateSeuratObject(counts = sce@assays$RNA@counts[,sce$celltype==cluster.type], project = "sample", min.cells = paramSeuMinCells)
  
  
  #discard genes that are not detected in any cell
  #Ismini note: I think is already implemented with paramSeuMinCells
  seurat.sub.nonzero = seurat[rowSums(seurat)!=0]
  
  #discard cells that none gene is detected ### NO NEED, it has been taken care of when the seurat object was created, with min. features argument
  #or with thw subsetting when checking the assay for QC
  seurat.sub.nonzero = seurat.sub.nonzero[,which(colSums(seurat.sub.nonzero)!=0)]
  
  ###insert normalized data in the seurat object
  #after revoving zero-rows and columns
  seurat.sub.nonzero <-  NormalizeData(seurat.sub.nonzero)
  
  #initialize a list to keep correlation analysis results
  cor.res <- list()
  
  #capitalize genes names in the seurat object
  #rownames(x=seurat.sub.nonzero$RNA) <- toupper(rownames(x=seurat.sub.nonzero$RNA))
  seurat.sub.nonzero$RNA@counts@Dimnames[1][[1]] <- toupper(dimnames(seurat.sub.nonzero)[1][[1]])
  seurat.sub.nonzero$RNA@data@Dimnames[1][[1]] <- toupper(dimnames(seurat.sub.nonzero)[1][[1]])
  
  #set target genes
  genes.incl = strsplit(targetGenes, ",")
  genes.incl = genes.incl[[1]]
  genes.incl = toupper(genes.incl)
  
  message("\nAnalysis will run with ", length(genes.incl), " gene(s): ", genes.incl, " against cluster ", cluster.type, "\n")
  
  
  #message(cat(genes.incl, sep = " ", "\n"))
  
  
  #initialize a list to keep genes excluded from analysis
  genes.xcl = c()
  #setdiff(genes.incl, intersect(names(cor.res), genes.incl))
  
  #create a directory to store results of the form ~/directory-running-iRNA/results/timestamp
  #check if the directory exists and if not create it
  #Not for shiny: res.dir = ifelse(!dir.exists(file.path(getwd(), "results")), dir.create(file.path(getwd(), "results")), FALSE)
  #create a sub directory, name it with a timestamp
  #Not for shiny:st=format(Sys.time(), "%Y-%m-%d_%H:%M")
  #Not for shiny:dir.create(file.path(getwd(), "results", st))
  #Not for shiny:the.path = file.path(getwd(), "results", st)
  #Not for shiny:message("\nCreating a directory to store results\n")
  #Not for shiny:message(" ",the.path )
  
  #Not for shiny:start.cor.time = Sys.time()
  
  ####################### Correlation & Wilcoxon tests ################################3
  for( x in 1:length(genes.incl)){
    
    selectedGene         = genes.incl[x]
    message("\n ------Working on ", selectedGene, "------")
    
    refRowIndex = which(row.names(seurat.sub.nonzero)==selectedGene)
    
    #in case a target gene is not included in the assay (or after the subsetting of the matrix)
    if(length(refRowIndex)==0){ 
      message("\n ERROR: the selected gene ", selectedGene, " is not included in the RNA-seq assay\n")
      #Not for shiny:message("\n ERROR: the selected gene ", selectedGene, " is not included in the RNA-seq assay\n Proceeding to the next gene")
      next
    }
    
    ######subsetting data##############
    
    #the row of the selected gene is in table div
    div = seurat.sub.nonzero[refRowIndex,]
    #div = div[["RNA"]]@data DON"T DO THAT!!!!! it skips rownames and colnames. Stick with the seurat object!!!!
    
    #Ea table holds only cells where the selected gene is NOT expressed
    
    #we must check if the selected gene is expressed in ALL cells, otherwise Ea gives an error
    
    #this variable is needed to go to the next iteration if becomes TRUE
    skip_to_next <- FALSE
    
    check.it <- tryCatch( div[,which(div[["RNA"]]@data[1,]==0)], error = function(e) { 
      
      e
      genes.xcl <<- c(genes.xcl, selectedGene)
      skip_to_next <<- TRUE 
      
      e$message <- message(e, selectedGene, " is expressed in ALL cells")
    })
    
    if(skip_to_next) { next }
    
    #if everything is fine continue with Ea
    Ea = div[,which(div[["RNA"]]@data[1,]==0)]
    
    #Eb holds the cells where the selected gene is expressed
    Eb = div[, which(div[["RNA"]]@data[1,]>0)]
    
    #table F: a subset of the original data  where the selected gene is not zero
    F = seurat.sub.nonzero[,colnames(Eb)]
    #exclude the row of the selected gene
    F = F[-refRowIndex,]
    
    #compute the percentage of cells (in the whole dataset) where the selected gene is expressed
    sel.gene.percent = ncol(Eb)/ncol(div)
    #Check for correlation between the selected gene and all other genes in all cells in F
    #use the transpose matrices
    Eb.t = as.matrix(t(Eb[["RNA"]]@data))
    
    F.t = as.matrix(t(F[["RNA"]]@data))
    
    #make a function to supress warnings
    cor.fun <- function(){
      #initialize a matrix to keep the results of the correlation test and the wilcoxon test
      corr.matrix = matrix(nrow = nrow(F), ncol = 5)
      rownames(corr.matrix) <- rownames(F)
      colnames(corr.matrix) <- c("r", "p_r", "p_w", "percentage in Ea", "percentage in Eb")
      
      for(j in 1:nrow(F)){
        
        #run the correlation test
        corr.output =  corr.test(as.numeric(Eb.t), as.numeric(F.t[,j]), method = paramCorrMethod, 
                                 adjust = "fdr", alpha = alpha.level, ci = FALSE )
        #append results to corr.matrix
        #corr.matrix[j,1] = colnames(F.t)[j]
        corr.matrix[j,1] = corr.output$r
        corr.matrix[j,2] = corr.output$p
        
      }
      corr.matrix
    }
    #alpha.level = 0.05
    #c.r = 0.3
    message("\n Running correlation...")
    corr.matrix = suppressWarnings(cor.fun())
    message("\n Correlation test run completed")
    #keep only those genes where p values are less than alpha.level = 0.05 c.r = 0.1
    corr.sign <- corr.matrix[which(corr.matrix[,2]<alpha.level & abs(as.numeric(corr.matrix[,1]))>c.r), ,drop = FALSE ]
    
    if(length(corr.sign)==0){
      message("None correlated gene returned at a level of significance below ", alpha.level, " and correlation above ", c.r)
      next
    }
    ####check whether the distibution of expression of each correlated gene in cells 
    ####where the selected gene is expressed or not is the same 
    
    ###run a Wilcoxon test###
    
    message("\n Running a Wilcoxon test...")
    for(i in 1:nrow(corr.sign)){
      
      #define the subsets of the data to include expression values for the correlated genes 
      
      #in cells where the selected gene is not expressed
      x.set = as.numeric(seurat.sub.nonzero[["RNA"]]@data[rownames(corr.sign)[i], colnames(Ea)])
      
      #the percentage of cells in Ea the selected gene is expressed
      percent.Ea = length(x.set[which(x.set>0)])/length(x.set)
      #append it to corr.sign
      corr.sign[i,4] = percent.Ea
      
      #where the selected gene is expressed
      y.set = as.numeric(seurat.sub.nonzero[["RNA"]]@data[rownames(corr.sign)[i], colnames(Eb)])
      
      #the percentage of cells in Eb the selected gene is expressed
      percent.Eb= length(y.set[which(y.set>0)])/length(y.set)
      #append it to corr.sign
      corr.sign[i,5] = percent.Eb
      
      #run the test
      p.value.wtest = wilcox.test(x.set,y.set, conf.int = FALSE)$p.value
      
      #append the p value to the corr.sign matrix
      corr.sign[i,3] = p.value.wtest
      
    }
    
    
    message("\n Wilcoxon test run completed")
    
    #again keep only those genes where p values are less than alpha.level
    corr.sign.genes <- corr.sign[which(corr.sign[,3]<alpha.level), , drop = FALSE]
    
    if(length(corr.sign.genes)==0){
      message(paste0("\n None correlated gene returned with p value < ", alpha.level))
      next
    }
    #order the list by the r coefficient
    corr.ordered <- corr.sign.genes[order(abs(as.numeric(corr.sign.genes[,1])),decreasing = TRUE), ,drop = FALSE ]
    
    #keep genes with r> c.r
    corr.ordered <- corr.ordered[which(abs(as.numeric(corr.ordered[,1]))>c.r), , drop = FALSE]
    
    #Not for shiny:message("\n Writing results to files")
    
    #Writing results to files
    #Not for shiny:outputPath          = paste0(the.path,'/', selectedGene, '/')
    
    #check if the directory exists and if not create it
    #Not for shiny:ifelse(!dir.exists(outputPath), dir.create(outputPath), FALSE)
    
    #Not for shiny:f.path = paste0(outputPath, paste0(selectedGene,"_", '.csv'))
    
    #append parameters of analysis to file
    #Not for shiny:cat(paste0("Correlation Method: ",paramCorrMethod,"\n"), 
    #Not for shiny:     paste0("correlation level: ", c.r, "\n"),
    #Not for shiny:    paste0("level of significance: ", alpha.level, "\n"),
    #Not for shiny:   paste0("Gene Sparcity Threshold: ", paramGeneSparcityThr, "\n"), 
    #Not for shiny:   paste0("percentage of cells where the selected gene is expressed: ", sel.gene.percent, "\n"),
    #Not for shiny:   file=f.path)
    
    #write the ordered matrix of correlated genes with r coefficients and p values to the designated directory
    #Not for shiny:suppressWarnings(write.table(corr.ordered, f.path, sep=",", append=TRUE, col.names=NA))
    #write.csv(corr.ordered, f.path)
    
    #append to cor.res
    cor.res[selectedGene] <- list(corr.ordered)
    cor.res[paste0(selectedGene, "_pos.cor")] <-list(corr.ordered[which(corr.ordered[,"r"]>0), , drop =FALSE])
    cor.res[paste0(selectedGene, "_neg.cor")] <-list(corr.ordered[which(corr.ordered[,"r"]<0), , drop =FALSE])
    #print(cor.res)
    
    #Not for shiny:f.path = paste0(outputPath, paste0(selectedGene, "_pos.cor"), "_", '.csv')
    #append parameters of analysis to file
    #Not for shiny:cat(paste0("Correlation Method: ",paramCorrMethod,"\n"), 
    #Not for shiny:    paste0("correlation level: ", c.r, "\n"),
    #Not for shiny:   paste0("level of significance: ", alpha.level, "\n"),
    #Not for shiny:   paste0("Gene Sparcity Threshold: ", paramGeneSparcityThr, "\n"), 
    #Not for shiny:   paste0("percentage of cells where the selected gene is expressed: ", sel.gene.percent, "\n"),
    #Not for shiny:   file=f.path)
    #write the ordered matrix of correlated genes with r coefficients and p values to the designated directory
    #Not for shiny:suppressWarnings(write.table( cor.res[paste0(selectedGene, "_pos.cor")], f.path, sep=",", append=TRUE, col.names=NA))
    
    #Not for shiny:f.path = paste0(outputPath, paste0(selectedGene, "_neg.cor"), "_", '.csv')
    #append parameters of analysis to file
    #Not for shiny:cat(paste0("Correlation Method: ",paramCorrMethod,"\n"), 
    #Not for shiny:    paste0("correlation level: ", c.r, "\n"),
    #Not for shiny:   paste0("level of significance: ", alpha.level, "\n"),
    #Not for shiny:   paste0("Gene Sparcity Threshold: ", paramGeneSparcityThr, "\n"), 
    #Not for shiny:   paste0("percentage of cells where the selected gene is expressed: ", sel.gene.percent, "\n"),
    #Not for shiny:   file=f.path)
    #write the ordered matrix of correlated genes with r coefficients and p values to the designated directory
    #Not for shiny:suppressWarnings(write.table( cor.res[paste0(selectedGene, "_neg.cor")], f.path, sep=",", append=TRUE, col.names=NA))
  }
  #Not for shiny:end.cor.time = Sys.time()
  message(paste0("\n ", "Selected gene is expressed among: ", sel.gene.percent, " of cells present in the set provided"))
  message("\n Correlation analysis completed")
  return(cor.res)
  
  
  
}