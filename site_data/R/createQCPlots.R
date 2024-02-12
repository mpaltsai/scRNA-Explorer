#functions to draw QC plots

#number of UMIs per cell
createUMIsPlot <- function(my_metrics, minNUMIs, maxNUMIs){
  as.data.frame(my_metrics) %>% 
    ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
    geom_density() + 
    scale_x_log10() + 
    ylab("log10 cell density") +
    geom_vline(xintercept = minNUMIs)+
    geom_vline(xintercept = maxNUMIs)
}

#nGenes histogram
createGenesHistPlot <- function(my_metrics, minNGenes, maxNGenes){
  as.data.frame(my_metrics) %>% 
    ggplot(aes(color=sample, x=nGene, fill= sample)) + 
    geom_density() + 
    scale_x_log10() + 
    ylab("log10 cell density") +
    geom_vline(xintercept = minNGenes)+
    geom_vline(xintercept = maxNGenes)
}

#nGenes boxplot
createGenesBoxPlot <- function(my_metrics){
  as.data.frame(my_metrics) %>% 
    ggplot(aes(x=sample, y=nGene, fill=sample)) + 
    geom_boxplot() + 
    ggtitle("NCells vs NGenes")
}

#UMIs vs nGenes with mitoRatio
createUMIsGenesPlot<- function(my_metrics, minNGenes, maxNGenes){
  as.data.frame(my_metrics) %>% 
    ggplot(aes(x=nUMI, y=nGene)) + 
    geom_point(aes(color=mitoRatio)) +  
    stat_smooth(method=lm) +
    scale_x_log10() + 
    scale_y_log10() + 
    geom_hline(yintercept = minNGenes) +
    geom_hline(yintercept = maxNGenes) +
    facet_wrap(~sample)
}

#UMIs vs nGenes with hemoRatio
createUMIsGenesPlotRBC <- function(my_metrics, minNGenes, maxNGenes){
  as.data.frame(my_metrics) %>% 
    ggplot(aes(x=nUMI, y=nGene)) + 
    geom_point(aes(color=rbcRatio)) + 
    stat_smooth(method=lm) +
    scale_x_log10() + 
    scale_y_log10() + 
    geom_hline(yintercept = minNGenes) +
    geom_hline(yintercept = maxNGenes) +
    facet_wrap(~sample)
}

#plot mitoRatio
createmitoRatioPlot <- function(my_metrics, minMitoRatio, maxMitoRatio){
  as.data.frame(my_metrics) %>% 
    ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
    geom_density() + 
    #scale_x_log10() + 
    geom_vline(xintercept = minMitoRatio)+
    geom_vline(xintercept = maxMitoRatio)
}

#plot hemoRatio
createrbcRatioPlot <- function(my_metrics, minRbcRatio, maxRbcRatio){
  as.data.frame(my_metrics) %>% 
    ggplot(aes(color=sample, x=rbcRatio, fill=sample)) + 
    geom_density() + 
    #scale_x_log10() + 
    geom_vline(xintercept = minRbcRatio)+
    geom_vline(xintercept = maxRbcRatio)
}

#Genes vs UMIs
createnoveltyPlot <- function(my_metrics, minlog10GenesPerUMI, maxlog10GenesPerUMI){
  as.data.frame(my_metrics) %>%
    ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
    geom_density() + 
    geom_vline(xintercept = minlog10GenesPerUMI)+
    geom_vline(xintercept = maxlog10GenesPerUMI)
}

#plot scaterPlot
createlibraryPlot <- function(my_se){
  plotScater(my_se, nfeatures = 300, exprs_values = "counts")
}

#plot highest expressed features
createhighExprPlot <- function(my_se){
  plotHighestExprs(my_se, exprs_values = "counts")
}