filter_cells <- function(raw_data, metrics, minUMIsPerCell, maxUMIsPerCell, minGenesPerCell, maxGenesPerCell, noveltyPerCell, minMitoPerCell, maxMitoPerCell, minRbcPerCell, maxRbcPerCell){
  filt_cells <-    metrics %>%
    dplyr::filter(nUMI > minUMIsPerCell , 
                  nUMI < maxUMIsPerCell,
                  nGene > minGenesPerCell,
                  nGene < maxGenesPerCell,
                  log10GenesPerUMI > noveltyPerCell,
                  mitoRatio > minMitoPerCell,
                  mitoRatio < maxMitoPerCell,
                  rbcRatio > minRbcPerCell,
                  rbcRatio < maxRbcPerCell
    ) %>% 
    pull(cells)
  #my_se_c <- calc_se()$se
  my_se_c <- raw_data
  #return filtered data
  if(dim(my_se_c[ ,filt_cells])[[2]] == 0){
    return(c(se = my_se_c[ ,filt_cells],
             mes= "There are no cells left after filtering. You have to change your filtering thresholds"))
  }else{
    return(c(se = my_se_c[ ,filt_cells],
             mes= "Data filtered"))
  }
  
}