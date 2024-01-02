#load libraries
#x <- c("tidyverse", "SingleCellExperiment", "Matrix", "AnnotationHub", "ensembldb", "ggplot2", "scales", "Seurat", "scuttle", "scater", "mltools")
#invisible(lapply(x, suppressPackageStartupMessages(library), character.only = TRUE))

x <- c("dplyr", "SingleCellExperiment", "Matrix", "AnnotationHub", "ensembldb", "ggplot2", "scales", "Seurat", "scuttle", "scater", "mltools","DT", "gprofiler2")
invisible(lapply(x, suppressPackageStartupMessages(library), character.only = TRUE))
