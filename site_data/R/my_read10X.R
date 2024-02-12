my_read10X <- function (mtxPath, genesPath, barcodesPath = NULL, DGEList = FALSE) 
{
  
  mtx <- mtxPath
  genes <- genesPath
  if (!is.null(barcodesPath)) 
    barcodes <- barcodesPath
  m <- substring(readLines(mtx, n = 8), 1, 1)
  n.mtx.comment.lines <- which(m != "%")[1] - 1L
  N <- scan(mtx, skip = n.mtx.comment.lines, what = 0L, sep = " ", 
            nmax = 3, quiet = TRUE)
  ngenes <- N[1]
  ncells <- N[2]
  nmtx <- N[3]
  Genes <- read.table(genes, header = FALSE, comment.char = "", 
                      sep = "\t", row.names = 1, colClasses = "character")
  if (nrow(Genes) != ngenes) 
    stop("Number of feature IDs doesn't agree with header information in mtx file")
  names(Genes)[1] <- "Symbol"
  if (ncol(Genes) > 1L) 
    names(Genes)[2] <- "Type"
  m <- read.table(mtx, skip = n.mtx.comment.lines + 1L, header = FALSE, 
                  comment.char = "", sep = " ", colClasses = "integer", 
                  nrows = nmtx)
  y <- matrix(0L, ngenes, ncells)
  i <- m[, 1] + (m[, 2] - 1L) * ngenes
  y[i] <- m[, 3]
  dimnames(y) <- list(Gene = row.names(Genes), Cell = 1:ncells)
  if (is.null(barcodesPath)) {
    Samples <- NULL
  }
  else {
    Barcodes <- scan(barcodes, what = "", quiet = TRUE)
    if (length(Barcodes) != ncells) 
      stop("Number of barcodes doesn't agree with header information in mtx file")
    Samples <- data.frame(Barcode = Barcodes)
  }
  if (DGEList) {
    DGEList(count = y, genes = Genes, samples = Samples)
  }
  else {
    counts <- as(y, "dgCMatrix")
    list(counts = counts, samples = Samples, genes = Genes)
    
  }
}