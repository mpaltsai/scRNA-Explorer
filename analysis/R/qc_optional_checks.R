# Optional QC evidence. These functions never filter or regress out cells/signals.
.qc_counts <- function(object) {
  if (inherits(object,"SingleCellExperiment")) return(SummarizedExperiment::assay(object,"counts"))
  if (inherits(object,"Seurat")) {
    layers <- SeuratObject::Layers(object[["RNA"]],search="^counts($|\\.)")
    if (length(layers)!=1L) stop("Choose an object with one RNA counts layer.")
    return(SeuratObject::LayerData(object,assay="RNA",layer=layers))
  }
  stop("Supply an existing SingleCellExperiment or Seurat object.")
}

seurat_cell_cycle_reference <- function(organism="Homo sapiens") {
  if (!identical(organism,"Homo sapiens")) stop("Supply an explicit species-matched reference; no automatic ortholog conversion.")
  if (!requireNamespace("Seurat",quietly=TRUE)) stop("Existing Seurat installation required.")
  e <- new.env(parent=baseenv())
  utils::data(list="cc.genes.updated.2019",package="Seurat",envir=e)
  list(name="Seurat cc.genes.updated.2019",species="Homo sapiens",
    s_genes=e$cc.genes.updated.2019$s.genes,g2m_genes=e$cc.genes.updated.2019$g2m.genes)
}

cell_cycle_qc_step <- function(object,parameters) {
  if (!requireNamespace("Seurat",quietly=TRUE)) stop("Existing Seurat installation required.")
  reference <- parameters$reference
  if (is.null(reference)) reference <- seurat_cell_cycle_reference(parameters$organism)
  if (is.null(parameters$organism) || !identical(reference$species,parameters$organism)) {
    stop("Cell-cycle reference species must match the explicitly selected organism.")
  }
  if (!is.character(reference$name) || length(reference$name)!=1L || !nzchar(reference$name)) stop("Name the cell-cycle reference.")
  counts <- .qc_counts(object)
  if (is.null(rownames(counts)) || anyNA(rownames(counts)) || any(!nzchar(rownames(counts))) || anyDuplicated(rownames(counts))) {
    stop("Cell-cycle scoring requires unique nonempty original gene identifiers.")
  }
  requested <- list(S=unique(reference$s_genes),G2M=unique(reference$g2m_genes))
  if (!all(vapply(requested,function(x) is.character(x) && length(x)>0L && !anyNA(x) && all(nzchar(x)),logical(1)))) {
    stop("Supply nonempty S and G2M gene sets.")
  }
  matched <- lapply(requested,intersect,y=rownames(counts))
  if (any(lengths(matched)==0L)) stop("A cell-cycle gene set has no matching features; scores are not assigned zero.")
  scale_factor <- if (is.null(parameters$scale_factor)) 10000 else parameters$scale_factor
  nbin <- if (is.null(parameters$nbin)) 24L else parameters$nbin
  ctrl <- if (is.null(parameters$ctrl)) min(lengths(matched)) else parameters$ctrl
  method_seed <- parameters$method_seed
  if (is.null(method_seed) || length(method_seed)!=1L || !is.numeric(method_seed) || is.na(method_seed) ||
      method_seed<0 || method_seed>.Machine$integer.max || method_seed!=floor(method_seed)) stop("Supply a nonnegative integer method_seed.")
  if (any(Matrix::colSums(counts)<=0)) stop("Review zero-count cells before scoring; none are removed automatically.")
  # Opaque, reversible names avoid Seurat's underscore replacement/collisions.
  # This changes only dimnames on the temporary scorer's sparse counts.
  scoring_names <- if (inherits(object,"Seurat")) rownames(counts) else paste0("qcFeature",seq_len(nrow(counts)))
  scoring_features <- lapply(matched,function(x) scoring_names[match(x,rownames(counts))])
  if (inherits(object,"Seurat")) temporary <- object else {
    scoring_counts <- counts
    rownames(scoring_counts) <- scoring_names
    temporary <- SeuratObject::CreateSeuratObject(counts=scoring_counts)
  }
  SeuratObject::DefaultAssay(temporary) <- "RNA"
  if (!identical(rownames(temporary),scoring_names) || !identical(colnames(temporary),colnames(object))) {
    stop("Seurat changed temporary feature/cell identifiers.")
  }
  message("Cell-cycle diagnostics: LogNormalize -> CellCycleScoring; no regression.")
  temporary <- Seurat::NormalizeData(temporary,assay="RNA",normalization.method="LogNormalize",
    scale.factor=scale_factor,verbose=FALSE)
  temporary <- Seurat::CellCycleScoring(temporary,s.features=scoring_features$S,g2m.features=scoring_features$G2M,
    ctrl=ctrl,nbin=nbin,set.ident=FALSE,seed=method_seed,assay="RNA")
  object$cell_cycle_S_score <- temporary$S.Score
  object$cell_cycle_G2M_score <- temporary$G2M.Score
  object$cell_cycle_phase <- as.character(temporary$Phase)
  list(object=object,diagnostics=list(status="scored",reference=reference,
    coverage=data.frame(program=names(requested),requested=as.integer(lengths(requested)),matched=as.integer(lengths(matched))),
    matched_genes=matched,missing_genes=Map(setdiff,requested,matched),
    scoring_gene_map=data.frame(original_gene=unique(unlist(matched)),
      temporary_gene=scoring_names[match(unique(unlist(matched)),rownames(counts))]),
    normalization="Temporary LogNormalize for diagnostic scoring only",scale_factor=scale_factor,
    nbin=nbin,ctrl=ctrl,method_seed=method_seed,regression_applied=FALSE,excluded_cells=0L,
    phase_counts=as.data.frame(table(phase=object$cell_cycle_phase)),
    interpretation="Program scores and predicted phases; not confirmed phases or proof of a downstream confounding effect."))
}

ambient_rna_readiness_step <- function(object,parameters) {
  reasons <- character()
  if (!identical(parameters$droplet_compatible,TRUE)) {
    reasons <- c(reasons,"SoupX droplet/protocol compatibility has not been established.")
  }
  if (is.null(parameters$raw_background_path)) {
    reasons <- c(reasons,"Matching raw background-droplet counts are not supplied; called-cell counts alone are insufficient.")
  } else if (length(parameters$raw_background_path)!=1L || !file.exists(parameters$raw_background_path)) {
    reasons <- c(reasons,"The declared raw background input is unavailable.")
  }
  columns <- if (inherits(object,"SingleCellExperiment")) names(SummarizedExperiment::colData(object)) else names(object[[]])
  if (is.null(parameters$capture_column) || !parameters$capture_column %in% columns) {
    reasons <- c(reasons,"Capture mapping is not supplied.")
  } else if (anyNA(object[[parameters$capture_column]])) reasons <- c(reasons,"Capture identities remain unresolved.")
  if (!requireNamespace("SoupX",quietly=TRUE)) reasons <- c(reasons,"SoupX is not installed; no installation is attempted automatically.")
  if (!length(reasons)) reasons <- "Declared inputs still require barcode/feature pairing and preliminary-cluster validation before SoupX estimation."
  list(object=object,diagnostics=list(method="SoupX",status="not_assessed",reasons=reasons,
    technology=parameters$technology,excluded_cells=0L,counts_corrected=FALSE,
    interpretation="Contamination is unknown, not zero. No correction has been performed."))
}

doublet_qc_step <- function(object,parameters) {
  columns <- if (inherits(object,"SingleCellExperiment")) names(SummarizedExperiment::colData(object)) else names(object[[]])
  reasons <- character()
  capture_column <- parameters$capture_column
  captures <- if (is.null(capture_column) || !capture_column %in% columns) NULL else as.character(object[[capture_column]])
  if (is.null(captures) || anyNA(captures) || any(!nzchar(captures))) reasons <- c(reasons,"Capture identities are unresolved; donors/files are not substituted for captures.")
  if (is.null(parameters$capture_evidence) || !nzchar(parameters$capture_evidence)) reasons <- c(reasons,"Capture provenance/confirmation is required.")
  rates <- parameters$rate_table
  if (is.null(rates)) reasons <- c(reasons,"Expected doublet rate, uncertainty and their source require review; no 10x rate is assumed.")
  if (!requireNamespace("scDblFinder",quietly=TRUE)) reasons <- c(reasons,"scDblFinder is not available in this R session's library paths.")
  if (length(reasons)) return(list(object=object,diagnostics=list(method="scDblFinder",status="not_assessed",
    reasons=reasons,excluded_cells=0L,technology=parameters$technology)))
  required <- c("capture","expected_rate","uncertainty","source")
  if (!is.data.frame(rates) || !all(required %in% names(rates)) || anyNA(rates[,required]) ||
      anyDuplicated(rates$capture) || !setequal(as.character(rates$capture),unique(captures)) ||
      !is.numeric(rates$expected_rate) || !is.numeric(rates$uncertainty) ||
      any(!is.finite(rates$expected_rate) | rates$expected_rate<=0 | rates$expected_rate>=1) ||
      any(!is.finite(rates$uncertainty) | rates$uncertainty<0 | rates$uncertainty>1) ||
      any(!nzchar(as.character(rates$source)))) stop("Supply a complete capture-specific rate_table: capture, expected_rate (0..1), uncertainty, source.")
  method_seed <- parameters$method_seed
  if (is.null(method_seed) || length(method_seed)!=1L || !is.numeric(method_seed) || is.na(method_seed) ||
      method_seed<0 || method_seed>.Machine$integer.max || method_seed!=floor(method_seed)) stop("Supply a nonnegative integer method_seed.")
  counts <- .qc_counts(object)
  if (any(Matrix::colSums(counts)<=0)) stop("Review zero-count cells first; no automatic cell exclusions.")
  scores <- rep(NA_real_,ncol(object)); classes <- rep(NA_character_,ncol(object)); summaries <- list()
  for (capture in unique(captures)) {
    idx <- which(captures==capture)
    rate <- rates[as.character(rates$capture)==capture,,drop=FALSE]
    # Temporary count-only SCE prevents stale reductions/normalized assays from
    # influencing a QC rerun. Counts stay sparse and no checkpoint is written.
    temporary <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=counts[,idx,drop=FALSE]))
    message("scDblFinder capture: ",capture," (",length(idx)," cells)")
    temporary <- scDblFinder::scDblFinder(temporary,clusters=FALSE,dbr=rate$expected_rate,
      dbr.sd=rate$uncertainty,BPPARAM=BiocParallel::SerialParam(RNGseed=method_seed),verbose=TRUE)
    if (!identical(colnames(temporary),colnames(object)[idx])) stop("Doublet method changed cell identities/order.")
    scores[idx] <- temporary$scDblFinder.score
    classes[idx] <- as.character(temporary$scDblFinder.class)
    summaries[[capture]] <- data.frame(capture=capture,cells=length(idx),predicted_doublets=sum(classes[idx]=="doublet"))
  }
  object$scDblFinder.score <- scores
  object$scDblFinder.class <- classes
  list(object=object,diagnostics=list(method="scDblFinder",status="scored",capture_column=capture_column,
    capture_evidence=parameters$capture_evidence,technology=parameters$technology,rate_table=rates,
    method_seed=method_seed,clusters=FALSE,excluded_cells=0L,per_capture=do.call(rbind,summaries),
    interpretation="Predicted doublet evidence only; high counts alone do not establish a doublet. No cells removed."))
}

decontx_qc_step <- function(object,parameters) {
  columns <- if (inherits(object,"SingleCellExperiment")) names(SummarizedExperiment::colData(object)) else names(object[[]])
  batch_column <- parameters$batch_column
  if (is.null(batch_column) && !identical(parameters$exploratory_pooled,TRUE)) {
    stop("Supply batch_column or explicitly enable exploratory_pooled; capture/batch identities are not inferred.")
  }
  batch <- if (is.null(batch_column)) NULL else {
    if (!batch_column %in% columns) stop("Batch metadata column is unavailable.")
    x <- as.character(object[[batch_column]])
    if (anyNA(x) || any(!nzchar(x))) stop("Batch metadata remain unresolved.")
    x
  }
  z <- if (is.null(parameters$clusters_column)) NULL else {
    if (!parameters$clusters_column %in% columns) stop("Preliminary cluster column is unavailable.")
    x <- as.character(object[[parameters$clusters_column]])
    if (anyNA(x) || any(!nzchar(x))) stop("Preliminary clusters must be nonempty.")
    x
  }
  if (!requireNamespace("decontX",quietly=TRUE)) stop("decontX must already be installed; no automatic installation.")
  seed <- parameters$method_seed
  if (is.null(seed) || length(seed)!=1L || !is.numeric(seed) || is.na(seed) ||
      seed<0 || seed>.Machine$integer.max || seed!=floor(seed)) stop("Supply a nonnegative integer method_seed.")
  defaults <- list(maxIter=500L,delta=c(10,10),estimateDelta=TRUE,convergence=0.001,varGenes=5000L,dbscanEps=1)
  settings <- parameters$model_args
  if (is.null(settings)) settings <- list()
  if (!is.list(settings) || (length(settings) && (is.null(names(settings)) || anyDuplicated(names(settings)) ||
      !all(names(settings) %in% names(defaults))))) stop("Unsupported DecontX model_args.")
  defaults[names(settings)] <- settings
  counts <- .qc_counts(object)
  sizes <- Matrix::colSums(counts)
  if (any(!is.finite(sizes) | sizes<=0)) stop("Review zero/invalid-count cells before fitting; none are removed automatically.")
  temporary <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=counts))
  warnings <- character()
  message("DecontX exploratory estimation on original counts; no correction will be applied.")
  temporary <- withCallingHandlers(do.call(decontX::decontX,c(list(x=temporary,z=z,batch=batch,
    assayName="counts",seed=seed,verbose=TRUE),defaults)),warning=function(w) {
      warnings <<- c(warnings,conditionMessage(w));invokeRestart("muffleWarning")
    })
  if (!identical(colnames(temporary),colnames(object)) ||
      !identical(SummarizedExperiment::assay(temporary,"counts"),counts)) stop("DecontX changed original count/cell identities.")
  fractions <- temporary$decontX_contamination
  clusters <- as.character(temporary$decontX_clusters)
  if (length(fractions)!=ncol(object) || any(!is.finite(fractions) | fractions<0 | fractions>1)) {
    stop("DecontX returned incomplete/invalid contamination estimates.")
  }
  object$ambient_decontX_fraction <- as.numeric(fractions)
  object$ambient_decontX_preliminary_cluster <- clusters
  fit_diagnostics <- lapply(S4Vectors::metadata(temporary)$decontX$estimates,function(x) {
    list(reported_iteration=x$iteration,iteration_limit=defaults$maxIter,
      reached_iteration_limit=x$iteration>=defaults$maxIter,log_likelihood=x$logLikelihood)
  })
  if (any(vapply(fit_diagnostics,function(x) isTRUE(x$reached_iteration_limit),logical(1)))) {
    warnings <- c(warnings,"DecontX reached an iteration limit; review fit convergence before interpreting estimates.")
  }
  list(object=object,diagnostics=list(method="DecontX",status="estimated",batch_column=batch_column,
    exploratory_pooled=is.null(batch),clusters_column=parameters$clusters_column,
    preliminary_clustering=if (is.null(z)) "DecontX heuristic clustering; not final analysis clusters" else "User-supplied preliminary labels",
    model_args=defaults,method_seed=seed,excluded_cells=0L,counts_corrected=FALSE,
    corrected_counts_retained=FALSE,warnings=unique(warnings),fit_diagnostics=fit_diagnostics,
    estimated_fraction_summary=c(min=min(fractions),median=stats::median(fractions),max=max(fractions)),
    interpretation=paste("Cluster-dependent model estimates, not independent contamination truth.",
      if (is.null(batch)) "Pooled exploration: capture mapping unresolved; not a capture-specific validated estimate." else "Fit by declared batches.",
      "Original counts remain selected; corrected counts were not applied or checkpointed.")))
}
