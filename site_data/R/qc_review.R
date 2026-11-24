# The QC baseline is immutable. Review helpers return tables and selection
# masks; only Apply subsets the existing SCE. No expression checkpoint is saved.
qc_rule_spec <- function() data.frame(
  input=c("minUMIsPerCell","maxUMIsPerCell","minGenesPerCell","maxGenesPerCell",
    "noveltyPerCell","minMitoPerCell","maxMitoPerCell","minRbcPerCell","maxRbcPerCell"),
  metric=c("nUMI","nUMI","nGene","nGene","log10GenesPerUMI","mitoRatio","mitoRatio","rbcRatio","rbcRatio"),
  lower=c(TRUE,FALSE,TRUE,FALSE,TRUE,TRUE,FALSE,TRUE,FALSE),
  label=c("Low UMIs","High UMIs","Few genes","Many genes","Low complexity",
    "Low mitochondrial fraction","High mitochondrial fraction","Low hemoglobin fraction","High hemoglobin fraction"),
  stringsAsFactors=FALSE)

qc_review_metrics <- function(object) {
  metrics <- baseline_qc_metrics(object)
  if (!identical(metrics$cells,colnames(object)) || anyDuplicated(metrics$cells)) stop("QC metrics must match unique baseline cell IDs in order.")
  metrics$cell_id <- as.character(object$qc_cell_id)
  counts <- SummarizedExperiment::assay(object,"counts")
  annotations <- SummarizedExperiment::rowData(object)
  recipe <- S4Vectors::metadata(object)$scrna_explorer_qc$config
  # Numerators use positively annotated members; every denominator is the
  # original total UMI count. Unresolved-feature counts are reported alongside,
  # so a partial annotated share cannot be mistaken for complete classification.
  for (group in .qc_ribosomal_fields) {
    flag <- annotations[[paste0("qc_",group)]]
    available <- length(flag)==nrow(object) && any(!is.na(flag))
    if (group %in% .qc_ribosomal_fields[1:3] && !identical(recipe$ribosomal_reference$status,"available")) available <- FALSE
    known <- if (available) as.numeric(Matrix::colSums(counts[which(flag %in% TRUE),,drop=FALSE]))/metrics$nUMI else rep(NA_real_,ncol(object))
    unknown <- if (length(flag)==nrow(object)) as.numeric(Matrix::colSums(counts[which(is.na(flag)),,drop=FALSE]))/metrics$nUMI else rep(NA_real_,ncol(object))
    metrics[[paste0("ribo_",group)]] <- known
    metrics[[paste0("ribo_unknown_",group)]] <- unknown
  }
  status <- as.character(annotations$qc_annotation_status)
  for (group in c("matched","ambiguous","partial","unmapped")) metrics[[paste0("annotation_",group)]] <-
    if (length(status)==nrow(object)) as.numeric(Matrix::colSums(counts[which(status==group),,drop=FALSE]))/metrics$nUMI else rep(NA_real_,ncol(object))
  metrics
}

qc_default_thresholds <- function(metrics) {
  if (!nrow(metrics)) stop("This sample has no baseline cells.")
  list(minUMIsPerCell=100,maxUMIsPerCell=max(metrics$nUMI),minGenesPerCell=0,maxGenesPerCell=max(metrics$nGene),
    noveltyPerCell=0,minMitoPerCell=0,maxMitoPerCell=1,minRbcPerCell=0,maxRbcPerCell=1)
}

qc_validate_thresholds <- function(thresholds) {
  spec <- qc_rule_spec()
  if (!is.list(thresholds) || anyDuplicated(names(thresholds)) || !setequal(names(thresholds),spec$input)) stop("Supply all nine named QC thresholds.")
  values <- vapply(thresholds[spec$input],function(x) {
    if (!is.numeric(x) || length(x)!=1L || !is.finite(x)) stop("Every threshold must be one finite number.")
    as.numeric(x)
  },numeric(1))
  if (any(values<0) || any(values[5:9]>1)) stop("Counts must be nonnegative; complexity and fractions must lie between 0 and 1.")
  if (values[1]<100) stop("The QC baseline requires UMIs >= 100; its UMI minimum cannot be lowered.")
  if (any(values[1:4]!=floor(values[1:4]))) stop("UMI and gene thresholds must be integers.")
  if (any(values[c(1,3,6,8)]>values[c(2,4,7,9)])) stop("Each minimum must be less than or equal to its maximum.")
  as.list(values)
}

qc_baseline_signature <- function(object) {
  # Compute once per baseline revision, not on each slider movement. The hash
  # binds decisions to sparse counts, cell/feature IDs, metrics and annotations.
  digest::digest(list(counts=SummarizedExperiment::assay(object,"counts"),
    row_data=SummarizedExperiment::rowData(object),col_data=SummarizedExperiment::colData(object),
    qc=S4Vectors::metadata(object)$scrna_explorer_qc),algo="sha256")
}

.qc_rule_failures <- function(metrics,thresholds) {
  spec <- qc_rule_spec()
  result <- vapply(seq_len(nrow(spec)),function(i) {
    x <- metrics[[spec$metric[i]]]; cutoff <- thresholds[[spec$input[i]]]
    if (spec$lower[i]) x<cutoff else x>cutoff
  },logical(nrow(metrics)))
  dim(result) <- c(nrow(metrics),nrow(spec)); colnames(result) <- spec$input
  result
}

qc_filter_proposal <- function(metrics,thresholds,baseline_hash) {
  thresholds <- qc_validate_thresholds(thresholds)
  spec <- qc_rule_spec()
  unavailable <- unique(spec$metric[vapply(spec$metric,function(n) !is.numeric(metrics[[n]]) || any(!is.finite(metrics[[n]])),logical(1))])
  if (!nrow(metrics) || length(unavailable)) stop("QC evidence unavailable: ",paste(unavailable,collapse=", "),
    ". Inspect annotation/metric coverage before previewing the existing nine-rule filter.")
  failures <- .qc_rule_failures(metrics,thresholds)
  failed <- rowSums(failures)
  reasons <- apply(failures,1,function(v) paste(spec$label[v],collapse="; "))
  cells <- data.frame(cell_id=metrics$cell_id,cells=metrics$cells,sample=metrics$sample,
    decision=ifelse(failed==0L,"retain","proposed exclusion"),failed_rules=failed,reasons=reasons,stringsAsFactors=FALSE)
  rule_counts <- data.frame(rule=spec$label,failed=as.integer(colSums(failures)),
    failed_only_this_rule=as.integer(colSums(failures & (failed==1L))),stringsAsFactors=FALSE)
  combinations <- table(reasons[failed>0])
  overlap <- data.frame(combination=as.character(names(combinations)),cells=as.integer(combinations),stringsAsFactors=FALSE)
  list(baseline_hash=baseline_hash,thresholds=thresholds,
    proposal_id=digest::digest(list(baseline_hash=baseline_hash,thresholds=thresholds),algo="sha256"),
    keep=failed==0L,per_cell=cells,rule_counts=rule_counts,
    overlap=overlap[order(overlap$cells,decreasing=TRUE),,drop=FALSE],
    summary=data.frame(baseline_cells=nrow(metrics),retained=sum(failed==0L),excluded=sum(failed>0L)))
}

qc_apply_filter <- function(object,proposal,reason,history=list()) {
  if (length(reason)!=1L || is.na(reason) || !nzchar(trimws(reason))) stop("Provide a reason before applying filtering.")
  if (!identical(proposal$baseline_hash,qc_baseline_signature(object)) ||
      !identical(proposal$per_cell$cells,colnames(object))) stop("Preview is not for this baseline; preview again.")
  # Recheck the decision from the actual baseline. Never trust a supplied mask
  # or a stored preview from another method/threshold revision.
  verified <- qc_filter_proposal(qc_review_metrics(object),proposal$thresholds,proposal$baseline_hash)
  if (!identical(verified$proposal_id,proposal$proposal_id) || !identical(verified$keep,proposal$keep)) stop("Preview selection differs from the verified rules.")
  previous <- if (length(history)) history[[length(history)]]$revision else NA_character_
  record <- list(sample=unique(verified$per_cell$sample),baseline_hash=proposal$baseline_hash,
    proposal_id=proposal$proposal_id,previous_revision=previous,thresholds=proposal$thresholds,
    reason=trimws(reason),applied_at=format(Sys.time(),tz="UTC",usetz=TRUE),
    per_cell=verified$per_cell,rule_counts=verified$rule_counts,overlap=verified$overlap,summary=verified$summary,
    qc_config=S4Vectors::metadata(object)$scrna_explorer_qc$config,
    method="Inclusive bounds on the nine existing QC rules; no optional diagnostic exclusions",
    # Hash function text/formals, not closure environments containing app state.
    implementation_hash=digest::digest(lapply(list(qc_rule_spec,qc_review_metrics,qc_validate_thresholds,
      .qc_rule_failures,qc_filter_proposal,qc_apply_filter),function(f) list(formals=formals(f),body=deparse(body(f)))),algo="sha256"))
  record$revision <- digest::digest(list(record,sequence=length(history)+1L),algo="sha256")
  history <- c(history,list(record))
  filtered <- object[,verified$keep,drop=FALSE]
  md <- S4Vectors::metadata(filtered)
  md$scrna_explorer_qc$filtering <- list(revision=record$revision,history=history)
  S4Vectors::metadata(filtered) <- md
  list(object=filtered,decision=record,history=history)
}

# Hold all other cutoffs fixed while varying one cutoff. This is sensitivity to
# settings, not validation of cell quality and not a target retention rate.
qc_threshold_sensitivity <- function(metrics,thresholds,points=31L) {
  spec <- qc_rule_spec(); thresholds <- qc_validate_thresholds(thresholds)
  do.call(rbind,lapply(seq_len(nrow(spec)),function(i) {
    values <- metrics[[spec$metric[i]]]
    candidates <- unique(c(seq(min(values),max(values),length.out=points),thresholds[[spec$input[i]]]))
    if (i<=4) candidates <- unique(round(candidates))
    if (i==1) candidates <- candidates[candidates>=100]
    # Keep min/max relationships feasible while varying one rule.
    other <- switch(as.character(i),`1`=2L,`2`=1L,`3`=4L,`4`=3L,`6`=7L,`7`=6L,`8`=9L,`9`=8L,NA_integer_)
    if (!is.na(other)) candidates <- candidates[if(spec$lower[i]) candidates<=thresholds[[spec$input[other]]] else candidates>=thresholds[[spec$input[other]]]]
    data.frame(rule=spec$label[i],threshold=candidates,
      retained=vapply(candidates,function(v) {t<-thresholds;t[[spec$input[i]]]<-v;sum(rowSums(.qc_rule_failures(metrics,t))==0L)},integer(1)),
      current=thresholds[[spec$input[i]]],baseline_cells=nrow(metrics))
  }))
}

qc_require_applied <- function(baselines,filtered) {
  if (is.null(baselines) || !length(baselines)) stop("Create the QC baselines first.")
  missing <- names(baselines)[!vapply(names(baselines),function(n) !is.null(filtered[[n]]) && ncol(filtered[[n]])>0L,logical(1))]
  if (length(missing)) stop("Apply filtering with at least one retained cell in every sample. Pending/empty: ",paste(missing,collapse=", "))
  filtered[names(baselines)]
}
