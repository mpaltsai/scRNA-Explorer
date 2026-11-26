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

qc_optional_filter_defaults <- function() list(doublets="none",doublet_score=NA_real_,ambient=FALSE,ambient_fraction=NA_real_)

qc_validate_optional_filters <- function(rules=qc_optional_filter_defaults()) {
  if(!is.list(rules) || anyDuplicated(names(rules)) || !setequal(names(rules),names(qc_optional_filter_defaults())))
    stop("Supply the named optional exclusion settings.")
  if(length(rules$doublets)!=1L || !rules$doublets %in% c("none","prediction","score") ||
    !is.logical(rules$ambient) || length(rules$ambient)!=1L || is.na(rules$ambient))
    stop("Choose valid doublet and ambient exclusion options.")
  fraction<-function(x,label) {
    if(!is.numeric(x) || length(x)!=1L || !is.finite(x) || x<0 || x>1)
      stop("Enter a maximum retained ",label," between 0 and 1; no cutoff is chosen automatically.")
    as.numeric(x)
  }
  list(doublets=rules$doublets,
    doublet_score=if(rules$doublets=="score") fraction(rules$doublet_score,"doublet score") else NA_real_,
    ambient=rules$ambient,
    ambient_fraction=if(rules$ambient) fraction(rules$ambient_fraction,"ambient fraction") else NA_real_)
}

qc_optional_filters_enabled <- function(rules) rules$doublets!="none" || isTRUE(rules$ambient)

qc_optional_filter_state <- function(metrics,rules,results,baseline_hash) {
  rules<-qc_validate_optional_filters(rules)
  failures<-matrix(FALSE,nrow=nrow(metrics),ncol=0L);labels<-character();evidence<-list();values<-list()
  add<-function(method,field,valid,fail,label) {
    r<-results[[method]];status<-if(method=="doublets") "scored" else "estimated"
    if(is.null(r) || !identical(r$status,status) || !identical(r$baseline_hash,baseline_hash) ||
      length(r$revision)!=1L || is.na(r$revision) || !nzchar(r$revision))
      stop("Optional exclusion requires current ",method," results for this baseline. Run/review the check and preview again; missing or stale evidence cannot exclude cells.")
    v<-r$values
    if(!is.data.frame(v) || !all(c("cell_id",field) %in% names(v)) || anyNA(v$cell_id) ||
      anyDuplicated(v$cell_id) || nrow(v)!=nrow(metrics) || !setequal(v$cell_id,metrics$cell_id))
      stop("Optional exclusion cell identifiers do not match the complete baseline.")
    x<-v[[field]][match(metrics$cell_id,v$cell_id)]
    if(!valid(x)) stop("Optional exclusion has missing/invalid ",method," evidence; no unavailable value is treated as zero or singlet.")
    failures<<-cbind(failures,fail(x));colnames(failures)[ncol(failures)]<<-paste0("optional_",method)
    labels<<-c(labels,label)
    values[[paste0(field,"_used")]]<<-x
    fields<-c("method","status","revision","baseline_hash","parameter_hash","analysis_context","software","config","core_md5","created_at")
    ref<-stats::setNames(lapply(fields,function(field) r[[field]]),fields)
    ref$values_hash<-digest::digest(data.frame(cell_id=metrics$cell_id,value=x),algo="sha256")
    evidence[[method]]<<-ref
  }
  valid_fraction<-function(x) is.numeric(x) && length(x)==nrow(metrics) && all(is.finite(x)&x>=0&x<=1)
  if(rules$doublets=="prediction") add("doublets","qc_diag_doublet_class",
    function(x) !anyNA(x) && all(x %in% c("singlet","doublet")),function(x) x=="doublet","Predicted doublet")
  if(rules$doublets=="score") add("doublets","qc_diag_doublet_score",valid_fraction,
    function(x) x>rules$doublet_score,"High doublet score")
  if(rules$ambient) add("ambient","qc_diag_ambient_fraction",valid_fraction,
    function(x) x>rules$ambient_fraction,"High estimated ambient fraction")
  list(rules=rules,failures=failures,labels=labels,evidence=evidence,values=values)
}

qc_filter_proposal <- function(metrics,thresholds,baseline_hash,
    optional_rules=qc_optional_filter_defaults(),optional_results=list()) {
  thresholds <- qc_validate_thresholds(thresholds)
  spec <- qc_rule_spec()
  unavailable <- unique(spec$metric[vapply(spec$metric,function(n) !is.numeric(metrics[[n]]) || any(!is.finite(metrics[[n]])),logical(1))])
  if (!nrow(metrics) || length(unavailable)) stop("QC evidence unavailable: ",paste(unavailable,collapse=", "),
    ". Inspect annotation/metric coverage before previewing the existing nine-rule filter.")
  failures <- .qc_rule_failures(metrics,thresholds)
  optional<-qc_optional_filter_state(metrics,optional_rules,optional_results,baseline_hash)
  failures<-cbind(failures,optional$failures);labels<-c(spec$label,optional$labels)
  failed <- rowSums(failures)
  reasons <- apply(failures,1,function(v) paste(labels[v],collapse="; "))
  cells <- data.frame(cell_id=metrics$cell_id,cells=metrics$cells,sample=metrics$sample,
    decision=ifelse(failed==0L,"retain","proposed exclusion"),failed_rules=failed,reasons=reasons,stringsAsFactors=FALSE)
  # Record only the per-cell fields actually used, not another expression
  # object. These decision-time values survive later model reruns/exports.
  for(field in names(optional$values)) cells[[field]]<-optional$values[[field]]
  rule_counts <- data.frame(rule=labels,failed=as.integer(colSums(failures)),
    failed_only_this_rule=as.integer(colSums(failures & (failed==1L))),stringsAsFactors=FALSE)
  combinations <- table(reasons[failed>0])
  overlap <- data.frame(combination=as.character(names(combinations)),cells=as.integer(combinations),stringsAsFactors=FALSE)
  identity<-list(baseline_hash=baseline_hash,thresholds=thresholds)
  if(qc_optional_filters_enabled(optional$rules)) identity$optional<-optional[c("rules","evidence")]
  list(baseline_hash=baseline_hash,thresholds=thresholds,
    optional_rules=optional$rules,optional_evidence=optional$evidence,optional_keep=rowSums(optional$failures)==0L,
    proposal_id=digest::digest(identity,algo="sha256"),
    keep=failed==0L,per_cell=cells,rule_counts=rule_counts,
    overlap=overlap[order(overlap$cells,decreasing=TRUE),,drop=FALSE],
    summary=data.frame(baseline_cells=nrow(metrics),retained=sum(failed==0L),excluded=sum(failed>0L)))
}

qc_apply_filter <- function(object,proposal,reason,history=list(),optional_evidence=NULL,optional_results=list()) {
  if (length(reason)!=1L || is.na(reason) || !nzchar(trimws(reason))) stop("Provide a reason before applying filtering.")
  if (!identical(proposal$baseline_hash,qc_baseline_signature(object)) ||
      !identical(proposal$per_cell$cells,colnames(object))) stop("Preview is not for this baseline; preview again.")
  # Recheck the decision from the actual baseline. Never trust a supplied mask
  # or a stored preview from another method/threshold revision.
  rules<-if(is.null(proposal$optional_rules)) qc_optional_filter_defaults() else proposal$optional_rules
  verified <- qc_filter_proposal(qc_review_metrics(object),proposal$thresholds,proposal$baseline_hash,rules,optional_results)
  if (!identical(verified$proposal_id,proposal$proposal_id) || !identical(verified$keep,proposal$keep)) stop("Preview selection differs from the verified rules.")
  previous <- if (length(history)) history[[length(history)]]$revision else NA_character_
  record <- list(sample=unique(verified$per_cell$sample),baseline_hash=proposal$baseline_hash,
    proposal_id=proposal$proposal_id,previous_revision=previous,thresholds=proposal$thresholds,
    reason=trimws(reason),applied_at=format(Sys.time(),tz="UTC",usetz=TRUE),
    per_cell=verified$per_cell,rule_counts=verified$rule_counts,overlap=verified$overlap,summary=verified$summary,
    qc_config=S4Vectors::metadata(object)$scrna_explorer_qc$config,
    optional_rules=verified$optional_rules,optional_evidence_used=verified$optional_evidence,optional_keep=verified$optional_keep,
    method=if(qc_optional_filters_enabled(verified$optional_rules)) "Inclusive QC bounds plus explicitly selected doublet/ambient exclusions" else "Inclusive bounds on the nine existing QC rules; no optional diagnostic exclusions",
    # Hash function text/formals, not closure environments containing app state.
    implementation_hash=digest::digest(lapply(list(qc_rule_spec,qc_review_metrics,qc_validate_thresholds,
      .qc_rule_failures,qc_optional_filter_defaults,qc_validate_optional_filters,qc_optional_filters_enabled,qc_optional_filter_state,
      qc_filter_proposal,qc_apply_filter),function(f) list(formals=formals(f),body=deparse(body(f)))),algo="sha256"))
  if(length(optional_evidence)) record$optional_evidence<-optional_evidence
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
qc_threshold_sensitivity <- function(metrics,thresholds,points=31L,optional_keep=rep(TRUE,nrow(metrics))) {
  spec <- qc_rule_spec(); thresholds <- qc_validate_thresholds(thresholds)
  if(!is.logical(optional_keep) || length(optional_keep)!=nrow(metrics) || anyNA(optional_keep)) stop("Optional sensitivity mask must match the baseline.")
  do.call(rbind,lapply(seq_len(nrow(spec)),function(i) {
    values <- metrics[[spec$metric[i]]]
    candidates <- unique(c(seq(min(values),max(values),length.out=points),thresholds[[spec$input[i]]]))
    if (i<=4) candidates <- unique(round(candidates))
    if (i==1) candidates <- candidates[candidates>=100]
    # Keep min/max relationships feasible while varying one rule.
    other <- switch(as.character(i),`1`=2L,`2`=1L,`3`=4L,`4`=3L,`6`=7L,`7`=6L,`8`=9L,`9`=8L,NA_integer_)
    if (!is.na(other)) candidates <- candidates[if(spec$lower[i]) candidates<=thresholds[[spec$input[other]]] else candidates>=thresholds[[spec$input[other]]]]
    data.frame(rule=spec$label[i],threshold=candidates,
      retained=vapply(candidates,function(v) {t<-thresholds;t[[spec$input[i]]]<-v;sum(rowSums(.qc_rule_failures(metrics,t))==0L & optional_keep)},integer(1)),
      current=thresholds[[spec$input[i]]],baseline_cells=nrow(metrics))
  }))
}

qc_require_applied <- function(baselines,filtered) {
  if (is.null(baselines) || !length(baselines)) stop("Create the QC baselines first.")
  missing <- names(baselines)[!vapply(names(baselines),function(n) !is.null(filtered[[n]]) && ncol(filtered[[n]])>0L,logical(1))]
  if (length(missing)) stop("Apply filtering with at least one retained cell in every sample. Pending/empty: ",paste(missing,collapse=", "))
  filtered[names(baselines)]
}
