# Optional models run on the fixed baseline and return compact evidence only.
# This bridge reuses the standalone algorithms; it never persists their temporary
# SCE/Seurat scorers or corrected-count assays.
qc_optional_library <- function(app_root=".") {
  configured<-getOption("scrnaExplorer.qcLibrary",Sys.getenv("SCRNA_EXPLORER_QC_LIBRARY",""))
  if(nzchar(configured)) return(configured)
  local<-file.path(app_root,"..","analysis","local","R-library")
  if(dir.exists(local)) return(normalizePath(local))
  # Reuse the already-approved library from the common repository when running
  # in an isolated Git worktree. No settings or package files are changed.
  if(nzchar(Sys.which("git"))) {
    common<-suppressWarnings(system2("git",c("-C",shQuote(normalizePath(app_root)),"rev-parse","--path-format=absolute","--git-common-dir"),stdout=TRUE,stderr=FALSE))
    if(length(common)==1L) {
      candidate<-file.path(dirname(common),"analysis","local","R-library")
      if(dir.exists(candidate)) return(normalizePath(candidate))
    }
  }
  ""
}

qc_optional_runtime <- function(library_path,seed,execute) {
  libraries<-.libPaths();kind<-RNGkind()
  had_seed<-exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  old_seed<-if(had_seed) get(".Random.seed",envir=.GlobalEnv) else NULL
  on.exit({
    .libPaths(libraries);do.call(RNGkind,as.list(kind))
    if(had_seed) assign(".Random.seed",old_seed,envir=.GlobalEnv)
    else if(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)) rm(".Random.seed",envir=.GlobalEnv)
  },add=TRUE)
  if(length(seed)!=1L || !is.numeric(seed) || !is.finite(seed) || seed<0 || seed> .Machine$integer.max || seed!=floor(seed)) stop("Choose a nonnegative integer seed.")
  if(!is.null(library_path) && nzchar(library_path)) {
    if(!dir.exists(library_path)) stop("The selected existing R library is unavailable; no installation attempted.")
    .libPaths(c(normalizePath(library_path),libraries))
  }
  set.seed(as.integer(seed));execute()
}

.qc_optional_fields <- list(doublets=c("scDblFinder.score","scDblFinder.class"),
  ambient=c("ambient_decontX_fraction","ambient_decontX_preliminary_cluster"),
  cell_cycle=c("cell_cycle_S_score","cell_cycle_G2M_score","cell_cycle_phase"))
.qc_optional_columns <- list(doublets=c("qc_diag_doublet_score","qc_diag_doublet_class"),
  ambient=c("qc_diag_ambient_fraction","qc_diag_ambient_population"),
  cell_cycle=c("qc_diag_S_score","qc_diag_G2M_score","qc_diag_phase"))

qc_optional_core <- function(app_root=".") {
  path<-normalizePath(file.path(app_root,"R","qc_optional_checks.R"),mustWork=TRUE)
  e<-new.env(parent=baseenv());sys.source(path,envir=e)
  list(functions=e,path=path,md5=unname(tools::md5sum(path)))
}

qc_capture_assignment <- function(object,config,cohort=list()) {
  mode<-config$mode
  values<-NULL
  if(identical(mode,"single")) {
    if(length(config$id)!=1L || is.na(config$id) || !nzchar(trimws(config$id))) stop("Name the explicitly declared capture.")
    values<-rep(trimws(config$id),ncol(object))
  } else if(identical(mode,"metadata")) {
    column<-config$column
    if(length(column)!=1L || !column %in% names(SummarizedExperiment::colData(object))) stop("Select an existing capture metadata column.")
    values<-as.character(object[[column]])
    if(length(values)!=ncol(object) || anyNA(values) || any(!nzchar(trimws(values)))) stop("Capture identity is unresolved; missing metadata are not capture labels.")
    # A shared capture must be modelled together. Do not score an input subset
    # as a complete capture simply because it is stored in a separate file.
    for(other in cohort) if(inherits(other,"SingleCellExperiment") && column %in% names(SummarizedExperiment::colData(other)) &&
      any(as.character(other[[column]]) %in% values)) stop("A declared capture spans multiple input samples. Assemble/review the complete capture before doublet scoring.")
  }
  if(is.null(values) || length(values)!=ncol(object) || anyNA(values) || any(!nzchar(trimws(values)))) stop("Capture identity is unresolved; files, donors and GEO accessions are not substituted.")
  for(other in cohort) if(is.list(other) && !inherits(other,"SingleCellExperiment") &&
    any(values %in% other$declared_captures)) stop("A capture is declared in multiple input samples. Review/assemble the complete capture before scoring.")
  if(!isTRUE(config$complete) || length(config$evidence)!=1L || is.na(config$evidence) || !nzchar(trimws(config$evidence)))
    stop("Confirm complete capture input and record the evidence for this assignment.")
  values
}

qc_pinned_cycle_annotation <- function(object) {
  cfg<-S4Vectors::metadata(object)$scrna_explorer_qc$config
  p<-cfg$reference$path
  if(length(p)!=1L || !file.exists(p) || !identical(unname(tools::md5sum(p)),cfg$reference$md5))
    stop("The baseline's pinned gene annotation file is missing or changed. Restore it or prepare annotations again.")
  if(grepl("[.]rds$",p,ignore.case=TRUE)) table<-readRDS(p) else {
    e<-new.env(parent=baseenv());load(p,envir=e);table<-e$annotations
  }
  if(!identical(digest::digest(table,algo="sha256"),cfg$reference_table_hash)) stop("Pinned annotation table hash differs from the baseline.")
  list(table=table,species=cfg$organism,feature_id_type=cfg$feature_id_type,provenance=cfg$reference,
    feature_map=data.frame(feature=rownames(object),lookup_id=as.character(SummarizedExperiment::rowData(object)$qc_source_feature_id)))
}

qc_optional_prepare <- function(object,method,config,core,cohort=list()) {
  if(!is.null(config$invalid_reason)) stop(config$invalid_reason)
  params<-config$parameters;params$method_seed<-config$seed
  if(method=="doublets") {
    captures<-qc_capture_assignment(object,config$capture,cohort)
    if(!identical(config$rates_mode,"csv")) {
      rates<-data.frame(capture=unique(captures),expected_rate=config$expected_rate,
        uncertainty=config$uncertainty,source=config$rate_source,stringsAsFactors=FALSE)
    } else rates<-config$rate_table
    needed<-c("capture","expected_rate","uncertainty","source")
    if(!is.data.frame(rates) || !all(needed %in% names(rates)) || anyNA(rates[,needed]) ||
      anyDuplicated(rates$capture) || !setequal(as.character(rates$capture),unique(captures)) ||
      !is.numeric(rates$expected_rate) || !is.numeric(rates$uncertainty) ||
      any(!is.finite(rates$expected_rate)|rates$expected_rate<=0|rates$expected_rate>=1) ||
      any(!is.finite(rates$uncertainty)|rates$uncertainty<0|rates$uncertainty>1) || any(!nzchar(trimws(rates$source))))
      stop("Review a complete per-capture rate table: capture, expected_rate (0–1), uncertainty and source. No default 10X rate is assumed.")
    object$qc_diagnostic_capture<-captures
    params<-c(params,list(capture_column="qc_diagnostic_capture",capture_evidence=config$capture$evidence,
      rate_table=rates,technology=config$technology))
  } else if(method=="ambient") {
    captures<-tryCatch(qc_capture_assignment(object,config$capture),error=function(e) NULL)
    if(identical(config$ambient_method,"SoupX")) {
      if(!is.null(captures)) object$qc_diagnostic_capture<-captures
      return(list(object=object,parameters=list(technology=config$technology,
        droplet_compatible=isTRUE(config$droplet_compatible),raw_background_path=NULL,
        capture_column=if(is.null(captures)) NULL else "qc_diagnostic_capture")))
    }
    if(is.null(captures) && !isTRUE(config$exploratory_pooled)) stop("Capture/batch identity is unresolved. Explicitly enable exploratory pooled DecontX or provide a reviewed capture assignment.")
    if(!is.null(captures)) object$qc_diagnostic_capture<-captures
    params$batch_column<-if(is.null(captures)) NULL else "qc_diagnostic_capture"
    params$exploratory_pooled<-is.null(captures)
    params$clusters_column<-if(nzchar(config$clusters_column)) config$clusters_column else NULL
    a<-params$model_args
    for(n in c("maxIter","varGenes")) if(length(a[[n]])!=1L || !is.finite(a[[n]]) || a[[n]]<1 || a[[n]]!=floor(a[[n]])) stop("DecontX iteration/variable-gene settings must be positive integers.")
    if(length(a$delta)!=2L || any(!is.finite(a$delta)|a$delta<=0) || length(a$convergence)!=1L ||
      !is.finite(a$convergence) || a$convergence<=0 || length(a$dbscanEps)!=1L || !is.finite(a$dbscanEps) || a$dbscanEps<=0)
      stop("Review positive DecontX prior, convergence and clustering settings.")
  } else {
    params$organism<-S4Vectors::metadata(object)$scrna_explorer_qc$config$organism
    if(is.null(params$reference)) params$reference<-core$functions$seurat_cell_cycle_reference(params$organism)
    if(!identical(params$reference$species,params$organism)) stop("Cell-cycle markers must explicitly match the baseline species; no automatic ortholog conversion.")
    params$annotation<-qc_pinned_cycle_annotation(object)
    params$strip_ensembl_version<-isTRUE(S4Vectors::metadata(object)$scrna_explorer_qc$config$strip_ensembl_version)
    mapping<-core$functions$.qc_cell_cycle_mapping(rownames(object),list(S=params$reference$s_genes,G2M=params$reference$g2m_genes),
      params$annotation,params$organism,params$strip_ensembl_version)
    if(length(config$minimum_markers)!=1L || !is.finite(config$minimum_markers) || config$minimum_markers<1 || config$minimum_markers!=floor(config$minimum_markers)) stop("Choose a positive integer minimum marker count.")
    if(any(lengths(mapping$matched)<config$minimum_markers)) stop("Too few unambiguously mapped cell-cycle markers. Review marker coverage; unavailable scores are not zero.")
    if(length(params$nbin)!=1L || !is.finite(params$nbin) || params$nbin<2 || params$nbin!=floor(params$nbin) ||
      length(params$scale_factor)!=1L || !is.finite(params$scale_factor) || params$scale_factor<=0 ||
      (!is.null(params$ctrl) && (length(params$ctrl)!=1L || !is.finite(params$ctrl) || params$ctrl<1 || params$ctrl!=floor(params$ctrl)))) stop("Review the cell-cycle bin, control-gene and normalization settings.")
  }
  list(object=object,parameters=params)
}

qc_run_optional <- function(object,method,config,baseline_hash,app_root=".",cohort=list(),execute=NULL) {
  if(!method %in% names(.qc_optional_fields)) stop("Choose a supported optional QC check.")
  if(!identical(qc_baseline_signature(object),baseline_hash)) stop("The optional-check input differs from its declared baseline.")
  fields<-.qc_optional_fields[[method]];columns<-.qc_optional_columns[[method]]
  values<-data.frame(cell_id=as.character(object$qc_cell_id),cells=colnames(object),stringsAsFactors=FALSE)
  for(i in seq_along(columns)) values[[columns[i]]]<-if(grepl("class|population|phase",columns[i])) NA_character_ else NA_real_
  result<-list(method=method,baseline_hash=baseline_hash,config=config,parameter_hash=digest::digest(config,algo="sha256"),
    created_at=format(Sys.time(),tz="UTC",usetz=TRUE),status="not_assessed",reasons=character(),warnings=character(),values=values)
  result<-tryCatch(qc_optional_runtime(config$library_path,config$seed,function() {
    core<-qc_optional_core(app_root);out<-result
    out$core_md5<-core$md5
    readiness<-method=="ambient" && identical(config$ambient_method,"SoupX")
    package<-if(readiness) "SoupX" else c(doublets="scDblFinder",ambient="decontX",cell_cycle="Seurat")[[method]]
    available<-requireNamespace(package,quietly=TRUE)
    if(!available && !readiness) {out$reasons<-paste(package,"is unavailable in the selected existing libraries. No package installation attempted.");return(out)}
    prepared<-tryCatch(qc_optional_prepare(object,method,config,core,cohort),error=function(e) e)
    if(inherits(prepared,"error")) {out$reasons<-conditionMessage(prepared);return(out)}
    run<-if(is.null(execute)) core$functions[[if(readiness) "ambient_rna_readiness_step" else c(doublets="doublet_qc_step",ambient="decontx_qc_step",cell_cycle="cell_cycle_qc_step")[[method]]]] else execute
    captured<-character()
    fitted<-tryCatch(withCallingHandlers(run(prepared$object,prepared$parameters),warning=function(w) {
      captured<<-c(captured,conditionMessage(w));invokeRestart("muffleWarning")
    }),error=function(e) e)
    out$software<-list(R=R.version.string,package=package,version=if(available) as.character(getNamespaceVersion(package)) else NA_character_,
      path=if(available) getNamespaceInfo(asNamespace(package),"path") else NA_character_)
    out$warnings<-captured
    if(inherits(fitted,"error")) {out$status<-"failed";out$reasons<-conditionMessage(fitted);return(out)}
    if(!identical(SummarizedExperiment::assay(fitted$object,"counts"),SummarizedExperiment::assay(object,"counts")) ||
      !identical(colnames(fitted$object),colnames(object)) || !identical(rownames(fitted$object),rownames(object)) ||
      isTRUE(fitted$diagnostics$counts_corrected) || isTRUE(fitted$diagnostics$excluded_cells>0)) stop("Optional method changed original counts/cell identities or applied exclusions; result rejected.")
    out$diagnostics<-fitted$diagnostics
    out$status<-fitted$diagnostics$status
    out$warnings<-unique(c(out$warnings,fitted$diagnostics$warnings))
    if(out$status %in% c("scored","estimated")) {
      for(i in seq_along(fields)) {
        v<-fitted$object[[fields[i]]]
        if(length(v)!=ncol(object)) stop("Optional method returned incomplete per-cell evidence.")
        out$values[[columns[i]]]<-if(is.character(values[[columns[i]]])) as.character(v) else as.numeric(v)
      }
      numeric_fields<-columns[!grepl("class|population|phase",columns)]
      if(any(!is.finite(as.matrix(out$values[numeric_fields])))) stop("Optional method returned unavailable/nonfinite scores; result rejected.")
      if(method=="doublets" && (any(!out$values$qc_diag_doublet_class %in% c("singlet","doublet")) || any(out$values$qc_diag_doublet_score<0|out$values$qc_diag_doublet_score>1))) stop("Invalid doublet predictions.")
      if(method=="ambient" && any(out$values$qc_diag_ambient_fraction<0|out$values$qc_diag_ambient_fraction>1)) stop("Invalid ambient estimates.")
      if(method=="ambient") {
        group<-if(is.null(prepared$parameters$batch_column)) rep("pooled",ncol(object)) else as.character(prepared$object[[prepared$parameters$batch_column]])
        populations<-vapply(split(out$values$qc_diag_ambient_population,group),function(v) length(unique(v[!is.na(v)&nzchar(v)])),integer(1))
        if(any(populations<2L)) {
          out$status<-"not_assessed"
          out$reasons<-paste("Insufficient preliminary population structure without background droplets in fit(s):",paste(names(populations)[populations<2L],collapse=", "))
          for(field in columns) out$values[[field]][]<-NA
          out$diagnostics$status<-"not_assessed";out$diagnostics$estimated_fraction_summary<-NULL
        }
      }
      if(method=="cell_cycle" && (anyNA(out$values$qc_diag_phase) || any(!out$values$qc_diag_phase %in% c("G1","S","G2M","Undecided")))) stop("Invalid or unavailable predicted cell-cycle phases.")
    } else out$reasons<-fitted$diagnostics$reasons
    if(method=="doublets") out$values$qc_diag_capture<-as.character(prepared$object$qc_diagnostic_capture)
    if(method=="ambient" && out$status=="estimated") out$values$qc_diag_ambient_fit_group<-
      if(is.null(prepared$parameters$batch_column)) rep("Pooled exploration (capture unresolved)",ncol(object)) else as.character(prepared$object[[prepared$parameters$batch_column]])
    out
  }),error=function(e) {out<-result;out$status<-"failed";out$reasons<-conditionMessage(e);out})
  result$revision<-digest::digest(result,algo="sha256")
  result
}

qc_cycle_marker_review <- function(object,config,app_root=".") qc_optional_runtime(config$library_path,config$seed,function() {
  if(!is.null(config$invalid_reason)) stop(config$invalid_reason)
  core<-qc_optional_core(app_root)
  species<-S4Vectors::metadata(object)$scrna_explorer_qc$config$organism
  reference<-config$parameters$reference
  if(is.null(reference)) reference<-core$functions$seurat_cell_cycle_reference(species)
  if(!identical(reference$species,species)) stop("Supply explicitly species-matched cell-cycle markers.")
  annotation<-qc_pinned_cycle_annotation(object)
  mapping<-core$functions$.qc_cell_cycle_mapping(rownames(object),list(S=reference$s_genes,G2M=reference$g2m_genes),
    annotation,species,isTRUE(S4Vectors::metadata(object)$scrna_explorer_qc$config$strip_ensembl_version))
  list(reference=reference,annotation=mapping[c("provenance","table_hash","feature_map_hash")],mapping=mapping$mapping,
    coverage=data.frame(program=c("S","G2M"),requested=c(length(reference$s_genes),length(reference$g2m_genes)),matched=as.integer(lengths(mapping$matched))))
})

qc_attach_optional <- function(object,results,baseline_hash) {
  if(!length(results)) return(object)
  md<-S4Vectors::metadata(object);profile<-md$scrna_explorer_qc
  own<-profile$optional_diagnostics
  prior<-profile$previous_qc
  while(is.null(own) && is.list(prior) && identical(prior$producer,"scRNA-Explorer")) {own<-prior$optional_diagnostics;prior<-prior$previous_qc}
  if(!is.null(own) && (!identical(own$producer,"scRNA-Explorer") || !identical(own$schema_version,1L)))
    stop("Optional-QC metadata ownership/version is unrecognised; no fields were overwritten.")
  existing<-grep("^qc_diag_",names(SummarizedExperiment::colData(object)),value=TRUE)
  if(length(setdiff(existing,own$fields))) stop("Imported metadata uses reserved qc_diag_* fields; none were overwritten.")
  for(n in existing) object[[n]]<-NULL
  for(r in results) {
    if(!identical(r$baseline_hash,baseline_hash)) stop("Optional evidence belongs to another baseline.")
    i<-match(as.character(object$qc_cell_id),r$values$cell_id)
    if(anyNA(i) || anyDuplicated(r$values$cell_id)) stop("Optional QC cell identifiers do not match the object.")
    for(n in setdiff(names(r$values),c("cell_id","cells"))) object[[n]]<-r$values[[n]][i]
  }
  md$scrna_explorer_qc$optional_diagnostics<-list(producer="scRNA-Explorer",schema_version=1L,
    baseline_hash=baseline_hash,fields=grep("^qc_diag_",names(SummarizedExperiment::colData(object)),value=TRUE),results=results)
  S4Vectors::metadata(object)<-md;object
}
