# Annotation/metric preparation is separate from baseline construction. Prepared
# tables contain no expression matrix; the baseline SCE is made only on request.
annotation_qc_server <- function(input,output,session,records,app_root=".",reference_services=list()) {
  # Test seams replace catalogue/resource access, never the QC implementation.
  references <- do.call(annotation_reference_server,c(list(input=input,output=output,session=session,app_root=app_root),reference_services))
  prepared <- shiny::reactiveVal(NULL)
  baseline <- shiny::reactiveVal(NULL)
  status <- shiny::reactiveVal("Read data, choose the species and identifier types, then add annotations.")
  baseline_status <- shiny::reactiveVal("")
  busy <- shiny::reactiveVal(FALSE)
  revision <- 0L
  previous_records <- previous_options <- NULL
  configured <- FALSE
  suggestions <- shiny::reactive({
    x <- records();if (is.null(x)) return(NULL)
    species <- input$typeOrganism;reference <- NULL
    if (length(species)==1L && species %in% c("Homo sapiens","Mus musculus"))
      reference <- tryCatch(load_gene_annotation_cache(species,app_root,include_ribosomal=FALSE,
        reference_id=references$configuration()$reference_id,cache_dir=references$cache_dir),error=function(e) NULL)
    lapply(x,suggest_feature_identifiers,reference=reference)
  })
  version_defaults <- shiny::reactive({
    x <- records();if (is.null(x)) return(NULL)
    formats <- lapply(x,suggest_feature_identifiers)
    list(id=paste0("annotationStripVersion_",substr(digest::digest(vapply(formats,`[[`,character(1),"token"),algo="sha256"),1,16)),
      value=any(vapply(formats,`[[`,logical(1),"strip_version_suggested")))
  })
  options <- shiny::reactive({
    x <- records(); if (is.null(x)) return(NULL)
    s <- suggestions()
    decisions <- lapply(seq_along(x),function(i) resolve_feature_identifier_choice(s[[i]],input[[paste0("annotationKey_",i)]]))
    keys <- stats::setNames(vapply(decisions,`[[`,character(1),"key"),names(x))
    modes <- stats::setNames(vapply(decisions,`[[`,character(1),"mode"),names(x))
    defaults <- version_defaults()
    # Read this control only. Enumerating names(input) would subscribe to every
    # new app input and clear prepared QC when unrelated panels are used.
    suffix_value <- input[[defaults$id]]
    strip <- if (!is.null(suffix_value)) isTRUE(suffix_value) else
      if (!is.null(input$annotationStripVersion)) isTRUE(input$annotationStripVersion) else defaults$value
    c(list(organism=input$typeOrganism,keys=keys,identifier_modes=modes,
      strip_ensembl_version=strip,suffix_mode=if (identical(strip,defaults$value)) "suggested_default" else "user_override"),references$configuration())
  })
  output$annotationIdControls <- shiny::renderUI({
    x <- records(); if (is.null(x)) return(NULL)
    s <- suggestions()
    shiny::tagList(lapply(seq_along(x),function(i) {
      choices <- c("Automatic detection"="auto","Gene symbols"=paste0("gene_name:",s[[i]]$token),
        "Gene IDs (Ensembl/reference IDs)"=paste0("gene_id:",s[[i]]$token))
      selected <- shiny::isolate(input[[paste0("annotationKey_",i)]])
      if (length(selected)==1L && selected %in% c("gene_id","gene_name")) selected <- paste(selected,s[[i]]$token,sep=":")
      if (length(selected)!=1L || !selected %in% choices) selected <- "auto"
      shiny::tagList(shiny::selectInput(paste0("annotationKey_",i),paste("Feature identifiers:",names(x)[i]),
        choices=choices,selected=selected),shiny::helpText(s[[i]]$reason))
    }))
  })
  output$annotationVersionControl <- shiny::renderUI({
    defaults <- version_defaults();if (is.null(defaults)) return(NULL)
    shiny::checkboxInput(defaults$id,"Ignore Ensembl version suffixes for lookup only (allows annotation matching when Ensembl IDs contain a version suffix)",defaults$value)
  })
  output$annotationIdentifierSummary <- shiny::renderTable({
    s <- suggestions();if (is.null(s)) return(NULL)
    cfg <- options()
    type_label <- function(key,missing) if (key=="gene_id") "Gene IDs" else if (key=="gene_name") "Gene symbols" else missing
    do.call(rbind,lapply(names(s),function(n) data.frame(sample=n,suggested_type=type_label(s[[n]]$suggested_type,"Needs review"),
      selected_type=type_label(cfg$keys[[n]],"Choose explicitly"),mode=cfg$identifier_modes[[n]],
      features=s[[n]]$features,ensembl_ids=s[[n]]$ensembl_features,versioned_ids=s[[n]]$versioned_features,
      exact_symbol_matches=s[[n]]$symbol_matches,suffix_lookup=if(cfg$keys[[n]]=="gene_id") {
        if(cfg$strip_ensembl_version) "Enabled" else "Disabled"
      } else if(cfg$keys[[n]]=="gene_name") "Not applicable" else "Choose identifier type")))
  })
  output$paramsAddMetadata <- shiny::renderText(status())
  output$initDataFilt <- shiny::renderText(baseline_status())
  output$annotationBusy <- shiny::renderUI({
    if (busy()) shiny::tags$div(role="status",shiny::icon("spinner",class="fa-spin")," Preparing annotations and QC metrics…")
  })
  output$annotationSummary <- shiny::renderTable({
    p <- prepared(); if (is.null(p)) return(NULL)
    do.call(rbind,lapply(p,`[[`,"summary"))
  },rownames=FALSE)
  output$ribosomalSummary <- shiny::renderTable({
    p <- prepared(); if (is.null(p)) return(NULL)
    do.call(rbind,lapply(p,`[[`,"ribosomal_summary"))
  },rownames=FALSE)
  output$annotationProvenance <- shiny::renderTable({
    p <- prepared(); if (is.null(p)) return(NULL)
    r <- p[[1]]$config$reference
    value <- function(field) {x<-r[[field]];if(is.null(x) || !length(x) || is.na(x)) "Unknown" else x}
    data.frame(species=p[[1]]$config$organism,cache=basename(r$path),md5=r$md5,
      resource_id=value("resource_id"),snapshot=value("snapshot_date"),
      annotation_release=value("annotation_release"),genome_build=value("genome_build"),
      source=r$source,reference_id=value("reference_id"),data_genome_build=p[[1]]$config$data_genome_build,
      supplement=if (is.null(r$supplement)) "None" else paste(r$supplement$source,
        "gene models",r$supplement$annotation_release,r$supplement$mart_version,"retrieved",r$supplement$retrieved_at),
      ribosomal_source=if (identical(p[[1]]$config$ribosomal_reference$status,"available")) {
        paste(p[[1]]$config$ribosomal_reference$package,p[[1]]$config$ribosomal_reference$package_version)
      } else "Protein reference unavailable; RNA biotypes only")
  },rownames=FALSE)
  output$annotationIssues <- shiny::renderTable({
    p <- prepared(); if (is.null(p)) return(NULL)
    rows <- lapply(names(p),function(n) {
      d <- p[[n]]$feature_annotation
      d <- d[d$annotation_status!="matched" | d$ribosomal_annotation_status=="unresolved",,drop=FALSE]
      cbind(sample=rep(n,nrow(d)),d[,c("source_feature_id","annotation_status","annotation_candidates","ribosomal_annotation_status"),drop=FALSE])
    })
    utils::head(do.call(rbind,rows),20L)
  },rownames=FALSE)
  output$annotationMappingDownload <- shiny::downloadHandler(filename=function() "gene-annotation-mapping.csv",content=function(path) {
    p <- shiny::req(prepared()); x <- shiny::req(records())
    d <- do.call(rbind,lapply(names(p),function(n) data.frame(sample=n,
      internal_feature_id=rownames(x[[n]]$counts),p[[n]]$feature_annotation,check.names=FALSE)))
    utils::write.csv(d,path,row.names=FALSE)
  })
  output$qcMetricsDownload <- shiny::downloadHandler(filename=function() "qc-metrics-before-initial-filter.csv",content=function(path) {
    p <- shiny::req(prepared())
    utils::write.csv(do.call(rbind,lapply(p,`[[`,"metrics")),path,row.names=FALSE)
  })
  output$baselineSummary <- shiny::renderTable({
    b <- baseline(); if (is.null(b)) return(NULL)
    do.call(rbind,lapply(names(b),function(n) {
      info <- S4Vectors::metadata(b[[n]])$scrna_explorer_qc$initial_filter
      data.frame(sample=n,input_cells=info$input_cells,baseline_cells=info$retained_cells,
        excluded_below_100=length(info$excluded_cells))
    }))
  },rownames=FALSE)
  shiny::observeEvent(list(records(),options()),{
    x <- records();cfg <- options()
    # Initial browser bindings can report the same automatic defaults again.
    # Preserve prepared results unless the actual import/configuration changes.
    if (configured && identical(x,previous_records) && identical(cfg,previous_options)) return()
    configured <<- TRUE;previous_records <<- x;previous_options <<- cfg
    revision <<- revision+1L
    prepared(NULL); baseline(NULL); busy(FALSE); baseline_status("")
    status(if (is.null(records())) "Read data before adding annotations." else "Annotation settings not yet applied. Choose species and identifier types, then add annotations.")
  },ignoreNULL=FALSE,priority=50)
  shiny::observeEvent(input$addMetadata,{
    revision <<- revision+1L
    token <- revision
    x <- records(); cfg <- options()
    prepared(NULL); baseline(NULL); baseline_status("")
    if (is.null(x)) {status("Read and review data first.");return()}
    if (length(cfg$organism)!=1L || !cfg$organism %in% c("Homo sapiens","Mus musculus") || any(!cfg$keys %in% c("gene_id","gene_name"))) {
      status("Select a species and an identifier type for every sample.");return()
    }
    busy(TRUE); status("Preparing cached gene annotations and QC metrics…")
    session$onFlushed(function() shiny::isolate({
      if (session$isClosed() || token!=revision) return(invisible(NULL))
      tryCatch({
        reference <- load_gene_annotation_cache(cfg$organism,app_root,reference_id=cfg$reference_id,
          cache_dir=references$cache_dir)
        known_build <- reference$provenance$genome_build
        if (nzchar(cfg$data_genome_build) && !is.na(known_build) &&
            !identical(known_build,cfg$data_genome_build)) stop("Selected reference build differs from the declared counts build. Choose a compatible reference.")
        if (cfg$reference_id!="bundled" && !nzchar(cfg$data_genome_build))
          stop("Declare the counts genome build before using a downloaded reference.")
        result <- stats::setNames(lapply(names(x),function(n)
          prepare_sample_qc(x[[n]],reference,cfg$keys[[n]],cfg$strip_ensembl_version)),names(x))
        for (n in names(result)) {
          result[[n]]$config$data_genome_build <- if (nzchar(cfg$data_genome_build)) cfg$data_genome_build else "Unknown"
          result[[n]]$config$identifier_decision <- c(suggestions()[[n]],list(mode=cfg$identifier_modes[[n]],selected_type=cfg$keys[[n]],suffix_mode=cfg$suffix_mode))
        }
        if (token!=revision) return(invisible(NULL))
        prepared(result)
        missing_hemoglobin <- any(vapply(result,function(p) identical(p$summary$hemoglobin_status,"missing_reference_definition"),logical(1)))
        status(paste("Annotations and QC metrics prepared. Review mapping coverage and reference provenance, then run Minimal filtering (UMIs >= 100). Original feature IDs and counts are preserved; unmapped/ambiguous features were not removed.",
          if (missing_hemoglobin) "Hemoglobin QC is unavailable: this reference lacks the description-based gene definition. QC filtering remains blocked until that evidence is supplied; unavailable values are not zero contamination." else ""))
      },error=function(e) {if(token==revision) status(paste("Annotation stopped:",conditionMessage(e)))})
      if (token==revision) busy(FALSE)
    }),once=TRUE)
  },ignoreInit=TRUE)
  shiny::observeEvent(input$initFiltering,{
    x <- records(); p <- prepared()
    if (is.null(x) || is.null(p)) {baseline_status("Prepare and review gene annotations first.");return()}
    tryCatch({
      b <- stats::setNames(lapply(names(x),function(n) make_baseline_sce(x[[n]],p[[n]])),names(x))
      baseline(b)
      empty <- names(b)[vapply(b,ncol,integer(1))==0L]
      baseline_status(paste("Minimal filtering completed: retained UMIs >= 100.",
        if(length(empty)) paste("No baseline cells in:",paste(empty,collapse=", ")) else ""))
    },error=function(e) {baseline(NULL);baseline_status(paste("Baseline construction stopped:",conditionMessage(e)))})
  },ignoreInit=TRUE,priority=10)
  list(prepared=shiny::reactive(prepared()),baseline=shiny::reactive(baseline()),
    status=shiny::reactive(status()),baseline_status=shiny::reactive(baseline_status()),
    identifier_suggestions=suggestions,identifier_options=options,version_control_id=shiny::reactive(version_defaults()$id))
}
